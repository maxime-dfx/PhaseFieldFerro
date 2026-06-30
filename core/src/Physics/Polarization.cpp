#include "Physics/Mechanics.h"
#include "Physics/Electrostatics.h"
#include "Physics/Polarization.h"
#include "Core/ShapeFunctions.h"
#include "Utils/Logger.h"
#include "Physics/Fracture.h"
#include "Physics/Math.h"
#include "Core/Quadrature.h"
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

Polarization::Polarization(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager) : config(config), mesh(mesh), bc_manager(bc_manager) {
    size_t n_nodes = mesh.get_num_nodes();
    Px_prev_iter.setZero(n_nodes);
    Py_prev_iter.setZero(n_nodes);
    Px_current.resize(n_nodes);
    Px_current.setZero();
    Py_current.resize(n_nodes);
    Py_current.setZero();
    Px_n.setZero(n_nodes);
    Py_n.setZero(n_nodes);

    if (config.get_initial_polarization() == PolarizationInitializationType::UNIFORM) {
        Px_0 = config.get_Px_0();
        Py_0 = config.get_Py_0();
        for (size_t i = 0; i < Px_current.size(); ++i) {
            Px_current[i] = Px_0;
            Py_current[i] = Py_0;
        }
    } else if (config.get_initial_polarization() == PolarizationInitializationType::RANDOM) {
        for (size_t i = 0; i < Px_current.size(); ++i) {
            Px_current[i] = static_cast<double>(rand()) / RAND_MAX;
            Py_current[i] = static_cast<double>(rand()) / RAND_MAX;
        }
    } else {
        Logger::error("Type d'initialisation de polarisation non reconnu.");
    }

    // Securite : Px_n/Py_n synchronises avec l'etat initial, au cas ou
    // freeze_time_step() ne serait pas encore appele avant le tout premier pas.
    Px_n = Px_current;
    Py_n = Py_current;
}

// ------------------------------------------------------------
// Conditions aux limites : methode par penalite (preserve la symetrie,
// compatible avec SimplicialLDLT). Remplace l'ancienne elimination de ligne
// qui cassait la symetrie et faussait la factorisation.
// ------------------------------------------------------------
void Polarization::apply_boundary_conditions(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const std::vector<NodeBC>& bcs) {
    double max_diag = 0.0;
    for (int k = 0; k < A.outerSize(); ++k)
        max_diag = std::max(max_diag, std::abs(A.coeff(k, k)));
    const double penalty = max_diag * 1e7;

    for (size_t i = 0; i < bcs.size(); ++i) {
        if (bcs[i].type == BCType::DIRICHLET) {
            int dof = static_cast<int>(i);
            A.coeffRef(dof, dof) += penalty;
            b(dof) += penalty * bcs[i].value;
        }
        // NEUMANN : rien a faire, condition naturelle (flux nul)
    }

    A.makeCompressed();
}

void Polarization::update_polarization_component(double time, const Fracture& fracture, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math, int component) {
    double dt = config.get_dt();
    const int n_dof = static_cast<int>(Px_current.size());

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.get_elements().size() * 16);

    Eigen::VectorXd b;
    b.setZero(n_dof);

    // Logger::debug("[Polarization] update_polarization_component comp=" + std::to_string(component) +
    //               " mu_p=" + std::to_string(math.mu_p) +
    //               " a0=" + std::to_string(math.a0) +
    //               " eta_k=" + std::to_string(math.eta_k) +
    //               " dt=" + std::to_string(dt), config.debug_enabled());
    Eigen::VectorXd diag_check = Eigen::VectorXd::Zero(n_dof);

    const auto& elements = mesh.get_elements();

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);

        Eigen::MatrixXd local_A = Eigen::MatrixXd::Zero(n_nodes, n_nodes);
        Eigen::VectorXd local_b = Eigen::VectorXd::Zero(n_nodes);
        Eigen::VectorXd Px_local = Eigen::VectorXd::Zero(n_nodes);
        Eigen::VectorXd Py_local = Eigen::VectorXd::Zero(n_nodes);
        Eigen::VectorXd v_local  = Eigen::VectorXd::Zero(n_nodes);
        // Etat fige P_n (debut du pas de temps physique), utilise UNIQUEMENT
        // dans le terme de masse implicite. Px_local/Py_local ci-dessus restent
        // le dernier itere de Picard et continuent d'alimenter les termes non
        // lineaires (Pi_gp, dW, dchi) -- c'est volontaire (schema IMEX).
        Eigen::VectorXd Px_local_n = Eigen::VectorXd::Zero(n_nodes);
        Eigen::VectorXd Py_local_n = Eigen::VectorXd::Zero(n_nodes);

        const auto& indices = elem.get_node_indices();
        for (int i = 0; i < n_nodes; ++i) {
            Px_local(i) = Px_current[indices[i]];
            Py_local(i) = Py_current[indices[i]];
            v_local(i)  = fracture.get_v()[indices[i]];
            Px_local_n(i) = Px_n[indices[i]];
            Py_local_n(i) = Py_n[indices[i]];
        }

        // Lambda commune d'assemblage pour un point d'integration (T3 ou Q4)
        auto assemble_gp = [&](const Eigen::RowVectorXd& N,
                                const Eigen::MatrixXd& grad_N,
                                double dV,
                                const GaussPoint2D& gp_for_coupling) {
            double v_gp = N.dot(v_local);
            Eigen::Vector2d Pi_gp(N.dot(Px_local), N.dot(Py_local));
            Eigen::Matrix2d eps_gp = mechanics.get_strain_at_gp(elem, gp_for_coupling, coords);

            double E_gp = 0.0, dW = 0.0, dchi = 0.0;
            const Eigen::VectorXd& P_comp_local = (component == 0) ? Px_local : Py_local;
            // Terme source explicite du couplage temporel : doit utiliser P_n
            // (etat fige au debut du pas de temps), PAS le dernier itere Picard.
            const Eigen::VectorXd& P_comp_local_n = (component == 0) ? Px_local_n : Py_local_n;

            if (component == 0) {
                E_gp = electrostatics.get_Ex_at_gp(elem, gp_for_coupling);
                dW   = math.dW_dp1(Pi_gp, eps_gp);
                dchi = math.dchi_dp1(Pi_gp);
            } else {
                E_gp = electrostatics.get_Ey_at_gp(elem, gp_for_coupling);
                dW   = math.dW_dp2(Pi_gp, eps_gp);
                dchi = math.dchi_dp2(Pi_gp);
            }

            if (!std::isfinite(E_gp) || !std::isfinite(dW) || !std::isfinite(dchi)) {
                Logger::debug("[Polarization] elem_idx=" + std::to_string(elem_idx) +
                              " gp=(" + std::to_string(gp_for_coupling.xi) + "," + std::to_string(gp_for_coupling.eta) + ")" +
                              " v_gp=" + std::to_string(v_gp) +
                              " Pi_gp=(" + std::to_string(Pi_gp[0]) + "," + std::to_string(Pi_gp[1]) + ")" +
                              " eps_gp=[" + std::to_string(eps_gp(0,0)) + "," + std::to_string(eps_gp(0,1)) + ";" +
                                           std::to_string(eps_gp(1,0)) + "," + std::to_string(eps_gp(1,1)) + "]", config.debug_enabled());
                Logger::debug("[Polarization] elem_idx=" + std::to_string(elem_idx) +
                              " gp=(" + std::to_string(gp_for_coupling.xi) + "," + std::to_string(gp_for_coupling.eta) + ")" +
                              " v_gp=" + std::to_string(v_gp) +
                              " Pi_gp=(" + std::to_string(Pi_gp[0]) + "," + std::to_string(Pi_gp[1]) + ")" +
                              " eps_gp=[" + std::to_string(eps_gp(0,0)) + "," + std::to_string(eps_gp(0,1)) + ";" +
                                           std::to_string(eps_gp(1,0)) + "," + std::to_string(eps_gp(1,1)) + "]", config.debug_enabled());
                return;
            }

            double penalite_fracture = (v_gp * v_gp) + math.eta_k;
            double force_thermo = 0.0;

            if (config.get_fracture_mode() == CrackBCType::PERMEABLE) {
                force_thermo = - (penalite_fracture * (dW + dchi) - E_gp);
            } else if (config.get_fracture_mode() == CrackBCType::IMPERMEABLE) {
                force_thermo = - (penalite_fracture * (dW + dchi - E_gp));
            }

            local_A.noalias() += (math.mu_p / dt) * (N.transpose() * N * dV)
                               + (math.a0 * penalite_fracture) * (grad_N.transpose() * grad_N * dV);

            local_b.noalias() += (math.mu_p / dt) * (N.transpose() * N * dV) * P_comp_local_n
                               + N.transpose() * force_thermo * dV;
        };

        if (n_nodes == 3) {
            // --- TRIANGLE T3 : un seul point d'integration (CST) ---
            double xi = 1.0 / 3.0, eta = 1.0 / 3.0;
            auto N_std = ShapeFunctions::get_shape_functions_tri(xi, eta);
            auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);

            if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) {
                Logger::debug("[Polarization] elem_idx=" + std::to_string(elem_idx) +
                              " gp=(" + std::to_string(xi) + "," + std::to_string(eta) + ")" +
                              " detJ=" + std::to_string(detJ), config.debug_enabled());
                continue;
            }

            double dV = std::abs(detJ) / 2.0;

            Eigen::RowVectorXd N(3);
            N << N_std[0], N_std[1], N_std[2];

            Eigen::MatrixXd grad_N(2, 3);
            for (int i = 0; i < 3; ++i) {
                grad_N(0, i) = dN_xy[i][0];
                grad_N(1, i) = dN_xy[i][1];
            }

            GaussPoint2D gp_fake{xi, eta, 1.0};
            assemble_gp(N, grad_N, dV, gp_fake);

        } else if (n_nodes == 4) {
            // --- QUAD Q4 : integration 2x2 ---
            const auto& gauss_points = Quadrature::get_gauss_2x2();
            for (const auto& gp : gauss_points) {
                auto N_std = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
                auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
                auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);

                if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) {
                    Logger::debug("[Polarization] elem_idx=" + std::to_string(elem_idx) +
                                  " gp=(" + std::to_string(gp.xi) + "," + std::to_string(gp.eta) + ")" +
                                  " detJ=" + std::to_string(detJ), config.debug_enabled());
                    continue;
                }

                double dV = gp.weight * std::abs(detJ);

                Eigen::RowVectorXd N(4);
                N << N_std[0], N_std[1], N_std[2], N_std[3];

                Eigen::MatrixXd grad_N(2, 4);
                for (int i = 0; i < 4; ++i) {
                    grad_N(0, i) = dN_xy[0][i];
                    grad_N(1, i) = dN_xy[1][i];
                }

                assemble_gp(N, grad_N, dV, gp);
            }
        } else {
            Logger::debug("[Polarization] Type d'element non supporte: n_nodes=" + std::to_string(n_nodes) +
                          " idx=" + std::to_string(elem_idx) + "\n", config.debug_enabled());
            continue;
        }

        // Assemblage global -> triplets sparse, + accumulation de la diagonale
        for (int i = 0; i < n_nodes; ++i) {
            for (int j = 0; j < n_nodes; ++j) {
                double val = local_A(i, j);
                if (val != 0.0) {
                    triplets.emplace_back(indices[i], indices[j], val);
                    if (i == j) {
                        diag_check(indices[i]) += val;
                    }
                }
            }
            b(indices[i]) += local_b(i);
        }
    }

    // --- Garde-fou : DOF sans contribution diagonale (noeud isole ou zone
    // totalement degradee par la fracture) -> on impose une diagonale unite
    // pour eviter une matrice singuliere, sans casser la symetrie. ---
    const auto& bcs = (component == 0) ? bc_manager.get_px_bcs() : bc_manager.get_py_bcs();
    for (int k = 0; k < n_dof; ++k) {
        bool is_dirichlet = (static_cast<size_t>(k) < bcs.size()) && (bcs[k].type == BCType::DIRICHLET);
        if (!is_dirichlet && std::abs(diag_check(k)) < 1e-12) {
            triplets.emplace_back(k, k, 1.0);
            b(k) = 0.0;
            Logger::debug("[Polarization] DOF flottant detecte idx=" + std::to_string(k) +
                          " (comp=" + std::to_string(component) + ") -> diagonale forcee a 1.0\n", config.debug_enabled());
        }
    }

    // Construction de la matrice sparse globale
    Eigen::SparseMatrix<double> A_sparse(n_dof, n_dof);
    A_sparse.setFromTriplets(triplets.begin(), triplets.end());

    // --- APPLICATION DES CL DIRICHLET ET RÉSOLUTION ---
    apply_boundary_conditions(A_sparse, b, bcs);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A_sparse);

    if (solver.info() != Eigen::Success) {
        Logger::debug("[Polarization] Echec factorisation (comp=" + std::to_string(component) +
                      ") - matrice singuliere ou mal conditionnee\n", config.debug_enabled());
        return;
    }

    Eigen::VectorXd solution = solver.solve(b);
    Logger::debug("[Polarization] comp=" + std::to_string(component) +
                  " solution min=" + std::to_string(solution.minCoeff()) +
                  " max=" + std::to_string(solution.maxCoeff()) + "\n", config.debug_enabled());

    if (solver.info() != Eigen::Success || !solution.allFinite()) {
        Logger::debug("[Polarization] Echec resolution ou solution non finie (comp=" + std::to_string(component) + ")\n", config.debug_enabled());
        return;
    }

    if (component == 0) {
        Px_current = solution;
    } else {
        Py_current = solution;
    }
}