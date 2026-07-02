#include "Physics/Mechanics.h"
#include "Core/Quadrature.h"
#include "Physics/Math.h"
#include "Physics/Polarization.h"
#include "Physics/Fracture.h"
#include "Core/ShapeFunctions.h"
#include <cmath>
#include <iostream>

Mechanics::Mechanics(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager)
    : config(config), mesh(mesh), bc_manager(bc_manager) {
    size_t n_nodes = mesh.get_num_nodes();
    ux_prev_iter.setZero(n_nodes);
    uy_prev_iter.setZero(n_nodes);
    ux_current.resize(n_nodes);
    ux_current.setZero();
    uy_current.resize(n_nodes);
    uy_current.setZero();
    if (config.get_initial_mechanics() == MechanicsInitializationType::UNIFORM) {
        ux_0 = config.get_ux_0();
        uy_0 = config.get_uy_0();
        for (size_t i = 0; i < ux_current.size(); ++i) {
            ux_current[i] = ux_0;
            uy_current[i] = uy_0;
        }
    } else if (config.get_initial_mechanics() == MechanicsInitializationType::RANDOM) {
        for (size_t i = 0; i < ux_current.size(); ++i) {
            ux_current[i] = static_cast<double>(rand()) / RAND_MAX;
            uy_current[i] = static_cast<double>(rand()) / RAND_MAX;
        }
    } else {
        Logger::error("Type d'initialisation de mechanics non reconnu.");
    }
}

void Mechanics::update_u(double time, const Polarization& polarization, const Fracture& fracture, const Math& math) {
    size_t n_nodes = mesh.get_num_nodes();
    size_t system_size = 2 * n_nodes;

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.get_elements().size() * 64);

    F_global.setZero(system_size);

    const auto& elements = mesh.get_elements();
    Logger::debug("[DEBUG][Mechanics] c1=" + std::to_string(math.c1) + " c2=" + std::to_string(math.c2) + " c3=" + std::to_string(math.c3) +
                  " b1=" + std::to_string(math.b1) + " b2=" + std::to_string(math.b2) + " b3=" + std::to_string(math.b3) +
                  " eta_k=" + std::to_string(math.eta_k) + "\n", config.debug_enabled());
    // Matrice de rigidite (Loi de Hooke)
    Eigen::Matrix3d C;
    C << math.c1, math.c2, 0,
        math.c2, math.c1, 0,
        0, 0, 0.5 * math.c3;

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_elem_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);

        Eigen::MatrixXd K_local = Eigen::MatrixXd::Zero(2 * n_elem_nodes, 2 * n_elem_nodes);
        Eigen::VectorXd F_local = Eigen::VectorXd::Zero(2 * n_elem_nodes);

        const auto& indices = elem.get_node_indices();

        Eigen::VectorXd Px_local(n_elem_nodes), Py_local(n_elem_nodes), v_local(n_elem_nodes);
        for (int i = 0; i < n_elem_nodes; ++i) {
            Px_local(i) = polarization.get_Px()[indices[i]];
            Py_local(i) = polarization.get_Py()[indices[i]];
            v_local(i)  = fracture.get_v()[indices[i]];
        }

        // Lambda commune pour assembler la contribution d'un point d'integration,
        // quel que soit le type d'element (evite la duplication entre T3 et Q4)
        auto assemble_gp = [&](const Eigen::RowVectorXd& N,
                                const Eigen::MatrixXd& grad_N,
                                double dV) {
            int nn = N.size();
            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2 * nn);
            for (int i = 0; i < nn; ++i) {
                B(0, 2 * i)     = grad_N(0, i); // dNi/dx
                B(1, 2 * i + 1) = grad_N(1, i); // dNi/dy
                B(2, 2 * i)     = grad_N(1, i); // dNi/dy
                B(2, 2 * i + 1) = grad_N(0, i); // dNi/dx
            }

            double p1 = N.dot(Px_local);
            double p2 = N.dot(Py_local);
            double v_gp = N.dot(v_local);

            if (!std::isfinite(p1) || !std::isfinite(p2) || !std::isfinite(v_gp)) {
                Logger::debug("[Mechanics] Champ non fini p1=" + std::to_string(p1) + " p2=" + std::to_string(p2) +
                              " v_gp=" + std::to_string(v_gp) + " -> contribution ignoree\n", config.debug_enabled());
                return;
            }

            double penalite_fracture = (v_gp * v_gp) + math.eta_k;

            Eigen::Vector3d sigma_0;
            sigma_0(0) = - (math.b1 / 2.0) * p1 * p1 - (math.b2 / 2.0) * p2 * p2;
            sigma_0(1) = - (math.b2 / 2.0) * p1 * p1 - (math.b1 / 2.0) * p2 * p2;
            sigma_0(2) = - math.b3 * p1 * p2;

            K_local.noalias() += B.transpose() * (penalite_fracture * C) * B * dV;
            F_local.noalias() += B.transpose() * (penalite_fracture * sigma_0) * dV;
        };

        if (n_elem_nodes == 3) {
            // --- TRIANGLE T3 : un seul point d'integration (CST, gradients constants) ---
            double xi = 1.0 / 3.0, eta = 1.0 / 3.0;
            auto N_std = ShapeFunctions::get_shape_functions_tri(xi, eta);
            auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);

            if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) {
                Logger::debug("[Mechanics] Triangle degenere idx=" + std::to_string(elem_idx) +
                              " detJ=" + std::to_string(detJ) + " -> contribution ignoree\n", config.debug_enabled());
                continue;
            }

            double dV = std::abs(detJ) / 2.0; // aire du triangle

            Eigen::RowVectorXd N(3);
            N << N_std[0], N_std[1], N_std[2];

            Eigen::MatrixXd grad_N(2, 3);
            for (int i = 0; i < 3; ++i) {
                grad_N(0, i) = dN_xy[i][0];
                grad_N(1, i) = dN_xy[i][1];
            }

            assemble_gp(N, grad_N, dV);

        } else if (n_elem_nodes == 4) {
            // --- QUAD Q4 : integration 2x2 ---
            const auto& gauss_points = Quadrature::get_gauss_2x2();
            for (const auto& gp : gauss_points) {
                auto N_std = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
                auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
                auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);

                if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) {
                    Logger::debug("[Mechanics] Element degenere idx=" + std::to_string(elem_idx) +
                                  " detJ=" + std::to_string(detJ) + " -> contribution ignoree\n", config.debug_enabled());
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

                assemble_gp(N, grad_N, dV);
            }
        } else {
            Logger::debug("[Mechanics] Type d'element non supporte: n_nodes=" + std::to_string(n_elem_nodes) +
                          " idx=" + std::to_string(elem_idx) + "\n", config.debug_enabled());
            continue;
        }

        // Assemblage global (format Triplet) -- generique en n_elem_nodes
        for (int i = 0; i < n_elem_nodes; ++i) {
            for (int j = 0; j < n_elem_nodes; ++j) {
                triplets.emplace_back(2 * indices[i],     2 * indices[j],     K_local(2 * i,     2 * j));
                triplets.emplace_back(2 * indices[i],     2 * indices[j] + 1, K_local(2 * i,     2 * j + 1));
                triplets.emplace_back(2 * indices[i] + 1, 2 * indices[j],     K_local(2 * i + 1, 2 * j));
                triplets.emplace_back(2 * indices[i] + 1, 2 * indices[j] + 1, K_local(2 * i + 1, 2 * j + 1));
            }
            F_global(2 * indices[i])     += F_local(2 * i);
            F_global(2 * indices[i] + 1) += F_local(2 * i + 1);
        }
    }

    K_global.resize(system_size, system_size);
    K_global.setFromTriplets(triplets.begin(), triplets.end());

    const auto& bcs_ux = bc_manager.get_ux_bcs();
    const auto& bcs_uy = bc_manager.get_uy_bcs();
    apply_boundary_conditions(K_global, F_global, bcs_ux, bcs_uy);

    // --- Solveur persistant : analyzePattern() une seule fois ---
    if (!pattern_analyzed_) {
        solver_.analyzePattern(K_global);
        pattern_analyzed_ = true;
    }
    solver_.factorize(K_global);

    if (solver_.info() != Eigen::Success) {
        Logger::debug("[Mechanics] Echec factorisation - matrice singuliere ou mal conditionnee\n", config.debug_enabled());
        return;
    }

    Eigen::VectorXd U_new = solver_.solve(F_global);
    Logger::debug("[DEBUG][Mechanics] solver.info=" + std::to_string(solver_.info()) + " U min=" + std::to_string(U_new.minCoeff()) + " max=" + std::to_string(U_new.maxCoeff()) + "\n", config.debug_enabled());
    if (solver_.info() != Eigen::Success || !U_new.allFinite()) {
        Logger::debug("[Mechanics] Echec resolution ou solution non finie\n", config.debug_enabled());
        return;
    }

    for (size_t i = 0; i < n_nodes; ++i) {
        ux_current(i) = U_new(2 * i);
        uy_current(i) = U_new(2 * i + 1);
    }
}

Eigen::Matrix2d Mechanics::get_strain_at_gp(const Element& elem, const GaussPoint2D& gp, const std::vector<std::array<double, 2>>& coords) const {
    const auto& indices = elem.get_node_indices();
    int n_nodes = elem.get_num_nodes();

    double eps_11 = 0.0, eps_22 = 0.0, eps_12 = 0.0;

    if (n_nodes == 3) {
        // Triangle : gradients constants, le "gp" passe n'est pas utilise (CST)
        auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);
        if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) {
            Logger::debug("[Mechanics] get_strain_at_gp: triangle degenere detJ=" + std::to_string(detJ) + "\n", config.debug_enabled());
            return Eigen::Matrix2d::Zero();
        }
        for (int i = 0; i < 3; ++i) {
            double dN_dx = dN_xy[i][0];
            double dN_dy = dN_xy[i][1];
            eps_11 += dN_dx * ux_current[indices[i]];
            eps_22 += dN_dy * uy_current[indices[i]];
            eps_12 += dN_dy * ux_current[indices[i]] + dN_dx * uy_current[indices[i]];
        }
    } else if (n_nodes == 4) {
        auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
        auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);
        if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) {
            Logger::debug("[Mechanics] get_strain_at_gp: element degenere detJ=" + std::to_string(detJ) + "\n", config.debug_enabled());
            return Eigen::Matrix2d::Zero();
        }
        for (int i = 0; i < n_nodes; ++i) {
            double dN_dx = dN_xy[0][i];
            double dN_dy = dN_xy[1][i];
            eps_11 += dN_dx * ux_current[indices[i]];
            eps_22 += dN_dy * uy_current[indices[i]];
            eps_12 += dN_dy * ux_current[indices[i]] + dN_dx * uy_current[indices[i]];
        }
    } else {
        Logger::debug("[Mechanics] get_strain_at_gp: type d'element non supporte n_nodes=" + std::to_string(n_nodes) + "\n", config.debug_enabled());
        return Eigen::Matrix2d::Zero();
    }

    eps_12 *= 0.5;

    Eigen::Matrix2d strain;
    strain << eps_11, eps_12,
              eps_12, eps_22;
    return strain;
}

void Mechanics::apply_boundary_conditions(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F,
                                          const std::vector<NodeBC>& bcs_x, const std::vector<NodeBC>& bcs_y) {
    double max_diag = 0.0;
    for (int k = 0; k < K.outerSize(); ++k)
        max_diag = std::max(max_diag, std::abs(K.coeff(k, k)));
    const double penalty = max_diag * 1e7;

    Logger::debug("[DEBUG][BC phi] n_nodes=" + std::to_string(mesh.get_num_nodes()) +
                  " bcs.size()=" + std::to_string(bcs_x.size()) + "\n", config.debug_enabled());
    int n_dirichlet = 0;
    for (size_t i = 0; i < bcs_x.size(); ++i)
        if (bcs_x[i].type == BCType::DIRICHLET) n_dirichlet++;
    Logger::debug("[DEBUG][BC phi] n_dirichlet=" + std::to_string(n_dirichlet) + "\n", config.debug_enabled());

    for (size_t i = 0; i < bcs_x.size(); ++i) {
        if (bcs_x[i].type == BCType::DIRICHLET) {
            int dof = 2 * static_cast<int>(i);
            K.coeffRef(dof, dof) += penalty;
            F(dof) += penalty * bcs_x[i].value;
        }
    }

    for (size_t i = 0; i < bcs_y.size(); ++i) {
        if (bcs_y[i].type == BCType::DIRICHLET) {
            int dof = 2 * static_cast<int>(i) + 1;
            K.coeffRef(dof, dof) += penalty;
            F(dof) += penalty * bcs_y[i].value;
        }
    }

    int nd_x = 0, nd_y = 0;
    for (auto& bc : bcs_x) if (bc.type == BCType::DIRICHLET) nd_x++;
    for (auto& bc : bcs_y) if (bc.type == BCType::DIRICHLET) nd_y++;
    Logger::debug("[DEBUG][Mechanics BC] n_dirichlet_ux=" + std::to_string(nd_x) + " n_dirichlet_uy=" + std::to_string(nd_y) + "\n", config.debug_enabled());
    K.makeCompressed();
}