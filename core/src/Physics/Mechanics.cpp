#include "Physics/Mechanics.h"
#include "Core/Quadrature.h"
#include "Physics/Math.h"
#include "Physics/Polarization.h"
#include "Physics/Fracture.h"
#include "Core/ShapeFunctions.h"
<<<<<<< HEAD
#include "Utils/Logger.h"
#include "Core/ElementIntegrator.h"
=======
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
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
<<<<<<< HEAD

=======
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
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

<<<<<<< HEAD
// =========================================================================
// 1. LE CHEF D'ORCHESTRE ( update_u )
// =========================================================================
=======
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
void Mechanics::update_u(double time, const Polarization& polarization, const Fracture& fracture, const Math& math) {
    size_t n_nodes = mesh.get_num_nodes();
    size_t system_size = 2 * n_nodes;

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.get_elements().size() * 64);
<<<<<<< HEAD
    F_global.setZero(system_size);

    // 1. Assemblage global (Matrice K et vecteur Force F_global)
    assemble_mechanical_system(polarization, fracture, math, triplets);

    Eigen::SparseMatrix<double> K_global_sparse(system_size, system_size);
    K_global_sparse.setFromTriplets(triplets.begin(), triplets.end());

    // 2. Application des Conditions aux Limites
    apply_boundary_conditions(K_global_sparse, F_global, bc_manager.get_ux_bcs(), bc_manager.get_uy_bcs());

    // 3. Résolution du système
    // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    // solver.compute(K_global_sparse);
    // if (solver.info() != Eigen::Success) {
    //     Logger::debug("[Mechanics] Echec factorisation - matrice singuliere ou mal conditionnee\n", config.debug_enabled());
    //     return;
    // }
    // Eigen::VectorXd U_new = solver.solve(F_global);
    Eigen::VectorXd U_guess(system_size);
    for (size_t i = 0; i < n_nodes; ++i) {
        U_guess(2 * i)     = ux_current(i);
        U_guess(2 * i + 1) = uy_current(i);
    }

    // 2. Configuration du Solveur Itératif (Gradient Conjugué + Incomplete Cholesky)
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> solver;
    solver.setTolerance(1e-10); // Tolérance stricte pour préserver la physique
    solver.setMaxIterations(1000);
    solver.compute(K_global_sparse);

    // 3. Résolution avec Démarrage à chaud
    Eigen::VectorXd U_new = solver.solveWithGuess(F_global, U_guess);

    // 4. Fallback de sécurité (Si la matrice est trop mal conditionnée)
    if (solver.info() != Eigen::Success) {
        Logger::debug("[Mechanics] CG a echoue (iters: " + std::to_string(solver.iterations()) + "). Bascule sur LDLT Direct.", config.debug_enabled());
        
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> direct_solver;
        direct_solver.compute(K_global_sparse);
        U_new = direct_solver.solve(F_global);
    } else {
        Logger::debug("[Mechanics] CG converge en " + std::to_string(solver.iterations()) + " iterations.", config.debug_enabled());
    }
    Logger::debug("[DEBUG][Mechanics] solver.info=" + std::to_string(solver.info()) + 
                  " U min=" + std::to_string(U_new.minCoeff()) + 
                  " max=" + std::to_string(U_new.maxCoeff()) + "\n", config.debug_enabled());

    // 4. Extraction et séparation des composantes (ux, uy)
    if (solver.info() == Eigen::Success && U_new.allFinite()) {
        map_global_vector_to_components(U_new);
    } else {
        Logger::debug("[Mechanics] Echec resolution ou solution non finie\n", config.debug_enabled());
    }
}

// =========================================================================
// 2. L'ASSEMBLEUR GLOBAL
// =========================================================================
void Mechanics::assemble_mechanical_system(const Polarization& polarization, const Fracture& fracture, const Math& math, std::vector<Eigen::Triplet<double>>& triplets) {
    const auto& elements = mesh.get_elements();
    Eigen::Matrix3d C = math.get_elastic_matrix(); 

    // 1. Allocations hors boucle
    const int MAX_NODES = 4;
    Eigen::MatrixXd K_local(2 * MAX_NODES, 2 * MAX_NODES);
    Eigen::VectorXd F_local(2 * MAX_NODES);
    Eigen::VectorXd Px_local(MAX_NODES), Py_local(MAX_NODES), v_local(MAX_NODES);
    
    // Tampons géométriques pour l'intégrateur
    Eigen::RowVectorXd N(MAX_NODES);
    Eigen::MatrixXd grad_N(2, MAX_NODES);
    Eigen::MatrixXd B(3, 2 * MAX_NODES);
=======

    F_global.setZero(system_size);

    const auto& elements = mesh.get_elements();
    Logger::debug("[DEBUG][Mechanics] c1=" + std::to_string(math.c1) + " c2=" + std::to_string(math.c2) + " c3=" + std::to_string(math.c3) +
                  " b1=" + std::to_string(math.b1) + " b2=" + std::to_string(math.b2) + " b3=" + std::to_string(math.b3) +
                  " eta_k=" + std::to_string(math.eta_k) + "\n", config.debug_enabled());
    // Matrice de rigidité (Loi de Hooke)
    Eigen::Matrix3d C;
    C << 2.0 * math.c1, math.c2, 0,
         math.c2, 2.0 * math.c1, 0,
         0, 0, math.c3;
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_elem_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
<<<<<<< HEAD
        const auto& indices = elem.get_node_indices();

        K_local.topLeftCorner(2 * n_elem_nodes, 2 * n_elem_nodes).setZero();
        F_local.head(2 * n_elem_nodes).setZero();
        
=======

        Eigen::MatrixXd K_local = Eigen::MatrixXd::Zero(2 * n_elem_nodes, 2 * n_elem_nodes);
        Eigen::VectorXd F_local = Eigen::VectorXd::Zero(2 * n_elem_nodes);

        const auto& indices = elem.get_node_indices();

        Eigen::VectorXd Px_local(n_elem_nodes), Py_local(n_elem_nodes), v_local(n_elem_nodes);
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
        for (int i = 0; i < n_elem_nodes; ++i) {
            Px_local(i) = polarization.get_Px()[indices[i]];
            Py_local(i) = polarization.get_Py()[indices[i]];
            v_local(i)  = fracture.get_v()[indices[i]];
        }

<<<<<<< HEAD
        // ===================================================================
        // 2. MAGIE DU DRY : L'Intégrateur fait tout le travail géométrique
        // ===================================================================
        ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
            
            // LA PHYSIQUE PURE COMMENCE ICI
            B.leftCols(2 * n).setZero();
            for (int i = 0; i < n; ++i) {
                B(0, 2 * i)     = grad_N(0, i); 
                B(1, 2 * i + 1) = grad_N(1, i); 
                B(2, 2 * i)     = grad_N(1, i); 
                B(2, 2 * i + 1) = grad_N(0, i); 
            }

            Eigen::Vector2d Pi_gp(N.head(n).dot(Px_local.head(n)), N.head(n).dot(Py_local.head(n)));
            double v_gp = N.head(n).dot(v_local.head(n));

            if (!std::isfinite(Pi_gp(0)) || !std::isfinite(Pi_gp(1)) || !std::isfinite(v_gp)) return;

            double penalite_fracture = (v_gp * v_gp) + math.eta_k;
            Eigen::Vector3d sigma_0 = math.compute_sigma_0(Pi_gp);

            auto B_active = B.leftCols(2 * n);
            
            K_local.topLeftCorner(2 * n, 2 * n).noalias() += B_active.transpose() * (penalite_fracture * C) * B_active * dV;
            F_local.head(2 * n).noalias() += B_active.transpose() * (penalite_fracture * sigma_0) * dV;
        });

        // 3. Transmission des contributions locales au système global
=======
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
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
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
<<<<<<< HEAD
}

// =========================================================================
// 3. SOUS-FONCTIONS UTILES (Post-traitement et CL)
// =========================================================================
void Mechanics::map_global_vector_to_components(const Eigen::VectorXd& U_new) {
    for (size_t i = 0; i < mesh.get_num_nodes(); ++i) {
=======

    K_global.resize(system_size, system_size);
    K_global.setFromTriplets(triplets.begin(), triplets.end());

    const auto& bcs_ux = bc_manager.get_ux_bcs();
    const auto& bcs_uy = bc_manager.get_uy_bcs();
    apply_boundary_conditions(K_global, F_global, bcs_ux, bcs_uy);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K_global);

    if (solver.info() != Eigen::Success) {
        Logger::debug("[Mechanics] Echec factorisation - matrice singuliere ou mal conditionnee\n", config.debug_enabled());
        return;
    }

    Eigen::VectorXd U_new = solver.solve(F_global);
    Logger::debug("[DEBUG][Mechanics] solver.info=" + std::to_string(solver.info()) + " U min=" + std::to_string(U_new.minCoeff()) + " max=" + std::to_string(U_new.maxCoeff()) + "\n", config.debug_enabled());
    if (solver.info() != Eigen::Success || !U_new.allFinite()) {
        Logger::debug("[Mechanics] Echec resolution ou solution non finie\n", config.debug_enabled());
        return;
    }

    for (size_t i = 0; i < n_nodes; ++i) {
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
        ux_current(i) = U_new(2 * i);
        uy_current(i) = U_new(2 * i + 1);
    }
}

Eigen::Matrix2d Mechanics::get_strain_at_gp(const Element& elem, const GaussPoint2D& gp, const std::vector<std::array<double, 2>>& coords) const {
    const auto& indices = elem.get_node_indices();
    int n_nodes = elem.get_num_nodes();

    double eps_11 = 0.0, eps_22 = 0.0, eps_12 = 0.0;

    if (n_nodes == 3) {
<<<<<<< HEAD
        auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);
        if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) return Eigen::Matrix2d::Zero();
        
=======
        // Triangle : gradients constants, le "gp" passe n'est pas utilise (CST)
        auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);
        if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) {
            Logger::debug("[Mechanics] get_strain_at_gp: triangle degenere detJ=" + std::to_string(detJ) + "\n", config.debug_enabled());
            return Eigen::Matrix2d::Zero();
        }
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
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
<<<<<<< HEAD
        if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) return Eigen::Matrix2d::Zero();
        
=======
        if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) {
            Logger::debug("[Mechanics] get_strain_at_gp: element degenere detJ=" + std::to_string(detJ) + "\n", config.debug_enabled());
            return Eigen::Matrix2d::Zero();
        }
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
        for (int i = 0; i < n_nodes; ++i) {
            double dN_dx = dN_xy[0][i];
            double dN_dy = dN_xy[1][i];
            eps_11 += dN_dx * ux_current[indices[i]];
            eps_22 += dN_dy * uy_current[indices[i]];
            eps_12 += dN_dy * ux_current[indices[i]] + dN_dx * uy_current[indices[i]];
        }
<<<<<<< HEAD
    }

    eps_12 *= 0.5; 
=======
    } else {
        Logger::debug("[Mechanics] get_strain_at_gp: type d'element non supporte n_nodes=" + std::to_string(n_nodes) + "\n", config.debug_enabled());
        return Eigen::Matrix2d::Zero();
    }

    eps_12 *= 0.5;
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da

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

<<<<<<< HEAD
=======
    Logger::debug("[DEBUG][BC phi] n_nodes=" + std::to_string(mesh.get_num_nodes()) +
                  " bcs.size()=" + std::to_string(bcs_x.size()) + "\n", config.debug_enabled());
    int n_dirichlet = 0;
    for (size_t i = 0; i < bcs_x.size(); ++i)
        if (bcs_x[i].type == BCType::DIRICHLET) n_dirichlet++;
    Logger::debug("[DEBUG][BC phi] n_dirichlet=" + std::to_string(n_dirichlet) + "\n", config.debug_enabled());

>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
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

<<<<<<< HEAD
=======
    int nd_x = 0, nd_y = 0;
    for (auto& bc : bcs_x) if (bc.type == BCType::DIRICHLET) nd_x++;
    for (auto& bc : bcs_y) if (bc.type == BCType::DIRICHLET) nd_y++;
    Logger::debug("[DEBUG][Mechanics BC] n_dirichlet_ux=" + std::to_string(nd_x) + " n_dirichlet_uy=" + std::to_string(nd_y) + "\n", config.debug_enabled());
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
    K.makeCompressed();
}