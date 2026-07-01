#include "Physics/Mechanics.h"
#include "Core/Quadrature.h"
#include "Physics/Math.h"
#include "Physics/Polarization.h"
#include "Physics/Fracture.h"
#include "Core/ShapeFunctions.h"
#include "Utils/Logger.h"
#include "Core/ElementIntegrator.h"
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

// =========================================================================
// 1. LE CHEF D'ORCHESTRE ( update_u )
// =========================================================================
void Mechanics::update_u(double time, const Polarization& polarization, const Fracture& fracture, const Math& math) {
    size_t n_nodes = mesh.get_num_nodes();
    size_t system_size = 2 * n_nodes;

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.get_elements().size() * 64);
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

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_elem_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        const auto& indices = elem.get_node_indices();

        K_local.topLeftCorner(2 * n_elem_nodes, 2 * n_elem_nodes).setZero();
        F_local.head(2 * n_elem_nodes).setZero();
        
        for (int i = 0; i < n_elem_nodes; ++i) {
            Px_local(i) = polarization.get_Px()[indices[i]];
            Py_local(i) = polarization.get_Py()[indices[i]];
            v_local(i)  = fracture.get_v()[indices[i]];
        }

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
}

// =========================================================================
// 3. SOUS-FONCTIONS UTILES (Post-traitement et CL)
// =========================================================================
void Mechanics::map_global_vector_to_components(const Eigen::VectorXd& U_new) {
    for (size_t i = 0; i < mesh.get_num_nodes(); ++i) {
        ux_current(i) = U_new(2 * i);
        uy_current(i) = U_new(2 * i + 1);
    }
}

Eigen::Matrix2d Mechanics::get_strain_at_gp(const Element& elem, const GaussPoint2D& gp, const std::vector<std::array<double, 2>>& coords) const {
    const auto& indices = elem.get_node_indices();
    int n_nodes = elem.get_num_nodes();

    double eps_11 = 0.0, eps_22 = 0.0, eps_12 = 0.0;

    if (n_nodes == 3) {
        auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);
        if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) return Eigen::Matrix2d::Zero();
        
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
        if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) return Eigen::Matrix2d::Zero();
        
        for (int i = 0; i < n_nodes; ++i) {
            double dN_dx = dN_xy[0][i];
            double dN_dy = dN_xy[1][i];
            eps_11 += dN_dx * ux_current[indices[i]];
            eps_22 += dN_dy * uy_current[indices[i]];
            eps_12 += dN_dy * ux_current[indices[i]] + dN_dx * uy_current[indices[i]];
        }
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

    K.makeCompressed();
}