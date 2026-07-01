#include "Physics/Fracture.h"
#include "Core/ShapeFunctions.h"
#include "Core/Quadrature.h"
#include "Physics/Mechanics.h"
#include "Physics/Polarization.h"
#include "Physics/Electrostatics.h"
#include "Physics/Math.h"
#include "Core/ElementIntegrator.h"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

Fracture::Fracture(const Datafile& config, const Mesh& mesh) : config(config), mesh(mesh) {
    size_t n_nodes = mesh.get_num_nodes();
    v_prev_iter.setOnes(n_nodes);
    v_current.resize(n_nodes); 
    v_current.setOnes();
    v_n.setOnes(n_nodes);
}

// void Fracture::update_v(double dt, const Polarization& polarization, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math) {
//     size_t n_nodes = mesh.get_num_nodes();
//     std::vector<Eigen::Triplet<double>> triplets;
//     Eigen::VectorXd F_global = Eigen::VectorXd::Zero(n_nodes);

//     // 1. Assemblage Éléments Finis
//     assemble_system(dt, polarization, mechanics, electrostatics, math, triplets, F_global);

//     Eigen::SparseMatrix<double> K_global(n_nodes, n_nodes);
//     K_global.setFromTriplets(triplets.begin(), triplets.end());

//     // 2. Application des contraintes physiques d'irréversibilité
//     apply_irreversibility(K_global, F_global);

//     // 3. Résolution du système linéaire
//     Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
//     solver.compute(K_global);
//     Eigen::VectorXd v_new = solver.solve(F_global);

//     // 4. Post-traitement et bornage mathématique
//     enforce_physical_bounds(v_new);
// }

void Fracture::update_v(double dt, const Polarization& polarization, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math) {
    size_t n_nodes = mesh.get_num_nodes();
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::VectorXd F_global = Eigen::VectorXd::Zero(n_nodes);

    // 1. Assemblage Éléments Finis
    assemble_system(dt, polarization, mechanics, electrostatics, math, triplets, F_global);

    Eigen::SparseMatrix<double> K_global(n_nodes, n_nodes);
    K_global.setFromTriplets(triplets.begin(), triplets.end());

    // 2. Application des contraintes physiques d'irréversibilité
    // (Ceci renforce massivement la diagonale de la matrice)
    apply_irreversibility(K_global, F_global);

    // =====================================================================
    // 3. RÉSOLUTION ITERATIVE AVEC WARM START ET FALLBACK
    // =====================================================================
    
    // Démarrage à chaud : on part de la solution de la sous-itération précédente
    Eigen::VectorXd v_guess = v_current;
    Eigen::VectorXd v_new;

    // Solveur itératif (CG + Preconditionneur Diagonal)
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> solver;
    solver.setTolerance(1e-10); // Précision stricte pour ne pas abîmer le champ de phase
    solver.setMaxIterations(1000);
    
    solver.compute(K_global);
    v_new = solver.solveWithGuess(F_global, v_guess);

    // Fallback de sécurité si le CG n'a pas convergé
    if (solver.info() != Eigen::Success) {
        Logger::debug("[Fracture] CG a echoue (iters: " + std::to_string(solver.iterations()) + "). Bascule sur LDLT.", config.debug_enabled());
        
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> direct_solver;
        direct_solver.compute(K_global);
        
        if (direct_solver.info() != Eigen::Success) {
            Logger::debug("[Fracture] Echec factorisation LDLT - Matrice singuliere\n", config.debug_enabled());
            return;
        }
        v_new = direct_solver.solve(F_global);
    } else {
        Logger::debug("[Fracture] CG converge en " + std::to_string(solver.iterations()) + " iterations.", config.debug_enabled());
    }

    // Vérification mathématique
    if (!v_new.allFinite()) {
        Logger::debug("[Fracture] Echec resolution : solution non finie\n", config.debug_enabled());
        return;
    }

    // 4. Post-traitement et bornage mathématique (v dans [0, 1])
    enforce_physical_bounds(v_new);
}

void Fracture::assemble_system(double dt, const Polarization& polarization, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math, std::vector<Eigen::Triplet<double>>& triplets, Eigen::VectorXd& F_global) {
    const auto& elements = mesh.get_elements();
    
    // 1. Allocations hors boucle
    const int MAX_NODES = 4;
    Eigen::MatrixXd K_local(MAX_NODES, MAX_NODES);
    Eigen::VectorXd F_local(MAX_NODES);
    Eigen::VectorXd Px_local(MAX_NODES), Py_local(MAX_NODES), v_local_n(MAX_NODES);
    
    // Tampons géométriques
    Eigen::RowVectorXd N(MAX_NODES);
    Eigen::MatrixXd grad_N(2, MAX_NODES);

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_elem_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        const auto& indices = elem.get_node_indices();

        K_local.topLeftCorner(n_elem_nodes, n_elem_nodes).setZero();
        F_local.head(n_elem_nodes).setZero();

        for (int i = 0; i < n_elem_nodes; ++i) {
            Px_local(i)     = polarization.get_Px()[indices[i]];
            Py_local(i)     = polarization.get_Py()[indices[i]];
            v_local_n(i)    = v_n[indices[i]];
        }

        // ===================================================================
        // 2. MAGIE DU DRY : L'Intégrateur fait tout le travail géométrique
        // ===================================================================
        ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
            Eigen::Vector2d Pi_gp(N.head(n).dot(Px_local.head(n)), N.head(n).dot(Py_local.head(n)));
            auto grad_N_active = grad_N.leftCols(n);

            Eigen::Matrix2d Pij_gp;
            Pij_gp(0, 0) = grad_N_active.row(0).dot(Px_local.head(n));
            Pij_gp(0, 1) = grad_N_active.row(1).dot(Px_local.head(n));
            Pij_gp(1, 0) = grad_N_active.row(0).dot(Py_local.head(n));
            Pij_gp(1, 1) = grad_N_active.row(1).dot(Py_local.head(n));

            Eigen::Matrix2d eps_gp = mechanics.get_strain_at_gp(elem, gp, coords);

            bool is_impermeable = (config.get_fracture_mode() == CrackBCType::IMPERMEABLE);
            Eigen::Vector2d E_gp = Eigen::Vector2d::Zero();
            if (is_impermeable) {
                E_gp = Eigen::Vector2d(electrostatics.get_Ex_at_gp(elem, gp), electrostatics.get_Ey_at_gp(elem, gp));
            }

            double H_drive = math.compute_H_drive(Pij_gp, Pi_gp, eps_gp, E_gp, is_impermeable);
            if (!std::isfinite(H_drive)) H_drive = 0.0;

            double mass_coeff = (config.get_mu_v() / dt) + (config.get_Gc() / (2.0 * config.get_kappa())) + 2.0 * H_drive;
            double diff_coeff = 2.0 * config.get_Gc() * config.get_kappa();
            double rhs_coeff  = (config.get_mu_v() / dt) * N.head(n).dot(v_local_n.head(n)) + (config.get_Gc() / (2.0 * config.get_kappa()));

            auto N_active = N.head(n);
            K_local.topLeftCorner(n, n).noalias() += (mass_coeff * N_active.transpose() * N_active + diff_coeff * grad_N_active.transpose() * grad_N_active) * dV;
            F_local.head(n).noalias() += (rhs_coeff * N_active.transpose()) * dV;
        });

        // 3. Transmission locale -> globale
        for (int i = 0; i < n_elem_nodes; ++i) {
            for (int j = 0; j < n_elem_nodes; ++j) {
                triplets.emplace_back(indices[i], indices[j], K_local(i, j));
            }
            F_global(indices[i]) += F_local(i);
        }
    }
}

void Fracture::apply_irreversibility(Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) {
    double alpha = 2e-2; 
    for (size_t i = 0; i < mesh.get_num_nodes(); ++i) {
        if (v_prev_iter(i) <= alpha) {
            K_global.coeffRef(i, i) += 1e15;
            F_global(i) += 1e15 * 0.0; 
        }
    }
}

void Fracture::enforce_physical_bounds(const Eigen::VectorXd& v_new) {
    for (size_t i = 0; i < mesh.get_num_nodes(); ++i) {
        v_current(i) = std::max(0.0, std::min(v_new(i), v_prev_iter(i))); 
    }
}