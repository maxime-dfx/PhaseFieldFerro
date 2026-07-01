#include "Physics/Fracture.h"
#include "Core/ShapeFunctions.h"
#include "Core/Quadrature.h"
#include "Physics/Mechanics.h"
#include "Physics/Polarization.h"
#include "Physics/Electrostatics.h"
#include "Physics/Math.h"
<<<<<<< HEAD
#include "Core/ElementIntegrator.h"
=======
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
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

<<<<<<< HEAD
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
=======
void Fracture::update_v(double dt, const Polarization& polarization, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math) {
    size_t n_nodes = mesh.get_num_nodes();
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.get_elements().size() * 16);

    Eigen::VectorXd F_global = Eigen::VectorXd::Zero(n_nodes);
    
    double Gc = config.get_Gc();
    double kappa = config.get_kappa();
    
    // On suppose que la mobilité d'endommagement mu_v est définie (ex: dans config ou en dur comme dans le papier)
    double mu_v = math.mu_v; 
    Logger::debug("[DEBUG][Fracture] mu_v=" + std::to_string(mu_v) + " Gc=" + std::to_string(Gc) + " kappa=" + std::to_string(kappa) +
                  " eta_k=" + std::to_string(math.eta_k) + "\n", config.debug_enabled());
    
    const auto& elements = mesh.get_elements();
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_elem_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
<<<<<<< HEAD
        const auto& indices = elem.get_node_indices();

        K_local.topLeftCorner(n_elem_nodes, n_elem_nodes).setZero();
        F_local.head(n_elem_nodes).setZero();
=======

        Eigen::MatrixXd K_local = Eigen::MatrixXd::Zero(n_elem_nodes, n_elem_nodes);
        Eigen::VectorXd F_local = Eigen::VectorXd::Zero(n_elem_nodes);

        const auto& indices = elem.get_node_indices();
        
        Eigen::VectorXd Px_local(n_elem_nodes);
        Eigen::VectorXd Py_local(n_elem_nodes);
        // Etat fige v_n (debut du pas de temps physique), utilise UNIQUEMENT
        // dans le terme de masse implicite (rhs_coeff).
        Eigen::VectorXd v_local_n(n_elem_nodes);
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da

        for (int i = 0; i < n_elem_nodes; ++i) {
            Px_local(i)     = polarization.get_Px()[indices[i]];
            Py_local(i)     = polarization.get_Py()[indices[i]];
            v_local_n(i)    = v_n[indices[i]];
        }

<<<<<<< HEAD
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
=======
        auto assemble_gp = [&](const Eigen::RowVectorXd& N, const Eigen::MatrixXd& grad_N, double dV, const GaussPoint2D& gp) {
            
            // 1. Évaluation de la polarisation et de son gradient au point de Gauss
            Eigen::Vector2d Pi_gp(N.dot(Px_local), N.dot(Py_local));
            Eigen::Matrix2d Pij_gp;
            Pij_gp(0, 0) = grad_N.row(0).dot(Px_local); // dPx/dx
            Pij_gp(0, 1) = grad_N.row(1).dot(Px_local); // dPx/dy
            Pij_gp(1, 0) = grad_N.row(0).dot(Py_local); // dPy/dx
            Pij_gp(1, 1) = grad_N.row(1).dot(Py_local); // dPy/dy

            // 2. Récupération des déformations
            Eigen::Matrix2d eps_gp = mechanics.get_strain_at_gp(elem, gp, coords);
            
            // 3. Calcul de la force motrice (Driving force H_drive) grâce à la classe Math
            double U = math.U_energy(Pij_gp);
            double W = math.W_energy(Pi_gp, eps_gp);
            double chi = math.chi_energy(Pi_gp); 
            double H_drive = U + W + chi;        

            if (config.get_fracture_mode() == CrackBCType::IMPERMEABLE) {
                double Ex = electrostatics.get_Ex_at_gp(elem, gp);
                double Ey = electrostatics.get_Ey_at_gp(elem, gp);
                Eigen::Vector2d E_gp(Ex, Ey);
                
                // Correction thermodynamique pour fissure imperméable
                H_drive += -0.5 * math.eps0 * E_gp.squaredNorm() - E_gp.dot(Pi_gp);
            }

            if (!std::isfinite(H_drive)) H_drive = 0.0; // Sécurité numérique

            // 4. Assemblage des matrices (Équation de Ginzburg-Landau standard)
            double mass_coeff = (mu_v / dt) + (Gc / (2.0 * kappa)) + 2.0 * H_drive;
            double diff_coeff = 2.0 * Gc * kappa;
            double rhs_coeff  = (mu_v / dt) * N.dot(v_local_n) + (Gc / (2.0 * kappa));

            K_local.noalias() += (mass_coeff * N.transpose() * N + diff_coeff * grad_N.transpose() * grad_N) * dV;
            F_local.noalias() += (rhs_coeff * N.transpose()) * dV;
        };

        // --- Intégration numérique ---
        if (n_elem_nodes == 3) {
            double xi = 1.0 / 3.0, eta = 1.0 / 3.0;
            auto N_std = ShapeFunctions::get_shape_functions_tri(xi, eta);
            auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);
            double dV = std::abs(detJ) / 2.0;
            
            Eigen::RowVectorXd N(3); N << N_std[0], N_std[1], N_std[2];
            Eigen::MatrixXd grad_N(2, 3);
            
            for (int i = 0; i < 3; ++i) {
                grad_N(0, i) = dN_xy[i][0];
                grad_N(1, i) = dN_xy[i][1];
            }

            GaussPoint2D gp_fake{xi, eta, 1.0};
            assemble_gp(N, grad_N, dV, gp_fake);

        } else if (n_elem_nodes == 4) {
            const auto& gauss_points = Quadrature::get_gauss_2x2();
            for (const auto& gp : gauss_points) {
                auto N_std = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
                auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
                auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);
                
                double dV = gp.weight * std::abs(detJ);
                Eigen::RowVectorXd N(4); N << N_std[0], N_std[1], N_std[2], N_std[3];
                Eigen::MatrixXd grad_N(2, 4);
                for(int i=0; i<4; i++) { grad_N(0,i) = dN_xy[0][i]; grad_N(1,i) = dN_xy[1][i]; }

                assemble_gp(N, grad_N, dV, gp);
            }
        }

>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
        for (int i = 0; i < n_elem_nodes; ++i) {
            for (int j = 0; j < n_elem_nodes; ++j) {
                triplets.emplace_back(indices[i], indices[j], K_local(i, j));
            }
            F_global(indices[i]) += F_local(i);
        }
    }
<<<<<<< HEAD
}

void Fracture::apply_irreversibility(Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) {
    double alpha = 2e-2; 
    for (size_t i = 0; i < mesh.get_num_nodes(); ++i) {
=======

    Eigen::SparseMatrix<double> K_global(n_nodes, n_nodes);
    K_global.setFromTriplets(triplets.begin(), triplets.end());

    // 5. Irréversibilité de la fracture (Condition unilatérale)
    double alpha = 2e-2; // Seuil pour blocage complet
    for (size_t i = 0; i < n_nodes; ++i) {
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
        if (v_prev_iter(i) <= alpha) {
            K_global.coeffRef(i, i) += 1e15;
            F_global(i) += 1e15 * 0.0; 
        }
    }
<<<<<<< HEAD
}

void Fracture::enforce_physical_bounds(const Eigen::VectorXd& v_new) {
    for (size_t i = 0; i < mesh.get_num_nodes(); ++i) {
=======

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K_global);
    Eigen::VectorXd v_new = solver.solve(F_global);
    Logger::debug("[DEBUG][Fracture] v_new min=" + std::to_string(v_new.minCoeff()) +
                  " max=" + std::to_string(v_new.maxCoeff()) +
                  " v_current min=" + std::to_string(v_current.minCoeff()) +
                  " max=" + std::to_string(v_current.maxCoeff()) + "\n", config.debug_enabled());
    // Borner la solution (Le paramètre d'endommagement ne peut que décroître et rester positif)
    for (size_t i = 0; i < n_nodes; ++i) {
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
        v_current(i) = std::max(0.0, std::min(v_new(i), v_prev_iter(i))); 
    }
}