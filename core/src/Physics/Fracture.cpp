#include "Physics/Fracture.h"
#include "Core/ShapeFunctions.h"
#include "Core/Quadrature.h"
#include "Physics/Mechanics.h"
#include "Physics/Polarization.h"
#include "Physics/Electrostatics.h"
#include "Physics/Math.h"
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

void Fracture::update_v(double dt_relax, const Polarization& polarization, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math) {
    size_t n_nodes = mesh.get_num_nodes();
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.get_elements().size() * 16);

    Eigen::VectorXd F_global = Eigen::VectorXd::Zero(n_nodes);
    
    double Gc = config.get_Gc();
    double kappa = config.get_kappa();
    
    // On suppose que la mobilite d'endommagement mu_v est definie (ex: dans config ou en dur comme dans le papier)
    double mu_v = math.mu_v; 
    Logger::debug("[DEBUG][Fracture] mu_v=" + std::to_string(mu_v) + " Gc=" + std::to_string(Gc) + " kappa=" + std::to_string(kappa) +
                  " eta_k=" + std::to_string(math.eta_k) + "\n", config.debug_enabled());
    
    const auto& elements = mesh.get_elements();

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_elem_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);

        Eigen::MatrixXd K_local = Eigen::MatrixXd::Zero(n_elem_nodes, n_elem_nodes);
        Eigen::VectorXd F_local = Eigen::VectorXd::Zero(n_elem_nodes);

        const auto& indices = elem.get_node_indices();
        
        Eigen::VectorXd Px_local(n_elem_nodes);
        Eigen::VectorXd Py_local(n_elem_nodes);
        // Etat fige v_n (debut du pas de temps physique), utilise UNIQUEMENT
        // dans le terme de masse implicite (rhs_coeff).
        Eigen::VectorXd v_local_n(n_elem_nodes);

        for (int i = 0; i < n_elem_nodes; ++i) {
            Px_local(i)     = polarization.get_Px()[indices[i]];
            Py_local(i)     = polarization.get_Py()[indices[i]];
            v_local_n(i)    = v_n[indices[i]];
        }

        auto assemble_gp = [&](const Eigen::RowVectorXd& N, const Eigen::MatrixXd& grad_N, double dV, const GaussPoint2D& gp) {
            
            // 1. Evaluation de la polarisation et de son gradient au point de Gauss
            Eigen::Vector2d Pi_gp(N.dot(Px_local), N.dot(Py_local));
            Eigen::Matrix2d Pij_gp;
            Pij_gp(0, 0) = grad_N.row(0).dot(Px_local); // dPx/dx
            Pij_gp(0, 1) = grad_N.row(1).dot(Px_local); // dPx/dy
            Pij_gp(1, 0) = grad_N.row(0).dot(Py_local); // dPy/dx
            Pij_gp(1, 1) = grad_N.row(1).dot(Py_local); // dPy/dy

            // 2. Recuperation des deformations
            Eigen::Matrix2d eps_gp = mechanics.get_strain_at_gp(elem, gp, coords);
            
            // 3. Calcul de la force motrice (Driving force H_drive) grace a la classe Math
            double U = math.U_energy(Pij_gp);
            double W = math.W_energy(Pi_gp, eps_gp);
            double H_drive = math.compute_H_drive(Pij_gp, Pi_gp, eps_gp, Eigen::Vector2d(0.0, 0.0), false);   

            if (config.get_fracture_mode() == CrackBCType::IMPERMEABLE) {
                double Ex = electrostatics.get_Ex_at_gp(elem, gp);
                double Ey = electrostatics.get_Ey_at_gp(elem, gp);
                Eigen::Vector2d E_gp(Ex, Ey);
                H_drive = math.compute_H_drive(Pij_gp, Pi_gp, eps_gp, E_gp, true);
            }

            if (!std::isfinite(H_drive)) H_drive = 0.0; // Securite numerique

            // 4. Assemblage des matrices (Equation de Ginzburg-Landau standard)
            double mass_coeff = (mu_v / dt_relax) + (Gc / (2.0 * kappa)) + 2.0 * H_drive;
            double diff_coeff = 2.0 * Gc * kappa;
            double rhs_coeff  = (mu_v / dt_relax) * N.dot(v_local_n) + (Gc / (2.0 * kappa));

            K_local.noalias() += (mass_coeff * N.transpose() * N + diff_coeff * grad_N.transpose() * grad_N) * dV;
            F_local.noalias() += (rhs_coeff * N.transpose()) * dV;
        };

        // --- Integration numerique ---
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

        for (int i = 0; i < n_elem_nodes; ++i) {
            for (int j = 0; j < n_elem_nodes; ++j) {
                triplets.emplace_back(indices[i], indices[j], K_local(i, j));
            }
            F_global(indices[i]) += F_local(i);
        }
    }

    Eigen::SparseMatrix<double> K_global(n_nodes, n_nodes);
    K_global.setFromTriplets(triplets.begin(), triplets.end());

    // 5. Irreversibilite de la fracture (Condition unilaterale)
    double alpha = 2e-2; // Seuil pour blocage complet
    for (size_t i = 0; i < n_nodes; ++i) {
        if (v_prev_iter(i) <= alpha) {
            K_global.coeffRef(i, i) += 1e15;
            F_global(i) += 1e15 * 0.0; 
        }
    }

    // --- Solveur persistant : analyzePattern() une seule fois.
    // NB: coeffRef(i,i) ci-dessus ne modifie qu'une entree deja presente dans
    // le pattern (diagonale toujours peuplee par mass_coeff*N^T*N), donc le
    // pattern de sparsite reste identique d'un appel a l'autre. ---
    if (!pattern_analyzed_) {
        solver_.analyzePattern(K_global);
        pattern_analyzed_ = true;
    }
    solver_.factorize(K_global);

    if (solver_.info() != Eigen::Success) {
        Logger::debug("[Fracture] Echec factorisation - matrice singuliere ou mal conditionnee\n", config.debug_enabled());
        return;
    }

    Eigen::VectorXd v_new = solver_.solve(F_global);
    Logger::debug("[DEBUG][Fracture] v_new min=" + std::to_string(v_new.minCoeff()) +
                  " max=" + std::to_string(v_new.maxCoeff()) +
                  " v_current min=" + std::to_string(v_current.minCoeff()) +
                  " max=" + std::to_string(v_current.maxCoeff()) + "\n", config.debug_enabled());

    if (solver_.info() != Eigen::Success || !v_new.allFinite()) {
        Logger::debug("[Fracture] Echec resolution ou solution non finie\n", config.debug_enabled());
        return;
    }

    // Borner la solution (Le parametre d'endommagement ne peut que decroitre et rester positif)
    for (size_t i = 0; i < n_nodes; ++i) {
        v_current(i) = std::max(0.0, std::min(v_new(i), v_prev_iter(i))); 
    }
}