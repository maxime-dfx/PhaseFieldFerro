#include "Physics/Electrostatics.h"
#include "Physics/Polarization.h"
#include "Physics/Fracture.h"
#include "Physics/Mechanics.h"
#include "Physics/Math.h"
#include "Core/ShapeFunctions.h"
#include "Core/Quadrature.h"
#include "Utils/Logger.h"
#include "Core/ElementIntegrator.h"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

Electrostatics::Electrostatics(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager) 
    : config(config), mesh(mesh), bc_manager(bc_manager) {
    size_t n_nodes = mesh.get_num_nodes();
    phi_current.resize(n_nodes);
    phi_current.setZero();
    Ex_current.resize(n_nodes); Ex_current.setZero();
    Ey_current.resize(n_nodes); Ey_current.setZero();
    Ex_prev_iter.resize(n_nodes); Ex_prev_iter.setZero();
    Ey_prev_iter.resize(n_nodes); Ey_prev_iter.setZero();
}

void Electrostatics::update_electric_potential(double time, const Polarization& polarization, const Fracture& fracture, const Math& math) {
    size_t n_nodes = mesh.get_num_nodes();
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.get_elements().size() * 16);
    Eigen::VectorXd F_global = Eigen::VectorXd::Zero(n_nodes);
    
    const auto& elements = mesh.get_elements();
    
    // Utilisation directe des membres de ta classe Math
    double eps0 = math.eps0;
    double eta_k = math.eta_k;
    Logger::debug("[DEBUG][Electrostatics] eps0=" + std::to_string(eps0) + " eta_k=" + std::to_string(eta_k) + "\n", config.debug_enabled());
    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_nodes_elem = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);

        Eigen::MatrixXd K_local = Eigen::MatrixXd::Zero(n_nodes_elem, n_nodes_elem);
        Eigen::VectorXd F_local = Eigen::VectorXd::Zero(n_nodes_elem);
        const auto& indices = elem.get_node_indices();

        Eigen::VectorXd Px_loc(n_nodes_elem), Py_loc(n_nodes_elem), v_loc(n_nodes_elem);
        for(int i=0; i<n_nodes_elem; ++i) {
            Px_loc(i) = polarization.get_Px()[indices[i]];
            Py_loc(i) = polarization.get_Py()[indices[i]];
            v_loc(i)  = fracture.get_v()[indices[i]];
        }

        auto assemble_gp = [&](const Eigen::RowVectorXd& N, const Eigen::MatrixXd& grad_N, double dV) {
            double v_gp = N.dot(v_loc);
            double px_gp = N.dot(Px_loc);
            double py_gp = N.dot(Py_loc);

            double phase_factor = 1.0;
            if (config.get_fracture_mode() == CrackBCType::IMPERMEABLE) {
                phase_factor = (v_gp * v_gp) + eta_k;
            }

            double eps_eff = eps0 * phase_factor;
            Eigen::Vector2d P_eff(px_gp * phase_factor, py_gp * phase_factor);

            // Matrice de rigidite dielectrique et vecteur force surfacique due a la polarisation
            K_local.noalias() += grad_N.transpose() * eps_eff * grad_N * dV;
            F_local.noalias() += grad_N.transpose() * P_eff * dV;
        };

        if (n_nodes_elem == 3) {
            double xi = 1.0 / 3.0, eta = 1.0 / 3.0;
            auto N_std = ShapeFunctions::get_shape_functions_tri(xi, eta);
            auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);

            if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) {
                Logger::debug("[DEBUG][Electrostatics] Triangle degenere idx=" + std::to_string(elem_idx) +
                              " detJ=" + std::to_string(detJ) + " -> contribution ignoree\n", config.debug_enabled());
                continue;
            }

            double dV = std::abs(detJ) / 2.0;
            Eigen::RowVectorXd N(3); N << N_std[0], N_std[1], N_std[2];
            Eigen::MatrixXd grad_N(2, 3);
            
            for (int i = 0; i < 3; ++i) {
                grad_N(0, i) = dN_xy[i][0];
                grad_N(1, i) = dN_xy[i][1];
            }

            assemble_gp(N, grad_N, dV);

        } else if (n_nodes_elem == 4) {
            const auto& gauss_points = Quadrature::get_gauss_2x2();
            for (const auto& gp : gauss_points) {
                auto N_std = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
                auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
                auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);
                double dV = gp.weight * std::abs(detJ);
                Eigen::RowVectorXd N(4); N << N_std[0], N_std[1], N_std[2], N_std[3];
                Eigen::MatrixXd grad_N(2, 4);
                for(int i=0; i<4; i++) { grad_N(0,i) = dN_xy[0][i]; grad_N(1,i) = dN_xy[1][i]; }

                assemble_gp(N, grad_N, dV);
            }
        }
        
        for (int i = 0; i < n_nodes_elem; ++i) {
            for (int j = 0; j < n_nodes_elem; ++j) {
                triplets.emplace_back(indices[i], indices[j], K_local(i, j));
            }
            F_global(indices[i]) += F_local(i);
        }
    }

    // --- Construction de la matrice sparse globale (MANQUAIT precedemment) ---
    Eigen::SparseMatrix<double> K_global(n_nodes, n_nodes);
    K_global.setFromTriplets(triplets.begin(), triplets.end());

    // --- Application des CL Dirichlet sur le potentiel ---
    const auto& bcs_phi = bc_manager.get_phi_bcs();
    apply_boundary_conditions(K_global, F_global, bcs_phi);

    // --- Resolution (solveur persistant : analyzePattern() une seule fois) ---
    if (!pattern_analyzed_) {
        solver_.analyzePattern(K_global);
        pattern_analyzed_ = true;
    }
    solver_.factorize(K_global);

    if (solver_.info() != Eigen::Success) {
        Logger::debug("[Electrostatics] Echec factorisation - matrice singuliere ou mal conditionnee\n", config.debug_enabled());
        return;
    }

    Eigen::VectorXd phi_new = solver_.solve(F_global);
    Logger::debug("[DEBUG][Electrostatics] phi_new min=" + std::to_string(phi_new.minCoeff()) +
                  " max=" + std::to_string(phi_new.maxCoeff()) + "\n", config.debug_enabled());

    if (solver_.info() != Eigen::Success || !phi_new.allFinite()) {
        Logger::debug("[Electrostatics] Echec resolution ou solution non finie\n", config.debug_enabled());
        return;
    }

    phi_current = phi_new;

    // Recalcule le champ electrique nodal (Ex, Ey) a partir du nouveau potentiel
    compute_nodal_electric_field();
}

// =========================================================================
// 3. METHODES AUXILIAIRES ET POST-TRAITEMENT
// =========================================================================
void Electrostatics::apply_boundary_conditions(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, const std::vector<NodeBC>& bcs) {
    double max_diag = 0.0;
    for (int k = 0; k < K.outerSize(); ++k)
        max_diag = std::max(max_diag, std::abs(K.coeff(k, k)));
    const double penalty = max_diag * 1e7;

    for (size_t i = 0; i < bcs.size(); ++i) {
        if (bcs[i].type == BCType::DIRICHLET) {
            int dof = static_cast<int>(i); 
            K.coeffRef(dof, dof) += penalty;
            F(dof) += penalty * bcs[i].value;
        }
    }
    K.makeCompressed();
}

void Electrostatics::compute_nodal_electric_field() {
    size_t n_nodes = mesh.get_num_nodes();
    Eigen::VectorXd Fx = Eigen::VectorXd::Zero(n_nodes);
    Eigen::VectorXd Fy = Eigen::VectorXd::Zero(n_nodes);
    Eigen::VectorXd M_lumped = Eigen::VectorXd::Zero(n_nodes);

    const auto& elements = mesh.get_elements();
    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_nodes_elem = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        const auto& indices = elem.get_node_indices();
        
        Eigen::VectorXd phi_loc(n_nodes_elem);
        for(int i=0; i<n_nodes_elem; ++i) phi_loc(i) = phi_current(indices[i]);

        if (n_nodes_elem == 3) {
            auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);
            double area = std::abs(detJ) / 2.0;
            double mass_contrib = area / 3.0;

            double Ex_elem = 0.0, Ey_elem = 0.0;
            for(int i=0; i<3; ++i) {
                Ex_elem -= dN_xy[i][0] * phi_loc(i);
                Ey_elem -= dN_xy[i][1] * phi_loc(i);
            }

            for(int i=0; i<3; ++i) {
                Fx(indices[i]) += Ex_elem * mass_contrib;
                Fy(indices[i]) += Ey_elem * mass_contrib;
                M_lumped(indices[i]) += mass_contrib;
            }
        } 
        // Ajouter la logique Q4 si necessaire...
    }

    for(size_t i=0; i<n_nodes; ++i) {
        if(M_lumped(i) > 1e-12) {
            Ex_current(i) = Fx(i) / M_lumped(i);
            Ey_current(i) = Fy(i) / M_lumped(i);
        }
    }
}

// Interpolation classique aux points de Gauss depuis les variables nodales
double Electrostatics::get_Ex_at_gp(const Element& elem, const GaussPoint2D& gp) const {
    const auto& indices = elem.get_node_indices();
    int n_nodes = elem.get_num_nodes();
    double Ex = 0.0;
    
    if (n_nodes == 3) {
        auto N = ShapeFunctions::get_shape_functions_tri(1.0/3.0, 1.0/3.0);
        for(int i=0; i<3; ++i) Ex += N[i] * Ex_current[indices[i]];
    } else if (n_nodes == 4) {
        auto N = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
        for(int i=0; i<4; ++i) Ex += N[i] * Ex_current[indices[i]];
    }
    return Ex;
}

double Electrostatics::get_Ey_at_gp(const Element& elem, const GaussPoint2D& gp) const {
    const auto& indices = elem.get_node_indices();
    int n_nodes = elem.get_num_nodes();
    double Ey = 0.0;
    
    if (n_nodes == 3) {
        auto N = ShapeFunctions::get_shape_functions_tri(1.0/3.0, 1.0/3.0);
        for(int i=0; i<3; ++i) Ey += N[i] * Ey_current[indices[i]];
    } else if (n_nodes == 4) {
        auto N = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
        for(int i=0; i<4; ++i) Ey += N[i] * Ey_current[indices[i]];
    }
    return Ey;
}