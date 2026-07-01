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
}

// =========================================================================
// 1. LE CHEF D'ORCHESTRE ( update_electric_potential )
// =========================================================================
void Electrostatics::update_electric_potential(double time, const Polarization& polarization, const Fracture& fracture, const Math& math) {
    size_t n_nodes = mesh.get_num_nodes();
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.get_elements().size() * 16);
    Eigen::VectorXd F_global = Eigen::VectorXd::Zero(n_nodes);
    
    // 1. Assemblage du laplacien diélectrique
    assemble_electrostatic_system(polarization, fracture, math, triplets, F_global);

    Eigen::SparseMatrix<double> K_global(n_nodes, n_nodes);
    K_global.setFromTriplets(triplets.begin(), triplets.end());

    // 2. Application des Conditions aux Limites de Dirichlet
    apply_boundary_conditions(K_global, F_global, bc_manager.get_phi_bcs());

    // 3. Résolution du système (potentiel phi)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K_global);

    if (solver.info() != Eigen::Success) {
        Logger::debug("[Electrostatics] Echec factorisation - matrice singuliere\n", config.debug_enabled());
        return;
    }

    phi_current = solver.solve(F_global);
    Logger::debug("[DEBUG][Electrostatics] phi_current min=" + std::to_string(phi_current.minCoeff()) +
                  " max=" + std::to_string(phi_current.maxCoeff()) + "\n", config.debug_enabled());

    // 4. Post-traitement : Calcul du champ électrique dérivé (E = -grad phi)
    if (solver.info() == Eigen::Success && phi_current.allFinite()) {
        compute_nodal_electric_field(); 
    }
}

// =========================================================================
// 2. L'ASSEMBLEUR GLOBAL
// =========================================================================
void Electrostatics::assemble_electrostatic_system(const Polarization& polarization, const Fracture& fracture, const Math& math, std::vector<Eigen::Triplet<double>>& triplets, Eigen::VectorXd& F_global) {
    const auto& elements = mesh.get_elements();
    double eps0 = math.eps0;
    double eta_k = math.eta_k;

    // 1. Allocations hors boucle
    const int MAX_NODES = 4;
    Eigen::MatrixXd K_local(MAX_NODES, MAX_NODES);
    Eigen::VectorXd F_local(MAX_NODES);
    Eigen::VectorXd Px_loc(MAX_NODES), Py_loc(MAX_NODES), v_loc(MAX_NODES);
    
    // Tampons géométriques
    Eigen::RowVectorXd N(MAX_NODES);
    Eigen::MatrixXd grad_N(2, MAX_NODES);

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_nodes_elem = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        const auto& indices = elem.get_node_indices();

        K_local.topLeftCorner(n_nodes_elem, n_nodes_elem).setZero();
        F_local.head(n_nodes_elem).setZero();

        for(int i = 0; i < n_nodes_elem; ++i) {
            Px_loc(i) = polarization.get_Px()[indices[i]];
            Py_loc(i) = polarization.get_Py()[indices[i]];
            v_loc(i)  = fracture.get_v()[indices[i]];
        }

        // ===================================================================
        // 2. MAGIE DU DRY : L'Intégrateur fait tout le travail géométrique
        // ===================================================================
        ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
            double v_gp = N.head(n).dot(v_loc.head(n));
            double px_gp = N.head(n).dot(Px_loc.head(n));
            double py_gp = N.head(n).dot(Py_loc.head(n));

            double phase_factor = 1.0;
            if (config.get_fracture_mode() == CrackBCType::IMPERMEABLE) {
                phase_factor = (v_gp * v_gp) + eta_k;
            }

            double eps_eff = eps0 * phase_factor;
            Eigen::Vector2d P_eff(px_gp * phase_factor, py_gp * phase_factor);

            auto grad_N_active = grad_N.leftCols(n); 
            
            K_local.topLeftCorner(n, n).noalias() += grad_N_active.transpose() * eps_eff * grad_N_active * dV;
            F_local.head(n).noalias() += grad_N_active.transpose() * P_eff * dV;
        });
        
        // 3. Transmission locale -> globale
        for (int i = 0; i < n_nodes_elem; ++i) {
            for (int j = 0; j < n_nodes_elem; ++j) {
                triplets.emplace_back(indices[i], indices[j], K_local(i, j));
            }
            F_global(indices[i]) += F_local(i);
        }
    }
}

// =========================================================================
// 3. MÉTHODES AUXILIAIRES ET POST-TRAITEMENT
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
            // Logique T3 (Constante sur l'élément)
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
        } else if (n_nodes_elem == 4) {
            // Logique Q4 Ajoutée : Intégration sur les points de Gauss
            const auto& gauss_points = Quadrature::get_gauss_2x2();
            for (const auto& gp : gauss_points) {
                auto N_std = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
                auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
                auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);
                
                double dV = gp.weight * std::abs(detJ);
                double Ex_gp = 0.0, Ey_gp = 0.0;
                
                for(int i = 0; i < 4; ++i) {
                    Ex_gp -= dN_xy[0][i] * phi_loc(i);
                    Ey_gp -= dN_xy[1][i] * phi_loc(i);
                }

                // Distribution aux noeuds (Lumped Mass via fonctions de forme)
                for(int i = 0; i < 4; ++i) {
                    Fx(indices[i]) += Ex_gp * N_std[i] * dV;
                    Fy(indices[i]) += Ey_gp * N_std[i] * dV;
                    M_lumped(indices[i]) += N_std[i] * dV;
                }
            }
        }
    }

    // Projection finale pour obtenir E
    for(size_t i = 0; i < n_nodes; ++i) {
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