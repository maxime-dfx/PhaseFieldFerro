#include "Physics/Mechanics.h"
#include "Physics/Electrostatics.h"
#include "Physics/Polarization.h"
#include "Core/ShapeFunctions.h"
#include "Utils/Logger.h"
#include "Physics/Fracture.h"
#include "Physics/Math.h"
#include "Core/Quadrature.h"
#include <cmath>
#include <Eigen/Dense>

Polarization::Polarization(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager) : config(config), mesh(mesh), bc_manager(bc_manager) {
    size_t n_nodes = mesh.get_num_nodes();
    Px_prev_iter.setZero(n_nodes);
    Py_prev_iter.setZero(n_nodes);
    Px_current.resize(n_nodes);
    Px_current.setZero();
    Py_current.resize(n_nodes);
    Py_current.setZero();

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
}

void Polarization::apply_boundary_conditions(Eigen::MatrixXd& A, Eigen::VectorXd& b, const std::vector<NodeBC>& bcs) {
    for (size_t i = 0; i < bcs.size(); ++i) {
        if (bcs[i].type == BCType::DIRICHLET) {
            A.row(i).setZero();
            A(i, i) = 1.0;
            b(i) = bcs[i].value;
        }
        // Note: Si type == NEUMANN, vous devriez théoriquement ajouter 
        // une intégrale de surface au vecteur 'b', mais si votre 
        // formulation traite Neumann comme une condition naturelle (0 flux),
        // vous n'avez rien à faire ici.
    }
}

void Polarization::update_polarization_component(double time, const Fracture& fracture, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math, int component) {
    double dt = config.get_dt();
    
    // Initialisation des matrices globales
    A.setZero(Px_current.size(), Px_current.size());
    b.setZero(Px_current.size());

    const auto& elements = mesh.get_elements();
    
    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        
        local_A.setZero(n_nodes, n_nodes);
        local_b.setZero(n_nodes);
        Px_local.setZero(n_nodes);
        Py_local.setZero(n_nodes);
        v_local.setZero(n_nodes);

        const auto& indices = elem.get_node_indices();
        for (int i = 0; i < n_nodes; ++i) {
            Px_local(i) = Px_current[indices[i]];
            Py_local(i) = Py_current[indices[i]];
            v_local(i) = fracture.get_v()[indices[i]];
        }

        const auto& gauss_points = Quadrature::get_gauss_2x2();

        for (const auto& gp : gauss_points) {
            auto N_std = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
            auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
            auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);

            double dV = gp.weight * detJ;

            Eigen::RowVectorXd N(4);
            N << N_std[0], N_std[1], N_std[2], N_std[3];

            Eigen::MatrixXd grad_N(2, 4);
            for(int i = 0; i < 4; ++i) {
                grad_N(0, i) = dN_xy[0][i];
                grad_N(1, i) = dN_xy[1][i];
            }

            double v_gp = N.dot(v_local);
            Eigen::Vector2d Pi_gp(N.dot(Px_local), N.dot(Py_local));
            Eigen::Matrix2d eps_gp = mechanics.get_strain_at_gp(elem, gp, coords);

            // --- DIFFÉRENCIATION X vs Y ---
            double E_gp = 0.0;
            double dW = 0.0;
            double dchi = 0.0;
            
            // Pointe vers le vecteur local approprié pour b
            const Eigen::VectorXd& P_comp_local = (component == 0) ? Px_local : Py_local;

            if (component == 0) { // Composante X
                E_gp = electrostatics.get_Ex_at_gp(elem, gp);
                dW = math.dW_dp1(Pi_gp, eps_gp);
                dchi = math.dchi_dp1(Pi_gp);
            } else {              // Composante Y
                E_gp = electrostatics.get_Ey_at_gp(elem, gp);
                dW = math.dW_dp2(Pi_gp, eps_gp);   // Correction : dp2 au lieu de dp1
                dchi = math.dchi_dp2(Pi_gp);       // Correction : dp2 au lieu de dp1
            }

            double penalite_fracture = (v_gp * v_gp) + math.eta_k;
            double force_thermo = 0.0;
            
            // Note : remplacez fracture_mode par votre vrai appel de configuration
            if (config.get_fracture_mode() == CrackBCType::PERMEABLE) {
                force_thermo = - (penalite_fracture * dW + dchi - E_gp);
            } else if (config.get_fracture_mode() == CrackBCType::IMPERMEABLE) {
                force_thermo = - (penalite_fracture * dW + dchi - (penalite_fracture * E_gp)); 
            }

            // Assemblage local
            local_A.noalias() += (math.mu_p / dt) * (N.transpose() * N * dV) 
                               + (math.a0 * penalite_fracture) * (grad_N.transpose() * grad_N * dV);
                               
            local_b.noalias() += (math.mu_p / dt) * (N.transpose() * N * dV) * P_comp_local 
                               + N.transpose() * force_thermo * dV;
        }

        // Assemblage global
        for (int i = 0; i < n_nodes; ++i) {
            for (int j = 0; j < n_nodes; ++j) {
                A(indices[i], indices[j]) += local_A(i, j);
            }
            b(indices[i]) += local_b(i);
        }
    }
    if (component == 0) {
        apply_boundary_conditions(A, b, bc_manager.get_px_bcs());
        
        // DEBUG: Vérifier la matrice
        std::cout << "Condition number estimate: " << A.norm() * A.inverse().norm() << std::endl;
        std::cout << "Matrix A - min coeff: " << A.minCoeff() << ", max coeff: " << A.maxCoeff() << std::endl;
        std::cout << "Vector b - min coeff: " << b.minCoeff() << ", max coeff: " << b.maxCoeff() << std::endl;
        
        auto solution = A.colPivHouseholderQr().solve(b);
        
        // Vérifier la solution
        std::cout << "Solution - min: " << solution.minCoeff() << ", max: " << solution.maxCoeff() << std::endl;
        std::cout << "Solution has NaN: " << solution.hasNaN() << ", has Inf: " << solution.hasInf() << std::endl;
        
        Px_current = solution;
    }
    // --- APPLICATION DES CL ET RÉSOLUTION ---
    if (component == 0) {
        apply_boundary_conditions(A, b, bc_manager.get_px_bcs());
        Px_current = A.colPivHouseholderQr().solve(b);
    } else {
        apply_boundary_conditions(A, b, bc_manager.get_py_bcs());
        Py_current = A.colPivHouseholderQr().solve(b);
    }
}