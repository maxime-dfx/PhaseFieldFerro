#include "Physics/Mechanics.h"
#include "Core/Quadrature.h"
#include "Physics/Math.h"
#include "Physics/Polarization.h"
#include "Physics/Fracture.h"
#include "Core/ShapeFunctions.h"

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
    
    // Utilisation d'un format creux (Triplet) pour l'assemblage
    std::vector<Eigen::Triplet<double>> triplets;
    F_global.setZero(system_size);

    const auto& elements = mesh.get_elements();
    
    // Matrice de rigidité (Loi de Hooke)
    Eigen::Matrix3d C;
    C << math.c1, math.c2, 0,
         math.c2, math.c1, 0,
         0, 0, math.c3 / 2.0;

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_elem_nodes = elem.get_num_nodes(); // 4 pour un quad
        auto coords = mesh.get_element_coords(elem_idx);
        
        Eigen::MatrixXd K_local = Eigen::MatrixXd::Zero(2 * n_elem_nodes, 2 * n_elem_nodes);
        Eigen::VectorXd F_local = Eigen::VectorXd::Zero(2 * n_elem_nodes);

        const auto& indices = elem.get_node_indices();
        
        // Récupération des champs locaux
        Eigen::VectorXd Px_local(n_elem_nodes), Py_local(n_elem_nodes), v_local(n_elem_nodes);
        for (int i = 0; i < n_elem_nodes; ++i) {
            Px_local(i) = polarization.get_Px()[indices[i]];
            Py_local(i) = polarization.get_Py()[indices[i]];
            v_local(i) = fracture.get_v()[indices[i]];
        }

        const auto& gauss_points = Quadrature::get_gauss_2x2();

        for (const auto& gp : gauss_points) {
            auto N_std = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
            auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
            auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);

            double dV = gp.weight * detJ;

            // Construction de la matrice B (3 x 8)
            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2 * n_elem_nodes);
            for (int i = 0; i < n_elem_nodes; ++i) {
                B(0, 2 * i)     = dN_xy[0][i]; // dNi/dx
                B(1, 2 * i + 1) = dN_xy[1][i]; // dNi/dy
                B(2, 2 * i)     = dN_xy[1][i]; // dNi/dy
                B(2, 2 * i + 1) = dN_xy[0][i]; // dNi/dx
            }

            // Interpolation aux points de Gauss
            double p1 = 0, p2 = 0, v_gp = 0;
            for (int i = 0; i < n_elem_nodes; ++i) {
                p1 += N_std[i] * Px_local(i);
                p2 += N_std[i] * Py_local(i);
                v_gp += N_std[i] * v_local(i);
            }

            double penalite_fracture = (v_gp * v_gp) + math.eta_k; // (v^2 + eta_k)

            // Contraintes propres dues à la polarisation (dérivée de W_coup par rapport à eps)
            Eigen::Vector3d sigma_0;
            sigma_0(0) = - (math.b1 / 2.0) * p1 * p1 - (math.b2 / 2.0) * p2 * p2; // sigma_11^0
            sigma_0(1) = - (math.b2 / 2.0) * p1 * p1 - (math.b1 / 2.0) * p2 * p2; // sigma_22^0
            sigma_0(2) = - math.b3 * p1 * p2;                                     // sigma_12^0

            // Assemblage local
            K_local.noalias() += B.transpose() * (penalite_fracture * C) * B * dV;
            F_local.noalias() += B.transpose() * (penalite_fracture * sigma_0) * dV; 
        }

        // Assemblage global (format Triplet)
        for (int i = 0; i < n_elem_nodes; ++i) {
            for (int j = 0; j < n_elem_nodes; ++j) {
                triplets.push_back(Eigen::Triplet<double>(2 * indices[i], 2 * indices[j], K_local(2 * i, 2 * j)));
                triplets.push_back(Eigen::Triplet<double>(2 * indices[i], 2 * indices[j] + 1, K_local(2 * i, 2 * j + 1)));
                triplets.push_back(Eigen::Triplet<double>(2 * indices[i] + 1, 2 * indices[j], K_local(2 * i + 1, 2 * j)));
                triplets.push_back(Eigen::Triplet<double>(2 * indices[i] + 1, 2 * indices[j] + 1, K_local(2 * i + 1, 2 * j + 1)));
            }
            F_global(2 * indices[i])     += F_local(2 * i);
            F_global(2 * indices[i] + 1) += F_local(2 * i + 1);
        }
    }

    K_global.resize(system_size, system_size);
    K_global.setFromTriplets(triplets.begin(), triplets.end());

    // Appliquer les CL (à adapter selon la logique de votre projet)
    const auto& bcs_ux = bc_manager.get_ux_bcs();
    const auto& bcs_uy = bc_manager.get_uy_bcs();
    apply_boundary_conditions(K_global, F_global, bcs_ux, bcs_uy);

    // Résolution avec un solveur creux direct (Très rapide pour K symétrique)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K_global);
    Eigen::VectorXd U_new = solver.solve(F_global);

    // Mise à jour des vecteurs ux et uy
    for(size_t i = 0; i < n_nodes; ++i) {
        ux_current(i) = U_new(2 * i);
        uy_current(i) = U_new(2 * i + 1);
    }
}

Eigen::Matrix2d Mechanics::get_strain_at_gp(const Element& elem, const GaussPoint2D& gp, const std::vector<std::array<double, 2>>& coords) const {
    const auto& indices = elem.get_node_indices();
    int n_nodes = elem.get_num_nodes();

    auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
    auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);

    double eps_11 = 0.0, eps_22 = 0.0, eps_12 = 0.0;

    for (int i = 0; i < n_nodes; ++i) {
        double dN_dx = dN_xy[0][i];
        double dN_dy = dN_xy[1][i];
        
        eps_11 += dN_dx * ux_current[indices[i]];
        eps_22 += dN_dy * uy_current[indices[i]];
        
        eps_12 += dN_dy * ux_current[indices[i]] + dN_dx * uy_current[indices[i]];
    }
    eps_12 *= 0.5;

    Eigen::Matrix2d strain;
    strain << eps_11, eps_12,
              eps_12, eps_22;
    return strain;
}

void Mechanics::apply_boundary_conditions(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, 
                                          const std::vector<NodeBC>& bcs_x, const std::vector<NodeBC>& bcs_y) {
    const double penalty = 1e15; 

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