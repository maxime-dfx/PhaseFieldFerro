#include "Physics/MechanicsAssembler.h"
#include "Core/ElementIntegrator.h"
#include <cmath>

void MechanicsAssembler::assemble_system(
    const Mesh& mesh, const Polarization& polarization, const Fracture& fracture, 
    const Math& math, Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) 
{
    const auto& elements = mesh.get_elements();
    Eigen::Matrix3d C = math.get_elastic_matrix(); 
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(elements.size() * 64);

    const int MAX_NODES = 4;
    Eigen::MatrixXd K_local(2 * MAX_NODES, 2 * MAX_NODES);
    Eigen::VectorXd F_local(2 * MAX_NODES);
    Eigen::VectorXd Px_local(MAX_NODES), Py_local(MAX_NODES), v_local(MAX_NODES);
    
    Eigen::RowVectorXd N(MAX_NODES);
    Eigen::MatrixXd grad_N(2, MAX_NODES);
    Eigen::MatrixXd B(3, 2 * MAX_NODES);

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        const auto& indices = elem.get_node_indices();

        K_local.topLeftCorner(2 * n_nodes, 2 * n_nodes).setZero();
        F_local.head(2 * n_nodes).setZero();
        
        for (int i = 0; i < n_nodes; ++i) {
            Px_local(i) = polarization.get_Px()[indices[i]];
            Py_local(i) = polarization.get_Py()[indices[i]];
            v_local(i)  = fracture.get_v()[indices[i]];
        }

        ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
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
            F_local.head(2 * n).noalias() -= B_active.transpose() * (penalite_fracture * sigma_0) * dV;
        });

        for (int i = 0; i < n_nodes; ++i) {
            for (int j = 0; j < n_nodes; ++j) {
                triplets.emplace_back(2 * indices[i],     2 * indices[j],     K_local(2 * i,     2 * j));
                triplets.emplace_back(2 * indices[i],     2 * indices[j] + 1, K_local(2 * i,     2 * j + 1));
                triplets.emplace_back(2 * indices[i] + 1, 2 * indices[j],     K_local(2 * i + 1, 2 * j));
                triplets.emplace_back(2 * indices[i] + 1, 2 * indices[j] + 1, K_local(2 * i + 1, 2 * j + 1));
            }
            F_global(2 * indices[i])     += F_local(2 * i);
            F_global(2 * indices[i] + 1) += F_local(2 * i + 1);
        }
    }
    K_global.setFromTriplets(triplets.begin(), triplets.end());
}

void MechanicsAssembler::apply_boundary_conditions(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, 
                                                   const std::vector<NodeBC>& bcs_x, const std::vector<NodeBC>& bcs_y) {
    
    double max_diag = K.diagonal().cwiseAbs().maxCoeff();

    // --- LE BOUCLIER ANTI-ZERO ---
    if (max_diag < 1e-12) {
        // std::cout << "\n[ALERTE ROUGE] La matrice de rigidite K est VIDE (max_diag = 0) !\n";
        // std::cout << "-> L'assemblage n'a rien calcule (Integrateur inactif ou Matrice C nulle).\n";
        max_diag = 1.0; // On force une valeur pour que la pénalité fonctionne quand meme !
    }

    const double penalty = max_diag * 1e7;
    int forced_x = 0, forced_y = 0;

    for (size_t i = 0; i < bcs_x.size(); ++i) {
        if (bcs_x[i].type == BCType::DIRICHLET) {
            int dof = 2 * static_cast<int>(i);
            K.coeffRef(dof, dof) += penalty;
            F(dof) += penalty * bcs_x[i].value;
            forced_x++;
        }
    }
    for (size_t i = 0; i < bcs_y.size(); ++i) {
        if (bcs_y[i].type == BCType::DIRICHLET) {
            int dof = 2 * static_cast<int>(i) + 1;
            K.coeffRef(dof, dof) += penalty;
            F(dof) += penalty * bcs_y[i].value;
            forced_y++;
        }
    }

    K.makeCompressed();
}