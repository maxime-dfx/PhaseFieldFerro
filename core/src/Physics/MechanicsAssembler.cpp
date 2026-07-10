#include "Physics/MechanicsAssembler.h"
#include "Core/ElementIntegrator.h"
#include <omp.h>

void MechanicsAssembler::assemble_system(
    const Mesh& mesh, const Polarization& polarization, const Fracture& fracture, 
    const Math& math, const std::vector<NodeBC>& bcs_x, const std::vector<NodeBC>& bcs_y,
    Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) 
{
    const auto& elements = mesh.get_elements();
    const int n_dof = static_cast<int>(mesh.get_num_nodes());
    const int system_size = 2 * n_dof;

    Eigen::Matrix3d C = math.get_elastic_matrix(); 

    int num_threads = omp_get_max_threads();
    std::vector<std::vector<Eigen::Triplet<double>>> thread_triplets(num_threads);
    std::vector<Eigen::VectorXd> thread_F(num_threads, Eigen::VectorXd::Zero(system_size));
    std::vector<Eigen::VectorXd> thread_diag(num_threads, Eigen::VectorXd::Zero(system_size));

    for (int t = 0; t < num_threads; ++t) {
        thread_triplets[t].reserve((elements.size() * 64) / num_threads + 64);
    }

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        constexpr int MAX_NODES = 4;
        constexpr int MAX_ELEMENT_DOF = 2 * MAX_NODES;
        Eigen::Matrix<double, MAX_ELEMENT_DOF, MAX_ELEMENT_DOF> K_local;
        Eigen::Matrix<double, MAX_ELEMENT_DOF, 1> F_local;
        Eigen::Matrix<double, MAX_NODES, 1> Px_local, Py_local, v_local;
        Eigen::Matrix<double, 3, MAX_ELEMENT_DOF> B;

        // Allocation dynamique UNIQUE par thread pour satisfaire le template
        Eigen::RowVectorXd N(MAX_NODES);
        Eigen::MatrixXd grad_N(2, MAX_NODES);

        #pragma omp for schedule(static)
        for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
            const auto& elem = elements[elem_idx];
            const int n_nodes = elem.get_num_nodes();
            auto coords = mesh.get_element_coords(elem_idx);
            const auto& indices = elem.get_node_indices();
            const int n_local_dof = 2 * n_nodes;

            K_local.setZero();
            F_local.setZero();
            
            for (int i = 0; i < n_nodes; ++i) {
                Px_local(i) = polarization.get_Px()[indices[i]];
                Py_local(i) = polarization.get_Py()[indices[i]];
                v_local(i)  = fracture.get_v()[indices[i]];
            }

            ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
                (void)gp;
                B.setZero();
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
                double factor = penalite_fracture * dV;
                K_local.topLeftCorner(n_local_dof, n_local_dof).noalias() += (B_active.transpose() * C * B_active) * factor;
                F_local.head(n_local_dof).noalias()                       -= B_active.transpose() * sigma_0 * factor;
            });

            for (int i = 0; i < n_nodes; ++i) {
                int g_i_x = 2 * indices[i], g_i_y = 2 * indices[i] + 1;
                int l_i_x = 2 * i, l_i_y = 2 * i + 1;

                for (int j = 0; j < n_nodes; ++j) {
                    int g_j_x = 2 * indices[j], g_j_y = 2 * indices[j] + 1;
                    int l_j_x = 2 * j, l_j_y = 2 * j + 1;

                    thread_triplets[tid].emplace_back(g_i_x, g_j_x, K_local(l_i_x, l_j_x));
                    thread_triplets[tid].emplace_back(g_i_x, g_j_y, K_local(l_i_x, l_j_y));
                    thread_triplets[tid].emplace_back(g_i_y, g_j_x, K_local(l_i_y, l_j_x));
                    thread_triplets[tid].emplace_back(g_i_y, g_j_y, K_local(l_i_y, l_j_y));
                }
                
                thread_F[tid](g_i_x) += F_local(l_i_x);
                thread_F[tid](g_i_y) += F_local(l_i_y);
                thread_diag[tid](g_i_x) += K_local(l_i_x, l_i_x);
                thread_diag[tid](g_i_y) += K_local(l_i_y, l_i_y);
            }
        }
    }

    size_t total_triplets = 0;
    for (int t = 0; t < num_threads; ++t) total_triplets += thread_triplets[t].size();
    total_triplets += bcs_x.size() + bcs_y.size();

    std::vector<Eigen::Triplet<double>> global_triplets;
    global_triplets.reserve(total_triplets);

    F_global.setZero(system_size);
    Eigen::VectorXd diag_global = Eigen::VectorXd::Zero(system_size);

    for (int t = 0; t < num_threads; ++t) {
        global_triplets.insert(global_triplets.end(), thread_triplets[t].begin(), thread_triplets[t].end());
        F_global += thread_F[t];
        diag_global += thread_diag[t];
    }

    double max_diag = diag_global.cwiseAbs().maxCoeff();
    if (max_diag < 1e-12) max_diag = 1.0;
    const double penalty = max_diag * 1e4;

    for (size_t i = 0; i < bcs_x.size(); ++i) {
        if (bcs_x[i].type == BCType::DIRICHLET) {
            int dof = 2 * static_cast<int>(i);
            global_triplets.emplace_back(dof, dof, penalty);
            F_global(dof) += penalty * bcs_x[i].value;
        }
    }
    for (size_t i = 0; i < bcs_y.size(); ++i) {
        if (bcs_y[i].type == BCType::DIRICHLET) {
            int dof = 2 * static_cast<int>(i) + 1;
            global_triplets.emplace_back(dof, dof, penalty);
            F_global(dof) += penalty * bcs_y[i].value;
        }
    }

    K_global.setFromTriplets(global_triplets.begin(), global_triplets.end());
}