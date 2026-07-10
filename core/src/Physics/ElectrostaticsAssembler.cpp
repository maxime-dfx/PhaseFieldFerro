#include "Physics/ElectrostaticsAssembler.h"
#include "Core/ElementIntegrator.h"
#include <omp.h>

void ElectrostaticsAssembler::assemble_system(
    const Mesh& mesh, const Polarization& polarization, const Fracture& fracture, 
    const Math& math, const Datafile& config, const std::vector<NodeBC>& bcs,
    Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) 
{
    (void)math;
    const auto& elements = mesh.get_elements();
    const int n_dof = static_cast<int>(mesh.get_num_nodes());

    double eps0 = config.material.eps0;
    double eta_k = config.material.eta_k;
    bool is_impermeable = (config.fracture_mode == CrackBCType::IMPERMEABLE);

    int num_threads = omp_get_max_threads();
    std::vector<std::vector<Eigen::Triplet<double>>> thread_triplets(num_threads);
    std::vector<Eigen::VectorXd> thread_F(num_threads, Eigen::VectorXd::Zero(n_dof));
    std::vector<Eigen::VectorXd> thread_diag(num_threads, Eigen::VectorXd::Zero(n_dof));

    for (int t = 0; t < num_threads; ++t) {
        thread_triplets[t].reserve((elements.size() * 16) / num_threads + 16);
    }

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        constexpr int MAX_NODES = 4;
        Eigen::Matrix<double, MAX_NODES, MAX_NODES> K_local;
        Eigen::Matrix<double, MAX_NODES, 1> F_local, Px_local, Py_local, v_local;

        // Allocation dynamique UNIQUE par thread pour satisfaire le template
        Eigen::RowVectorXd N(MAX_NODES);
        Eigen::MatrixXd grad_N(2, MAX_NODES);

        #pragma omp for schedule(static)
        for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
            const auto& elem = elements[elem_idx];
            const int n_nodes = elem.get_num_nodes();
            auto coords = mesh.get_element_coords(elem_idx);
            const auto& indices = elem.get_node_indices();

            K_local.setZero();
            F_local.setZero();

            for (int i = 0; i < n_nodes; ++i) {
                Px_local(i) = polarization.get_Px()[indices[i]];
                Py_local(i) = polarization.get_Py()[indices[i]];
                v_local(i)  = fracture.get_v()[indices[i]];
            }

            ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
                (void)gp;
                double v_gp = N.head(n).dot(v_local.head(n));
                double px_gp = N.head(n).dot(Px_local.head(n));
                double py_gp = N.head(n).dot(Py_local.head(n));

                double phase_factor = is_impermeable ? (v_gp * v_gp) + eta_k : 1.0;
                double eps_eff = eps0 * phase_factor;
                Eigen::Vector2d P_eff(px_gp * phase_factor, py_gp * phase_factor);

                auto grad_N_active = grad_N.leftCols(n); 
                
                K_local.topLeftCorner(n, n).noalias() += (grad_N_active.transpose() * grad_N_active) * (eps_eff * dV);
                F_local.head(n).noalias()             += grad_N_active.transpose() * P_eff * dV;
            });
            
            for (int i = 0; i < n_nodes; ++i) {
                for (int j = 0; j < n_nodes; ++j) {
                    thread_triplets[tid].emplace_back(indices[i], indices[j], K_local(i, j));
                    if (i == j) thread_diag[tid](indices[i]) += K_local(i, i);
                }
                thread_F[tid](indices[i]) += F_local(i);
            }
        }
    }

    size_t total_triplets = 0;
    for (int t = 0; t < num_threads; ++t) total_triplets += thread_triplets[t].size();
    total_triplets += bcs.size();

    std::vector<Eigen::Triplet<double>> global_triplets;
    global_triplets.reserve(total_triplets);

    F_global.setZero(n_dof);
    Eigen::VectorXd diag_global = Eigen::VectorXd::Zero(n_dof);

    for (int t = 0; t < num_threads; ++t) {
        global_triplets.insert(global_triplets.end(), thread_triplets[t].begin(), thread_triplets[t].end());
        F_global += thread_F[t];
        diag_global += thread_diag[t];
    }

    double max_diag = diag_global.cwiseAbs().maxCoeff();
    if (max_diag < 1e-12) max_diag = 1.0;
    const double penalty = max_diag * 1e7;

    for (size_t i = 0; i < bcs.size(); ++i) {
        if (bcs[i].type == BCType::DIRICHLET) {
            global_triplets.emplace_back(static_cast<int>(i), static_cast<int>(i), penalty);
            F_global(i) += penalty * bcs[i].value;
        }
    }
    K_global.setFromTriplets(global_triplets.begin(), global_triplets.end());
}