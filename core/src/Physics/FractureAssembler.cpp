#include "Physics/FractureAssembler.h"
#include "Core/ElementIntegrator.h"
#include <omp.h>
#include <algorithm>

void FractureAssembler::assemble_system(
    double dt, const Eigen::VectorXd& v_n,
    const Mesh& mesh, const Polarization& polarization, 
    const Mechanics& mechanics, const Electrostatics& electrostatics, 
    const Math& math, const Datafile& config,
    Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) 
{
    const auto& elements = mesh.get_elements();
    const int n_dof = static_cast<int>(mesh.get_num_nodes());

    bool is_impermeable = (config.fracture_mode == CrackBCType::IMPERMEABLE);
    double mu_v   = config.material.mu_v;
    double Gc     = config.material.Gc;
    double kappa  = config.material.kappa;

    int num_threads = omp_get_max_threads();
    std::vector<std::vector<Eigen::Triplet<double>>> thread_triplets(num_threads);
    std::vector<Eigen::VectorXd> thread_F(num_threads, Eigen::VectorXd::Zero(n_dof));

    for (int t = 0; t < num_threads; ++t) {
        thread_triplets[t].reserve((elements.size() * 16) / num_threads + 16);
    }

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        constexpr int MAX_NODES = 4;
        Eigen::Matrix<double, MAX_NODES, MAX_NODES> K_local;
        Eigen::Matrix<double, MAX_NODES, 1> F_local, Px_local, Py_local, v_local_n;

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
                Px_local(i)  = polarization.get_Px()[indices[i]];
                Py_local(i)  = polarization.get_Py()[indices[i]];
                v_local_n(i) = v_n[indices[i]];
            }

            ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
                Eigen::Vector2d Pi_gp(N.head(n).dot(Px_local.head(n)), N.head(n).dot(Py_local.head(n)));
                auto grad_N_active = grad_N.leftCols(n);

                Eigen::Matrix2d Pij_gp;
                Pij_gp(0, 0) = grad_N_active.row(0).dot(Px_local.head(n));
                Pij_gp(0, 1) = grad_N_active.row(1).dot(Px_local.head(n));
                Pij_gp(1, 0) = grad_N_active.row(0).dot(Py_local.head(n));
                Pij_gp(1, 1) = grad_N_active.row(1).dot(Py_local.head(n));

                Eigen::Matrix2d eps_gp = mechanics.get_strain_at_gp(elem, gp, coords);
                Eigen::Vector2d E_gp = Eigen::Vector2d::Zero();
                if (is_impermeable) {
                    E_gp = Eigen::Vector2d(electrostatics.get_Ex_at_gp(elem, gp), electrostatics.get_Ey_at_gp(elem, gp));
                }

                double H_drive = std::max(0.0, math.compute_H_drive(Pij_gp, Pi_gp, eps_gp, E_gp, is_impermeable));

                double mass_coeff = (mu_v / dt) + (Gc / (2.0 * kappa)) + 2.0 * H_drive;
                double diff_coeff = 2.0 * Gc * kappa;
                double rhs_coeff  = (mu_v / dt) * N.head(n).dot(v_local_n.head(n)) + (Gc / (2.0 * kappa));

                auto N_active = N.head(n);
                K_local.topLeftCorner(n, n).noalias() += (N_active.transpose() * N_active) * (mass_coeff * dV);
                K_local.topLeftCorner(n, n).noalias() += (grad_N_active.transpose() * grad_N_active) * (diff_coeff * dV);
                F_local.head(n).noalias()             += N_active.transpose() * (rhs_coeff * dV);
            });

            for (int i = 0; i < n_nodes; ++i) {
                for (int j = 0; j < n_nodes; ++j) {
                    thread_triplets[tid].emplace_back(indices[i], indices[j], K_local(i, j));
                }
                thread_F[tid](indices[i]) += F_local(i);
            }
        }
    }

    size_t total_triplets = 0;
    for (int t = 0; t < num_threads; ++t) total_triplets += thread_triplets[t].size();
    total_triplets += n_dof;

    std::vector<Eigen::Triplet<double>> global_triplets;
    global_triplets.reserve(total_triplets);
    F_global.setZero(n_dof);

    for (int t = 0; t < num_threads; ++t) {
        global_triplets.insert(global_triplets.end(), thread_triplets[t].begin(), thread_triplets[t].end());
        F_global += thread_F[t];
    }

    const double alpha = 2e-2; 
    for (int i = 0; i < n_dof; ++i) { 
        if (v_n(i) <= alpha) {
            global_triplets.emplace_back(i, i, 1e15);
            F_global(i) = 0.0; 
        }
    }

    K_global.setFromTriplets(global_triplets.begin(), global_triplets.end());
}