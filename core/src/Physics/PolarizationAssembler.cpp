#include "Physics/PolarizationAssembler.h"
#include "Core/ElementIntegrator.h"
#include "Utils/Logger.h"
#include <omp.h>
#include <cmath>

void PolarizationAssembler::assemble_system(
    double dt,
    const Eigen::VectorXd& Px_current, const Eigen::VectorXd& Py_current,
    const Eigen::VectorXd& Px_n, const Eigen::VectorXd& Py_n,
    const Mesh& mesh, const Fracture& fracture, const Mechanics& mechanics,
    const Electrostatics& electrostatics, const Math& math, const Datafile& config,
    const BoundaryManager& bc_manager,
    Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global)
{
    const auto& elements = mesh.get_elements();
    const int n_dof = mesh.get_num_nodes();
    const int system_size = 2 * n_dof;

    bool is_impermeable = (config.fracture_mode == CrackBCType::IMPERMEABLE);

    double alpha1 = config.material.alpha_1, alpha11 = config.material.alpha_11, alpha12 = config.material.alpha_12;
    double alpha111 = config.material.alpha_111, alpha112 = config.material.alpha_112;
    double alpha1111 = config.material.alpha_1111, alpha1112 = config.material.alpha_1112, alpha1122 = config.material.alpha_1122;

    int num_threads = omp_get_max_threads();
    std::vector<std::vector<Eigen::Triplet<double>>> thread_triplets(num_threads);
    std::vector<Eigen::VectorXd> thread_F(num_threads, Eigen::VectorXd::Zero(system_size));
    std::vector<Eigen::VectorXd> thread_diag_px(num_threads, Eigen::VectorXd::Zero(n_dof));
    std::vector<Eigen::VectorXd> thread_diag_py(num_threads, Eigen::VectorXd::Zero(n_dof));

    for (int t = 0; t < num_threads; ++t) thread_triplets[t].reserve((elements.size() * 64) / num_threads + 64);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        constexpr int MAX_NODES = 4;
        constexpr int MAX_ELEMENT_DOF = 2 * MAX_NODES;
        Eigen::Matrix<double, MAX_ELEMENT_DOF, MAX_ELEMENT_DOF> K_local;
        Eigen::Matrix<double, MAX_ELEMENT_DOF, 1> F_local;
        Eigen::Matrix<double, MAX_NODES, 1> Px_local, Py_local, v_local, Px_local_n, Py_local_n;

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
                Px_local(i) = Px_current[indices[i]]; Py_local(i) = Py_current[indices[i]];
                v_local(i) = fracture.get_v()[indices[i]];
                Px_local_n(i) = Px_n[indices[i]]; Py_local_n(i) = Py_n[indices[i]];
            }

            ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
                double v_gp = N.head(n).dot(v_local.head(n));
                Eigen::Vector2d Pi_gp(N.head(n).dot(Px_local.head(n)), N.head(n).dot(Py_local.head(n)));
                Eigen::Matrix2d eps_gp = mechanics.get_strain_at_gp(elem, gp, coords);
                double Ex_gp = electrostatics.get_Ex_at_gp(elem, gp);
                double Ey_gp = electrostatics.get_Ey_at_gp(elem, gp);

                double penalite_fracture = (v_gp * v_gp) + config.material.eta_k;

                double force_px = math.compute_polarization_force(0, Pi_gp, eps_gp, Ex_gp, penalite_fracture, is_impermeable);
                double force_py = math.compute_polarization_force(1, Pi_gp, eps_gp, Ey_gp, penalite_fracture, is_impermeable);

                double p1 = Pi_gp(0), p2 = Pi_gp(1);
                double p1_2 = p1 * p1, p1_4 = p1_2 * p1_2, p1_6 = p1_4 * p1_2;
                double p2_2 = p2 * p2, p2_4 = p2_2 * p2_2, p2_6 = p2_4 * p2_2;

                double d2chi_dp12 = 2.0 * alpha1 + 12.0 * alpha11 * p1_2 + 2.0 * alpha12 * p2_2 
                                  + 30.0 * alpha111 * p1_4 + 12.0 * alpha112 * p1_2 * p2_2 + 2.0 * alpha112 * p2_4 
                                  + 56.0 * alpha1111 * p1_6 + 30.0 * alpha1112 * p1_4 * p2_2 + 2.0 * alpha1112 * p2_6 
                                  + 12.0 * alpha1122 * p1_2 * p2_4;

                double d2chi_dp22 = 2.0 * alpha1 + 12.0 * alpha11 * p2_2 + 2.0 * alpha12 * p1_2 
                                  + 30.0 * alpha111 * p2_4 + 12.0 * alpha112 * p2_2 * p1_2 + 2.0 * alpha112 * p1_4 
                                  + 56.0 * alpha1111 * p2_6 + 30.0 * alpha1112 * p2_4 * p1_2 + 2.0 * alpha1112 * p1_6 
                                  + 12.0 * alpha1122 * p2_2 * p1_4;

                double d2chi_dp1dp2 = 4.0 * alpha12 * p1 * p2 + 8.0 * alpha112 * p1 * p2 * (p1_2 + p2_2) 
                                    + 12.0 * alpha1112 * p1 * p2 * (p1_4 + p2_4) + 16.0 * alpha1122 * p1_2 * p1 * p2_2 * p2;

                double d2W_dp12 = -config.material.b1 * eps_gp(0,0) - config.material.b2 * eps_gp(1,1);
                double d2W_dp22 = -config.material.b1 * eps_gp(1,1) - config.material.b2 * eps_gp(0,0);
                double d2W_dp1dp2 = -config.material.b3 * (eps_gp(1,0) + eps_gp(0,1));

                double J_11 = penalite_fracture * d2W_dp12 + d2chi_dp12;
                double J_22 = penalite_fracture * d2W_dp22 + d2chi_dp22;
                double J_12 = penalite_fracture * d2W_dp1dp2 + d2chi_dp1dp2;

                double Px_n_gp = N.head(n).dot(Px_local_n.head(n));
                double Py_n_gp = N.head(n).dot(Py_local_n.head(n));

                for (int i = 0; i < n; ++i) {
                    F_local(2 * i)     += ( (config.material.mu_p / dt) * Px_n_gp + J_11 * p1 + J_12 * p2 - force_px ) * N[i] * dV;
                    F_local(2 * i + 1) += ( (config.material.mu_p / dt) * Py_n_gp + J_12 * p1 + J_22 * p2 - force_py ) * N[i] * dV;

                    for (int j = 0; j < n; ++j) {
                        double mass = (config.material.mu_p / dt) * N[i] * N[j] * dV;
                        double stiff = (config.material.a0 * penalite_fracture) * (grad_N(0, i) * grad_N(0, j) + grad_N(1, i) * grad_N(1, j)) * dV;

                        K_local(2 * i, 2 * j)         += mass + stiff + J_11 * N[i] * N[j] * dV;
                        K_local(2 * i, 2 * j + 1)     += J_12 * N[i] * N[j] * dV;
                        K_local(2 * i + 1, 2 * j)     += J_12 * N[i] * N[j] * dV;
                        K_local(2 * i + 1, 2 * j + 1) += mass + stiff + J_22 * N[i] * N[j] * dV;
                    }
                }
            });

            for (int i = 0; i < n_nodes; ++i) {
                int g_i_px = 2 * indices[i], g_i_py = 2 * indices[i] + 1;
                int l_i_px = 2 * i, l_i_py = 2 * i + 1;

                for (int j = 0; j < n_nodes; ++j) {
                    int g_j_px = 2 * indices[j], g_j_py = 2 * indices[j] + 1;
                    int l_j_px = 2 * j, l_j_py = 2 * j + 1;

                    thread_triplets[tid].emplace_back(g_i_px, g_j_px, K_local(l_i_px, l_j_px));
                    thread_triplets[tid].emplace_back(g_i_px, g_j_py, K_local(l_i_px, l_j_py));
                    thread_triplets[tid].emplace_back(g_i_py, g_j_px, K_local(l_i_py, l_j_px));
                    thread_triplets[tid].emplace_back(g_i_py, g_j_py, K_local(l_i_py, l_j_py));
                }
                
                thread_F[tid](g_i_px) += F_local(l_i_px);
                thread_F[tid](g_i_py) += F_local(l_i_py);
                thread_diag_px[tid](indices[i]) += K_local(l_i_px, l_i_px);
                thread_diag_py[tid](indices[i]) += K_local(l_i_py, l_i_py);
            }
        }
    }

    const auto& bcs_px = bc_manager.get_px_bcs();
    const auto& bcs_py = bc_manager.get_py_bcs();

    size_t total_triplets = 0;
    for (int t = 0; t < num_threads; ++t) total_triplets += thread_triplets[t].size();
    total_triplets += 2 * n_dof + bcs_px.size() + bcs_py.size();

    std::vector<Eigen::Triplet<double>> global_triplets;
    global_triplets.reserve(total_triplets);

    F_global.setZero(system_size);
    Eigen::VectorXd diag_global_px = Eigen::VectorXd::Zero(n_dof);
    Eigen::VectorXd diag_global_py = Eigen::VectorXd::Zero(n_dof);

    for (int t = 0; t < num_threads; ++t) {
        global_triplets.insert(global_triplets.end(), thread_triplets[t].begin(), thread_triplets[t].end());
        F_global += thread_F[t];
        diag_global_px += thread_diag_px[t];
        diag_global_py += thread_diag_py[t];
    }

    bool debug_enabled = config.simulation.debug_enabled;
    for (int k = 0; k < n_dof; ++k) {
        bool is_dirichlet_px = (static_cast<size_t>(k) < bcs_px.size()) && (bcs_px[k].type == BCType::DIRICHLET);
        if (!is_dirichlet_px && std::abs(diag_global_px(k)) < 1e-12) {
            global_triplets.emplace_back(2 * k, 2 * k, 1.0);
            F_global(2 * k) = 0.0;
            Logger::debug("[Polarization] DOF Px flottant force a 1.0", debug_enabled);
        }

        bool is_dirichlet_py = (static_cast<size_t>(k) < bcs_py.size()) && (bcs_py[k].type == BCType::DIRICHLET);
        if (!is_dirichlet_py && std::abs(diag_global_py(k)) < 1e-12) {
            global_triplets.emplace_back(2 * k + 1, 2 * k + 1, 1.0);
            F_global(2 * k + 1) = 0.0;
            Logger::debug("[Polarization] DOF Py flottant force a 1.0", debug_enabled);
        }
    }

    double max_diag = std::max(diag_global_px.cwiseAbs().maxCoeff(), diag_global_py.cwiseAbs().maxCoeff());
    if (max_diag < 1e-12) max_diag = 1.0;
    const double penalty = max_diag * 1e7;

    for (size_t i = 0; i < bcs_px.size(); ++i) {
        if (bcs_px[i].type == BCType::DIRICHLET) {
            int dof = 2 * static_cast<int>(i);
            global_triplets.emplace_back(dof, dof, penalty);
            F_global(dof) += penalty * bcs_px[i].value;
        }
    }
    for (size_t i = 0; i < bcs_py.size(); ++i) {
        if (bcs_py[i].type == BCType::DIRICHLET) {
            int dof = 2 * static_cast<int>(i) + 1;
            global_triplets.emplace_back(dof, dof, penalty);
            F_global(dof) += penalty * bcs_py[i].value;
        }
    }

    K_global.setFromTriplets(global_triplets.begin(), global_triplets.end());
}