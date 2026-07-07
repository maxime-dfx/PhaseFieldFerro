#include "Physics/FractureAssembler.h"
#include "Core/ElementIntegrator.h"
#include "Physics/Math.h"
#include "Physics/Mechanics.h"
#include "Physics/Polarization.h"
#include "Physics/Electrostatics.h"
#include <cmath>

void FractureAssembler::assemble_system(
    double dt, const Eigen::VectorXd& v_n,
    const Mesh& mesh, const Polarization& polarization, 
    const Mechanics& mechanics, const Electrostatics& electrostatics, 
    const Math& math, const Datafile& config,
    Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) 
{
    const auto& elements = mesh.get_elements();
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(elements.size() * 16);

    const int MAX_NODES = 4;
    Eigen::MatrixXd K_local(MAX_NODES, MAX_NODES);
    Eigen::VectorXd F_local(MAX_NODES);
    Eigen::VectorXd Px_local(MAX_NODES), Py_local(MAX_NODES), v_local_n(MAX_NODES);
    Eigen::RowVectorXd N(MAX_NODES);
    Eigen::MatrixXd grad_N(2, MAX_NODES);

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        const auto& indices = elem.get_node_indices();

        K_local.topLeftCorner(n_nodes, n_nodes).setZero();
        F_local.head(n_nodes).setZero();

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
            bool is_impermeable = (config.fracture_mode == CrackBCType::IMPERMEABLE);
            
            Eigen::Vector2d E_gp = Eigen::Vector2d::Zero();
            if (is_impermeable) {
                E_gp = Eigen::Vector2d(electrostatics.get_Ex_at_gp(elem, gp), electrostatics.get_Ey_at_gp(elem, gp));
            }

            double H_drive_current = math.compute_H_drive(Pij_gp, Pi_gp, eps_gp, E_gp, is_impermeable);
            double H_drive = std::max(0.0, H_drive_current);

            double mass_coeff = (config.material.mu_v / dt) + (config.material.Gc / (2.0 * config.material.kappa)) + 2.0 * H_drive;
            double diff_coeff = 2.0 * config.material.Gc * config.material.kappa;
            double rhs_coeff  = (config.material.mu_v / dt) * N.head(n).dot(v_local_n.head(n)) + (config.material.Gc / (2.0 * config.material.kappa));

            auto N_active = N.head(n);
            K_local.topLeftCorner(n, n).noalias() += (mass_coeff * N_active.transpose() * N_active + diff_coeff * grad_N_active.transpose() * grad_N_active) * dV;
            F_local.head(n).noalias() += (rhs_coeff * N_active.transpose()) * dV;
        });

        for (int i = 0; i < n_nodes; ++i) {
            for (int j = 0; j < n_nodes; ++j) {
                triplets.emplace_back(indices[i], indices[j], K_local(i, j));
            }
            F_global(indices[i]) += F_local(i);
        }
    }
    K_global.setFromTriplets(triplets.begin(), triplets.end());
}

void FractureAssembler::apply_irreversibility(const Mesh& mesh, const Eigen::VectorXd& v_prev_iter, Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) {
    double alpha = 2e-2; 
    for (size_t i = 0; i < mesh.get_num_nodes(); ++i) {
        if (v_prev_iter(i) <= alpha) {
            K_global.coeffRef(i, i) += 1e15;
            F_global(i) += 1e15 * 0.0; 
        }
    }
    K_global.makeCompressed();
}