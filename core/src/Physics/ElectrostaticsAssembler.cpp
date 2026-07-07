#include "Physics/ElectrostaticsAssembler.h"
#include "Core/ElementIntegrator.h"
#include "Physics/Math.h"
#include "Physics/Polarization.h"
#include "Physics/Fracture.h"
#include <cmath>

void ElectrostaticsAssembler::assemble_system(
    const Mesh& mesh, const Polarization& polarization, const Fracture& fracture, 
    const Math& math, const Datafile& config, 
    Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) 
{
    (void)math; // Suppression du warning "unused parameter"

    const auto& elements = mesh.get_elements();
    double eps0 = config.material.eps0;
    double eta_k = config.material.eta_k;

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(elements.size() * 16);

    const int MAX_NODES = 4;
    Eigen::MatrixXd K_local(MAX_NODES, MAX_NODES);
    Eigen::VectorXd F_local(MAX_NODES);
    Eigen::VectorXd Px_loc(MAX_NODES), Py_loc(MAX_NODES), v_loc(MAX_NODES);
    Eigen::RowVectorXd N(MAX_NODES);
    Eigen::MatrixXd grad_N(2, MAX_NODES);

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        const auto& indices = elem.get_node_indices();

        K_local.topLeftCorner(n_nodes, n_nodes).setZero();
        F_local.head(n_nodes).setZero();

        for(int i = 0; i < n_nodes; ++i) {
            Px_loc(i) = polarization.get_Px()[indices[i]];
            Py_loc(i) = polarization.get_Py()[indices[i]];
            v_loc(i)  = fracture.get_v()[indices[i]];
        }

        ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
            (void)gp; // Suppression du warning "unused parameter"

            double v_gp = N.head(n).dot(v_loc.head(n));
            double px_gp = N.head(n).dot(Px_loc.head(n));
            double py_gp = N.head(n).dot(Py_loc.head(n));

            double phase_factor = 1.0;
            if (config.fracture_mode == CrackBCType::IMPERMEABLE) {
                phase_factor = (v_gp * v_gp) + eta_k;
            }

            double eps_eff = eps0 * phase_factor;
            Eigen::Vector2d P_eff(px_gp * phase_factor, py_gp * phase_factor);

            auto grad_N_active = grad_N.leftCols(n); 
            K_local.topLeftCorner(n, n).noalias() += grad_N_active.transpose() * eps_eff * grad_N_active * dV;
            F_local.head(n).noalias() += grad_N_active.transpose() * P_eff * dV;
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

void ElectrostaticsAssembler::apply_boundary_conditions(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, const std::vector<NodeBC>& bcs) {
    double max_diag = 0.0;
    for (int k = 0; k < K.outerSize(); ++k) {
        max_diag = std::max(max_diag, std::abs(K.coeff(k, k)));
    }
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