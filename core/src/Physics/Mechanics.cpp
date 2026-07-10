#include "Physics/Mechanics.h"
#include "Physics/MechanicsAssembler.h"
#include "Solvers/LinearSolver.h"
#include "Core/ShapeFunctions.h"
#include "Core/ElementIntegrator.h"
#include "Core/Types.h"
#include "Utils/Logger.h"
#include <omp.h>

Mechanics::Mechanics(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager)
    : config(config), mesh(mesh), bc_manager(bc_manager) 
{
    size_t n_nodes = mesh.get_num_nodes();
    ux_prev_iter.setZero(n_nodes);
    uy_prev_iter.setZero(n_nodes);
    ux_current.setZero(n_nodes);
    uy_current.setZero(n_nodes);

    if (config.mechanics.type == InitializationType::UNIFORM) {
        ux_current.setConstant(config.mechanics.val_x_0);
        uy_current.setConstant(config.mechanics.val_y_0);
    } else if (config.mechanics.type == InitializationType::RANDOM) {
        for (size_t i = 0; i < n_nodes; ++i) {
            ux_current[i] = static_cast<double>(rand()) / RAND_MAX;
            uy_current[i] = static_cast<double>(rand()) / RAND_MAX;
        }
    }
}

void Mechanics::update_u(double time, const Polarization& polarization, const Fracture& fracture, const Math& math) {
    (void)time;
    size_t n_nodes = mesh.get_num_nodes();
    size_t system_size = 2 * n_nodes;

    Eigen::SparseMatrix<double> K_global(system_size, system_size);
    Eigen::VectorXd F_global = Eigen::VectorXd::Zero(system_size);

    MechanicsAssembler::assemble_system(mesh, polarization, fracture, math, bc_manager.get_ux_bcs(), bc_manager.get_uy_bcs(), K_global, F_global);

    Eigen::VectorXd U_guess(system_size);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_nodes; ++i) {
        U_guess(2 * i)     = ux_current(i);
        U_guess(2 * i + 1) = uy_current(i);
    }

    Eigen::VectorXd U_new = LinearSolver::solve_with_guess(K_global, F_global, U_guess, config.simulation.debug_enabled);

    if (U_new.allFinite()) {
        map_global_vector_to_components(U_new);
    } else {
        Logger::error("[Mechanics] Echec resolution ou solution non finie");
    }
}

void Mechanics::map_global_vector_to_components(const Eigen::VectorXd& U_new) {
    size_t n_nodes = mesh.get_num_nodes();
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_nodes; ++i) {
        ux_current(i) = U_new(2 * i);
        uy_current(i) = U_new(2 * i + 1);
    }
}

double Mechanics::calculate_error() const { return (ux_current - ux_prev_iter).norm() + (uy_current - uy_prev_iter).norm(); }
void Mechanics::save_previous_iteration() { ux_prev_iter = ux_current; uy_prev_iter = uy_current; }
void Mechanics::save_previous_state() { ux_backup = ux_current; uy_backup = uy_current; }
void Mechanics::restore_previous_state() { ux_current = ux_backup; uy_current = uy_backup; }
void Mechanics::update_history() {}

Eigen::Matrix2d Mechanics::get_strain_at_gp(const Element& elem, const GaussPoint2D& gp, const std::vector<std::array<double, 2>>& coords) const {
    const auto& indices = elem.get_node_indices();
    int n_nodes = elem.get_num_nodes();
    double eps_11 = 0.0, eps_22 = 0.0, eps_12 = 0.0;

    if (n_nodes == 3) {
        auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);
        if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) return Eigen::Matrix2d::Zero();
        
        for (int i = 0; i < 3; ++i) {
            eps_11 += dN_xy[i][0] * ux_current[indices[i]];
            eps_22 += dN_xy[i][1] * uy_current[indices[i]];
            eps_12 += dN_xy[i][1] * ux_current[indices[i]] + dN_xy[i][0] * uy_current[indices[i]];
        }
    } else if (n_nodes == 4) {
        auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
        auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);
        if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) return Eigen::Matrix2d::Zero();
        
        for (int i = 0; i < n_nodes; ++i) {
            eps_11 += dN_xy[0][i] * ux_current[indices[i]];
            eps_22 += dN_xy[1][i] * uy_current[indices[i]];
            eps_12 += dN_xy[1][i] * ux_current[indices[i]] + dN_xy[0][i] * uy_current[indices[i]];
        }
    }

    eps_12 *= 0.5; 
    Eigen::Matrix2d strain;
    strain << eps_11, eps_12, eps_12, eps_22;
    return strain;
}