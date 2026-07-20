#include "Physics/Polarization.h"
#include "Physics/PolarizationAssembler.h"
#include "Solvers/LinearSolver.h"
#include "Utils/Logger.h"
#include <omp.h>
#include <algorithm>
#include <cstdlib>

Polarization::Polarization(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager)
    : config(config), mesh(mesh), bc_manager(bc_manager)
{
    size_t n_nodes = mesh.get_num_nodes();
    Px_prev_iter.setZero(n_nodes); Py_prev_iter.setZero(n_nodes);
    Px_current.setZero(n_nodes); Py_current.setZero(n_nodes);
    Px_n.setZero(n_nodes); Py_n.setZero(n_nodes);

    appliquer_conditions_initiales();
}

void Polarization::appliquer_conditions_initiales() {
    size_t n_nodes = mesh.get_num_nodes();
    if (config.polarization.type == InitializationType::UNIFORM) {
        Px_current.setConstant(config.polarization.val_x_0);
        Py_current.setConstant(config.polarization.val_y_0);
    } else if (config.polarization.type == InitializationType::RANDOM) {
        double amp_x = config.polarization.val_x_0; 
        double amp_y = config.polarization.val_y_0; 

        for (size_t i = 0; i < n_nodes; ++i) {
            Px_current[i] = ((static_cast<double>(rand()) / RAND_MAX) * 2.0 - 1.0) * amp_x;
            Py_current[i] = ((static_cast<double>(rand()) / RAND_MAX) * 2.0 - 1.0) * amp_y;
        }
    }
    Px_n = Px_current;
    Py_n = Py_current;
}

void Polarization::update_P(double time, double dt_relax,  const Fracture& fracture, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math)
{
    (void)time; 
    double dt = dt_relax;
    const int n_nodes = static_cast<int>(Px_current.size());
    const int system_size = 2 * n_nodes;

    Eigen::SparseMatrix<double> K_global(system_size, system_size);
    Eigen::VectorXd F_global = Eigen::VectorXd::Zero(system_size);

    PolarizationAssembler::assemble_system(
        dt, Px_current, Py_current, Px_current, Py_current,
        mesh, fracture, mechanics, electrostatics, math, config, bc_manager,
        K_global, F_global
    );

    Eigen::VectorXd P_guess(system_size);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n_nodes; ++i) {
        P_guess(2 * i)     = Px_current(i);
        P_guess(2 * i + 1) = Py_current(i);
    }

    Eigen::VectorXd P_new = LinearSolver::solve_with_guess(K_global, F_global, P_guess, config.simulation.debug_enabled);

    if (P_new.allFinite()) {
        map_global_vector_to_components(P_new);
    } else {
        Logger::error("[Polarization] Echec du solveur tangent monolithique : présence de NaN/Inf.");
    }
}

void Polarization::map_global_vector_to_components(const Eigen::VectorXd& P_new) {
    const int n_nodes = static_cast<int>(Px_current.size());
    const double omega = 1.0; 

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n_nodes; ++i) {
        Px_current(i) = omega * P_new(2 * i)     + (1.0 - omega) * Px_prev_iter(i);
        Py_current(i) = omega * P_new(2 * i + 1) + (1.0 - omega) * Py_prev_iter(i);
    }
}

double Polarization::calculate_error() const {
    return std::max((Px_current - Px_prev_iter).cwiseAbs().maxCoeff(), (Py_current - Py_prev_iter).cwiseAbs().maxCoeff());
}

void Polarization::save_previous_iteration() { Px_prev_iter = Px_current; Py_prev_iter = Py_current; }
void Polarization::save_previous_state() { Px_backup = Px_current; Py_backup = Py_current; }
void Polarization::restore_previous_state() { Px_current = Px_backup; Py_current = Py_backup; }
void Polarization::update_history() { Px_n = Px_current; Py_n = Py_current; }