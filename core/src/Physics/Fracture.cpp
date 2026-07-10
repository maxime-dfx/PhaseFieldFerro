#include "Physics/Fracture.h"
#include "Physics/FractureAssembler.h"
#include "Solvers/LinearSolver.h"
#include "Utils/Logger.h"
#include <omp.h>
#include <algorithm>

Fracture::Fracture(const Datafile& config, const Mesh& mesh) : config(config), mesh(mesh) {
    size_t n_nodes = mesh.get_num_nodes();
    v_prev_iter.setOnes(n_nodes);
    v_current.setOnes(n_nodes); 
    v_n.setOnes(n_nodes);
}

void Fracture::update_v(double dt, const Polarization& polarization, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math) {
    size_t n_nodes = mesh.get_num_nodes();
    Eigen::SparseMatrix<double> K_global(n_nodes, n_nodes);
    Eigen::VectorXd F_global = Eigen::VectorXd::Zero(n_nodes);

    FractureAssembler::assemble_system(dt, v_n, mesh, polarization, mechanics, electrostatics, math, config, K_global, F_global);

    Eigen::VectorXd v_new = LinearSolver::solve_with_guess(K_global, F_global, v_current, config.simulation.debug_enabled);

    if (v_new.allFinite()) {
        enforce_physical_bounds(v_new);
    } else {
        Logger::error("[Fracture] Echec resolution ou solution non finie");
    }
}

void Fracture::enforce_physical_bounds(const Eigen::VectorXd& v_new) {
    size_t n_nodes = mesh.get_num_nodes();
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_nodes; ++i) {
        v_current(i) = std::max(0.0, std::min(v_new(i), v_n(i))); 
    }
}

double Fracture::calculate_error() const { return (v_current - v_prev_iter).cwiseAbs().maxCoeff(); }
void Fracture::save_previous_iteration() { v_prev_iter = v_current; }
void Fracture::save_previous_state() { v_backup = v_current; }
void Fracture::restore_previous_state() { v_current = v_backup; }
void Fracture::update_history() { v_n = v_current; }