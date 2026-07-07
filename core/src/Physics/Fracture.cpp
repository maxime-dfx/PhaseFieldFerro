#include "Physics/Fracture.h"
#include "Physics/FractureAssembler.h"
#include "Solvers/LinearSolver.h"
#include "Utils/Logger.h"

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

    // 1. Assemblage
    FractureAssembler::assemble_system(dt, v_n, mesh, polarization, mechanics, electrostatics, math, config, K_global, F_global);

    // 2. Irréversibilité
    FractureAssembler::apply_irreversibility(mesh, v_prev_iter, K_global, F_global);

    // 3. Résolution (avec v_current comme guess initial)
    Eigen::VectorXd v_new = LinearSolver::solve_with_guess(K_global, F_global, v_current, config.simulation.debug_enabled);

    // 4. Post-traitement et bornage mathématique
    if (v_new.allFinite()) {
        enforce_physical_bounds(v_new);
    } else {
        Logger::error("[Fracture] Echec resolution ou solution non finie");
    }
}

void Fracture::enforce_physical_bounds(const Eigen::VectorXd& v_new) {
    for (size_t i = 0; i < mesh.get_num_nodes(); ++i) {
        v_current(i) = std::max(0.0, std::min(v_new(i), v_prev_iter(i))); 
    }
}