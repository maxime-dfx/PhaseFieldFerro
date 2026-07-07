#include "Physics/Polarization.h"
#include "Physics/PolarizationAssembler.h"
#include "Solvers/LinearSolver.h"
#include "Utils/Logger.h"

Polarization::Polarization(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager) 
    : config(config), mesh(mesh), bc_manager(bc_manager) 
{
    size_t n_nodes = mesh.get_num_nodes();
    Px_prev_iter.setZero(n_nodes);
    Py_prev_iter.setZero(n_nodes);
    Px_current.setZero(n_nodes);
    Py_current.setZero(n_nodes);
    Px_n.setZero(n_nodes);
    Py_n.setZero(n_nodes);

    if (config.polarization.type == InitializationType::UNIFORM) {
        Px_current.setConstant(config.polarization.val_x_0);
        Py_current.setConstant(config.polarization.val_y_0);
    } else if (config.polarization.type == InitializationType::RANDOM) {
        for (size_t i = 0; i < n_nodes; ++i) {
            Px_current[i] = static_cast<double>(rand()) / RAND_MAX;
            Py_current[i] = static_cast<double>(rand()) / RAND_MAX;
        }
    }

    Px_n = Px_current;
    Py_n = Py_current;
}

void Polarization::update_polarization_component(double time, const Fracture& fracture, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math, int component) {
    double dt = config.simulation.dt;
    const int n_dof = static_cast<int>(Px_current.size());

    Eigen::SparseMatrix<double> A_sparse(n_dof, n_dof);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n_dof);

    // 1. Assemblage du système (inclut le traitement des DOF flottants)
    PolarizationAssembler::assemble_system(
        dt, component, Px_current, Py_current, Px_n, Py_n, 
        mesh, fracture, mechanics, electrostatics, math, config, bc_manager, 
        A_sparse, b
    );

    // 2. Application des conditions de Dirichlet
    const auto& bcs = (component == 0) ? bc_manager.get_px_bcs() : bc_manager.get_py_bcs();
    PolarizationAssembler::apply_boundary_conditions(A_sparse, b, bcs);

    // 3. Résolution itérative avec Warm Start (via LinearSolver)
    Eigen::VectorXd P_guess = (component == 0) ? Px_current : Py_current;
    Eigen::VectorXd solution = LinearSolver::solve_with_guess(A_sparse, b, P_guess, config.simulation.debug_enabled);

    // 4. Mise à jour de l'état
    if (solution.allFinite()) {
        if (component == 0) Px_current = solution;
        else Py_current = solution;
    } else {
        Logger::error("[Polarization] Echec resolution ou solution non finie (comp=" + std::to_string(component) + ")");
    }
}