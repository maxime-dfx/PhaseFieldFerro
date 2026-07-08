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

// ==============================================================================
// NOUVELLE METHODE MONOLITHIQUE
// ==============================================================================
// Resout Px et Py dans un unique systeme lineaire de taille 2*n_dof, avec
// des DOFs interleaves (2*i = Px_i, 2*i+1 = Py_i), exactement comme
// Mechanics::update_u le fait deja pour ux/uy. Remplace les deux appels
// sequentiels update_Px() + update_Py() par une seule assemblage + une
// seule resolution -> gain de performance (un seul solve au lieu de deux),
// et coherence avec le reste du code (Mechanics).
//
// IMPORTANT : le couplage physique Px<->Py (terme b3, termes croises
// Landau) reste traite de facon semi-implicite, exactement comme avant
// (via les valeurs courantes Pi_gp dans PolarizationAssembler) - il n'y a
// PAS de blocs hors-diagonale Px-Py dans la matrice globale. Seule la
// RESOLUTION devient monolithique, pas le schema numerique.
// ==============================================================================
void Polarization::update_polarization(double time, const Fracture& fracture, const Mechanics& mechanics,
                                        const Electrostatics& electrostatics, const Math& math)
{
    double dt = config.simulation.dt;
    const int n_nodes = static_cast<int>(Px_current.size());
    const int system_size = 2 * n_nodes;

    Eigen::SparseMatrix<double> K_global(system_size, system_size);
    Eigen::VectorXd F_global = Eigen::VectorXd::Zero(system_size);

    // 1. Assemblage du systeme complet (Px et Py, DOFs interleaves)
    PolarizationAssembler::assemble_system_monolithic(
        dt, Px_current, Py_current, Px_n, Py_n,
        mesh, fracture, mechanics, electrostatics, math, config, bc_manager,
        K_global, F_global
    );

    // 2. Application des conditions de Dirichlet (Px et Py)
    PolarizationAssembler::apply_boundary_conditions_monolithic(
        K_global, F_global, bc_manager.get_px_bcs(), bc_manager.get_py_bcs()
    );

    // 3. Preparation du guess interleave (demarrage a chaud)
    Eigen::VectorXd P_guess(system_size);
    for (int i = 0; i < n_nodes; ++i) {
        P_guess(2 * i)     = Px_current(i);
        P_guess(2 * i + 1) = Py_current(i);
    }

    // 4. Resolution
    Eigen::VectorXd P_new = LinearSolver::solve_with_guess(K_global, F_global, P_guess, config.simulation.debug_enabled);

    // 5. Mise a jour de l'etat
    if (P_new.allFinite()) {
        map_global_vector_to_components(P_new);
    } else {
        Logger::error("[Polarization] Echec resolution monolithique ou solution non finie");
    }
}

void Polarization::map_global_vector_to_components(const Eigen::VectorXd& P_new) {
    const int n_nodes = static_cast<int>(Px_current.size());
    for (int i = 0; i < n_nodes; ++i) {
        Px_current(i) = P_new(2 * i);
        Py_current(i) = P_new(2 * i + 1);
    }
}

// ==============================================================================
// ANCIENNE METHODE PAR COMPOSANTE - conservee pour compatibilite ascendante
// (utilisee par update_Px()/update_Py() dans le header). A supprimer une
// fois la migration vers update_polarization() validee partout.
// ==============================================================================
void Polarization::update_polarization_component(double time, const Fracture& fracture, const Mechanics& mechanics,
                                                   const Electrostatics& electrostatics, const Math& math, int component)
{
    double dt = config.simulation.dt;
    const int n_dof = static_cast<int>(Px_current.size());

    Eigen::SparseMatrix<double> A_sparse(n_dof, n_dof);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n_dof);

    PolarizationAssembler::assemble_system(
        dt, component, Px_current, Py_current, Px_n, Py_n,
        mesh, fracture, mechanics, electrostatics, math, config, bc_manager,
        A_sparse, b
    );

    const auto& bcs = (component == 0) ? bc_manager.get_px_bcs() : bc_manager.get_py_bcs();
    PolarizationAssembler::apply_boundary_conditions(A_sparse, b, bcs);

    Eigen::VectorXd P_guess = (component == 0) ? Px_current : Py_current;
    Eigen::VectorXd solution = LinearSolver::solve_with_guess(A_sparse, b, P_guess, config.simulation.debug_enabled);

    if (solution.allFinite()) {
        if (component == 0) Px_current = solution;
        else Py_current = solution;
    } else {
        Logger::error("[Polarization] Echec resolution ou solution non finie (comp=" + std::to_string(component) + ")");
    }
}