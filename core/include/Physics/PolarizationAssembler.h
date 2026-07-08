#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include "Core/Mesh.h"
#include "IO/Datafile.h"
#include "Core/BoundaryManager.h"

class Fracture;
class Mechanics;
class Electrostatics;
class Math;

class PolarizationAssembler {
public:
    // ==========================================================================
    // NOUVELLES FONCTIONS MONOLITHIQUES
    // ==========================================================================
    // Assemble le systeme complet 2*n_dof x 2*n_dof (Px et Py interleaves,
    // DOF 2*i = Px_i, DOF 2*i+1 = Py_i), analogue a MechanicsAssembler::assemble_system.
    static void assemble_system_monolithic(
        double dt,
        const Eigen::VectorXd& Px_current, const Eigen::VectorXd& Py_current,
        const Eigen::VectorXd& Px_n, const Eigen::VectorXd& Py_n,
        const Mesh& mesh, const Fracture& fracture, const Mechanics& mechanics,
        const Electrostatics& electrostatics, const Math& math, const Datafile& config,
        const BoundaryManager& bc_manager,
        Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global
    );

    static void apply_boundary_conditions_monolithic(
        Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F,
        const std::vector<NodeBC>& bcs_px, const std::vector<NodeBC>& bcs_py
    );

    // ==========================================================================
    // ANCIENNES FONCTIONS PAR COMPOSANTE - conservees pour compatibilite
    // ascendante (utilisees par Polarization::update_Px/update_Py). A
    // supprimer une fois la migration validee.
    // ==========================================================================
    static void assemble_system(
        double dt, int component,
        const Eigen::VectorXd& Px_current, const Eigen::VectorXd& Py_current,
        const Eigen::VectorXd& Px_n, const Eigen::VectorXd& Py_n,
        const Mesh& mesh, const Fracture& fracture, const Mechanics& mechanics,
        const Electrostatics& electrostatics, const Math& math, const Datafile& config,
        const BoundaryManager& bc_manager,
        Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global
    );

    static void apply_boundary_conditions(
        Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, const std::vector<NodeBC>& bcs
    );

    static void handle_floating_dofs(
        std::vector<Eigen::Triplet<double>>& triplets, Eigen::VectorXd& b,
        const Eigen::VectorXd& diag_check, const std::vector<NodeBC>& bcs, bool debug_enabled
    );

private:
    // Version monolithique du traitement des DOFs flottants (verifie les
    // deux composantes Px et Py separement, comme deux blocs diagonaux
    // independants dans le systeme interleave).
    static void handle_floating_dofs_monolithic(
        std::vector<Eigen::Triplet<double>>& triplets, Eigen::VectorXd& b,
        const Eigen::VectorXd& diag_check_px, const Eigen::VectorXd& diag_check_py,
        const std::vector<NodeBC>& bcs_px, const std::vector<NodeBC>& bcs_py,
        bool debug_enabled
    );
};