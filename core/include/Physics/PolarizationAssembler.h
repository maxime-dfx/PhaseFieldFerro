#pragma once
#include <Eigen/Sparse>
#include <vector>
#include "Core/Mesh.h"
#include "Core/BoundaryManager.h"
#include "IO/Datafile.h"

class Math;
class Fracture;
class Mechanics;
class Electrostatics;

class PolarizationAssembler {
public:
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
        Eigen::SparseMatrix<double>& K, 
        Eigen::VectorXd& F, 
        const std::vector<NodeBC>& bcs
    );

private:
    static void handle_floating_dofs(
        std::vector<Eigen::Triplet<double>>& triplets, 
        Eigen::VectorXd& b, 
        const Eigen::VectorXd& diag_check, 
        const std::vector<NodeBC>& bcs,
        bool debug_enabled
    );
};