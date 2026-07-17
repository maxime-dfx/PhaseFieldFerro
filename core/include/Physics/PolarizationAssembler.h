#ifndef PHYSICS_POLARIZATION_ASSEMBLER_H
#define PHYSICS_POLARIZATION_ASSEMBLER_H

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>
#include <array>
#include "Physics/Fracture.h"
#include "Physics/Mechanics.h"
#include "Physics/Electrostatics.h"
#include "Physics/Math.h"
#include "IO/Datafile.h"
#include "Core/Mesh.h"
#include "Core/BoundaryManager.h"
#include "Core/ElementIntegrator.h"

class PolarizationAssembler {
public:
    static void assemble_system(
        double dt,
        const Eigen::VectorXd& Px_current, const Eigen::VectorXd& Py_current,
        const Eigen::VectorXd& Px_n, const Eigen::VectorXd& Py_n,
        const Mesh& mesh, const Fracture& fracture, const Mechanics& mechanics,
        const Electrostatics& electrostatics, const Math& math, const Datafile& config,
        const BoundaryManager& bc_manager,
        Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global);

private:
    // --- Les Travailleurs ---
    static void calculer_matrices_elementaires(
        const Element& elem, const std::vector<std::array<double, 2>>& coords,
        double dt, const Eigen::VectorXd& Px_current, const Eigen::VectorXd& Py_current,
        const Eigen::VectorXd& Px_n, const Eigen::VectorXd& Py_n,
        const Fracture& fracture, const Mechanics& mechanics,
        const Electrostatics& electrostatics, const Math& math, const Datafile& config,
        Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local);

    static void distribuer_local_vers_global(
        const std::vector<int>& indices, int n_nodes,
        const Eigen::Ref<const Eigen::MatrixXd>& K_local, const Eigen::Ref<const Eigen::VectorXd>& F_local,
        std::vector<Eigen::Triplet<double>>& thread_triplets, 
        Eigen::VectorXd& thread_F, Eigen::VectorXd& thread_diag_px, Eigen::VectorXd& thread_diag_py);

    static void traiter_dof_flottants(
        int n_dof, const std::vector<NodeBC>& bcs_px, const std::vector<NodeBC>& bcs_py,
        const Eigen::VectorXd& diag_global_px, const Eigen::VectorXd& diag_global_py,
        std::vector<Eigen::Triplet<double>>& global_triplets, Eigen::VectorXd& F_global, bool debug_enabled);

    static void appliquer_conditions_limites(
        const std::vector<NodeBC>& bcs_px, const std::vector<NodeBC>& bcs_py,
        double max_diag, std::vector<Eigen::Triplet<double>>& global_triplets, Eigen::VectorXd& F_global);
};

#endif // PHYSICS_POLARIZATION_ASSEMBLER_H