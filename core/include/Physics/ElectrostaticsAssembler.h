#ifndef PHYSICS_ELECTROSTATICS_ASSEMBLER_H
#define PHYSICS_ELECTROSTATICS_ASSEMBLER_H

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>
#include <array>
#include "Physics/Polarization.h"
#include "Physics/Fracture.h"
#include "Physics/MaterialModel.h"
#include "IO/Datafile.h"
#include "Mesh/Mesh.h"
#include "BC/BoundaryManager.h"
#include "FEM/ElementIntegrator.h"

class ElectrostaticsAssembler {
public:
    static void assemble_system(
        const Mesh& mesh, const Polarization& polarization, const Fracture& fracture, 
        const MaterialModel& material, const Datafile& config, const std::vector<NodeBC>& bcs,
        Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global);

private:
    // Nos fonctions spécialisées (Travailleurs)
    static void calculer_matrices_elementaires(
        const Element& elem, const std::vector<std::array<double, 2>>& coords,
        const Polarization& polarization, const Fracture& fracture, 
        const MaterialModel& material, const Datafile& config,
        Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local,
        Eigen::RowVectorXd& N_buffer, Eigen::MatrixXd& grad_N_buffer);

    static void distribuer_local_vers_global(
        const std::vector<int>& indices, int n_nodes,
        const Eigen::Ref<const Eigen::MatrixXd>& K_local, const Eigen::Ref<const Eigen::VectorXd>& F_local,
        std::vector<Eigen::Triplet<double>>& thread_triplets, 
        Eigen::VectorXd& thread_F, Eigen::VectorXd& thread_diag);
        
    static void appliquer_conditions_limites(
        const std::vector<NodeBC>& bcs, double max_diag, 
        std::vector<Eigen::Triplet<double>>& global_triplets, Eigen::VectorXd& F_global);
};

#endif // PHYSICS_ELECTROSTATICS_ASSEMBLER_H