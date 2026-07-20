#ifndef PHYSICS_FRACTURE_ASSEMBLER_H
#define PHYSICS_FRACTURE_ASSEMBLER_H

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>
#include <array>
#include "Physics/Polarization.h"
#include "Physics/Mechanics.h"
#include "Physics/Electrostatics.h"
#include "Physics/MaterialModel.h"
#include "IO/Datafile.h"
#include "Mesh/Mesh.h"
#include "FEM/ElementIntegrator.h"

class FractureAssembler {
public:
    static void assemble_system(
        double dt, const Eigen::VectorXd& v_n, const Mesh& mesh, 
        const Polarization& polarization, const Mechanics& mechanics, 
        const Electrostatics& electrostatics, const MaterialModel& material, 
        const Datafile& config, Eigen::SparseMatrix<double>& K_global, 
        Eigen::VectorXd& F_global);

private:
    // --- Les Travailleurs ---
    static void calculer_matrices_elementaires(
        const Element& elem, const std::vector<std::array<double, 2>>& coords,
        double dt, const Eigen::VectorXd& v_n,
        const Polarization& polarization, const Mechanics& mechanics,
        const Electrostatics& electrostatics, const MaterialModel& material, const Datafile& config,
        Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local,
        Eigen::RowVectorXd& N_buffer, Eigen::MatrixXd& grad_N_buffer);

    static void distribuer_local_vers_global(
        const std::vector<int>& indices, int n_nodes,
        const Eigen::Ref<const Eigen::MatrixXd>& K_local, const Eigen::Ref<const Eigen::VectorXd>& F_local,
        std::vector<Eigen::Triplet<double>>& thread_triplets, Eigen::VectorXd& thread_F, Eigen::VectorXd& thread_diag);

    static void appliquer_irreversibilite(
        int n_dof, const Eigen::VectorXd& v_n, double max_diag, 
        std::vector<Eigen::Triplet<double>>& global_triplets, Eigen::VectorXd& F_global);
};

#endif // PHYSICS_FRACTURE_ASSEMBLER_H