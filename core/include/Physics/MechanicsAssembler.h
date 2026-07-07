#pragma once
#include <Eigen/Sparse>
#include <vector>
#include "Core/Mesh.h"
#include "Core/BoundaryManager.h"
#include "Physics/Math.h"
#include "Physics/Polarization.h"
#include "Physics/Fracture.h"

class MechanicsAssembler {
public:
    // Assemble la matrice globale K et le vecteur force F
    static void assemble_system(
        const Mesh& mesh, 
        const Polarization& polarization, 
        const Fracture& fracture, 
        const Math& math,
        Eigen::SparseMatrix<double>& K_global,
        Eigen::VectorXd& F_global
    );

    // Applique les conditions aux limites (Méthode de pénalité)
    static void apply_boundary_conditions(
        Eigen::SparseMatrix<double>& K, 
        Eigen::VectorXd& F, 
        const std::vector<NodeBC>& bcs_x, 
        const std::vector<NodeBC>& bcs_y
    );
};