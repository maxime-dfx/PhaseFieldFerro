#ifndef PHYSICS_MECHANICS_ASSEMBLER_H
#define PHYSICS_MECHANICS_ASSEMBLER_H

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>
#include "Physics/Polarization.h"
#include "Physics/Fracture.h"
#include "Physics/Math.h"
#include "Core/Mesh.h"
#include "Core/BoundaryManager.h"

class MechanicsAssembler {
public:
    static void assemble_system(
        const Mesh& mesh, const Polarization& polarization, const Fracture& fracture, 
        const Math& math, const std::vector<NodeBC>& bcs_x, const std::vector<NodeBC>& bcs_y,
        Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global);
};

#endif // PHYSICS_MECHANICS_ASSEMBLER_H