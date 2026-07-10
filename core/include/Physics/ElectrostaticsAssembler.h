#ifndef PHYSICS_ELECTROSTATICS_ASSEMBLER_H
#define PHYSICS_ELECTROSTATICS_ASSEMBLER_H

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>
#include "Physics/Polarization.h"
#include "Physics/Fracture.h"
#include "Physics/Math.h"
#include "IO/Datafile.h"
#include "Core/Mesh.h"
#include "Core/BoundaryManager.h"

class ElectrostaticsAssembler {
public:
    static void assemble_system(
        const Mesh& mesh, const Polarization& polarization, const Fracture& fracture, 
        const Math& math, const Datafile& config, const std::vector<NodeBC>& bcs,
        Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global);
};

#endif // PHYSICS_ELECTROSTATICS_ASSEMBLER_H