#ifndef PHYSICS_POLARIZATION_ASSEMBLER_H
#define PHYSICS_POLARIZATION_ASSEMBLER_H

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>
#include "Physics/Fracture.h"
#include "Physics/Mechanics.h"
#include "Physics/Electrostatics.h"
#include "Physics/Math.h"
#include "IO/Datafile.h"
#include "Core/Mesh.h"
#include "Core/BoundaryManager.h"

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
};

#endif // PHYSICS_POLARIZATION_ASSEMBLER_H