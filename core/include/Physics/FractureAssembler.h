#ifndef PHYSICS_FRACTURE_ASSEMBLER_H
#define PHYSICS_FRACTURE_ASSEMBLER_H

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>
#include "Physics/Polarization.h"
#include "Physics/Mechanics.h"
#include "Physics/Electrostatics.h"
#include "Physics/Math.h"
#include "IO/Datafile.h"
#include "Core/Mesh.h"

class FractureAssembler {
public:
    static void assemble_system(
        double dt, const Eigen::VectorXd& v_n, const Mesh& mesh, 
        const Polarization& polarization, const Mechanics& mechanics, 
        const Electrostatics& electrostatics, const Math& math, 
        const Datafile& config, Eigen::SparseMatrix<double>& K_global, 
        Eigen::VectorXd& F_global);
};

#endif // PHYSICS_FRACTURE_ASSEMBLER_H