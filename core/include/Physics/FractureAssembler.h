#pragma once
#include <Eigen/Sparse>
#include <vector>
#include "Core/Mesh.h"
#include "IO/Datafile.h"

class Math;
class Polarization;
class Mechanics;
class Electrostatics;

class FractureAssembler {
public:
    static void assemble_system(
        double dt, const Eigen::VectorXd& v_n,
        const Mesh& mesh, const Polarization& polarization, 
        const Mechanics& mechanics, const Electrostatics& electrostatics, 
        const Math& math, const Datafile& config,
        Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global
    );

    static void apply_irreversibility(
        const Mesh& mesh, const Eigen::VectorXd& v_prev_iter,
        Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global
    );
};