#pragma once
#include <Eigen/Sparse>
#include <vector>
#include "Core/Mesh.h"
#include "Core/BoundaryManager.h"
#include "IO/Datafile.h"

class Math;
class Polarization;
class Fracture;

class ElectrostaticsAssembler {
public:
    // Assemble la matrice globale K et le vecteur force F pour l'électrostatique
    static void assemble_system(
        const Mesh& mesh, 
        const Polarization& polarization, 
        const Fracture& fracture, 
        const Math& math,
        const Datafile& config,
        Eigen::SparseMatrix<double>& K_global, 
        Eigen::VectorXd& F_global
    );

    // Applique les conditions aux limites (Méthode de pénalité)
    static void apply_boundary_conditions(
        Eigen::SparseMatrix<double>& K, 
        Eigen::VectorXd& F, 
        const std::vector<NodeBC>& bcs
    );
};