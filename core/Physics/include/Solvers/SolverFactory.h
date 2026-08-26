#pragma once
//
// SolverFactory.h
// -----------------
// Fusionne l'ancien SolverConfig.h (petite struct de configuration) avec la
// Factory elle-meme. ISolver reste defini dans Core/PhysicsConcepts.h (voir
// la note de layering dans ce fichier) ; ce header n'a besoin que de le
// consommer.

#include "Physics/include/Core/PhysicsConcepts.h"
#include "Physics/include/Solvers/internal/CGSolver.h"
#include "Physics/include/Solvers/internal/DirectSolver.h"
#include "Physics/include/Solvers/internal/PetscSolver.h"

#include "Utils/include/Logger.h"

#include <Eigen/IterativeLinearSolvers>

#include <memory>
#include <string>
#include <algorithm>

struct SolverConfig {
    std::string solver_type = "CG";
    std::string preconditioner_type = "DIAGONAL";
    double tolerance = 1e-8;
    int max_iterations = 1000;
};

class SolverFactory {
public:
    template <typename ModuleTag>
    static std::unique_ptr<ISolver> create(const SolverConfig& config) {
        std::string type = config.solver_type;
        std::string prec = config.preconditioner_type;
        std::transform(type.begin(), type.end(), type.begin(), ::toupper);
        std::transform(prec.begin(), prec.end(), prec.begin(), ::toupper);

        if (type == "DIRECT") {
            return std::make_unique<DirectSolver<ModuleTag>>();
        }

#ifdef USE_PETSC
        if (type == "PETSC" || type == "GMRES") {
            return std::make_unique<PetscSolver<ModuleTag>>(config.tolerance, config.max_iterations, prec);
        }
#endif

        // --- Famille Conjugate Gradient (Symetrique) ---
        if (type == "CG") {
            if (prec == "INCOMPLETE_CHOLESKY") return std::make_unique<CGSolver<Eigen::IncompleteCholesky<double>, ModuleTag>>(config.tolerance, config.max_iterations);
            if (prec == "ILUT") return std::make_unique<CGSolver<Eigen::IncompleteLUT<double>, ModuleTag>>(config.tolerance, config.max_iterations);
            return std::make_unique<CGSolver<Eigen::DiagonalPreconditioner<double>, ModuleTag>>(config.tolerance, config.max_iterations);
        }

        // --- Famille BiCGSTAB (Asymetrique) ---
        if (type == "BICGSTAB") {
            if (prec == "INCOMPLETE_CHOLESKY") return std::make_unique<BiCGSTABSolver<Eigen::IncompleteCholesky<double>, ModuleTag>>(config.tolerance, config.max_iterations);
            if (prec == "ILUT") return std::make_unique<BiCGSTABSolver<Eigen::IncompleteLUT<double>, ModuleTag>>(config.tolerance, config.max_iterations);
            return std::make_unique<BiCGSTABSolver<Eigen::DiagonalPreconditioner<double>, ModuleTag>>(config.tolerance, config.max_iterations);
        }

        Logger::warning("Solveur inconnu (", type, "), fallback sur CG + Diagonal");
        return std::make_unique<CGSolver<Eigen::DiagonalPreconditioner<double>, ModuleTag>>(config.tolerance, config.max_iterations);
    }
};
