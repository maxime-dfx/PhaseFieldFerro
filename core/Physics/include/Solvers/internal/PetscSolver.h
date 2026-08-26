#pragma once
//
// PetscSolver.h
// ---------------
// Solveur ISolver adossé a PETSc/Hypre (via PetscKspSolverWrapper), actif
// uniquement si USE_PETSC est defini au build.

#include "Physics/include/Core/PhysicsConcepts.h"

#ifdef USE_PETSC
#include "Utils/include/PetscKspWrapper.h"
#include "Utils/include/Profiling.h" // NECESSAIRE POUR LES MACROS TRACY
#include <string>
#include <utility>

// On ajoute le ModuleTag en parametre de template !
template <typename ModuleTag>
class PetscSolver : public ISolver {
private:
    PetscKspSolverWrapper m_wrapper;
    double m_tol;
    int m_max_iter;
    std::string m_prec;
public:
    // Le constructeur ne prend plus le label en argument libre, il l'extrait du tag
    PetscSolver(double tol, int max_iter, std::string prec)
        : m_wrapper(ModuleTag::ModuleName), m_tol(tol), m_max_iter(max_iter), m_prec(std::move(prec)) {}

    void analyze_pattern(Eigen::SparseMatrix<double>&) override {
        // Le wrapper PETSc ne separe pas l'analyse du pattern de la
        // configuration complete : la KSP est (re)configuree entierement
        // dans factorize(). Rien a faire ici, mais on garde la methode
        // pour respecter l'interface ISolver commune.
    }

    void factorize(Eigen::SparseMatrix<double>& mat) override {
        // Injection du nom et de la couleur a la compilation !
        PROFILE_ZONE_NC(ModuleTag::ZoneSolveCompute, ModuleTag::Color);
        m_wrapper.configure(mat, m_tol, m_max_iter, m_prec);
    }

    Eigen::VectorXd solve(const Eigen::VectorXd& rhs, const Eigen::VectorXd& guess) override {
        // Injection du nom et de la couleur a la compilation !
        PROFILE_ZONE_NC(ModuleTag::ZoneSolveApply, ModuleTag::Color);
        return m_wrapper.solve(rhs, guess);
    }

    bool has_failed() const override {
        return m_wrapper.reason() < 0;
    }

    int iterations() const override { return m_wrapper.iterations(); }
    double error() const override { return m_wrapper.residual(); }
};
#endif
