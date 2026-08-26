#pragma once
//
// CGSolver.h
// ------------
// Implementations ISolver de la famille des solveurs iteratifs Eigen a
// preconditionneur parametrable : Conjugate Gradient (symetrique) et
// BiCGSTAB (asymetrique). Regroupes ici car ce sont les deux seules
// familles iteratives "pures Eigen" du projet (a l'exclusion de PETSc).

#include "Physics/include/Core/PhysicsConcepts.h"
#include "Utils/include/Profiling.h"
#include <Eigen/IterativeLinearSolvers>

template <typename PreconditionerType, typename ModuleTag>
class CGSolver : public ISolver {
private:
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, PreconditionerType> m_cg;
    bool m_failed = false;

public:
    // Le constructeur ne prend plus que la tolerance et le max_iter
    CGSolver(double tol, int max_iter) {
        m_cg.setTolerance(tol);
        m_cg.setMaxIterations(max_iter);
    }

    void analyze_pattern(Eigen::SparseMatrix<double>& mat) override {
        // Le nom et la couleur sont figes a la compilation
        PROFILE_ZONE_NC(ModuleTag::ZoneSolveCompute, ModuleTag::Color);
        m_cg.analyzePattern(mat);
        m_failed = (m_cg.info() != Eigen::Success);
    }

    void factorize(Eigen::SparseMatrix<double>& mat) override {
        // Le nom et la couleur sont figes a la compilation
        PROFILE_ZONE_NC(ModuleTag::ZoneSolveCompute, ModuleTag::Color);
        m_cg.factorize(mat);
        m_failed = (m_cg.info() != Eigen::Success);
    }

    Eigen::VectorXd solve(const Eigen::VectorXd& rhs, const Eigen::VectorXd& guess) override {
        // Le nom et la couleur sont figes a la compilation
        PROFILE_ZONE_NC(ModuleTag::ZoneSolveApply, ModuleTag::Color);
        Eigen::VectorXd res = m_cg.solveWithGuess(rhs, guess);
        m_failed = (m_cg.info() != Eigen::Success);
        return res;
    }

    bool has_failed() const override { return m_failed; }
    int iterations() const override { return m_cg.iterations(); }
    double error() const override { return m_cg.error(); }
};

template <typename PreconditionerType, typename ModuleTag>
class BiCGSTABSolver : public ISolver {
private:
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, PreconditionerType> m_bicgstab;
    bool m_failed = false;
public:
    BiCGSTABSolver(double tol, int max_iter) {
        m_bicgstab.setTolerance(tol);
        m_bicgstab.setMaxIterations(max_iter);
    }
    void analyze_pattern(Eigen::SparseMatrix<double>& mat) override {
        PROFILE_ZONE_NC(ModuleTag::ZoneSolveCompute, ModuleTag::Color);
        m_bicgstab.analyzePattern(mat);
        m_failed = (m_bicgstab.info() != Eigen::Success);
    }

    void factorize(Eigen::SparseMatrix<double>& mat) override {
        PROFILE_ZONE_NC(ModuleTag::ZoneSolveCompute, ModuleTag::Color);
        m_bicgstab.factorize(mat);
        m_failed = (m_bicgstab.info() != Eigen::Success);
    }
    Eigen::VectorXd solve(const Eigen::VectorXd& rhs, const Eigen::VectorXd& guess) override {
        PROFILE_ZONE_NC(ModuleTag::ZoneSolveApply, ModuleTag::Color);
        Eigen::VectorXd res = m_bicgstab.solveWithGuess(rhs, guess);
        m_failed = (m_bicgstab.info() != Eigen::Success);
        return res;
    }
    bool has_failed() const override { return m_failed; }
    int iterations() const override { return m_bicgstab.iterations(); }
    double error() const override { return m_bicgstab.error(); }
};
