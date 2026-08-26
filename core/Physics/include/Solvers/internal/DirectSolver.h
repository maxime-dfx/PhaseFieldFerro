#pragma once
//
// DirectSolver.h
// ----------------
// Solveur direct : PARDISO (MKL) si USE_PARDISO est defini, sinon
// SimplicialLDLT d'Eigen.

#include "Physics/include/Core/PhysicsConcepts.h"
#include "Utils/include/Logger.h"
#include "Utils/include/Profiling.h"

#ifdef USE_PARDISO
#include <Eigen/PardisoSupport>
#include <mkl.h>
#endif

template <typename ModuleTag>
class DirectSolver : public ISolver {
private:
#ifdef USE_PARDISO
    Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> m_direct;
#else
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> m_direct;
#endif
    bool m_failed = false;

public:
    DirectSolver() {
#ifdef USE_PARDISO
        Logger::info("Using Pardiso direct solver");

        // iparm[1] : fill-reducing ordering. 2 = METIS.
        m_direct.pardisoParameterArray()[1] = 2;

        // iparm[23] : two-level factorization. Prerequis pour que le parallel
        // solve (iparm[24]) ait un effet reel — sans ca MKL ignore iparm[24]
        // silencieusement.
        m_direct.pardisoParameterArray()[23] = 1;

        // iparm[24] : active le solve forward/backward parallele.
        m_direct.pardisoParameterArray()[24] = 1;

        // iparm[26] (msglvl) : affiche les stats natives PARDISO (threads,
        // taille systeme, nnz, temps par phase) sur stdout a chaque appel.
        // A retirer une fois le diagnostic termine, ca spam la sortie.
        m_direct.pardisoParameterArray()[26] = 1;
#else
        Logger::info("Using Eigen direct solver");
#endif
    }

    void analyze_pattern(Eigen::SparseMatrix<double>& mat) override {
        PROFILE_ZONE_NC(ModuleTag::ZoneSolveCompute, ModuleTag::Color);
        m_direct.analyzePattern(mat);
        m_failed = (m_direct.info() != Eigen::Success);
    }

    void factorize(Eigen::SparseMatrix<double>& mat) override {
        PROFILE_ZONE_NC(ModuleTag::ZoneSolveCompute, ModuleTag::Color);
        m_direct.factorize(mat);
        m_failed = (m_direct.info() != Eigen::Success);
    }

    Eigen::VectorXd solve(const Eigen::VectorXd& rhs, const Eigen::VectorXd&) override {
        PROFILE_ZONE_NC(ModuleTag::ZoneSolveApply, ModuleTag::Color);
        Eigen::VectorXd res = m_direct.solve(rhs);
        m_failed = (m_direct.info() != Eigen::Success);
        return res;
    }

    bool has_failed() const override { return m_failed; }
    int iterations() const override { return 1; }
    double error() const override { return 0.0; }
};
