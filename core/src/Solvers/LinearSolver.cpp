#include "Solvers/LinearSolver.h"
#include "Utils/Logger.h"

Eigen::VectorXd LinearSolver::solve_with_guess(Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, const Eigen::VectorXd& guess, bool debug_enabled) {
    // 1. Configuration du Solveur Itératif (Gradient Conjugué + Incomplete Cholesky)
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> solver;
    solver.setTolerance(1e-10); 
    solver.setMaxIterations(300);
    solver.compute(A);

    // 2. Résolution avec Démarrage à chaud
    Eigen::VectorXd x = solver.solveWithGuess(b, guess);

    // 3. Fallback de sécurité (Si la matrice est trop mal conditionnée)
    if (solver.info() != Eigen::Success) {
        Logger::debug("[LinearSolver] CG a échoué (iters: " + std::to_string(solver.iterations()) + "). Bascule sur LDLT Direct.", debug_enabled);
        
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> direct_solver;
        direct_solver.compute(A);
        x = direct_solver.solve(b);
    } else {
        Logger::debug("[LinearSolver] CG converge en " + std::to_string(solver.iterations()) + " itérations.", debug_enabled);
    }
    
    return x;
}