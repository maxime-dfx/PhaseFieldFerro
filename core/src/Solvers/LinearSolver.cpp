#include "Solvers/LinearSolver.h"
#include <Eigen/SparseLU>
#include <iostream>

Eigen::VectorXd LinearSolver::solve_with_guess(
    Eigen::SparseMatrix<double>& A, 
    const Eigen::VectorXd& b, 
    const Eigen::VectorXd& guess, 
    bool debug_enabled
) {
    // Compression de la matrice (obligatoire pour SparseLU)
    A.makeCompressed();

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    
    // Analyse de la structure (sparsity pattern) et factorisation numérique
    solver.analyzePattern(A);
    solver.factorize(A);

    if (solver.info() != Eigen::Success) {
        std::cerr << "[ERREUR] La factorisation SparseLU a échoué (matrice singulière ?)\n";
        return guess; // Repli de sécurité
    }

    // Résolution directe (le "guess" initial n'est pas utilisé par les solveurs directs,
    // mais il est conservé dans la signature pour compatibilité avec votre architecture)
    Eigen::VectorXd x = solver.solve(b);

    if (solver.info() != Eigen::Success) {
        std::cerr << "[ERREUR] La résolution SparseLU a échoué.\n";
        return guess;
    }

    return x;
}