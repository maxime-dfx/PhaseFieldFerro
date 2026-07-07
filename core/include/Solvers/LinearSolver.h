#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

class LinearSolver {
public:
    // Résout A*x = b avec un "guess" initial. Utilise le Gradient Conjugué, avec un repli sur un solveur direct si échec.
    static Eigen::VectorXd solve_with_guess(
        Eigen::SparseMatrix<double>& A, 
        const Eigen::VectorXd& b, 
        const Eigen::VectorXd& guess, 
        bool debug_enabled = false
    );
};