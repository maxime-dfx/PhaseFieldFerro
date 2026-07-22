#include "Solvers/LinearSolver.h"
#include <Eigen/CholmodSupport>
#include <iostream>
// #include <tracy/Tracy.hpp> // Si tu utilises Tracy ici

Eigen::VectorXd LinearSolver::solve_with_guess(
    Eigen::SparseMatrix<double>& A, 
    const Eigen::VectorXd& b, 
    const Eigen::VectorXd& guess, 
    bool debug_enabled) 
{
    // ZoneScoped; // Pour Tracy

    // Le solveur et l'état d'analyse sont déclarés "static". 
    // Ils persistent en mémoire entre chaque appel de la fonction.
    static Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    static bool is_analyzed = false;

    if (debug_enabled) {
        std::cout << "[DEBUG] Resolution du systeme lineaire..." << std::endl;
    }

    // 1. Analyse symbolique : exécutee UNE SEULE FOIS lors du tout premier appel
    if (!is_analyzed) {
        solver.analyzePattern(A);
        if (solver.info() != Eigen::Success) {
            std::cerr << "[ERREUR] Echec de l'analyse symbolique CHOLMOD !" << std::endl;
        }
        is_analyzed = true;
    }

    // 2. Factorisation numérique : exécutée à chaque itération (très rapide)
    solver.factorize(A);
    if (solver.info() != Eigen::Success) {
        std::cerr << "[ERREUR] Echec de la factorisation CHOLMOD !" << std::endl;
    }

    // 3. Résolution du système (A * x = b)
    Eigen::VectorXd x = solver.solve(b);

    return x;
}