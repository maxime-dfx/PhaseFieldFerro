#ifndef UTILS_SOLVER_CONFIG_H
#define UTILS_SOLVER_CONFIG_H

// ============================================================================
// Bascule de compilation pour les solveurs des 4 modules physiques
// (Mechanics, Polarization, Fracture, Electrostatics).
//
// Les deux typedefs sont DECOUPLES (cf. bench Tracy jobs 345 vs 346, w102,
// 48 coeurs, même config.toml) :
//
//   DirectSolverType (Mechanics, système non-symétrique, CG volontairement
//   évité — cf. commentaire historique dans Mechanics.h sur la convergence
//   prématurée du CG face à la pénalité des CL) :
//       USE_PARDISO défini     -> Eigen::PardisoLDLT (MKL, multi-thread,
//                                  ~4x plus rapide que SimplicialLDLT sur le
//                                  bench : 443,95 ms -> 108,72 ms).
//       USE_PARDISO non défini -> Eigen::SimplicialLDLT (mono-thread, sans
//                                  dépendance MKL).
//
//   SymmetricSolverType (Polarization/Electrostatics/Fracture, systèmes
//   symétriques bien conditionnés) :
//       TOUJOURS Eigen::ConjugateGradient + préconditionneur diagonal,
//       QUEL QUE SOIT USE_PARDISO. Le CG+diag y est déjà quasi-instantané ;
//       le bench a montré qu'y substituer un direct PARDISO est ~100x plus
//       lent (Solve_Compute Polarization : 817,75 µs -> 83,98 ms), car la
//       factorisation directe complète est inutile sur des systèmes que le
//       CG résout en quelques itérations.
//
// Piloté depuis CMakeLists.txt via `option(USE_PARDISO ...)`. Ne jamais
// définir cette macro à la main dans le code — uniquement via le flag de
// build, pour que la bascule reste visible et traçable au niveau CMake.
// ============================================================================

#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

#ifdef USE_PARDISO
    #include <Eigen/PardisoSupport>

    using DirectSolverType = Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>;

    #define SOLVER_IS_DIRECT 1
#else
    using DirectSolverType = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;

    #define SOLVER_IS_DIRECT 0
#endif

// SymmetricSolverType : indépendant de USE_PARDISO, toujours CG.
//
// Préconditionneur : Diagonal (Jacobi), PAS IncompleteCholesky.
//
// IncompleteCholesky a été testé (bench job 393, w102) : il réduit bien le
// nombre d'itérations CG d'Electrostatics (~600 -> ~120-220 sur job 390 vs
// 393), mais son application par itération est bien plus coûteuse — la
// résolution triangulaire incomplète est intrinsèquement séquentielle et
// mal parallélisée par Eigen (contrairement au simple scaling diagonal du
// Jacobi, trivialement vectorisable). Résultat net : Solve_Apply
// Electrostatics est passé de 69 ms (Diagonal) à 202 ms (IncompleteCholesky)
// malgré 3x moins d'itérations. Le gain en itérations ne compense pas le
// surcoût par itération sur ce matériel (bi-socket Xeon E5-2680 v3, OpenMP).
// Diagonal reste donc le meilleur compromis mesuré à ce jour.
using SymmetricSolverType = Eigen::ConjugateGradient<
    Eigen::SparseMatrix<double>,
    Eigen::Lower | Eigen::Upper,
    Eigen::DiagonalPreconditioner<double>>;

#endif // UTILS_SOLVER_CONFIG_H