#pragma once
//
// SystemAssembler.h
// ------------------
// Couche d'assemblage CSR/Triplets. Fusionne l'ancien FlattenedBC.h ici
// (petite struct de donnees, ne justifiait pas son propre fichier).
//
// NOTE DE LAYERING : ce header est volontairement autonome vis-a-vis de
// Core/PhysicsConcepts.h (qui, lui, inclut SystemAssembler.h pour definir
// IElementEquation/IPhysicsDofMapper via ElementContext/FlattenedBC). Le
// type PhysicsState n'est donc que forward-declare ici : SystemAssembler
// ne fait que le transmettre par reference a EquationType::compute_element_matrices,
// jamais y accede directement.

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <array>

#include "Mesh/include/Mesh.h"
#include "Utils/include/Profiling.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(_MSC_VER)
    #define FORCE_INLINE inline __forceinline
#elif defined(__GNUC__) || defined(__clang__)
    #define FORCE_INLINE inline __attribute__((always_inline))
#else
    #define FORCE_INLINE inline
#endif

struct PhysicsState; // defini dans Core/PhysicsConcepts.h

// --- Conditions aux limites "aplaties", pretes a etre injectees dans le
//     systeme lineaire global par penalisation (voir apply_bcs ci-dessous).
struct FlattenedBC {
    int global_dof;
    double value;
    bool is_dirichlet;
};

// --- Contexte passe a chaque IElementEquation::compute_element_matrices.
struct ElementContext {
    const Element& elem;
    const std::array<std::array<double, 2>, 8>& coords;
    int elem_idx;
    int material_id;
    const Mesh& mesh;
};

class SystemAssembler {
public:
    // Transformation en template pour inlining parfait de l'equation
    template <typename EquationType>
    FORCE_INLINE static void assemble(
        int num_elements, int num_dofs_total, bool is_first_assembly,
        Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global,
        std::vector<long>& csr_mapping,
        const EquationType& equation, const std::vector<FlattenedBC>& active_bcs,
        const PhysicsState& state,
        const Mesh& mesh) // <-- Plus de module_name ici
    {
        if (is_first_assembly) {
            assemble_initial_triplets(num_elements, num_dofs_total, K_global, F_global,
                                       csr_mapping, equation, active_bcs, state, mesh);
        } else {
            assemble_fast_csr_bypass(num_elements, num_dofs_total, K_global, F_global,
                                      csr_mapping, equation, active_bcs, state, mesh);
        }
    }

private:
    static constexpr int MAX_LOCAL_DOFS = 16;

    // --- Utilitaires CSR ---
    static long find_csr_index(int row, int col, const int* outer_ptr, const int* inner_ptr) {
        const int* begin = inner_ptr + outer_ptr[col];
        const int* end = inner_ptr + outer_ptr[col + 1];
        const int* it = std::lower_bound(begin, end, row);
        return std::distance(inner_ptr, it);
    }

    static long find_diag_index(int dof, const int* outer_ptr, const int* inner_ptr) {
        return find_csr_index(dof, dof, outer_ptr, inner_ptr);
    }

    // --- Calcul local ---
    template <typename EquationType>
    FORCE_INLINE static void compute_local(int i, const Mesh& mesh, const EquationType& equation,
                               const PhysicsState& state, Eigen::MatrixXd& K_local,
                               Eigen::VectorXd& F_local, std::vector<int>& global_dofs) {
        const Element& elem = mesh.get_elements()[i];
        auto coords = mesh.get_element_coords(i);
        const int material_id = elem.ref_tag;

        ElementContext ctx{elem, coords, i, material_id, mesh};

        K_local.setZero();
        F_local.setZero();
        global_dofs.clear();

        equation.compute_element_matrices(ctx, state, K_local, F_local, global_dofs);
    }

    // --- Application des Conditions aux Limites ---
    static void apply_bcs(double max_diag, const Eigen::VectorXd& diag_global,
                           std::vector<Eigen::Triplet<double>>& triplets, Eigen::VectorXd& F,
                           double* val_ptr, const int* outer_ptr, const int* inner_ptr,
                           const std::vector<FlattenedBC>& active_bcs) {
        const double penalty = max_diag * 1e5;
        for (const auto& bc : active_bcs) {
            if (bc.is_dirichlet) {
                if (val_ptr == nullptr) triplets.emplace_back(bc.global_dof, bc.global_dof, penalty);
                else val_ptr[find_diag_index(bc.global_dof, outer_ptr, inner_ptr)] += penalty;

                F(bc.global_dof) += penalty * bc.value;
            } else if (std::abs(diag_global(bc.global_dof)) < 1e-12) {
                 if (val_ptr == nullptr) triplets.emplace_back(bc.global_dof, bc.global_dof, 1.0);
                 else val_ptr[find_diag_index(bc.global_dof, outer_ptr, inner_ptr)] = 1.0;

                 F(bc.global_dof) = 0.0;
            }
        }
    }

    // --- Algorithme 1 : Assemblage lent (Triplets) ---
    template <typename EquationType>
    static void assemble_initial_triplets(
        int num_elements, int num_dofs_total, Eigen::SparseMatrix<double>& K_global,
        Eigen::VectorXd& F_global, std::vector<long>& csr_mapping,
        const EquationType& equation, const std::vector<FlattenedBC>& active_bcs,
        const PhysicsState& state, const Mesh& mesh)
    {
        // Nom et couleur injectes statiquement a la compilation !
        PROFILE_ZONE_NC(EquationType::ZoneAssemblySlow, EquationType::Color);

        int num_threads = 1;
#ifdef _OPENMP
        num_threads = omp_get_max_threads();
#endif
        F_global = Eigen::VectorXd::Zero(num_dofs_total);
        Eigen::VectorXd diag_global = Eigen::VectorXd::Zero(num_dofs_total);

        std::vector<std::vector<Eigen::Triplet<double>>> thread_triplets(num_threads);
        std::vector<Eigen::VectorXd> thread_F(num_threads, Eigen::VectorXd::Zero(num_dofs_total));
        std::vector<Eigen::VectorXd> thread_diag(num_threads, Eigen::VectorXd::Zero(num_dofs_total));

        #pragma omp parallel
        {
#ifdef _OPENMP
            int tid = omp_get_thread_num();
#else
            int tid = 0;
#endif
            thread_triplets[tid].reserve((num_elements / num_threads + 1) * MAX_LOCAL_DOFS * MAX_LOCAL_DOFS);

            Eigen::MatrixXd K_local = Eigen::MatrixXd::Zero(MAX_LOCAL_DOFS, MAX_LOCAL_DOFS);
            Eigen::VectorXd F_local = Eigen::VectorXd::Zero(MAX_LOCAL_DOFS);
            std::vector<int> global_dofs;
            global_dofs.reserve(MAX_LOCAL_DOFS);

            #pragma omp for schedule(guided)
            for (int i = 0; i < num_elements; ++i) {
                compute_local(i, mesh, equation, state, K_local, F_local, global_dofs);
                const int n_local = static_cast<int>(global_dofs.size());

                for (int a = 0; a < n_local; ++a) {
                    const int gi = global_dofs[a];
                    thread_F[tid](gi) += F_local(a);
                    for (int b = 0; b < n_local; ++b) {
                        const int gj = global_dofs[b];
                        thread_triplets[tid].emplace_back(gi, gj, K_local(a, b));
                    }
                    thread_diag[tid](gi) += K_local(a, a);
                }
            }
        }

        // Fusion des threads
        size_t total_triplets = 0;
        for (const auto& v : thread_triplets) total_triplets += v.size();

        std::vector<Eigen::Triplet<double>> all_triplets;
        all_triplets.reserve(total_triplets);
        for (int t = 0; t < num_threads; ++t) {
            all_triplets.insert(all_triplets.end(), thread_triplets[t].begin(), thread_triplets[t].end());
            F_global += thread_F[t];
            diag_global += thread_diag[t];
        }

        const double max_diag = diag_global.size() > 0 ? diag_global.cwiseAbs().maxCoeff() : 1.0;
        apply_bcs(max_diag, diag_global, all_triplets, F_global, nullptr, nullptr, nullptr, active_bcs);

        K_global.resize(num_dofs_total, num_dofs_total);
        K_global.setFromTriplets(all_triplets.begin(), all_triplets.end());
        K_global.makeCompressed();

        // === CREATION DU CACHE CSR ===
        csr_mapping.assign(num_elements * MAX_LOCAL_DOFS * MAX_LOCAL_DOFS, -1);
        const int* outer_ptr = K_global.outerIndexPtr();
        const int* inner_ptr = K_global.innerIndexPtr();

        #pragma omp parallel
        {
            Eigen::MatrixXd K_local = Eigen::MatrixXd::Zero(MAX_LOCAL_DOFS, MAX_LOCAL_DOFS);
            Eigen::VectorXd F_local = Eigen::VectorXd::Zero(MAX_LOCAL_DOFS);
            std::vector<int> global_dofs;
            global_dofs.reserve(MAX_LOCAL_DOFS);

            #pragma omp for schedule(guided)
            for (int i = 0; i < num_elements; ++i) {
                compute_local(i, mesh, equation, state, K_local, F_local, global_dofs);
                const int n_local = static_cast<int>(global_dofs.size());

                for (int a = 0; a < n_local; ++a) {
                    for (int b = 0; b < n_local; ++b) {
                        long idx = find_csr_index(global_dofs[a], global_dofs[b], outer_ptr, inner_ptr);
                        csr_mapping[i * MAX_LOCAL_DOFS * MAX_LOCAL_DOFS + a * MAX_LOCAL_DOFS + b] = idx;
                    }
                }
            }
        }
    }

    template <typename EquationType>
    static void assemble_fast_csr_bypass(
        [[maybe_unused]] int num_elements, int num_dofs_total, Eigen::SparseMatrix<double>& K_global,
        Eigen::VectorXd& F_global, const std::vector<long>& csr_mapping,
        const EquationType& equation, const std::vector<FlattenedBC>& active_bcs,
        const PhysicsState& state, const Mesh& mesh)
    {
        // Nom et couleur injectes statiquement a la compilation !
        PROFILE_ZONE_NC(EquationType::ZoneAssemblyFast, EquationType::Color);

        F_global.setZero();
        double* val_ptr = K_global.valuePtr();
        const int* outer_ptr = K_global.outerIndexPtr();
        const int* inner_ptr = K_global.innerIndexPtr();

        std::fill(val_ptr, val_ptr + K_global.nonZeros(), 0.0);

        const auto& color_groups = mesh.get_color_groups();

        #pragma omp parallel
        {
            Eigen::MatrixXd K_local = Eigen::MatrixXd::Zero(MAX_LOCAL_DOFS, MAX_LOCAL_DOFS);
            Eigen::VectorXd F_local = Eigen::VectorXd::Zero(MAX_LOCAL_DOFS);
            std::vector<int> global_dofs;
            global_dofs.reserve(MAX_LOCAL_DOFS);

            for (size_t c = 0; c < color_groups.size(); ++c) {
                const auto& group = color_groups[c];

                #pragma omp for schedule(guided)
                for (size_t idx = 0; idx < group.size(); ++idx) {
                    int i = group[idx]; // Indice global de l'element

                    compute_local(i, mesh, equation, state, K_local, F_local, global_dofs);
                    const int n_local = static_cast<int>(global_dofs.size());

                    for (int a = 0; a < n_local; ++a) {
                        const int gi = global_dofs[a];

                        F_global(gi) += F_local(a);

                        for (int b = 0; b < n_local; ++b) {
                            const long ptr_idx = csr_mapping[i * MAX_LOCAL_DOFS * MAX_LOCAL_DOFS + a * MAX_LOCAL_DOFS + b];
                            val_ptr[ptr_idx] += K_local(a, b);
                        }
                    }
                }
            }
        }

        Eigen::VectorXd diag_global = Eigen::VectorXd::Zero(num_dofs_total);
        for (int dof = 0; dof < num_dofs_total; ++dof) {
            diag_global(dof) = val_ptr[find_diag_index(dof, outer_ptr, inner_ptr)];
        }

        const double max_diag = num_dofs_total > 0 ? diag_global.cwiseAbs().maxCoeff() : 1.0;
        std::vector<Eigen::Triplet<double>> unused;
        apply_bcs(max_diag, diag_global, unused, F_global, val_ptr, outer_ptr, inner_ptr, active_bcs);
    }
};
