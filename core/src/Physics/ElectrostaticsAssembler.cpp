#include "Physics/ElectrostaticsAssembler.h"
#include <iostream>
#include <tracy/Tracy.hpp>

#ifdef _OPENMP
#include <omp.h>

#endif

// =====================================================================
// 1. Fonction principale d'assemblage
// =====================================================================
void ElectrostaticsAssembler::assemble_system(
    const Mesh& mesh, const Polarization& polarization, const Fracture& fracture, 
    const MaterialModel& material, const Datafile& config, const std::vector<NodeBC>& bcs,
    Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) 
{
    ZoneScoped;
    int num_elements = mesh.get_num_elements();
    int num_dofs_total = mesh.get_num_nodes(); // 1 DDL par noeud (potentiel phi)

    // Initialisation
    F_global = Eigen::VectorXd::Zero(num_dofs_total);

    int num_threads = 1;
#ifdef _OPENMP
    num_threads = omp_get_max_threads();
#endif

    // Tableaux de stockage thread-locaux pour supprimer le #pragma omp critical
    std::vector<std::vector<Eigen::Triplet<double>>> thread_triplets_array(num_threads);
    std::vector<Eigen::VectorXd> thread_F_array(num_threads, Eigen::VectorXd::Zero(num_dofs_total));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& thread_triplets = thread_triplets_array[tid];
        auto& thread_F = thread_F_array[tid];

        // -------------------------------------------------------------
        // A. OPTIMISATION : ALLOCATION UNIQUE PAR THREAD (Zero Heap)
        // -------------------------------------------------------------
        constexpr int max_nodes = 8;
        constexpr int dofs_per_elem = max_nodes * 1;
        constexpr int max_triplets_per_elem = dofs_per_elem * dofs_per_elem; // 64
        
        int elem_per_thread = (num_elements / num_threads) + 1;
        thread_triplets.reserve(elem_per_thread * max_triplets_per_elem);

        Eigen::VectorXd thread_diag = Eigen::VectorXd::Zero(num_dofs_total); // Local dummy pour l'appel
        
        // Buffers locaux réutilisables (0 allocation dynamique dans la boucle)
        Eigen::MatrixXd K_local_buffer = Eigen::MatrixXd::Zero(dofs_per_elem, dofs_per_elem);
        Eigen::VectorXd F_local_buffer = Eigen::VectorXd::Zero(dofs_per_elem);
        Eigen::RowVectorXd N_buffer = Eigen::RowVectorXd::Zero(max_nodes);
        Eigen::MatrixXd grad_N_buffer = Eigen::MatrixXd::Zero(2, max_nodes);

        std::vector<int> indices_buffer; 
        indices_buffer.reserve(max_nodes);

        // -------------------------------------------------------------
        // B. BOUCLE SUR LES ÉLÉMENTS (100% Parallèle, Zéro Verrou)
        // -------------------------------------------------------------
        #pragma omp for schedule(guided)
        for (int i = 0; i < num_elements; ++i) {
            const auto& elem = mesh.get_elements()[i];
            int n_nodes = elem.get_num_nodes();
            int n_dof = n_nodes * 1;

            auto coords = mesh.get_element_coords(i);
            
            // Remplissage rapide du buffer d'indices sans allocation
            indices_buffer.assign(elem.begin(), elem.end());

            // Vues ("Ref") sur la portion utile de nos buffers
            Eigen::Ref<Eigen::MatrixXd> K_local = K_local_buffer.topLeftCorner(n_dof, n_dof);
            Eigen::Ref<Eigen::VectorXd> F_local = F_local_buffer.head(n_dof);
            K_local.setZero();
            F_local.setZero();

            // 1. Calcul des matrices élémentaires
            calculer_matrices_elementaires(
                elem, coords, polarization, fracture, material, config,
                K_local, F_local, N_buffer, grad_N_buffer
            );

            // 2. Distribution (Scatter)
            distribuer_local_vers_global(
                indices_buffer, n_nodes, K_local, F_local, 
                thread_triplets, thread_F, thread_diag
            );
        }
    } // Fin #pragma omp parallel (synchronisation implicite)

    // -------------------------------------------------------------
    // C. FUSION SÉQUENTIELLE HORS PARALLÉLISME (Ultra-rapide)
    // -------------------------------------------------------------
    size_t total_triplets = 0;
    for (int t = 0; t < num_threads; ++t) {
        total_triplets += thread_triplets_array[t].size();
    }

    std::vector<Eigen::Triplet<double>> global_triplets;
    global_triplets.reserve(total_triplets);

    for (int t = 0; t < num_threads; ++t) {
        global_triplets.insert(global_triplets.end(), 
                               thread_triplets_array[t].begin(), 
                               thread_triplets_array[t].end());
        F_global += thread_F_array[t];
    }

    // Application des conditions aux limites
    double max_diag = 1e12; // Ou calculer la diagonale max si nécessaire
    appliquer_conditions_limites(bcs, max_diag, global_triplets, F_global);

    // Assemblage final du solveur
    K_global.resize(num_dofs_total, num_dofs_total);
    K_global.setFromTriplets(global_triplets.begin(), global_triplets.end());
}

// =====================================================================
// 2. Calcul Physique de l'Élément
// =====================================================================
void ElectrostaticsAssembler::calculer_matrices_elementaires(
    const Element& elem, const std::array<std::array<double, 2>, 8>& coords,
    const Polarization& polarization, const Fracture& fracture, 
    const MaterialModel& material, const Datafile& config,
    Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local,
    Eigen::RowVectorXd& N_buffer, Eigen::MatrixXd& grad_N_buffer) 
{
    ZoneScoped;
    int n_nodes_elem = elem.get_num_nodes();
    std::array<int, 8> indices;
    for (int i = 0; i < n_nodes_elem; ++i) indices[i] = elem.get_node_index(i);

    bool is_impermeable = (config.fracture.mode == CrackBCType::IMPERMEABLE);
    double eta_k = config.material.eta_k;

    auto compute_physics = [&](int n_nodes, double dV, const GaussPoint2D& gp) {
        (void)gp;
        int n_dof = n_nodes * 1;

        // Matrice B (gradient des fonctions de forme)
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2, n_dof);
        for (int i = 0; i < n_nodes; ++i) {
            B(0, i) = grad_N_buffer(0, i);
            B(1, i) = grad_N_buffer(1, i);
        }

        // --- INTERPOLATION REELLE DES CHAMPS AU POINT DE GAUSS ---
        double v_gp = 0.0;
        Eigen::Vector2d P_gp = Eigen::Vector2d::Zero();
        for (int i = 0; i < n_nodes; ++i) {
            int idx = indices[i];
            double Ni = N_buffer(i);
            v_gp    += Ni * fracture.get_v()[idx];
            P_gp(0) += Ni * polarization.get_Px()[idx];
            P_gp(1) += Ni * polarization.get_Py()[idx];
        }

        // Proprietes effectives modulees par le champ de fracture (jump-set)
        double eps_eff = material.compute_effective_permittivity(v_gp, eta_k, is_impermeable);
        Eigen::Vector2d P_eff = material.compute_effective_polarization(P_gp, v_gp, eta_k, is_impermeable);

        // Assemblage : K = int(B^T * eps_eff * B) dV ,  F = int(B^T * P_eff) dV
        K_local.topLeftCorner(n_dof, n_dof).noalias() += (B.transpose() * B) * (eps_eff * dV);
        F_local.head(n_dof).noalias()                 += B.transpose() * P_eff * dV;
    };

    ElementIntegrator::integrate(elem, coords, N_buffer, grad_N_buffer, compute_physics);
}

// =====================================================================
// 3. Distribution dans le maillage global
// =====================================================================
void ElectrostaticsAssembler::distribuer_local_vers_global(
    const std::vector<int>& indices, int n_nodes,
    const Eigen::Ref<const Eigen::MatrixXd>& K_local, const Eigen::Ref<const Eigen::VectorXd>& F_local,
    std::vector<Eigen::Triplet<double>>& thread_triplets, Eigen::VectorXd& thread_F, Eigen::VectorXd& thread_diag) 
{
    for (int i = 0; i < n_nodes; ++i) {
        int global_dof_i = indices[i];
        
        thread_F(global_dof_i) += F_local(i);

        for (int j = 0; j < n_nodes; ++j) {
            int global_dof_j = indices[j];
            thread_triplets.emplace_back(global_dof_i, global_dof_j, K_local(i, j));
        }
    }
}

// =====================================================================
// 4. Application des Conditions aux Limites
// =====================================================================
void ElectrostaticsAssembler::appliquer_conditions_limites(
    const std::vector<NodeBC>& bcs, double max_diag, 
    std::vector<Eigen::Triplet<double>>& global_triplets, Eigen::VectorXd& F_global) 
{
    // Itération directe puisque l'indice du noeud correspond au DDL (1 DDL par noeud)
    for (size_t i = 0; i < bcs.size(); ++i) {
        if (bcs[i].type == BCType::DIRICHLET) {
            global_triplets.emplace_back(i, i, max_diag);
            F_global(i) = bcs[i].value * max_diag;
        } else if (bcs[i].type == BCType::NEUMANN) {
            F_global(i) += bcs[i].value;
        }
    }
}