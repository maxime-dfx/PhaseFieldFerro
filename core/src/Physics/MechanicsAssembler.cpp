#include "Physics/MechanicsAssembler.h"
#include <iostream>
#include <tracy/Tracy.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

// =====================================================================
// 1. Fonction principale d'assemblage
// =====================================================================
void MechanicsAssembler::assemble_system(
    const Mesh& mesh, const Polarization& polarization, const Fracture& fracture, 
    const MaterialModel& material, const std::vector<NodeBC>& bcs_x, const std::vector<NodeBC>& bcs_y,
    Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) 
{
    ZoneScoped;
    int num_elements = mesh.get_num_elements();
    int num_dofs_total = mesh.get_num_nodes() * 2; // 2 DDLs par noeud (ux, uy)

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

        // Calcul pour la réservation exacte des triplets (Max Q8 -> 8 noeuds -> 16 DDL -> 256 valeurs)
        constexpr int max_nodes = 8;
        constexpr int dofs_per_elem = max_nodes * 2;
        
        int elem_per_thread = (num_elements / num_threads) + 1;
        constexpr int max_triplets_per_elem = dofs_per_elem * dofs_per_elem;
        thread_triplets.reserve(elem_per_thread * max_triplets_per_elem);

        // Buffers locaux réutilisables (plus de création de vecteurs à la volée !)
        Eigen::MatrixXd K_local_buffer = Eigen::MatrixXd::Zero(dofs_per_elem, dofs_per_elem);
        Eigen::VectorXd F_local_buffer = Eigen::VectorXd::Zero(dofs_per_elem);
        Eigen::RowVectorXd N_buffer = Eigen::RowVectorXd::Zero(max_nodes);
        Eigen::MatrixXd grad_N_buffer = Eigen::MatrixXd::Zero(2, max_nodes);

        // Buffer pour les indices afin de respecter l'ancienne signature
        std::vector<int> indices_buffer; 
        indices_buffer.reserve(max_nodes);

        // -------------------------------------------------------------
        // B. BOUCLE SUR LES ÉLÉMENTS (100% Parallèle, Zéro Verrou)
        // -------------------------------------------------------------
        #pragma omp for schedule(guided)
        for (int i = 0; i < num_elements; ++i) {
            const auto& elem = mesh.get_elements()[i];
            int n_nodes = elem.get_num_nodes();
            int n_dof = n_nodes * 2;

            // Récupération des données géométriques
            auto coords = mesh.get_element_coords(i);
            
            // Remplissage rapide du buffer d'indices sans allocation
            indices_buffer.assign(elem.begin(), elem.end());

            // Vues ("Ref") sur la portion utile de nos buffers préalloués
            Eigen::Ref<Eigen::MatrixXd> K_local = K_local_buffer.topLeftCorner(n_dof, n_dof);
            Eigen::Ref<Eigen::VectorXd> F_local = F_local_buffer.head(n_dof);
            K_local.setZero();
            F_local.setZero();

            // 1. Calcul des matrices élémentaires
            calculer_matrices_elementaires(
                elem, coords, polarization, fracture, material, 
                K_local, F_local, N_buffer, grad_N_buffer
            );

            // 2. Distribution (Scatter) directe dans le vecteur local du thread (Zéro Conflit)
            for (int n_i = 0; n_i < n_nodes; ++n_i) {
                int idx_i = indices_buffer[n_i];
                int global_dof_x_i = 2 * idx_i;
                int global_dof_y_i = 2 * idx_i + 1;

                thread_F(global_dof_x_i) += F_local(2 * n_i);
                thread_F(global_dof_y_i) += F_local(2 * n_i + 1);

                for (int n_j = 0; n_j < n_nodes; ++n_j) {
                    int idx_j = indices_buffer[n_j];
                    int global_dof_x_j = 2 * idx_j;
                    int global_dof_y_j = 2 * idx_j + 1;

                    thread_triplets.emplace_back(global_dof_x_i, global_dof_x_j, K_local(2 * n_i, 2 * n_j));
                    thread_triplets.emplace_back(global_dof_x_i, global_dof_y_j, K_local(2 * n_i, 2 * n_j + 1));
                    thread_triplets.emplace_back(global_dof_y_i, global_dof_x_j, K_local(2 * n_i + 1, 2 * n_j));
                    thread_triplets.emplace_back(global_dof_y_i, global_dof_y_j, K_local(2 * n_i + 1, 2 * n_j + 1));
                }
            }
        }
    } // Fin de la région parallèle : la barrière implicite garantit que tous les threads ont fini.

    // =========================================================================
    // C. FUSION SÉQUENTIELLE HORS PARALLÉLISME (Ultra-rapide)
    // =========================================================================
    size_t total_triplets = 0;
    for (int t = 0; t < num_threads; ++t) {
        total_triplets += thread_triplets_array[t].size();
    }

    std::vector<Eigen::Triplet<double>> global_triplets;
    global_triplets.reserve(total_triplets); // Allocation globale unique exacte

    for (int t = 0; t < num_threads; ++t) {
        global_triplets.insert(global_triplets.end(), 
                               thread_triplets_array[t].begin(), 
                               thread_triplets_array[t].end());
        F_global += thread_F_array[t];
    }

    // Application des conditions aux limites
    double max_diag = 1e12; 
    appliquer_conditions_limites(bcs_x, bcs_y, max_diag, global_triplets, F_global);

    // Assemblage final du solveur
    K_global.resize(num_dofs_total, num_dofs_total);
    K_global.setFromTriplets(global_triplets.begin(), global_triplets.end());
}

// =====================================================================
// 2. Calcul Physique de l'Élément
// =====================================================================
void MechanicsAssembler::calculer_matrices_elementaires(
    const Element& elem, const std::array<std::array<double, 2>, 8>& coords,
    const Polarization& polarization, const Fracture& fracture, const MaterialModel& material,
    Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local,
    Eigen::RowVectorXd& N_buffer, Eigen::MatrixXd& grad_N_buffer) 
{
    ZoneScoped;
    int n_nodes_elem = elem.get_num_nodes();
    std::array<int, 8> indices;
    for (int i = 0; i < n_nodes_elem; ++i) indices[i] = elem.get_node_index(i);

    Eigen::Matrix3d C = material.get_elastic_matrix();

    auto compute_physics = [&](int n_nodes, double dV, const GaussPoint2D& gp) {
        (void)gp;
        int n_dof = n_nodes * 2;

        // Matrice B (Strain-Displacement)
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, n_dof);
        for (int i = 0; i < n_nodes; ++i) {
            B(0, 2*i)     = grad_N_buffer(0, i);
            B(1, 2*i + 1) = grad_N_buffer(1, i);
            B(2, 2*i)     = grad_N_buffer(1, i);
            B(2, 2*i + 1) = grad_N_buffer(0, i);
        }

        // --- INTERPOLATION REELLE DE v ET P AU POINT DE GAUSS ---
        double v_gp = 0.0;
        Eigen::Vector2d P_gp = Eigen::Vector2d::Zero();
        for (int i = 0; i < n_nodes; ++i) {
            int idx = indices[i];
            double Ni = N_buffer(i);
            v_gp    += Ni * fracture.get_v()[idx];
            P_gp(0) += Ni * polarization.get_Px()[idx];
            P_gp(1) += Ni * polarization.get_Py()[idx];
        }

        // Penalisation par la fracture : (v^2 + eta_k), cf. Eq. (8)/(13) du papier
        double degradation_factor = (v_gp * v_gp) + material.get_eta_k();

        // Rigidite degradee dans la zone fissuree
        K_local.topLeftCorner(n_dof, n_dof).noalias() += (B.transpose() * C * B) * (dV * degradation_factor);

        // Contrainte spontanee (couplage electrostrictif), egalement degradee
        Eigen::Vector3d sigma_0 = material.compute_sigma_0(P_gp);
        F_local.head(n_dof).noalias() -= B.transpose() * sigma_0 * (dV * degradation_factor);
    };

    ElementIntegrator::integrate(elem, coords, N_buffer, grad_N_buffer, compute_physics);
}

// =====================================================================
// 3. Distribution dans le maillage global
// =====================================================================
void MechanicsAssembler::distribuer_local_vers_global(
    const std::vector<int>& indices, int n_nodes,
    const Eigen::Ref<const Eigen::MatrixXd>& K_local, const Eigen::Ref<const Eigen::VectorXd>& F_local,
    std::vector<Eigen::Triplet<double>>& thread_triplets, Eigen::VectorXd& thread_F, Eigen::VectorXd& thread_diag) 
{
    for (int i = 0; i < n_nodes; ++i) {
        int idx_i = indices[i];
        int global_dof_x_i = 2 * idx_i;
        int global_dof_y_i = 2 * idx_i + 1;

        thread_F(global_dof_x_i) += F_local(2 * i);
        thread_F(global_dof_y_i) += F_local(2 * i + 1);

        for (int j = 0; j < n_nodes; ++j) {
            int idx_j = indices[j];
            int global_dof_x_j = 2 * idx_j;
            int global_dof_y_j = 2 * idx_j + 1;

            // Remplissage ultra-rapide (mémoire déjà réservée !)
            thread_triplets.emplace_back(global_dof_x_i, global_dof_x_j, K_local(2 * i, 2 * j));
            thread_triplets.emplace_back(global_dof_x_i, global_dof_y_j, K_local(2 * i, 2 * j + 1));
            thread_triplets.emplace_back(global_dof_y_i, global_dof_x_j, K_local(2 * i + 1, 2 * j));
            thread_triplets.emplace_back(global_dof_y_i, global_dof_y_j, K_local(2 * i + 1, 2 * j + 1));
        }
    }
}

// =====================================================================
// 4. Application des Conditions aux Limites
// =====================================================================
void MechanicsAssembler::appliquer_conditions_limites(
    const std::vector<NodeBC>& bcs_x, const std::vector<NodeBC>& bcs_y, double max_diag, 
    std::vector<Eigen::Triplet<double>>& global_triplets, Eigen::VectorXd& F_global) 
{
    // Dirichlet X
    for (size_t i = 0; i < bcs_x.size(); ++i) {
        if (bcs_x[i].type == BCType::DIRICHLET) {
            int dof = 2 * i;
            global_triplets.emplace_back(dof, dof, max_diag);
            F_global(dof) = bcs_x[i].value * max_diag;
        } else if (bcs_x[i].type == BCType::NEUMANN) {
            F_global(2 * i) += bcs_x[i].value;
        }
    }

    // Dirichlet Y
    for (size_t i = 0; i < bcs_y.size(); ++i) {
        if (bcs_y[i].type == BCType::DIRICHLET) {
            int dof = 2 * i + 1;
            global_triplets.emplace_back(dof, dof, max_diag);
            F_global(dof) = bcs_y[i].value * max_diag;
        } else if (bcs_y[i].type == BCType::NEUMANN) {
            F_global(2 * i + 1) += bcs_y[i].value;
        }
    }
}