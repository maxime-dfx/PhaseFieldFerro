#include "Physics/PolarizationAssembler.h"
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

// =====================================================================
// 1. Fonction principale d'assemblage
// =====================================================================
void PolarizationAssembler::assemble_system(
        double dt,
        const Eigen::VectorXd& Px_current, const Eigen::VectorXd& Py_current,
        const Eigen::VectorXd& Px_n, const Eigen::VectorXd& Py_n,
        const Mesh& mesh, const Fracture& fracture, const Mechanics& mechanics,
        const Electrostatics& electrostatics, const MaterialModel& material, const Datafile& config,
        const BoundaryManager& bc_manager,
        Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) 
{
    int num_elements = mesh.get_num_elements();
    int num_nodes_total = mesh.get_num_nodes();
    int num_dofs_total = num_nodes_total * 2; // 2 DDLs par noeud (Px, Py)

    // Initialisation globale
    F_global = Eigen::VectorXd::Zero(num_dofs_total);
    Eigen::VectorXd diag_global_px = Eigen::VectorXd::Zero(num_nodes_total);
    Eigen::VectorXd diag_global_py = Eigen::VectorXd::Zero(num_nodes_total);
    
    std::vector<Eigen::Triplet<double>> global_triplets;

    int num_threads = 1;
#ifdef _OPENMP
    num_threads = omp_get_max_threads();
#endif

    // Pré-allocation globale (Max Q8 -> 8 noeuds -> 16 DDL -> 256 triplets)
    global_triplets.reserve(num_elements * 256); 

    #pragma omp parallel
    {
        // -------------------------------------------------------------
        // A. OPTIMISATION : ALLOCATION UNIQUE PAR THREAD (Zero Heap)
        // -------------------------------------------------------------
        constexpr int max_nodes = 8;
        constexpr int dofs_per_elem = max_nodes * 2;
        constexpr int max_triplets_per_elem = dofs_per_elem * dofs_per_elem;
        
        std::vector<Eigen::Triplet<double>> thread_triplets;
        int elem_per_thread = (num_elements / num_threads) + 1;
        thread_triplets.reserve(elem_per_thread * max_triplets_per_elem);

        Eigen::VectorXd thread_F = Eigen::VectorXd::Zero(num_dofs_total);
        Eigen::VectorXd thread_diag_px = Eigen::VectorXd::Zero(num_nodes_total);
        Eigen::VectorXd thread_diag_py = Eigen::VectorXd::Zero(num_nodes_total);
        
        // Buffers locaux réutilisables (0 allocation dans la boucle)
        Eigen::MatrixXd K_local_buffer = Eigen::MatrixXd::Zero(dofs_per_elem, dofs_per_elem);
        Eigen::VectorXd F_local_buffer = Eigen::VectorXd::Zero(dofs_per_elem);
        Eigen::RowVectorXd N_buffer = Eigen::RowVectorXd::Zero(max_nodes);
        Eigen::MatrixXd grad_N_buffer = Eigen::MatrixXd::Zero(2, max_nodes);

        std::vector<int> indices_buffer; 
        indices_buffer.reserve(max_nodes);

        // -------------------------------------------------------------
        // B. BOUCLE SUR LES ÉLÉMENTS
        // -------------------------------------------------------------
        #pragma omp for schedule(guided)
        for (int i = 0; i < num_elements; ++i) {
            const auto& elem = mesh.get_elements()[i];
            int n_nodes = elem.get_num_nodes();
            int n_dof = n_nodes * 2;

            auto coords = mesh.get_element_coords(i);
            
            // Remplissage rapide du buffer d'indices
            indices_buffer.assign(elem.begin(), elem.end());

            // Vues ("Ref") sur la portion utile de nos buffers
            Eigen::Ref<Eigen::MatrixXd> K_local = K_local_buffer.topLeftCorner(n_dof, n_dof);
            Eigen::Ref<Eigen::VectorXd> F_local = F_local_buffer.head(n_dof);
            K_local.setZero();
            F_local.setZero();

            // 1. Calcul des matrices élémentaires
            calculer_matrices_elementaires(
                elem, coords, dt, Px_current, Py_current, Px_n, Py_n,
                fracture, mechanics, electrostatics, material, config,
                K_local, F_local, N_buffer, grad_N_buffer
            );

            // 2. Distribution (Scatter) + stockage des diagonales
            distribuer_local_vers_global(
                indices_buffer, n_nodes, K_local, F_local, 
                thread_triplets, thread_F, thread_diag_px, thread_diag_py
            );
        }

        // -------------------------------------------------------------
        // C. FUSION SÉCURISÉE DES DONNÉES THREAD-LOCALES
        // -------------------------------------------------------------
        #pragma omp critical
        {
            global_triplets.insert(global_triplets.end(), thread_triplets.begin(), thread_triplets.end());
            F_global += thread_F;
            diag_global_px += thread_diag_px;
            diag_global_py += thread_diag_py;
        }
    } // Fin #pragma omp parallel

    // Traitement des DDL flottants et Conditions aux Limites
    double max_diag = 1e12; 
    
    traiter_dof_flottants(num_nodes_total, bc_manager.get_px_bcs(), bc_manager.get_py_bcs(), 
                          diag_global_px, diag_global_py, global_triplets, F_global, config.simulation.debug_enabled);

    appliquer_conditions_limites(bc_manager.get_px_bcs(), bc_manager.get_py_bcs(), max_diag, global_triplets, F_global);

    // Assemblage final
    K_global.resize(num_dofs_total, num_dofs_total);
    K_global.setFromTriplets(global_triplets.begin(), global_triplets.end());
}

// =====================================================================
// 2. Calcul Physique de l'Élément
// =====================================================================
void PolarizationAssembler::calculer_matrices_elementaires(
    const Element& elem, const std::vector<std::array<double, 2>>& coords,
    double dt, const Eigen::VectorXd& Px_current, const Eigen::VectorXd& Py_current,
    const Eigen::VectorXd& Px_n, const Eigen::VectorXd& Py_n,
    const Fracture& fracture, const Mechanics& mechanics,
    const Electrostatics& electrostatics, const MaterialModel& material, const Datafile& config,
    Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local,
    Eigen::RowVectorXd& N_buffer, Eigen::MatrixXd& grad_N_buffer) 
{
    int n_nodes_elem = elem.get_num_nodes();
    std::array<int, 8> indices;
    for (int i = 0; i < n_nodes_elem; ++i) indices[i] = elem.get_node_index(i);

    bool is_impermeable = (config.fracture.mode == CrackBCType::IMPERMEABLE);
    const double mu_p = config.material.mu_p;
    const double a0   = config.material.a0;

    auto compute_physics = [&](int n_nodes, double dV, const GaussPoint2D& gp) {
        // --- INTERPOLATION REELLE DES CHAMPS AU POINT DE GAUSS ---
        double v_gp = 0.0;
        double Px_n_gp = 0.0, Py_n_gp = 0.0;
        Eigen::Vector2d P_gp = Eigen::Vector2d::Zero();
        for (int i = 0; i < n_nodes; ++i) {
            int idx = indices[i];
            double Ni = N_buffer(i);
            P_gp(0) += Ni * Px_current[idx];
            P_gp(1) += Ni * Py_current[idx];
            v_gp    += Ni * fracture.get_v()[idx];
            Px_n_gp += Ni * Px_n[idx];   // etat de reference temporel (m-1), cf. correctif precedent
            Py_n_gp += Ni * Py_n[idx];
        }

        // Deformation et champ electrique reels, couples via Mechanics/Electrostatics
        Eigen::Matrix2d strain = mechanics.get_strain_at_gp(elem, gp, coords);
        Eigen::Vector2d E_gp(electrostatics.get_Ex_at_gp(elem, gp), electrostatics.get_Ey_at_gp(elem, gp));

        double penalite = (v_gp * v_gp) + config.material.eta_k;
        GinzburgLandauTerms GL = material.compute_GL_terms(P_gp, strain, E_gp, penalite, is_impermeable);

        double inv_dt = mu_p / dt;

        for (int i = 0; i < n_nodes; ++i) {
            for (int j = 0; j < n_nodes; ++j) {
                double mass  = inv_dt * N_buffer(i) * N_buffer(j) * dV;
                double stiff = (a0 * penalite) *
                               (grad_N_buffer(0,i)*grad_N_buffer(0,j) + grad_N_buffer(1,i)*grad_N_buffer(1,j)) * dV;

                K_local(2*i,   2*j)   += mass + stiff + GL.J_11 * N_buffer(i) * N_buffer(j) * dV;
                K_local(2*i+1, 2*j+1) += mass + stiff + GL.J_22 * N_buffer(i) * N_buffer(j) * dV;
                K_local(2*i,   2*j+1) += GL.J_12 * N_buffer(i) * N_buffer(j) * dV;
                K_local(2*i+1, 2*j)   += GL.J_12 * N_buffer(i) * N_buffer(j) * dV;
            }

            F_local(2*i)   += (inv_dt * Px_n_gp + GL.J_11 * P_gp(0) + GL.J_12 * P_gp(1) - GL.force_px) * N_buffer(i) * dV;
            F_local(2*i+1) += (inv_dt * Py_n_gp + GL.J_12 * P_gp(0) + GL.J_22 * P_gp(1) - GL.force_py) * N_buffer(i) * dV;
        }
    };

    ElementIntegrator::integrate(elem, coords, N_buffer, grad_N_buffer, compute_physics);
}

// =====================================================================
// 3. Distribution dans le maillage global
// =====================================================================
void PolarizationAssembler::distribuer_local_vers_global(
    const std::vector<int>& indices, int n_nodes,
    const Eigen::Ref<const Eigen::MatrixXd>& K_local, const Eigen::Ref<const Eigen::VectorXd>& F_local,
    std::vector<Eigen::Triplet<double>>& thread_triplets, 
    Eigen::VectorXd& thread_F, Eigen::VectorXd& thread_diag_px, Eigen::VectorXd& thread_diag_py) 
{
    for (int i = 0; i < n_nodes; ++i) {
        int idx_i = indices[i];
        int global_dof_x_i = 2 * idx_i;
        int global_dof_y_i = 2 * idx_i + 1;

        thread_F(global_dof_x_i) += F_local(2 * i);
        thread_F(global_dof_y_i) += F_local(2 * i + 1);
        
        // Stockage de la diagonale pour le traitement des DDL flottants
        thread_diag_px(idx_i) += K_local(2 * i, 2 * i);
        thread_diag_py(idx_i) += K_local(2 * i + 1, 2 * i + 1);

        for (int j = 0; j < n_nodes; ++j) {
            int idx_j = indices[j];
            int global_dof_x_j = 2 * idx_j;
            int global_dof_y_j = 2 * idx_j + 1;

            thread_triplets.emplace_back(global_dof_x_i, global_dof_x_j, K_local(2 * i, 2 * j));
            thread_triplets.emplace_back(global_dof_x_i, global_dof_y_j, K_local(2 * i, 2 * j + 1));
            thread_triplets.emplace_back(global_dof_y_i, global_dof_x_j, K_local(2 * i + 1, 2 * j));
            thread_triplets.emplace_back(global_dof_y_i, global_dof_y_j, K_local(2 * i + 1, 2 * j + 1));
        }
    }
}

void PolarizationAssembler::traiter_dof_flottants(
    int n_dof, const std::vector<NodeBC>& bcs_px, const std::vector<NodeBC>& bcs_py,
    const Eigen::VectorXd& diag_global_px, const Eigen::VectorXd& diag_global_py,
    std::vector<Eigen::Triplet<double>>& global_triplets, Eigen::VectorXd& F_global, bool debug_enabled)
{
    for (int k = 0; k < n_dof; ++k) {
        bool is_dirichlet_px = (static_cast<size_t>(k) < bcs_px.size()) && (bcs_px[k].type == BCType::DIRICHLET);
        if (!is_dirichlet_px && std::abs(diag_global_px(k)) < 1e-12) {
            global_triplets.emplace_back(2 * k, 2 * k, 1.0);
            F_global(2 * k) = 0.0;
            Logger::debug("[Polarization] DOF Px flottant force a 1.0", debug_enabled);
        }

        bool is_dirichlet_py = (static_cast<size_t>(k) < bcs_py.size()) && (bcs_py[k].type == BCType::DIRICHLET);
        if (!is_dirichlet_py && std::abs(diag_global_py(k)) < 1e-12) {
            global_triplets.emplace_back(2 * k + 1, 2 * k + 1, 1.0);
            F_global(2 * k + 1) = 0.0;
            Logger::debug("[Polarization] DOF Py flottant force a 1.0", debug_enabled);
        }
    }
}

void PolarizationAssembler::appliquer_conditions_limites(
    const std::vector<NodeBC>& bcs_px, const std::vector<NodeBC>& bcs_py,
    double max_diag, std::vector<Eigen::Triplet<double>>& global_triplets, Eigen::VectorXd& F_global)
{
    const double penalty = max_diag * 1e5;
    for (size_t i = 0; i < bcs_px.size(); ++i) {
        if (bcs_px[i].type == BCType::DIRICHLET) {
            int dof = 2 * static_cast<int>(i);
            global_triplets.emplace_back(dof, dof, penalty);
            F_global(dof) += penalty * bcs_px[i].value;
        }
    }
    for (size_t i = 0; i < bcs_py.size(); ++i) {
        if (bcs_py[i].type == BCType::DIRICHLET) {
            int dof = 2 * static_cast<int>(i) + 1;
            global_triplets.emplace_back(dof, dof, penalty);
            F_global(dof) += penalty * bcs_py[i].value;
        }
    }
}