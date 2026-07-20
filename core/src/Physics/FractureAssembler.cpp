#include "Physics/FractureAssembler.h"
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

// =====================================================================
// 1. Fonction principale d'assemblage
// =====================================================================
void FractureAssembler::assemble_system(
    double dt, const Eigen::VectorXd& v_n, const Mesh& mesh, 
    const Polarization& polarization, const Mechanics& mechanics, 
    const Electrostatics& electrostatics, const MaterialModel& material, 
    const Datafile& config, Eigen::SparseMatrix<double>& K_global, 
    Eigen::VectorXd& F_global) 
{
    int num_elements = mesh.get_num_elements();
    int num_dofs_total = mesh.get_num_nodes(); // 1 DDL par noeud (champ v)

    F_global = Eigen::VectorXd::Zero(num_dofs_total);
    std::vector<Eigen::Triplet<double>> global_triplets;

    int num_threads = 1;
#ifdef _OPENMP
    num_threads = omp_get_max_threads();
#endif

    // Pré-allocation : 64 triplets max par élément (Q8 -> 8 noeuds -> 64 interactions)
    global_triplets.reserve(num_elements * 64); 

    #pragma omp parallel
    {
        // A. ALLOCATION UNIQUE PAR THREAD
        constexpr int max_nodes = 8;
        constexpr int dofs_per_elem = max_nodes;
        constexpr int max_triplets_per_elem = dofs_per_elem * dofs_per_elem;
        
        std::vector<Eigen::Triplet<double>> thread_triplets;
        int elem_per_thread = (num_elements / num_threads) + 1;
        thread_triplets.reserve(elem_per_thread * max_triplets_per_elem);

        Eigen::VectorXd thread_F = Eigen::VectorXd::Zero(num_dofs_total);
        Eigen::VectorXd thread_diag = Eigen::VectorXd::Zero(num_dofs_total);
        
        Eigen::MatrixXd K_local_buffer = Eigen::MatrixXd::Zero(dofs_per_elem, dofs_per_elem);
        Eigen::VectorXd F_local_buffer = Eigen::VectorXd::Zero(dofs_per_elem);
        Eigen::RowVectorXd N_buffer = Eigen::RowVectorXd::Zero(max_nodes);
        Eigen::MatrixXd grad_N_buffer = Eigen::MatrixXd::Zero(2, max_nodes);

        std::vector<int> indices_buffer; 
        indices_buffer.reserve(max_nodes);

        // B. BOUCLE SUR LES ÉLÉMENTS
        #pragma omp for schedule(guided)
        for (int i = 0; i < num_elements; ++i) {
            const auto& elem = mesh.get_elements()[i];
            int n_nodes = elem.get_num_nodes();
            auto coords = mesh.get_element_coords(i);
            indices_buffer.assign(elem.begin(), elem.end());

            Eigen::Ref<Eigen::MatrixXd> K_local = K_local_buffer.topLeftCorner(n_nodes, n_nodes);
            Eigen::Ref<Eigen::VectorXd> F_local = F_local_buffer.head(n_nodes);
            K_local.setZero();
            F_local.setZero();

            calculer_matrices_elementaires(
                elem, coords, dt, v_n, polarization, mechanics, electrostatics, 
                material, config, K_local, F_local, N_buffer, grad_N_buffer
            );

            distribuer_local_vers_global(
                indices_buffer, n_nodes, K_local, F_local, 
                thread_triplets, thread_F, thread_diag
            );
        }

        // C. FUSION CRITIQUE
        #pragma omp critical
        {
            global_triplets.insert(global_triplets.end(), thread_triplets.begin(), thread_triplets.end());
            F_global += thread_F;
        }
    }

    // Application de l'irréversibilité (souvent traitée post-assemblage dans les codes champ de phase)
    double max_diag = 1e12;
    appliquer_irreversibilite(num_dofs_total, v_n, max_diag, global_triplets, F_global);

    K_global.resize(num_dofs_total, num_dofs_total);
    K_global.setFromTriplets(global_triplets.begin(), global_triplets.end());
}

// =====================================================================
// 2. Calcul Physique
// =====================================================================
void FractureAssembler::calculer_matrices_elementaires(
    const Element& elem, const std::vector<std::array<double, 2>>& coords,
    double dt, const Eigen::VectorXd& v_n,
    const Polarization& polarization, const Mechanics& mechanics,
    const Electrostatics& electrostatics, const MaterialModel& material, const Datafile& config,
    Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local,
    Eigen::RowVectorXd& N_buffer, Eigen::MatrixXd& grad_N_buffer) 
{
    int n_nodes = elem.get_num_nodes();
    std::array<int, 8> indices;
    for (int i = 0; i < n_nodes; ++i) indices[i] = elem.get_node_index(i);

    bool is_impermeable = (config.fracture.mode == CrackBCType::IMPERMEABLE);
    const double mu_v  = config.material.mu_v;
    const double Gc    = config.material.Gc;
    const double kappa = config.material.kappa;

    // Polarisation locale : necessaire pour grad(P) (energie de paroi de domaine U)
    Eigen::VectorXd Px_local(n_nodes), Py_local(n_nodes), v_local_n(n_nodes);
    for (int i = 0; i < n_nodes; ++i) {
        Px_local(i)  = polarization.get_Px()[indices[i]];
        Py_local(i)  = polarization.get_Py()[indices[i]];
        v_local_n(i) = v_n[indices[i]];
    }

    auto compute_physics = [&](int num_n, double dV, const GaussPoint2D& gp) {
        // --- INTERPOLATION REELLE DE P, grad(P), strain ET E AU POINT DE GAUSS ---
        Eigen::Vector2d P_gp = Eigen::Vector2d::Zero();
        for (int i = 0; i < num_n; ++i) {
            P_gp(0) += N_buffer(i) * Px_local(i);
            P_gp(1) += N_buffer(i) * Py_local(i);
        }

        Eigen::Matrix2d grad_P_gp;
        grad_P_gp(0, 0) = grad_N_buffer.row(0).head(num_n).dot(Px_local.head(num_n));
        grad_P_gp(0, 1) = grad_N_buffer.row(1).head(num_n).dot(Px_local.head(num_n));
        grad_P_gp(1, 0) = grad_N_buffer.row(0).head(num_n).dot(Py_local.head(num_n));
        grad_P_gp(1, 1) = grad_N_buffer.row(1).head(num_n).dot(Py_local.head(num_n));

        Eigen::Matrix2d eps_gp = mechanics.get_strain_at_gp(elem, gp, coords);

        Eigen::Vector2d E_gp = Eigen::Vector2d::Zero();
        if (is_impermeable) {
            E_gp = Eigen::Vector2d(electrostatics.get_Ex_at_gp(elem, gp), electrostatics.get_Ey_at_gp(elem, gp));
        }

        // --- THERMODYNAMIQUE DE LA RUPTURE : force motrice H (fonction de eps, P, grad P, E) ---
        double H_drive = material.compute_H_drive(grad_P_gp, P_gp, eps_gp, E_gp, is_impermeable);

        // --- DYNAMIQUE DE ALLEN-CAHN (Eq. 16 du papier) ---
        double mass_coeff = (mu_v / dt) + (Gc / (2.0 * kappa)) + 2.0 * H_drive;
        double diff_coeff = 2.0 * Gc * kappa;

        double v_n_gp = N_buffer.head(num_n).dot(v_local_n.head(num_n));
        double rhs_coeff = (mu_v / dt) * v_n_gp + (Gc / (2.0 * kappa));

        for (int i = 0; i < num_n; ++i) {
            for (int j = 0; j < num_n; ++j) {
                double K_ij = N_buffer(i) * N_buffer(j) * mass_coeff * dV
                            + (grad_N_buffer(0,i)*grad_N_buffer(0,j) + grad_N_buffer(1,i)*grad_N_buffer(1,j)) * diff_coeff * dV;
                K_local(i, j) += K_ij;
            }
            F_local(i) += N_buffer(i) * rhs_coeff * dV;
        }
    };

    ElementIntegrator::integrate(elem, coords, N_buffer, grad_N_buffer, compute_physics);
}

// =====================================================================
// 3. Distribution
// =====================================================================
void FractureAssembler::distribuer_local_vers_global(
    const std::vector<int>& indices, int n_nodes,
    const Eigen::Ref<const Eigen::MatrixXd>& K_local, const Eigen::Ref<const Eigen::VectorXd>& F_local,
    std::vector<Eigen::Triplet<double>>& thread_triplets, Eigen::VectorXd& thread_F, Eigen::VectorXd& thread_diag) 
{
    for (int i = 0; i < n_nodes; ++i) {
        int idx_i = indices[i];
        thread_F(idx_i) += F_local(i);
        thread_diag(idx_i) += K_local(i, i);
        for (int j = 0; j < n_nodes; ++j) {
            thread_triplets.emplace_back(idx_i, indices[j], K_local(i, j));
        }
    }
}

// =====================================================================
// 4. Application de l'irréversibilité
// =====================================================================
void FractureAssembler::appliquer_irreversibilite(
    int n_dof, const Eigen::VectorXd& v_n, double max_diag,
    std::vector<Eigen::Triplet<double>>& global_triplets, Eigen::VectorXd& F_global)
{
    const double alpha = 2e-2;
    const double penalty = max_diag * 1e5;
   
    for (int i = 0; i < n_dof; ++i) {
        if (v_n(i) <= alpha) {
            global_triplets.emplace_back(i, i, penalty);
            F_global(i) = 0.0;
        }
    }
} 