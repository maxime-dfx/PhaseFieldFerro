#include "Physics/FractureAssembler.h"
#include "Core/ElementIntegrator.h"
#include <omp.h>
#include <algorithm>

void FractureAssembler::assemble_system(
    double dt, const Eigen::VectorXd& v_n,
    const Mesh& mesh, const Polarization& polarization, 
    const Mechanics& mechanics, const Electrostatics& electrostatics, 
    const Math& math, const Datafile& config,
    Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) 
{
    const auto& elements = mesh.get_elements();
    const int n_dof = static_cast<int>(mesh.get_num_nodes());
    int num_threads = omp_get_max_threads();

    // 1. Initialisation
    std::vector<std::vector<Eigen::Triplet<double>>> thread_triplets(num_threads);
    std::vector<Eigen::VectorXd> thread_F(num_threads, Eigen::VectorXd::Zero(n_dof));
    // NOUVEAU : Tableau pour stocker la diagonale par thread
    std::vector<Eigen::VectorXd> thread_diag(num_threads, Eigen::VectorXd::Zero(n_dof)); 

    for (int t = 0; t < num_threads; ++t) {
        thread_triplets[t].reserve((elements.size() * 16) / num_threads + 16);
    }

    // 2. Boucle Principale
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        constexpr int MAX_NODES = 4;
        Eigen::MatrixXd K_local(MAX_NODES, MAX_NODES);
        Eigen::VectorXd F_local(MAX_NODES);

        #pragma omp for schedule(static)
        for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
            const auto& elem = elements[elem_idx];
            auto coords = mesh.get_element_coords(elem_idx);
            
            calculer_matrices_elementaires(elem, coords, dt, v_n, polarization, mechanics, electrostatics, math, config, K_local, F_local);
            
            // MODIFIE : On passe thread_diag[tid] pour l'extraction de la diagonale
            distribuer_local_vers_global(elem.get_node_indices(), elem.get_num_nodes(), K_local, F_local, thread_triplets[tid], thread_F[tid], thread_diag[tid]);
        }
    }

    // 3. Fusion des données
    size_t total_triplets = n_dof;
    for (int t = 0; t < num_threads; ++t) total_triplets += thread_triplets[t].size();

    std::vector<Eigen::Triplet<double>> global_triplets;
    global_triplets.reserve(total_triplets);
    F_global.setZero(n_dof);
    // NOUVEAU : Vecteur pour la diagonale globale
    Eigen::VectorXd diag_global = Eigen::VectorXd::Zero(n_dof); 

    for (int t = 0; t < num_threads; ++t) {
        global_triplets.insert(global_triplets.end(), thread_triplets[t].begin(), thread_triplets[t].end());
        F_global += thread_F[t];
        diag_global += thread_diag[t]; // NOUVEAU : Accumulation
    }

    double max_diag = std::max(1.0, diag_global.cwiseAbs().maxCoeff());

    appliquer_irreversibilite(n_dof, v_n, max_diag, global_triplets, F_global);

    K_global.setFromTriplets(global_triplets.begin(), global_triplets.end());
}
// ==========================================
// --- LES TRAVAILLEURS (Implémentations) ---
// ==========================================

void FractureAssembler::calculer_matrices_elementaires(
    const Element& elem, const std::vector<std::array<double, 2>>& coords,
    double dt, const Eigen::VectorXd& v_n,
    const Polarization& polarization, const Mechanics& mechanics,
    const Electrostatics& electrostatics, const Math& math, const Datafile& config,
    Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local)
{
    const int n_nodes = elem.get_num_nodes();
    const auto& indices = elem.get_node_indices();
    
    bool is_impermeable = (config.fracture.mode == CrackBCType::IMPERMEABLE);
    double mu_v   = config.material.mu_v;
    double Gc     = config.material.Gc;
    double kappa  = config.material.kappa;

    K_local.setZero();
    F_local.setZero();

    Eigen::VectorXd Px_local(n_nodes), Py_local(n_nodes), v_local_n(n_nodes);
    for (int i = 0; i < n_nodes; ++i) {
        Px_local(i)  = polarization.get_Px()[indices[i]];
        Py_local(i)  = polarization.get_Py()[indices[i]];
        v_local_n(i) = v_n[indices[i]];
    }

    Eigen::RowVectorXd N(n_nodes);
    Eigen::MatrixXd grad_N(2, n_nodes);

    ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
        
        // --- Interpolation ---
        Eigen::Vector2d Pi_gp(N.head(n).dot(Px_local.head(n)), N.head(n).dot(Py_local.head(n)));
        auto grad_N_active = grad_N.leftCols(n);

        Eigen::Matrix2d grad_P_gp;
        grad_P_gp(0, 0) = grad_N_active.row(0).dot(Px_local.head(n));
        grad_P_gp(0, 1) = grad_N_active.row(1).dot(Px_local.head(n));
        grad_P_gp(1, 0) = grad_N_active.row(0).dot(Py_local.head(n));
        grad_P_gp(1, 1) = grad_N_active.row(1).dot(Py_local.head(n));

        Eigen::Matrix2d eps_gp = mechanics.get_strain_at_gp(elem, gp, coords);
        Eigen::Vector2d E_gp = Eigen::Vector2d::Zero();
        if (is_impermeable) {
            E_gp = Eigen::Vector2d(electrostatics.get_Ex_at_gp(elem, gp), electrostatics.get_Ey_at_gp(elem, gp));
        }

        // --- 1. Physique (Thermodynamique de la rupture) ---
        double H_drive = math.compute_H_drive(grad_P_gp, Pi_gp, eps_gp, E_gp, is_impermeable);

        // --- 2. Dynamique de la fracture (Équation d'Allen-Cahn) ---
        double mass_coeff = (mu_v / dt) + (Gc / (2.0 * kappa)) + 2.0 * H_drive;
        double diff_coeff = 2.0 * Gc * kappa;
        double rhs_coeff  = (mu_v / dt) * N.head(n).dot(v_local_n.head(n)) + (Gc / (2.0 * kappa));

        // --- 3. Assemblage ---
        auto N_active = N.head(n);
        K_local.topLeftCorner(n, n).noalias() += (N_active.transpose() * N_active) * (mass_coeff * dV);
        K_local.topLeftCorner(n, n).noalias() += (grad_N_active.transpose() * grad_N_active) * (diff_coeff * dV);
        F_local.head(n).noalias()             += N_active.transpose() * (rhs_coeff * dV);
    });
}

void FractureAssembler::distribuer_local_vers_global(
    const std::vector<int>& indices, int n_nodes,
    const Eigen::Ref<const Eigen::MatrixXd>& K_local, const Eigen::Ref<const Eigen::VectorXd>& F_local,
    std::vector<Eigen::Triplet<double>>& thread_triplets, Eigen::VectorXd& thread_F, Eigen::VectorXd& thread_diag) 
    {
    for (int i = 0; i < n_nodes; ++i) {
        for (int j = 0; j < n_nodes; ++j) {
            thread_triplets.emplace_back(indices[i], indices[j], K_local(i, j));
            if (i == j) thread_diag(indices[i]) += K_local(i, i); 
        }
        thread_F(indices[i]) += F_local(i);
    }
}

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