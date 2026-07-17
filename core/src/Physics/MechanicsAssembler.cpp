#include "Physics/MechanicsAssembler.h"
#include "Core/ElementIntegrator.h"
#include <omp.h>
#include <algorithm>

void MechanicsAssembler::assemble_system(
    const Mesh& mesh, const Polarization& polarization, const Fracture& fracture, 
    const Math& math, const std::vector<NodeBC>& bcs_x, const std::vector<NodeBC>& bcs_y,
    Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) 
{
    const auto& elements = mesh.get_elements();
    const int system_size = 2 * static_cast<int>(mesh.get_num_nodes());
    int num_threads = omp_get_max_threads();

    // 1. Initialisation des structures de données par thread
    std::vector<std::vector<Eigen::Triplet<double>>> thread_triplets(num_threads);
    std::vector<Eigen::VectorXd> thread_F(num_threads, Eigen::VectorXd::Zero(system_size));
    std::vector<Eigen::VectorXd> thread_diag(num_threads, Eigen::VectorXd::Zero(system_size));

    for (int t = 0; t < num_threads; ++t) {
        thread_triplets[t].reserve((elements.size() * 64) / num_threads + 64);
    }

    // 2. Boucle d'assemblage principal
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        constexpr int MAX_ELEMENT_DOF = 8; // Pour 4 noeuds max
        Eigen::MatrixXd K_local(MAX_ELEMENT_DOF, MAX_ELEMENT_DOF);
        Eigen::VectorXd F_local(MAX_ELEMENT_DOF);

        #pragma omp for schedule(static)
        for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
            const auto& elem = elements[elem_idx];
            auto coords = mesh.get_element_coords(elem_idx);
            
            // ÉTAPE A : Calcul mathématique pur
            calculer_matrices_elementaires(elem, coords, polarization, fracture, math, K_local, F_local);
            
            // ÉTAPE B : Distribution dans les vecteurs du thread
            distribuer_local_vers_global(elem.get_node_indices(), elem.get_num_nodes(), 
                                         K_local, F_local, thread_triplets[tid], thread_F[tid], thread_diag[tid]);
        }
    }

    // 3. Fusion des données threads vers le système global
    size_t total_triplets = bcs_x.size() + bcs_y.size();
    for (int t = 0; t < num_threads; ++t) total_triplets += thread_triplets[t].size();

    std::vector<Eigen::Triplet<double>> global_triplets;
    global_triplets.reserve(total_triplets);
    F_global.setZero(system_size);
    Eigen::VectorXd diag_global = Eigen::VectorXd::Zero(system_size);

    for (int t = 0; t < num_threads; ++t) {
        global_triplets.insert(global_triplets.end(), thread_triplets[t].begin(), thread_triplets[t].end());
        F_global += thread_F[t];
        diag_global += thread_diag[t];
    }

    // 4. Application des contraintes physiques (Pénalisation)
    double max_diag = std::max(1.0, diag_global.cwiseAbs().maxCoeff());
    appliquer_conditions_limites(bcs_x, bcs_y, max_diag, global_triplets, F_global);

    // 5. Construction finale
    K_global.setFromTriplets(global_triplets.begin(), global_triplets.end());
}

void MechanicsAssembler::calculer_matrices_elementaires(
    const Element& elem, const std::vector<std::array<double, 2>>& coords,
    const Polarization& polarization, const Fracture& fracture, const Math& math,
    Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local)
{
    const int n_nodes = elem.get_num_nodes();
    const int n_local_dof = 2 * n_nodes;
    const auto& indices = elem.get_node_indices();
    
    K_local.setZero();
    F_local.setZero();

    Eigen::VectorXd Px_local(n_nodes), Py_local(n_nodes), v_local(n_nodes);
    for (int i = 0; i < n_nodes; ++i) {
        Px_local(i) = polarization.get_Px()[indices[i]];
        Py_local(i) = polarization.get_Py()[indices[i]];
        v_local(i)  = fracture.get_v()[indices[i]];
    }

    Eigen::RowVectorXd N(n_nodes);
    Eigen::MatrixXd grad_N(2, n_nodes);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, n_local_dof);

    Eigen::Matrix3d C = math.get_elastic_matrix(); 

    ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
        (void)gp;
        B.setZero();
        for (int i = 0; i < n; ++i) {
            B(0, 2 * i)     = grad_N(0, i); 
            B(1, 2 * i + 1) = grad_N(1, i); 
            B(2, 2 * i)     = grad_N(1, i); 
            B(2, 2 * i + 1) = grad_N(0, i); 
        }

        Eigen::Vector2d Pi_gp(N.head(n).dot(Px_local.head(n)), N.head(n).dot(Py_local.head(n)));
        double v_gp = N.head(n).dot(v_local.head(n));

        if (!std::isfinite(Pi_gp(0)) || !std::isfinite(Pi_gp(1)) || !std::isfinite(v_gp)) return;

        // --- 1. Physique ---
        Eigen::Vector3d sigma_0 = math.compute_sigma_0(Pi_gp);
        double penalite_fracture = (v_gp * v_gp) + math.eta_k; // Ou config.material.eta_k selon votre accès
        
        // --- 2. Assemblage (Produits purs) ---
        auto B_active = B.leftCols(2 * n);
        double factor = penalite_fracture * dV;
        
        K_local.topLeftCorner(n_local_dof, n_local_dof).noalias() += (B_active.transpose() * C * B_active) * factor;
        F_local.head(n_local_dof).noalias()                       -= B_active.transpose() * sigma_0 * factor;
    });
}

void MechanicsAssembler::distribuer_local_vers_global(
    const std::vector<int>& indices, int n_nodes,
    const Eigen::Ref<const Eigen::MatrixXd>& K_local, const Eigen::Ref<const Eigen::VectorXd>& F_local,
    std::vector<Eigen::Triplet<double>>& thread_triplets, 
    Eigen::VectorXd& thread_F, Eigen::VectorXd& thread_diag)
{
    for (int i = 0; i < n_nodes; ++i) {
        // Indices globaux (g) et locaux (l) pour le noeud i
        int g_i_x = 2 * indices[i];
        int g_i_y = 2 * indices[i] + 1;
        int l_i_x = 2 * i;
        int l_i_y = 2 * i + 1;

        for (int j = 0; j < n_nodes; ++j) {
            // Indices globaux (g) et locaux (l) pour le noeud j
            int g_j_x = 2 * indices[j];
            int g_j_y = 2 * indices[j] + 1;
            int l_j_x = 2 * j;
            int l_j_y = 2 * j + 1;

            // Ajout aux triplets pour la matrice de rigidité
            thread_triplets.emplace_back(g_i_x, g_j_x, K_local(l_i_x, l_j_x));
            thread_triplets.emplace_back(g_i_x, g_j_y, K_local(l_i_x, l_j_y));
            thread_triplets.emplace_back(g_i_y, g_j_x, K_local(l_i_y, l_j_x));
            thread_triplets.emplace_back(g_i_y, g_j_y, K_local(l_i_y, l_j_y));
        }
        
        // Ajout au vecteur second membre et sauvegarde de la diagonale
        thread_F(g_i_x) += F_local(l_i_x);
        thread_F(g_i_y) += F_local(l_i_y);
        thread_diag(g_i_x) += K_local(l_i_x, l_i_x);
        thread_diag(g_i_y) += K_local(l_i_y, l_i_y);
    }
}

void MechanicsAssembler::appliquer_conditions_limites(
    const std::vector<NodeBC>& bcs_x, const std::vector<NodeBC>& bcs_y,
    double max_diag, std::vector<Eigen::Triplet<double>>& global_triplets, Eigen::VectorXd& F_global)
{
    // Calcul du coefficient de pénalité basé sur la rigidité maximale
    const double penalty = max_diag * 1e5;

    // Conditions de Dirichlet sur l'axe X
    for (size_t i = 0; i < bcs_x.size(); ++i) {
        if (bcs_x[i].type == BCType::DIRICHLET) {
            int dof = 2 * static_cast<int>(i);
            global_triplets.emplace_back(dof, dof, penalty);
            F_global(dof) += penalty * bcs_x[i].value;
        }
    }

    // Conditions de Dirichlet sur l'axe Y
    for (size_t i = 0; i < bcs_y.size(); ++i) {
        if (bcs_y[i].type == BCType::DIRICHLET) {
            int dof = 2 * static_cast<int>(i) + 1;
            global_triplets.emplace_back(dof, dof, penalty);
            F_global(dof) += penalty * bcs_y[i].value;
        }
    }
}