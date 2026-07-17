#include "Physics/ElectrostaticsAssembler.h"
#include "Core/ElementIntegrator.h"
#include <omp.h>
#include <algorithm>

void ElectrostaticsAssembler::assemble_system(
    const Mesh& mesh, const Polarization& polarization, const Fracture& fracture, 
    const Math& math, const Datafile& config, const std::vector<NodeBC>& bcs,
    Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global) 
{
    (void)math;
    const auto& elements = mesh.get_elements();
    const int n_dof = static_cast<int>(mesh.get_num_nodes());
    int num_threads = omp_get_max_threads();

    // 1. Initialisation des structures par thread
    std::vector<std::vector<Eigen::Triplet<double>>> thread_triplets(num_threads);
    std::vector<Eigen::VectorXd> thread_F(num_threads, Eigen::VectorXd::Zero(n_dof));
    std::vector<Eigen::VectorXd> thread_diag(num_threads, Eigen::VectorXd::Zero(n_dof));

    for (int t = 0; t < num_threads; ++t) {
        thread_triplets[t].reserve((elements.size() * 16) / num_threads + 16);
    }

    // 2. Boucle d'assemblage principal
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
            
            // ÉTAPE A : Physique
            calculer_matrices_elementaires(elem, coords, polarization, fracture, math, config, K_local, F_local); 
            
            // ÉTAPE B : Distribution
            distribuer_local_vers_global(elem.get_node_indices(), elem.get_num_nodes(), 
                                         K_local, F_local, thread_triplets[tid], thread_F[tid], thread_diag[tid]);
        }
    }

    // 3. Fusion des données threads
    size_t total_triplets = bcs.size();
    for (int t = 0; t < num_threads; ++t) total_triplets += thread_triplets[t].size();

    std::vector<Eigen::Triplet<double>> global_triplets;
    global_triplets.reserve(total_triplets);
    F_global.setZero(n_dof);
    Eigen::VectorXd diag_global = Eigen::VectorXd::Zero(n_dof);

    for (int t = 0; t < num_threads; ++t) {
        global_triplets.insert(global_triplets.end(), thread_triplets[t].begin(), thread_triplets[t].end());
        F_global += thread_F[t];
        diag_global += thread_diag[t];
    }

    // 4. Conditions limites
    double max_diag = std::max(1.0, diag_global.cwiseAbs().maxCoeff());
    appliquer_conditions_limites(bcs, max_diag, global_triplets, F_global);

    K_global.setFromTriplets(global_triplets.begin(), global_triplets.end());
}

// --- Les Travailleurs ---

void ElectrostaticsAssembler::calculer_matrices_elementaires(
    const Element& elem, const std::vector<std::array<double, 2>>& coords,
    const Polarization& polarization, const Fracture& fracture, 
    const Math& math, const Datafile& config,
    Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local)
{
    const int n_nodes = elem.get_num_nodes();
    const auto& indices = elem.get_node_indices();
    
    double eta_k = config.material.eta_k;
    bool is_impermeable = (config.fracture.mode == CrackBCType::IMPERMEABLE);

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

    ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
        (void)gp;
        
        // --- Interpolation ---
        double v_gp = N.head(n).dot(v_local.head(n));
        Eigen::Vector2d Pi_gp(N.head(n).dot(Px_local.head(n)), N.head(n).dot(Py_local.head(n)));

        // --- 1. Physique (Permittivité et Polarisation modifiées par la fracture) ---
        double eps_eff = math.compute_effective_permittivity(v_gp, eta_k, is_impermeable);
        Eigen::Vector2d P_eff = math.compute_effective_polarization(Pi_gp, v_gp, eta_k, is_impermeable);

        // --- 2. Assemblage ---
        auto grad_N_active = grad_N.leftCols(n); 
        
        K_local.topLeftCorner(n, n).noalias() += (grad_N_active.transpose() * grad_N_active) * (eps_eff * dV);
        F_local.head(n).noalias()             += grad_N_active.transpose() * P_eff * dV;
    });
}

void ElectrostaticsAssembler::distribuer_local_vers_global(
    const std::vector<int>& indices, int n_nodes,
    const Eigen::Ref<const Eigen::MatrixXd>& K_local, const Eigen::Ref<const Eigen::VectorXd>& F_local,
    std::vector<Eigen::Triplet<double>>& thread_triplets, 
    Eigen::VectorXd& thread_F, Eigen::VectorXd& thread_diag)
{
    for (int i = 0; i < n_nodes; ++i) {
        for (int j = 0; j < n_nodes; ++j) {
            thread_triplets.emplace_back(indices[i], indices[j], K_local(i, j));
            if (i == j) thread_diag(indices[i]) += K_local(i, i);
        }
        thread_F(indices[i]) += F_local(i);
    }
}

void ElectrostaticsAssembler::appliquer_conditions_limites(
    const std::vector<NodeBC>& bcs, double max_diag, 
    std::vector<Eigen::Triplet<double>>& global_triplets, Eigen::VectorXd& F_global)
{
    const double penalty = max_diag * 1e5;
    for (size_t i = 0; i < bcs.size(); ++i) {
        if (bcs[i].type == BCType::DIRICHLET) {
            int dof = static_cast<int>(i);
            global_triplets.emplace_back(dof, dof, penalty);
            F_global(dof) += penalty * bcs[i].value;
        }
    }
}