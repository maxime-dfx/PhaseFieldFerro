#include "Physics/PolarizationAssembler.h"
#include "Core/ElementIntegrator.h"
#include "Utils/Logger.h"
#include <omp.h>
#include <cmath>
#include <algorithm>

void PolarizationAssembler::assemble_system(
    double dt,
    const Eigen::VectorXd& Px_current, const Eigen::VectorXd& Py_current,
    const Eigen::VectorXd& Px_n, const Eigen::VectorXd& Py_n,
    const Mesh& mesh, const Fracture& fracture, const Mechanics& mechanics,
    const Electrostatics& electrostatics, const Math& math, const Datafile& config,
    const BoundaryManager& bc_manager,
    Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global)
{
    const auto& elements = mesh.get_elements();
    const int n_dof = mesh.get_num_nodes();
    const int system_size = 2 * n_dof;
    int num_threads = omp_get_max_threads();

    // 1. Initialisation
    std::vector<std::vector<Eigen::Triplet<double>>> thread_triplets(num_threads);
    std::vector<Eigen::VectorXd> thread_F(num_threads, Eigen::VectorXd::Zero(system_size));
    std::vector<Eigen::VectorXd> thread_diag_px(num_threads, Eigen::VectorXd::Zero(n_dof));
    std::vector<Eigen::VectorXd> thread_diag_py(num_threads, Eigen::VectorXd::Zero(n_dof));

    for (int t = 0; t < num_threads; ++t) thread_triplets[t].reserve((elements.size() * 64) / num_threads + 64);

    // 2. Boucle Principale
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        constexpr int MAX_ELEMENT_DOF = 8;
        Eigen::MatrixXd K_local(MAX_ELEMENT_DOF, MAX_ELEMENT_DOF);
        Eigen::VectorXd F_local(MAX_ELEMENT_DOF);

        #pragma omp for schedule(static)
        for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
            const auto& elem = elements[elem_idx];
            auto coords = mesh.get_element_coords(elem_idx);
            
            calculer_matrices_elementaires(elem, coords, dt, Px_current, Py_current, Px_n, Py_n, 
                                           fracture, mechanics, electrostatics, math, config, K_local, F_local);
            
            distribuer_local_vers_global(elem.get_node_indices(), elem.get_num_nodes(), 
                                         K_local, F_local, thread_triplets[tid], thread_F[tid], thread_diag_px[tid], thread_diag_py[tid]);
        }
    }

    // 3. Fusion des données
    const auto& bcs_px = bc_manager.get_px_bcs();
    const auto& bcs_py = bc_manager.get_py_bcs();
    size_t total_triplets = 2 * n_dof + bcs_px.size() + bcs_py.size();
    for (int t = 0; t < num_threads; ++t) total_triplets += thread_triplets[t].size();

    std::vector<Eigen::Triplet<double>> global_triplets;
    global_triplets.reserve(total_triplets);
    F_global.setZero(system_size);
    Eigen::VectorXd diag_global_px = Eigen::VectorXd::Zero(n_dof);
    Eigen::VectorXd diag_global_py = Eigen::VectorXd::Zero(n_dof);

    for (int t = 0; t < num_threads; ++t) {
        global_triplets.insert(global_triplets.end(), thread_triplets[t].begin(), thread_triplets[t].end());
        F_global += thread_F[t];
        diag_global_px += thread_diag_px[t];
        diag_global_py += thread_diag_py[t];
    }

    // 4. Plomberie : DOFs flottants et Conditions Limites
    traiter_dof_flottants(n_dof, bcs_px, bcs_py, diag_global_px, diag_global_py, global_triplets, F_global, config.simulation.debug_enabled);

    double max_diag = std::max(diag_global_px.cwiseAbs().maxCoeff(), diag_global_py.cwiseAbs().maxCoeff());
    appliquer_conditions_limites(bcs_px, bcs_py, std::max(1.0, max_diag), global_triplets, F_global);

    K_global.setFromTriplets(global_triplets.begin(), global_triplets.end());
}

// ==========================================
// --- LES TRAVAILLEURS (Implémentations) ---
// ==========================================

void PolarizationAssembler::calculer_matrices_elementaires(
    const Element& elem, const std::vector<std::array<double, 2>>& coords,
    double dt, const Eigen::VectorXd& Px_current, const Eigen::VectorXd& Py_current,
    const Eigen::VectorXd& Px_n, const Eigen::VectorXd& Py_n,
    const Fracture& fracture, const Mechanics& mechanics,
    const Electrostatics& electrostatics, const Math& math, const Datafile& config,
    Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local)
{
    const int n_nodes = elem.get_num_nodes();
    const auto& indices = elem.get_node_indices();
    bool is_impermeable = (config.fracture.mode == CrackBCType::IMPERMEABLE);

    K_local.setZero(); F_local.setZero();

    Eigen::VectorXd Px_local(n_nodes), Py_local(n_nodes), v_local(n_nodes), Px_local_n(n_nodes), Py_local_n(n_nodes);
    for (int i = 0; i < n_nodes; ++i) {
        Px_local(i) = Px_current[indices[i]]; Py_local(i) = Py_current[indices[i]];
        v_local(i) = fracture.get_v()[indices[i]];
        Px_local_n(i) = Px_n[indices[i]]; Py_local_n(i) = Py_n[indices[i]];
    }

    Eigen::RowVectorXd N(n_nodes);
    Eigen::MatrixXd grad_N(2, n_nodes);

    ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
        // 1. Interpolation locale
        double v_gp = N.head(n).dot(v_local.head(n));
        Eigen::Vector2d Pi_gp(N.head(n).dot(Px_local.head(n)), N.head(n).dot(Py_local.head(n)));
        
        Eigen::Matrix2d eps_gp = mechanics.get_strain_at_gp(elem, gp, coords);
        Eigen::Vector2d E_gp(electrostatics.get_Ex_at_gp(elem, gp), electrostatics.get_Ey_at_gp(elem, gp));

        // 2. APPEL À LA PHYSIQUE (Classe Math)
        double penalite = (v_gp * v_gp) + config.material.eta_k;
        auto GL = math.compute_GL_terms(Pi_gp, eps_gp, E_gp, penalite, is_impermeable);

        double Px_n_gp = N.head(n).dot(Px_local_n.head(n));
        double Py_n_gp = N.head(n).dot(Py_local_n.head(n));

        // 3. Assemblage local (Pur produit matriciel de la forme faible !)
        for (int i = 0; i < n; ++i) {
            F_local(2 * i)     += ((config.material.mu_p / dt) * Px_n_gp + GL.J_11 * Pi_gp(0) + GL.J_12 * Pi_gp(1) - GL.force_px) * N[i] * dV;
            F_local(2 * i + 1) += ((config.material.mu_p / dt) * Py_n_gp + GL.J_12 * Pi_gp(0) + GL.J_22 * Pi_gp(1) - GL.force_py) * N[i] * dV;

            for (int j = 0; j < n; ++j) {
                double mass = (config.material.mu_p / dt) * N[i] * N[j] * dV;
                double stiff = (config.material.a0 * penalite) * (grad_N(0, i) * grad_N(0, j) + grad_N(1, i) * grad_N(1, j)) * dV;

                K_local(2 * i, 2 * j)         += mass + stiff + GL.J_11 * N[i] * N[j] * dV;
                K_local(2 * i, 2 * j + 1)     += GL.J_12 * N[i] * N[j] * dV;
                K_local(2 * i + 1, 2 * j)     += GL.J_12 * N[i] * N[j] * dV;
                K_local(2 * i + 1, 2 * j + 1) += mass + stiff + GL.J_22 * N[i] * N[j] * dV;
            }
        }
    });
}

void PolarizationAssembler::distribuer_local_vers_global(
    const std::vector<int>& indices, int n_nodes,
    const Eigen::Ref<const Eigen::MatrixXd>& K_local, const Eigen::Ref<const Eigen::VectorXd>& F_local,
    std::vector<Eigen::Triplet<double>>& thread_triplets, 
    Eigen::VectorXd& thread_F, Eigen::VectorXd& thread_diag_px, Eigen::VectorXd& thread_diag_py)
{
    for (int i = 0; i < n_nodes; ++i) {
        int g_i_px = 2 * indices[i], g_i_py = 2 * indices[i] + 1;
        int l_i_px = 2 * i, l_i_py = 2 * i + 1;

        for (int j = 0; j < n_nodes; ++j) {
            int g_j_px = 2 * indices[j], g_j_py = 2 * indices[j] + 1;
            int l_j_px = 2 * j, l_j_py = 2 * j + 1;

            thread_triplets.emplace_back(g_i_px, g_j_px, K_local(l_i_px, l_j_px));
            thread_triplets.emplace_back(g_i_px, g_j_py, K_local(l_i_px, l_j_py));
            thread_triplets.emplace_back(g_i_py, g_j_px, K_local(l_i_py, l_j_px));
            thread_triplets.emplace_back(g_i_py, g_j_py, K_local(l_i_py, l_j_py));
        }
        
        thread_F(g_i_px) += F_local(l_i_px);
        thread_F(g_i_py) += F_local(l_i_py);
        thread_diag_px(indices[i]) += K_local(l_i_px, l_i_px);
        thread_diag_py(indices[i]) += K_local(l_i_py, l_i_py);
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