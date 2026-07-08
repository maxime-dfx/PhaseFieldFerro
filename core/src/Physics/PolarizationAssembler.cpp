#include "Physics/PolarizationAssembler.h"
#include "Core/ElementIntegrator.h"
#include "Physics/Math.h"
#include "Physics/Fracture.h"
#include "Physics/Mechanics.h"
#include "Physics/Electrostatics.h"
#include "Utils/Logger.h"
#include <cmath>

// ==============================================================================
// NOUVELLE FONCTION MONOLITHIQUE
// ==============================================================================
// Construit un systeme 2*n_dof x 2*n_dof avec DOFs interleaves :
//   DOF 2*i     = Px au noeud i
//   DOF 2*i + 1 = Py au noeud i
//
// La matrice globale a une structure BLOC-DIAGONALE par noeud (2x2 par
// paire de DOFs) uniquement dans la mesure ou le schema reste semi-implicite :
// les termes de couplage physique (dW_dp1 depend de p2 et inversement, via
// b3) sont evalues avec les valeurs COURANTES de Pi_gp (comme dans
// l'ancienne version par composante), donc ils apparaissent uniquement au
// second membre (force_thermo), pas comme blocs hors-diagonale dans K.
// Le systeme reste donc bloc-diagonal (2 blocs 1x1 independants pour
// Px et Py, structurellement identiques a l'ancienne version), mais
// assembles et resolus EN UNE SEULE FOIS.
// ==============================================================================
void PolarizationAssembler::assemble_system_monolithic(
    double dt,
    const Eigen::VectorXd& Px_current, const Eigen::VectorXd& Py_current,
    const Eigen::VectorXd& Px_n, const Eigen::VectorXd& Py_n,
    const Mesh& mesh, const Fracture& fracture, const Mechanics& mechanics,
    const Electrostatics& electrostatics, const Math& math, const Datafile& config,
    const BoundaryManager& bc_manager,
    Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global)
{
    const auto& elements = mesh.get_elements();
    int n_dof = mesh.get_num_nodes();
    int system_size = 2 * n_dof;

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(elements.size() * 32); // 2x plus que la version par composante (2 champs)

    Eigen::VectorXd b = Eigen::VectorXd::Zero(system_size);
    Eigen::VectorXd diag_check_px = Eigen::VectorXd::Zero(n_dof);
    Eigen::VectorXd diag_check_py = Eigen::VectorXd::Zero(n_dof);

    const int MAX_NODES = 4;
    // Blocs locaux SEPARES pour Px et Py (pas de couplage hors-diagonale
    // dans K, conformement au schema semi-implicite - voir commentaire
    // en tete de fichier)
    Eigen::MatrixXd local_A_px(MAX_NODES, MAX_NODES), local_A_py(MAX_NODES, MAX_NODES);
    Eigen::VectorXd local_b_px(MAX_NODES), local_b_py(MAX_NODES);
    Eigen::VectorXd Px_local(MAX_NODES), Py_local(MAX_NODES), v_local(MAX_NODES);
    Eigen::VectorXd Px_local_n(MAX_NODES), Py_local_n(MAX_NODES);
    Eigen::RowVectorXd N(MAX_NODES);
    Eigen::MatrixXd grad_N(2, MAX_NODES);

    bool is_impermeable = (config.fracture_mode == CrackBCType::IMPERMEABLE);

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        const auto& indices = elem.get_node_indices();

        local_A_px.topLeftCorner(n_nodes, n_nodes).setZero();
        local_A_py.topLeftCorner(n_nodes, n_nodes).setZero();
        local_b_px.head(n_nodes).setZero();
        local_b_py.head(n_nodes).setZero();

        for (int i = 0; i < n_nodes; ++i) {
            Px_local(i)   = Px_current[indices[i]];
            Py_local(i)   = Py_current[indices[i]];
            v_local(i)    = fracture.get_v()[indices[i]];
            Px_local_n(i) = Px_n[indices[i]];
            Py_local_n(i) = Py_n[indices[i]];
        }

        ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
            double v_gp = N.head(n).dot(v_local.head(n));
            Eigen::Vector2d Pi_gp(N.head(n).dot(Px_local.head(n)), N.head(n).dot(Py_local.head(n)));
            Eigen::Matrix2d eps_gp = mechanics.get_strain_at_gp(elem, gp, coords);

            double Ex_gp = electrostatics.get_Ex_at_gp(elem, gp);
            double Ey_gp = electrostatics.get_Ey_at_gp(elem, gp);

            double penalite_fracture = (v_gp * v_gp) + config.material.eta_k;

            // Force thermodynamique pour CHAQUE composante (Px et Py),
            // calculee avec le Pi_gp courant -> couplage semi-implicite
            // conserve exactement comme dans l'ancienne version.
            double force_px = math.compute_polarization_force(0, Pi_gp, eps_gp, Ex_gp, penalite_fracture, is_impermeable);
            double force_py = math.compute_polarization_force(1, Pi_gp, eps_gp, Ey_gp, penalite_fracture, is_impermeable);

            auto N_active = N.head(n);
            auto grad_N_active = grad_N.leftCols(n);

            Eigen::MatrixXd mass_term = (config.material.mu_p / dt) * (N_active.transpose() * N_active * dV);
            Eigen::MatrixXd stiff_term = (config.material.a0 * penalite_fracture) * (grad_N_active.transpose() * grad_N_active * dV);

            // Blocs Px et Py : IDENTIQUES en structure (meme mass_term/stiff_term),
            // seul le second membre differe (force_px vs force_py) - coherent
            // avec le fait que U(∇p) et le terme de mobilite mu_p sont
            // symetriques entre les deux composantes dans ce modele.
            local_A_px.topLeftCorner(n, n).noalias() += mass_term + stiff_term;
            local_A_py.topLeftCorner(n, n).noalias() += mass_term + stiff_term;

            local_b_px.head(n).noalias() += mass_term * Px_local_n.head(n) + N_active.transpose() * force_px * dV;
            local_b_py.head(n).noalias() += mass_term * Py_local_n.head(n) + N_active.transpose() * force_py * dV;
        });

        for (int i = 0; i < n_nodes; ++i) {
            for (int j = 0; j < n_nodes; ++j) {
                double val_px = local_A_px(i, j);
                double val_py = local_A_py(i, j);

                if (val_px != 0.0) {
                    // DOF Px au noeud indices[i] -> index global 2*indices[i]
                    triplets.emplace_back(2 * indices[i], 2 * indices[j], val_px);
                    if (i == j) diag_check_px(indices[i]) += val_px;
                }
                if (val_py != 0.0) {
                    // DOF Py au noeud indices[i] -> index global 2*indices[i]+1
                    triplets.emplace_back(2 * indices[i] + 1, 2 * indices[j] + 1, val_py);
                    if (i == j) diag_check_py(indices[i]) += val_py;
                }
            }
            b(2 * indices[i])     += local_b_px(i);
            b(2 * indices[i] + 1) += local_b_py(i);
        }
    }

    // Gestion des DOFs flottants (Px et Py separement)
    handle_floating_dofs_monolithic(triplets, b, diag_check_px, diag_check_py,
                                     bc_manager.get_px_bcs(), bc_manager.get_py_bcs(),
                                     config.simulation.debug_enabled);

    K_global.setFromTriplets(triplets.begin(), triplets.end());
    F_global = b;
}

void PolarizationAssembler::handle_floating_dofs_monolithic(
    std::vector<Eigen::Triplet<double>>& triplets, Eigen::VectorXd& b,
    const Eigen::VectorXd& diag_check_px, const Eigen::VectorXd& diag_check_py,
    const std::vector<NodeBC>& bcs_px, const std::vector<NodeBC>& bcs_py,
    bool debug_enabled)
{
    const int n_dof = static_cast<int>(diag_check_px.size());

    for (int k = 0; k < n_dof; ++k) {
        bool is_dirichlet_px = (static_cast<size_t>(k) < bcs_px.size()) && (bcs_px[k].type == BCType::DIRICHLET);
        if (!is_dirichlet_px && std::abs(diag_check_px(k)) < 1e-12) {
            triplets.emplace_back(2 * k, 2 * k, 1.0);
            b(2 * k) = 0.0;
            Logger::debug("[Polarization] DOF Px flottant force a 1.0 au noeud idx=" + std::to_string(k), debug_enabled);
        }

        bool is_dirichlet_py = (static_cast<size_t>(k) < bcs_py.size()) && (bcs_py[k].type == BCType::DIRICHLET);
        if (!is_dirichlet_py && std::abs(diag_check_py(k)) < 1e-12) {
            triplets.emplace_back(2 * k + 1, 2 * k + 1, 1.0);
            b(2 * k + 1) = 0.0;
            Logger::debug("[Polarization] DOF Py flottant force a 1.0 au noeud idx=" + std::to_string(k), debug_enabled);
        }
    }
}

void PolarizationAssembler::apply_boundary_conditions_monolithic(
    Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F,
    const std::vector<NodeBC>& bcs_px, const std::vector<NodeBC>& bcs_py)
{
    double max_diag = 0.0;
    for (int k = 0; k < K.outerSize(); ++k)
        max_diag = std::max(max_diag, std::abs(K.coeff(k, k)));
    if (max_diag < 1e-12) max_diag = 1.0; // garde-fou, comme dans MechanicsAssembler

    const double penalty = max_diag * 1e7;

    for (size_t i = 0; i < bcs_px.size(); ++i) {
        if (bcs_px[i].type == BCType::DIRICHLET) {
            int dof = 2 * static_cast<int>(i);
            K.coeffRef(dof, dof) += penalty;
            F(dof) += penalty * bcs_px[i].value;
        }
    }
    for (size_t i = 0; i < bcs_py.size(); ++i) {
        if (bcs_py[i].type == BCType::DIRICHLET) {
            int dof = 2 * static_cast<int>(i) + 1;
            K.coeffRef(dof, dof) += penalty;
            F(dof) += penalty * bcs_py[i].value;
        }
    }

    K.makeCompressed();
}

// ==============================================================================
// ANCIENNES FONCTIONS PAR COMPOSANTE - conservees telles quelles pour
// compatibilite ascendante (Polarization::update_Px/update_Py). A
// supprimer une fois la migration vers update_polarization() validee.
// ==============================================================================
void PolarizationAssembler::assemble_system(
    double dt, int component,
    const Eigen::VectorXd& Px_current, const Eigen::VectorXd& Py_current,
    const Eigen::VectorXd& Px_n, const Eigen::VectorXd& Py_n,
    const Mesh& mesh, const Fracture& fracture, const Mechanics& mechanics,
    const Electrostatics& electrostatics, const Math& math, const Datafile& config,
    const BoundaryManager& bc_manager,
    Eigen::SparseMatrix<double>& K_global, Eigen::VectorXd& F_global)
{
    const auto& elements = mesh.get_elements();
    int n_dof = mesh.get_num_nodes();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(elements.size() * 16);

    Eigen::VectorXd b = Eigen::VectorXd::Zero(n_dof);
    Eigen::VectorXd diag_check = Eigen::VectorXd::Zero(n_dof);

    const int MAX_NODES = 4;
    Eigen::MatrixXd local_A(MAX_NODES, MAX_NODES);
    Eigen::VectorXd local_b(MAX_NODES);
    Eigen::VectorXd Px_local(MAX_NODES), Py_local(MAX_NODES), v_local(MAX_NODES);
    Eigen::VectorXd Px_local_n(MAX_NODES), Py_local_n(MAX_NODES);
    Eigen::RowVectorXd N(MAX_NODES);
    Eigen::MatrixXd grad_N(2, MAX_NODES);

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        const auto& indices = elem.get_node_indices();

        local_A.topLeftCorner(n_nodes, n_nodes).setZero();
        local_b.head(n_nodes).setZero();

        for (int i = 0; i < n_nodes; ++i) {
            Px_local(i)   = Px_current[indices[i]];
            Py_local(i)   = Py_current[indices[i]];
            v_local(i)    = fracture.get_v()[indices[i]];
            Px_local_n(i) = Px_n[indices[i]];
            Py_local_n(i) = Py_n[indices[i]];
        }

        ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n, double dV, const GaussPoint2D& gp) {
            double v_gp = N.head(n).dot(v_local.head(n));
            Eigen::Vector2d Pi_gp(N.head(n).dot(Px_local.head(n)), N.head(n).dot(Py_local.head(n)));
            Eigen::Matrix2d eps_gp = mechanics.get_strain_at_gp(elem, gp, coords);

            double E_gp = (component == 0) ? electrostatics.get_Ex_at_gp(elem, gp)
                                           : electrostatics.get_Ey_at_gp(elem, gp);

            const Eigen::VectorXd& P_comp_local_n = (component == 0) ? Px_local_n : Py_local_n;

            double penalite_fracture = (v_gp * v_gp) + config.material.eta_k;
            bool is_impermeable = (config.fracture_mode == CrackBCType::IMPERMEABLE);

            double force_thermo = math.compute_polarization_force(component, Pi_gp, eps_gp, E_gp, penalite_fracture, is_impermeable);

            auto N_active = N.head(n);
            auto grad_N_active = grad_N.leftCols(n);

            local_A.topLeftCorner(n, n).noalias() += (config.material.mu_p / dt) * (N_active.transpose() * N_active * dV)
                                                   + (config.material.a0 * penalite_fracture) * (grad_N_active.transpose() * grad_N_active * dV);

            local_b.head(n).noalias() += (config.material.mu_p / dt) * (N_active.transpose() * N_active * dV) * P_comp_local_n.head(n)
                                       + N_active.transpose() * force_thermo * dV;
        });

        for (int i = 0; i < n_nodes; ++i) {
            for (int j = 0; j < n_nodes; ++j) {
                double val = local_A(i, j);
                if (val != 0.0) {
                    triplets.emplace_back(indices[i], indices[j], val);
                    if (i == j) {
                        diag_check(indices[i]) += val;
                    }
                }
            }
            b(indices[i]) += local_b(i);
        }
    }

    const auto& bcs = (component == 0) ? bc_manager.get_px_bcs() : bc_manager.get_py_bcs();
    handle_floating_dofs(triplets, b, diag_check, bcs, config.simulation.debug_enabled);

    K_global.setFromTriplets(triplets.begin(), triplets.end());
    F_global = b;
}

void PolarizationAssembler::handle_floating_dofs(std::vector<Eigen::Triplet<double>>& triplets, Eigen::VectorXd& b, const Eigen::VectorXd& diag_check, const std::vector<NodeBC>& bcs, bool debug_enabled) {
    const int n_dof = static_cast<int>(b.size());

    for (int k = 0; k < n_dof; ++k) {
        bool is_dirichlet = (static_cast<size_t>(k) < bcs.size()) && (bcs[k].type == BCType::DIRICHLET);

        if (!is_dirichlet && std::abs(diag_check(k)) < 1e-12) {
            triplets.emplace_back(k, k, 1.0);
            b(k) = 0.0;
            Logger::debug("[Polarization] DOF flottant force a 1.0 au noeud idx=" + std::to_string(k), debug_enabled);
        }
    }
}

void PolarizationAssembler::apply_boundary_conditions(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, const std::vector<NodeBC>& bcs) {
    double max_diag = 0.0;
    for (int k = 0; k < K.outerSize(); ++k)
        max_diag = std::max(max_diag, std::abs(K.coeff(k, k)));
    const double penalty = max_diag * 1e7;

    for (size_t i = 0; i < bcs.size(); ++i) {
        if (bcs[i].type == BCType::DIRICHLET) {
            int dof = static_cast<int>(i);
            K.coeffRef(dof, dof) += penalty;
            F(dof) += penalty * bcs[i].value;
        }
    }
    K.makeCompressed();
}