#include "Physics/Mechanics.h"
#include "Physics/Polarization.h"
#include "Physics/Fracture.h"
#include "Physics/MechanicsAssembler.h"
#include "Solvers/LinearSolver.h"
#include "Core/ShapeFunctions.h"
#include "Core/ElementIntegrator.h"
#include "Core/Types.h"
#include "Utils/Logger.h"
#include <omp.h>
#include <cmath>

Mechanics::Mechanics(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager)
    : config(config), mesh(mesh), bc_manager(bc_manager) 
{
    size_t n_nodes = mesh.get_num_nodes();
    ux_prev_iter.setZero(n_nodes);
    uy_prev_iter.setZero(n_nodes);
    ux_current.setZero(n_nodes);
    uy_current.setZero(n_nodes);
    sigma_xx_current.setZero(n_nodes);
    sigma_yy_current.setZero(n_nodes);
    sigma_xy_current.setZero(n_nodes);

    // ÉTAPE 1 : Appel du travailleur dédié aux conditions initiales
    appliquer_conditions_initiales();
}

void Mechanics::appliquer_conditions_initiales() {
    size_t n_nodes = mesh.get_num_nodes();
    if (config.mechanics.type == InitializationType::UNIFORM) {
        ux_current.setConstant(config.mechanics.val_x_0);
        uy_current.setConstant(config.mechanics.val_y_0);
    } else if (config.mechanics.type == InitializationType::RANDOM) {
        for (size_t i = 0; i < n_nodes; ++i) {
            ux_current[i] = static_cast<double>(rand()) / RAND_MAX;
            uy_current[i] = static_cast<double>(rand()) / RAND_MAX;
        }
    }
}

void Mechanics::update_u(double time, const Polarization& polarization, const Fracture& fracture, const Math& math) {
    (void)time;
    size_t n_nodes = mesh.get_num_nodes();
    size_t system_size = 2 * n_nodes;

    Eigen::SparseMatrix<double> K_global(system_size, system_size);
    Eigen::VectorXd F_global = Eigen::VectorXd::Zero(system_size);

    MechanicsAssembler::assemble_system(mesh, polarization, fracture, math, bc_manager.get_ux_bcs(), bc_manager.get_uy_bcs(), K_global, F_global);

    Eigen::VectorXd U_guess(system_size);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_nodes; ++i) {
        U_guess(2 * i)     = ux_current(i);
        U_guess(2 * i + 1) = uy_current(i);
    }

    Eigen::VectorXd U_new = LinearSolver::solve_with_guess(K_global, F_global, U_guess, config.simulation.debug_enabled);

    if (U_new.allFinite()) {
        map_global_vector_to_components(U_new);
    } else {
        Logger::error("[Mechanics] Echec resolution ou solution non finie");
    }
}

void Mechanics::map_global_vector_to_components(const Eigen::VectorXd& U_new) {
    size_t n_nodes = mesh.get_num_nodes();
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_nodes; ++i) {
        ux_current(i) = U_new(2 * i);
        uy_current(i) = U_new(2 * i + 1);
    }
}

double Mechanics::calculate_error() const { return (ux_current - ux_prev_iter).norm() + (uy_current - uy_prev_iter).norm(); }
void Mechanics::save_previous_iteration() { ux_prev_iter = ux_current; uy_prev_iter = uy_current; }
void Mechanics::save_previous_state() { ux_backup = ux_current; uy_backup = uy_current; }
void Mechanics::restore_previous_state() { ux_current = ux_backup; uy_current = uy_backup; }
void Mechanics::update_history() {}

// ÉTAPE 2 : Le routeur principal pour les déformations
Eigen::Matrix2d Mechanics::get_strain_at_gp(const Element& elem, const GaussPoint2D& gp, const std::vector<std::array<double, 2>>& coords) const {
    int n_nodes = elem.get_num_nodes();

    if (n_nodes == 3) {
        return calculer_deformation_triangle(elem, coords);
    } else if (n_nodes == 4) {
        return calculer_deformation_quadrangle(elem, gp, coords);
    }

    Logger::error("[Mechanics] Type d'element non supporte pour le calcul du tenseur des deformations");
    return Eigen::Matrix2d::Zero();
}
Eigen::Matrix2d Mechanics::get_stress_at_gp(const Element& elem, const GaussPoint2D& gp,
                                             const std::vector<std::array<double, 2>>& coords,
                                             const Polarization& polarization,
                                             const Fracture& fracture,
                                             const Math& math) const
{
    const auto& indices = elem.get_node_indices();
    const int n_nodes = elem.get_num_nodes();

    // 1. Déformation totale au point de Gauss (méthode existante)
    Eigen::Matrix2d strain = get_strain_at_gp(elem, gp, coords);

    // 2. Interpolation de P et v au point de Gauss
    //    (on reproduit ici la même logique d'interpolation que dans les Assemblers)
    Eigen::RowVectorXd N;
    if (n_nodes == 3) {
        auto Nvec = ShapeFunctions::get_shape_functions_tri(1.0/3.0, 1.0/3.0); // ou coords du gp si dispo
        N = Eigen::Map<Eigen::RowVectorXd>(Nvec.data(), static_cast<int>(Nvec.size()));
    } else {
        auto Narr = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
        N = Eigen::Map<Eigen::RowVectorXd>(Narr.data(), static_cast<int>(Narr.size()));
    }

    Eigen::Vector2d P_gp = Eigen::Vector2d::Zero();
    double v_gp = 1.0;
    for (int i = 0; i < n_nodes; ++i) {
        P_gp(0) += N(i) * polarization.get_Px()[indices[i]];
        P_gp(1) += N(i) * polarization.get_Py()[indices[i]];
        v_gp    += (i == 0 ? -1.0 : 0.0); // placeholder, corrigé ci-dessous
    }
    // interpolation propre de v (on écrase le placeholder ci-dessus)
    v_gp = 0.0;
    for (int i = 0; i < n_nodes; ++i) {
        v_gp += N(i) * fracture.get_v()[indices[i]];
    }

    // 3. Contrainte spontanée (couplage électrostrictif)
    Eigen::Vector3d sigma_0_voigt = math.compute_sigma_0(P_gp);

    // 4. Contrainte élastique pure en notation de Voigt : sigma = C : eps
    Eigen::Matrix3d C = math.get_elastic_matrix();
    Eigen::Vector3d eps_voigt;
    eps_voigt(0) = strain(0,0);
    eps_voigt(1) = strain(1,1);
    eps_voigt(2) = 2.0 * strain(0,1); // cisaillement ingénieur gamma_12 = 2*eps_12_tensoriel

    Eigen::Vector3d sigma_elastic_voigt = C * eps_voigt;

    // 5. Contrainte physique totale, modulée par le jump-set de fracture
    //    sigma = (v^2 + eta_k) * [C:eps - sigma_0(P)]   (cf. MechanicsAssembler)
    double penalite_fracture = (v_gp * v_gp) + math.eta_k;
    Eigen::Vector3d sigma_voigt = penalite_fracture * (sigma_elastic_voigt - sigma_0_voigt);

    // 6. Conversion Voigt -> tenseur 2x2
    Eigen::Matrix2d sigma;
    sigma(0,0) = sigma_voigt(0);
    sigma(1,1) = sigma_voigt(1);
    sigma(0,1) = sigma(1,0) = sigma_voigt(2); // sigma_12 = sigma_21, pas de facteur 2 en notation tensorielle

    return sigma;
}

void Mechanics::compute_stress_field(const Polarization& polarization, const Fracture& fracture, const Math& math) {
    const auto& elements = mesh.get_elements();
    const size_t n_nodes = mesh.get_num_nodes();

    int num_threads = omp_get_max_threads();
    std::vector<Eigen::VectorXd> thread_Fxx(num_threads, Eigen::VectorXd::Zero(n_nodes));
    std::vector<Eigen::VectorXd> thread_Fyy(num_threads, Eigen::VectorXd::Zero(n_nodes));
    std::vector<Eigen::VectorXd> thread_Fxy(num_threads, Eigen::VectorXd::Zero(n_nodes));
    std::vector<Eigen::VectorXd> thread_M(num_threads, Eigen::VectorXd::Zero(n_nodes));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        constexpr int MAX_NODES = 4;
        Eigen::RowVectorXd N_std(MAX_NODES);
        Eigen::MatrixXd grad_N(2, MAX_NODES);

        #pragma omp for schedule(static)
        for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
            const auto& elem = elements[elem_idx];
            auto coords = mesh.get_element_coords(elem_idx);
            const auto& indices = elem.get_node_indices();
            const int n = elem.get_num_nodes();

            ElementIntegrator::integrate(elem, coords, N_std, grad_N, [&](int num_n, double dV, const GaussPoint2D& gp) {
                Eigen::Matrix2d sigma = get_stress_at_gp(elem, gp, coords, polarization, fracture, math);

                for (int i = 0; i < num_n; ++i) {
                    thread_Fxx[tid](indices[i]) += sigma(0,0) * N_std[i] * dV;
                    thread_Fyy[tid](indices[i]) += sigma(1,1) * N_std[i] * dV;
                    thread_Fxy[tid](indices[i]) += sigma(0,1) * N_std[i] * dV;
                    thread_M[tid](indices[i])   += N_std[i] * dV;
                }
            });
        }
    }

    sigma_xx_current.setZero(n_nodes);
    sigma_yy_current.setZero(n_nodes);
    sigma_xy_current.setZero(n_nodes);
    Eigen::VectorXd M_lumped = Eigen::VectorXd::Zero(n_nodes);

    for (int t = 0; t < num_threads; ++t) {
        sigma_xx_current += thread_Fxx[t];
        sigma_yy_current += thread_Fyy[t];
        sigma_xy_current += thread_Fxy[t];
        M_lumped += thread_M[t];
    }

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_nodes; ++i) {
        if (M_lumped(i) > 1e-12) {
            sigma_xx_current(i) /= M_lumped(i);
            sigma_yy_current(i) /= M_lumped(i);
            sigma_xy_current(i) /= M_lumped(i);
        }
    }
}

// Travailleur A : Triangle
Eigen::Matrix2d Mechanics::calculer_deformation_triangle(const Element& elem, const std::vector<std::array<double, 2>>& coords) const {
    const auto& indices = elem.get_node_indices();
    double eps_11 = 0.0, eps_22 = 0.0, eps_12 = 0.0;

    auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);
    if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) return Eigen::Matrix2d::Zero();
    
    for (int i = 0; i < 3; ++i) {
        eps_11 += dN_xy[i][0] * ux_current[indices[i]];
        eps_22 += dN_xy[i][1] * uy_current[indices[i]];
        eps_12 += dN_xy[i][1] * ux_current[indices[i]] + dN_xy[i][0] * uy_current[indices[i]];
    }

    eps_12 *= 0.5; 
    Eigen::Matrix2d strain;
    strain << eps_11, eps_12, eps_12, eps_22;
    return strain;
}

// Travailleur B : Quadrangle
Eigen::Matrix2d Mechanics::calculer_deformation_quadrangle(const Element& elem, const GaussPoint2D& gp, const std::vector<std::array<double, 2>>& coords) const {
    const auto& indices = elem.get_node_indices();
    int n_nodes = elem.get_num_nodes();
    double eps_11 = 0.0, eps_22 = 0.0, eps_12 = 0.0;

    auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
    auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);
    if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) return Eigen::Matrix2d::Zero();
    
    for (int i = 0; i < n_nodes; ++i) {
        eps_11 += dN_xy[0][i] * ux_current[indices[i]];
        eps_22 += dN_xy[1][i] * uy_current[indices[i]];
        eps_12 += dN_xy[1][i] * ux_current[indices[i]] + dN_xy[0][i] * uy_current[indices[i]];
    }

    eps_12 *= 0.5; 
    Eigen::Matrix2d strain;
    strain << eps_11, eps_12, eps_12, eps_22;
    return strain;
}