#include "Physics/Electrostatics.h"
#include "Physics/ElectrostaticsAssembler.h"
#include "Solvers/LinearSolver.h"
#include "Core/ShapeFunctions.h"
#include "Core/ElementIntegrator.h"
#include "Utils/Logger.h"
#include <omp.h>
#include <cmath>
#include <cstdlib>

Electrostatics::Electrostatics(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager) 
    : config(config), mesh(mesh), bc_manager(bc_manager) 
{
    size_t n_nodes = mesh.get_num_nodes();
    phi_current.setZero(n_nodes);
    phi_prev_iter.setZero(n_nodes);
    phi_n.setZero(n_nodes);
    
    Ex_current.setZero(n_nodes);
    Ey_current.setZero(n_nodes);

    appliquer_conditions_initiales();
}

void Electrostatics::appliquer_conditions_initiales() {
    size_t n_nodes = mesh.get_num_nodes();
    if (config.electrostatics.type == InitializationType::UNIFORM) {
        Ex_current.setConstant(config.electrostatics.val_x_0);
        Ey_current.setConstant(config.electrostatics.val_y_0);
    } else if (config.electrostatics.type == InitializationType::RANDOM) {
        for (size_t i = 0; i < n_nodes; ++i) {
            Ex_current[i] = static_cast<double>(rand()) / RAND_MAX;
            Ey_current[i] = static_cast<double>(rand()) / RAND_MAX;
        }
    }
}

void Electrostatics::update_phi(double time, const Polarization& polarization, const Fracture& fracture, const Math& math) {
    (void)time; 
    size_t n_nodes = mesh.get_num_nodes();

    Eigen::SparseMatrix<double> K_global(n_nodes, n_nodes);
    Eigen::VectorXd F_global = Eigen::VectorXd::Zero(n_nodes);

    ElectrostaticsAssembler::assemble_system(mesh, polarization, fracture, math, config, bc_manager.get_phi_bcs(), K_global, F_global);

    Eigen::VectorXd phi_new = LinearSolver::solve_with_guess(K_global, F_global, phi_current, config.simulation.debug_enabled);

    if (phi_new.allFinite()) {
        phi_current = phi_new;
        compute_electric_field();
    } else {
        Logger::error("[Electrostatics] Echec resolution ou potentiel non fini");
    }
}

void Electrostatics::compute_electric_field() {
    // La méthode de projection L2 du champ électrique (inchangée car déjà très propre et spécifique)
    const auto& elements = mesh.get_elements();
    const size_t n_nodes = mesh.get_num_nodes();

    int num_threads = omp_get_max_threads();
    std::vector<Eigen::VectorXd> thread_Fx(num_threads, Eigen::VectorXd::Zero(n_nodes));
    std::vector<Eigen::VectorXd> thread_Fy(num_threads, Eigen::VectorXd::Zero(n_nodes));
    std::vector<Eigen::VectorXd> thread_M(num_threads, Eigen::VectorXd::Zero(n_nodes));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        constexpr int MAX_NODES = 4;
        
        Eigen::RowVectorXd N_std(MAX_NODES);
        Eigen::MatrixXd grad_N(2, MAX_NODES);
        Eigen::VectorXd phi_local(MAX_NODES);

        #pragma omp for schedule(static)
        for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
            const auto& elem = elements[elem_idx];
            const int n = elem.get_num_nodes();
            auto coords = mesh.get_element_coords(elem_idx);
            const auto& indices = elem.get_node_indices();

            for (int i = 0; i < n; ++i) phi_local(i) = phi_current(indices[i]);

            ElementIntegrator::integrate(elem, coords, N_std, grad_N, [&](int num_n, double dV, const GaussPoint2D& gp) {
                (void)gp;
                double Ex_gp = -grad_N.row(0).head(num_n).dot(phi_local.head(num_n));
                double Ey_gp = -grad_N.row(1).head(num_n).dot(phi_local.head(num_n));

                for (int i = 0; i < num_n; ++i) {
                    thread_Fx[tid](indices[i]) += Ex_gp * N_std[i] * dV;
                    thread_Fy[tid](indices[i]) += Ey_gp * N_std[i] * dV;
                    thread_M[tid](indices[i])  += N_std[i] * dV;
                }
            });
        }
    }

    Ex_current.setZero(n_nodes);
    Ey_current.setZero(n_nodes);
    Eigen::VectorXd M_lumped = Eigen::VectorXd::Zero(n_nodes);

    for (int t = 0; t < num_threads; ++t) {
        Ex_current += thread_Fx[t];
        Ey_current += thread_Fy[t];
        M_lumped += thread_M[t];
    }

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_nodes; ++i) {
        if (M_lumped(i) > 1e-12) {
            Ex_current(i) /= M_lumped(i);
            Ey_current(i) /= M_lumped(i);
        }
    }
}

// --- Le Routage Géométrique ---

double Electrostatics::get_Ex_at_gp(const Element& elem, const GaussPoint2D& gp) const {
    if (elem.get_num_nodes() == 3) return calculer_champ_triangle(elem, Ex_current);
    if (elem.get_num_nodes() == 4) return calculer_champ_quadrangle(elem, gp, Ex_current);
    
    Logger::error("[Electrostatics] Type d'element non supporte");
    return 0.0;
}

double Electrostatics::get_Ey_at_gp(const Element& elem, const GaussPoint2D& gp) const {
    if (elem.get_num_nodes() == 3) return calculer_champ_triangle(elem, Ey_current);
    if (elem.get_num_nodes() == 4) return calculer_champ_quadrangle(elem, gp, Ey_current);
    
    Logger::error("[Electrostatics] Type d'element non supporte");
    return 0.0;
}

double Electrostatics::calculer_champ_triangle(const Element& elem, const Eigen::VectorXd& champ_nodal) const {
    const auto& indices = elem.get_node_indices();
    double champ_val = 0.0;
    auto N = ShapeFunctions::get_shape_functions_tri(1.0/3.0, 1.0/3.0);
    for(int i = 0; i < 3; ++i) champ_val += N[i] * champ_nodal[indices[i]];
    return champ_val;
}

double Electrostatics::calculer_champ_quadrangle(const Element& elem, const GaussPoint2D& gp, const Eigen::VectorXd& champ_nodal) const {
    const auto& indices = elem.get_node_indices();
    double champ_val = 0.0;
    auto N = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
    for(int i = 0; i < 4; ++i) champ_val += N[i] * champ_nodal[indices[i]];
    return champ_val;
}

double Electrostatics::calculate_error() const { return (phi_current - phi_prev_iter).norm(); }
void Electrostatics::save_previous_iteration() { phi_prev_iter = phi_current; }
void Electrostatics::save_previous_state() { phi_backup = phi_current; }
void Electrostatics::restore_previous_state() { phi_current = phi_backup; }
void Electrostatics::update_history() { phi_n = phi_current; }