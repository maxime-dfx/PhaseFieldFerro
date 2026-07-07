#include "Physics/Electrostatics.h"
#include "Physics/ElectrostaticsAssembler.h"
#include "Solvers/LinearSolver.h"
#include "Core/ShapeFunctions.h"
#include "Core/Quadrature.h"
#include "Utils/Logger.h"

Electrostatics::Electrostatics(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager) 
    : config(config), mesh(mesh), bc_manager(bc_manager) 
{
    size_t n_nodes = mesh.get_num_nodes();
    phi_current.setZero(n_nodes);
    Ex_current.setZero(n_nodes);
    Ey_current.setZero(n_nodes);
    Ex_prev_iter.setZero(n_nodes);
    Ey_prev_iter.setZero(n_nodes);
    
    // Initialisation
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

void Electrostatics::update_electric_potential(double time, const Polarization& polarization, const Fracture& fracture, const Math& math) {
    size_t n_nodes = mesh.get_num_nodes();
    Eigen::SparseMatrix<double> K_global(n_nodes, n_nodes);
    Eigen::VectorXd F_global = Eigen::VectorXd::Zero(n_nodes);
    
    // 1. Assemblage
    ElectrostaticsAssembler::assemble_system(mesh, polarization, fracture, math, config, K_global, F_global);

    // 2. Application des Conditions aux Limites
    ElectrostaticsAssembler::apply_boundary_conditions(K_global, F_global, bc_manager.get_phi_bcs());

    // 3. Résolution avec démarrage à chaud
    phi_current = LinearSolver::solve_with_guess(K_global, F_global, phi_current, config.simulation.debug_enabled);

    // 4. Post-traitement
    if (phi_current.allFinite()) {
        compute_nodal_electric_field(); 
    } else {
        Logger::error("[Electrostatics] Echec de la resolution ou solution non finie.");
    }
}

void Electrostatics::compute_nodal_electric_field() {
    size_t n_nodes = mesh.get_num_nodes();
    Eigen::VectorXd Fx = Eigen::VectorXd::Zero(n_nodes);
    Eigen::VectorXd Fy = Eigen::VectorXd::Zero(n_nodes);
    Eigen::VectorXd M_lumped = Eigen::VectorXd::Zero(n_nodes);

    const auto& elements = mesh.get_elements();
    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_nodes_elem = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        const auto& indices = elem.get_node_indices();
        
        Eigen::VectorXd phi_loc(n_nodes_elem);
        for(int i=0; i<n_nodes_elem; ++i) phi_loc(i) = phi_current(indices[i]);

        if (n_nodes_elem == 3) {
            auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);
            double area = std::abs(detJ) / 2.0;
            double mass_contrib = area / 3.0;

            double Ex_elem = 0.0, Ey_elem = 0.0;
            for(int i=0; i<3; ++i) {
                Ex_elem -= dN_xy[i][0] * phi_loc(i);
                Ey_elem -= dN_xy[i][1] * phi_loc(i);
            }

            for(int i=0; i<3; ++i) {
                Fx(indices[i]) += Ex_elem * mass_contrib;
                Fy(indices[i]) += Ey_elem * mass_contrib;
                M_lumped(indices[i]) += mass_contrib;
            }
        } else if (n_nodes_elem == 4) {
            const auto& gauss_points = Quadrature::get_gauss_2x2();
            for (const auto& gp : gauss_points) {
                auto N_std = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
                auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
                auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);
                
                double dV = gp.weight * std::abs(detJ);
                double Ex_gp = 0.0, Ey_gp = 0.0;
                
                for(int i = 0; i < 4; ++i) {
                    Ex_gp -= dN_xy[0][i] * phi_loc(i);
                    Ey_gp -= dN_xy[1][i] * phi_loc(i);
                }

                for(int i = 0; i < 4; ++i) {
                    Fx(indices[i]) += Ex_gp * N_std[i] * dV;
                    Fy(indices[i]) += Ey_gp * N_std[i] * dV;
                    M_lumped(indices[i]) += N_std[i] * dV;
                }
            }
        }
    }

    for(size_t i = 0; i < n_nodes; ++i) {
        if(M_lumped(i) > 1e-12) {
            Ex_current(i) = Fx(i) / M_lumped(i);
            Ey_current(i) = Fy(i) / M_lumped(i);
        }
    }
}

double Electrostatics::get_Ex_at_gp(const Element& elem, const GaussPoint2D& gp) const {
    const auto& indices = elem.get_node_indices();
    int n_nodes = elem.get_num_nodes();
    double Ex = 0.0;
    
    if (n_nodes == 3) {
        auto N = ShapeFunctions::get_shape_functions_tri(1.0/3.0, 1.0/3.0);
        for(int i=0; i<3; ++i) Ex += N[i] * Ex_current[indices[i]];
    } else if (n_nodes == 4) {
        auto N = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
        for(int i=0; i<4; ++i) Ex += N[i] * Ex_current[indices[i]];
    }
    return Ex;
}

double Electrostatics::get_Ey_at_gp(const Element& elem, const GaussPoint2D& gp) const {
    const auto& indices = elem.get_node_indices();
    int n_nodes = elem.get_num_nodes();
    double Ey = 0.0;
    
    if (n_nodes == 3) {
        auto N = ShapeFunctions::get_shape_functions_tri(1.0/3.0, 1.0/3.0);
        for(int i=0; i<3; ++i) Ey += N[i] * Ey_current[indices[i]];
    } else if (n_nodes == 4) {
        auto N = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
        for(int i=0; i<4; ++i) Ey += N[i] * Ey_current[indices[i]];
    }
    return Ey;
}