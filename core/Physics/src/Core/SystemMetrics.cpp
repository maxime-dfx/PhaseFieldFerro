#include "Physics/include/Core/SystemMetrics.h"
#include "Materials/Core/MaterialConcepts.h"
#include "Materials/Microstructure/Polycrystal.h"
#include "Mesh/include/FEM.h"
#include "Mesh/include/IDofMap.h"
#include "Utils/include/Profiling.h"

SystemMetrics SystemMetricsCalculator::compute(const PhysicsState& state, bool solve_failed) const {
    PROFILE_ZONE_NC("SystemMetricsCalculator::compute", PROFILE_COLOR_SEQUENTIAL);
    SystemMetrics metrics;
    metrics.solve_failed = solve_failed;
    
    double U_sum = 0.0, W_sum = 0.0, chi_sum = 0.0, elec_sum = 0.0, surf_sum = 0.0;
    bool is_impermeable = (config.fracture.mode == CrackBCType::IMPERMEABLE);
    double Gc = config.material.Gc;
    double kappa = config.material.kappa;
    
    // Extraction securisee des champs depuis le dictionnaire magique
    const Eigen::VectorXd* Px = state.get_field("Px");
    const Eigen::VectorXd* Py = state.get_field("Py");
    const Eigen::VectorXd* ux = state.get_field("ux");
    const Eigen::VectorXd* uy = state.get_field("uy");
    const Eigen::VectorXd* phi = state.get_field("phi");
    const Eigen::VectorXd* v = state.get_field("v");
    const IDofMap* p_dof_map = state.p_dof_map;

    Eigen::RowVectorXd N(8);
    Eigen::MatrixXd grad_N(2, 8);
    
    for (int e = 0; e < mesh.get_num_elements(); ++e) {
        const Element& elem = mesh.get_elements()[e];
        auto coords = mesh.get_element_coords(e);
        int mat_id = elem.ref_tag;
        
        // On recupere le materiau oriente (polycristal) s'il existe, sinon le standard
        const MaterialModel& material = state.polycrystal 
                    ? state.polycrystal->get_material(e, mat_id, materials_manager)
                    : materials_manager.get_material(mat_id);

        ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n_nodes, double dV, const GaussPoint2D& gp) {
            (void)gp;
            Eigen::Vector2d P = Eigen::Vector2d::Zero();
            Eigen::Matrix2d grad_P = Eigen::Matrix2d::Zero();
            Eigen::Matrix2d strain = Eigen::Matrix2d::Zero();
            Eigen::Vector2d E = Eigen::Vector2d::Zero();
            Eigen::Vector2d grad_v = Eigen::Vector2d::Zero();
            double v_val = 0.0;
            
            for (int i = 0; i < n_nodes; ++i) {
                int n_id = mesh.get_node_index(e, i);
                int p_id = p_dof_map ? p_dof_map->dof_for(e, i) : n_id;
                
                // Si la physique est desactivee, le pointeur est nul, on utilise 0.0
                double px_val = Px ? (*Px)(p_id) : 0.0;
                double py_val = Py ? (*Py)(p_id) : 0.0;
                double ux_val = ux ? (*ux)(n_id) : 0.0;
                double uy_val = uy ? (*uy)(n_id) : 0.0;
                double phi_val = phi ? (*phi)(n_id) : 0.0;
                double v_val_node = v ? (*v)(n_id) : 0.0;

                P(0) += px_val * N(i); P(1) += py_val * N(i);
                grad_P(0,0) += px_val * grad_N(0, i); grad_P(0,1) += px_val * grad_N(1, i);
                grad_P(1,0) += py_val * grad_N(0, i); grad_P(1,1) += py_val * grad_N(1, i);
                
                strain(0,0) += ux_val * grad_N(0, i);
                strain(1,1) += uy_val * grad_N(1, i);
                strain(0,1) += 0.5 * (ux_val * grad_N(1, i) + uy_val * grad_N(0, i));
                strain(1,0) = strain(0,1);
                
                E(0) -= phi_val * grad_N(0, i);
                E(1) -= phi_val * grad_N(1, i);
                
                v_val += v_val_node * N(i);
                grad_v(0) += v_val_node * grad_N(0, i); grad_v(1) += v_val_node * grad_N(1, i);
            }
            
            double deg = (v_val * v_val + material.get_eta_k());
            U_sum += deg * material.U_energy(grad_P) * dV;
            W_sum += deg * material.W_energy(P, strain) * dV;
            
            if (is_impermeable) {
                chi_sum  += material.chi_energy(P) * dV;
                elec_sum += deg * (-P.dot(E) - 0.5 * material.get_eps0() * E.dot(E)) * dV;
            } else {
                chi_sum  += material.chi_energy(P) * dV;
                elec_sum += (-P.dot(E) - 0.5 * material.get_eps0() * E.dot(E)) * dV;
            }
            
            double surf = ((1.0 - v_val)*(1.0 - v_val)) / (4.0 * kappa) + kappa * grad_v.dot(grad_v);
            surf_sum += Gc * surf * dV;
        });
    }
    
    metrics.U_total = U_sum; metrics.W_total = W_sum; metrics.chi_total = chi_sum;
    metrics.electric_total = elec_sum; metrics.surface_energy = surf_sum;
    metrics.bulk_enthalpy = U_sum + W_sum + chi_sum + elec_sum;
    metrics.total_energy = metrics.bulk_enthalpy + metrics.surface_energy;
    
    if (v && v->size() > 0) {
        metrics.v_min = v->minCoeff();
        metrics.v_max = v->maxCoeff();
    } else {
        metrics.v_min = 1.0;
        metrics.v_max = 1.0;
    }
    
    return metrics;
}