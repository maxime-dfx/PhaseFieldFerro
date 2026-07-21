#pragma once
#include "Mesh/Mesh.h"
#include "FEM/ReferenceElement.h"
#include <Eigen/Dense>
#include <cmath>
#include <array>
#include <vector>

class ElementIntegrator {
private:
    // Implémentation templatisée qui sera inlinée à 100% par le compilateur
    template<typename RefElem, typename Func>
    static void integrate_impl(const std::array<std::array<double, 2>, 8>& coords,
                               Eigen::RowVectorXd& N_buffer,
                               Eigen::MatrixXd& grad_N_buffer,
                               Func&& compute_physics) 
    {
        constexpr int n_nodes = RefElem::get_num_nodes();
        
        std::array<double, 8> N_local;
        std::array<double, 8> dN_xi_local;
        std::array<double, 8> dN_eta_local;

        for (const auto& gp : RefElem::get_gauss_points()) {
            
            RefElem::compute_shape_functions(gp.xi, gp.eta, N_local.data());
            RefElem::compute_shape_gradients(gp.xi, gp.eta, dN_xi_local.data(), dN_eta_local.data());
            
            double dx_dxi = 0.0, dy_dxi = 0.0, dx_deta = 0.0, dy_deta = 0.0;
            for (int i = 0; i < n_nodes; ++i) {
                dx_dxi  += dN_xi_local[i] * coords[i][0];
                dy_dxi  += dN_xi_local[i] * coords[i][1];
                dx_deta += dN_eta_local[i] * coords[i][0];
                dy_deta += dN_eta_local[i] * coords[i][1];
            }
            double detJ = (dx_dxi * dy_deta) - (dx_deta * dy_dxi);

            if (std::isfinite(detJ) && std::abs(detJ) > 1e-12) {
                
                double dV = gp.weight * std::abs(detJ) * RefElem::get_weight_multiplier();
                
                for (int i = 0; i < n_nodes; ++i) {
                    N_buffer(i) = N_local[i];
                    grad_N_buffer(0, i) = ( dy_deta * dN_xi_local[i] - dy_dxi * dN_eta_local[i]) / detJ;
                    grad_N_buffer(1, i) = (-dx_deta * dN_xi_local[i] + dx_dxi * dN_eta_local[i]) / detJ;
                }
                
                compute_physics(n_nodes, dV, gp);
            }
        }
    }

public:
    // Le Dispatcher de branchement (pas de vtable, pas d'allocation)
    template<typename Func>
    static void integrate(const Element& elem, 
                          const std::array<std::array<double, 2>, 8>& coords,
                          Eigen::RowVectorXd& N_buffer,
                          Eigen::MatrixXd& grad_N_buffer,
                          Func&& compute_physics) 
    {
        switch (elem.get_num_nodes()) {
            case 3:
                integrate_impl<ElementT3>(coords, N_buffer, grad_N_buffer, std::forward<Func>(compute_physics));
                break;
            case 4:
                integrate_impl<ElementQ4>(coords, N_buffer, grad_N_buffer, std::forward<Func>(compute_physics));
                break;
            case 6:
                integrate_impl<ElementT6>(coords, N_buffer, grad_N_buffer, std::forward<Func>(compute_physics));
                break;
            case 8:
                integrate_impl<ElementQ8>(coords, N_buffer, grad_N_buffer, std::forward<Func>(compute_physics));
                break;
            default:
                throw std::runtime_error("Nombre de noeuds non supporte pour l'integration : " + std::to_string(elem.get_num_nodes()));
        }
    }
};