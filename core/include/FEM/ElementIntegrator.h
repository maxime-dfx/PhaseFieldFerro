#pragma once
#include "Mesh/Mesh.h"
#include "FEM/ReferenceElement.h" // <-- Remplacera ShapeFunctions.h et Quadrature.h
#include <Eigen/Dense>
#include <cmath>
#include <array>
#include <vector>

class ElementIntegrator {
public:
    template<typename Func>
    static void integrate(const Element& elem, 
                          const std::vector<std::array<double, 2>>& coords,
                          Eigen::RowVectorXd& N_buffer,
                          Eigen::MatrixXd& grad_N_buffer,
                          Func&& compute_physics) 
    {
        int n_nodes = elem.get_num_nodes();
        
        // 1. Récupération de l'interface polymorphe de l'élément (Résolution HORS boucle)
        const ReferenceElement& ref_elem = ElementFactory::get_element(n_nodes);
        
        // 2. Buffers locaux fixes (Max 8 noeuds, 0 allocation dynamique)
        std::array<double, 8> N_local;
        std::array<double, 8> dN_xi_local;
        std::array<double, 8> dN_eta_local;

        // 3. Boucle d'intégration universelle
        for (const auto& gp : ref_elem.get_gauss_points()) {
            
            // Délégation du calcul des fonctions de forme via l'interface
            ref_elem.compute_shape_functions(gp.xi, gp.eta, N_local.data());
            ref_elem.compute_shape_gradients(gp.xi, gp.eta, dN_xi_local.data(), dN_eta_local.data());

            // Calcul du Jacobien
            double dx_dxi = 0.0, dy_dxi = 0.0, dx_deta = 0.0, dy_deta = 0.0;
            for (int i = 0; i < n_nodes; ++i) {
                dx_dxi  += dN_xi_local[i] * coords[i][0];
                dy_dxi  += dN_xi_local[i] * coords[i][1];
                dx_deta += dN_eta_local[i] * coords[i][0];
                dy_deta += dN_eta_local[i] * coords[i][1];
            }
            double detJ = (dx_dxi * dy_deta) - (dx_deta * dy_dxi);

            // Validation et passage aux buffers Eigen
            if (std::isfinite(detJ) && std::abs(detJ) > 1e-12) {
                
                double dV = gp.weight * std::abs(detJ) * ref_elem.get_weight_multiplier();
                
                for (int i = 0; i < n_nodes; ++i) {
                    N_buffer(i) = N_local[i];
                    grad_N_buffer(0, i) = ( dy_deta * dN_xi_local[i] - dy_dxi * dN_eta_local[i]) / detJ;
                    grad_N_buffer(1, i) = (-dx_deta * dN_xi_local[i] + dx_dxi * dN_eta_local[i]) / detJ;
                }
                
                compute_physics(n_nodes, dV, gp);
            }
        }
    }
};