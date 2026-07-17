#pragma once
#include "Core/Mesh.h"
#include "Core/ShapeFunctions.h"
#include "Core/Quadrature.h"
#include <Eigen/Dense>
#include <cmath>

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
        std::vector<GaussPoint2D> gauss_points;

        // 1. Sélection des points de Gauss selon le type d'élément
        if (n_nodes == 3) {
            gauss_points = {{1.0/3.0, 1.0/3.0, 1.0}}; 
        } else if (n_nodes == 4) {
            gauss_points = Quadrature::get_gauss_2x2(); 
        } else if (n_nodes == 6) {
            // Supposant que tu as ajouté get_gauss_3_points_tri() dans Quadrature.h
            gauss_points = Quadrature::get_gauss_3_points_tri(); 
        } else if (n_nodes == 8) {
            // Supposant que tu as ajouté get_gauss_3x3() dans Quadrature.h
            gauss_points = Quadrature::get_gauss_3x3(); 
        }

        // 2. Boucle d'intégration numérique
        for (const auto& gp : gauss_points) {
            std::vector<double> N_std;
            std::array<std::vector<double>, 2> dN_local;

            // Chargement des fonctions de forme et dérivées locales
            if (n_nodes == 3) {
                N_std = ShapeFunctions::get_shape_functions_tri(gp.xi, gp.eta);
                dN_local = ShapeFunctions::get_shape_function_gradients_tri();
            } 
            else if (n_nodes == 4) {
                auto N_q4 = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
                N_std = {N_q4.begin(), N_q4.end()};
                
                auto dN_q4 = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
                dN_local[0] = {dN_q4[0].begin(), dN_q4[0].end()};
                dN_local[1] = {dN_q4[1].begin(), dN_q4[1].end()};
            } 
            else if (n_nodes == 6) {
                N_std = ShapeFunctions::get_shape_functions_t6(gp.xi, gp.eta);
                dN_local = ShapeFunctions::get_shape_function_gradients_t6(gp.xi, gp.eta);
            } 
            else if (n_nodes == 8) {
                N_std = ShapeFunctions::get_shape_functions_q8(gp.xi, gp.eta);
                dN_local = ShapeFunctions::get_shape_function_gradients_q8(gp.xi, gp.eta);
            }

            // Calcul du Jacobien
            double dx_dxi = 0.0, dy_dxi = 0.0, dx_deta = 0.0, dy_deta = 0.0;
            for (int i = 0; i < n_nodes; ++i) {
                dx_dxi  += dN_local[0][i] * coords[i][0];
                dy_dxi  += dN_local[0][i] * coords[i][1];
                dx_deta += dN_local[1][i] * coords[i][0];
                dy_deta += dN_local[1][i] * coords[i][1];
            }
            double detJ = (dx_dxi * dy_deta) - (dx_deta * dy_dxi);

            // Validation du Jacobien et passage aux buffers globaux
            if (std::isfinite(detJ) && std::abs(detJ) > 1e-12) {
                // Pour les triangles, l'aire dans l'espace de référence (0 à 1) nécessite un facteur 0.5
                double weight_multiplier = (n_nodes == 3 || n_nodes == 6) ? 0.5 : 1.0; 
                double dV = gp.weight * std::abs(detJ) * weight_multiplier;

                for (int i = 0; i < n_nodes; ++i) {
                    N_buffer(i) = N_std[i];
                    grad_N_buffer(0, i) = ( dy_deta * dN_local[0][i] - dy_dxi * dN_local[1][i]) / detJ;
                    grad_N_buffer(1, i) = (-dx_deta * dN_local[0][i] + dx_dxi * dN_local[1][i]) / detJ;
                }

                // Appel de la lambda / foncteur physique[cite: 2]
                compute_physics(n_nodes, dV, gp);
            }
        }
    }
};