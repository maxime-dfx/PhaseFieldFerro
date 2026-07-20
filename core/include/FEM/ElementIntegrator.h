#pragma once
#include "Mesh/Mesh.h"
#include "FEM/ShapeFunctions.h"
#include "FEM/Quadrature.h"
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
        const std::vector<GaussPoint2D>* gauss_points_ptr = nullptr;

        // 1. Selection (et mise en cache statique) des points de Gauss
        if (n_nodes == 3) {
            static const std::vector<GaussPoint2D> gp_t3 = {{1.0/3.0, 1.0/3.0, 1.0}};
            gauss_points_ptr = &gp_t3;
        } else if (n_nodes == 4) {
            static const std::vector<GaussPoint2D> gp_q4 = Quadrature::get_gauss_2x2();
            gauss_points_ptr = &gp_q4;
        } else if (n_nodes == 6) {
            static const std::vector<GaussPoint2D> gp_t6 = Quadrature::get_gauss_3_points_tri();
            gauss_points_ptr = &gp_t6;
        } else if (n_nodes == 8) {
            static const std::vector<GaussPoint2D> gp_q8 = Quadrature::get_gauss_3x3();
            gauss_points_ptr = &gp_q8;
        }
        
        if (!gauss_points_ptr) return;

        // 2. Buffers locaux fixes (0 allocation heap)
        std::array<double, 8> N_local;
        std::array<std::array<double, 8>, 2> dN_local;

        // 3. Boucle d'integration
        for (const auto& gp : *gauss_points_ptr) {
            
            // Chargement dans les buffers locaux
            if (n_nodes == 3) {
                auto N = ShapeFunctions::get_shape_functions_tri(gp.xi, gp.eta);
                auto dN = ShapeFunctions::get_shape_function_gradients_tri();
                for(int i=0; i<3; ++i) { N_local[i] = N[i]; dN_local[0][i] = dN[0][i]; dN_local[1][i] = dN[1][i]; }
            } 
            else if (n_nodes == 4) {
                auto N = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
                auto dN = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
                for(int i=0; i<4; ++i) { N_local[i] = N[i]; dN_local[0][i] = dN[0][i]; dN_local[1][i] = dN[1][i]; }
            } 
            else if (n_nodes == 6) {
                auto N = ShapeFunctions::get_shape_functions_t6(gp.xi, gp.eta);
                auto dN = ShapeFunctions::get_shape_function_gradients_t6(gp.xi, gp.eta);
                for(int i=0; i<6; ++i) { N_local[i] = N[i]; dN_local[0][i] = dN[0][i]; dN_local[1][i] = dN[1][i]; }
            } 
            else if (n_nodes == 8) {
                auto N = ShapeFunctions::get_shape_functions_q8(gp.xi, gp.eta);
                auto dN = ShapeFunctions::get_shape_function_gradients_q8(gp.xi, gp.eta);
                for(int i=0; i<8; ++i) { N_local[i] = N[i]; dN_local[0][i] = dN[0][i]; dN_local[1][i] = dN[1][i]; }
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

            // Validation et passage aux buffers Eigen
            if (std::isfinite(detJ) && std::abs(detJ) > 1e-12) {
                double weight_multiplier = (n_nodes == 3 || n_nodes == 6) ? 0.5 : 1.0; 
                double dV = gp.weight * std::abs(detJ) * weight_multiplier;

                for (int i = 0; i < n_nodes; ++i) {
                    N_buffer(i) = N_local[i];
                    grad_N_buffer(0, i) = ( dy_deta * dN_local[0][i] - dy_dxi * dN_local[1][i]) / detJ;
                    grad_N_buffer(1, i) = (-dx_deta * dN_local[0][i] + dx_dxi * dN_local[1][i]) / detJ;
                }

                compute_physics(n_nodes, dV, gp);
            }
        }
    }
};