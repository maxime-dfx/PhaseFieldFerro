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
                          Eigen::RowVectorXd& N_buffer,     // Tampons pré-alloués (Point 1)
                          Eigen::MatrixXd& grad_N_buffer,   // Tampons pré-alloués (Point 1)
                          Func&& compute_physics) 
    {
        int n = elem.get_num_nodes();

        if (n == 3) {
            double xi = 1.0 / 3.0, eta = 1.0 / 3.0;
            auto N_std = ShapeFunctions::get_shape_functions_tri(xi, eta);
            auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);

            if (std::isfinite(detJ) && std::abs(detJ) > 1e-12) {
                double dV = std::abs(detJ) / 2.0;
                for(int i = 0; i < 3; ++i) { 
                    N_buffer(i) = N_std[i]; 
                    grad_N_buffer(0, i) = dN_xy[i][0]; 
                    grad_N_buffer(1, i) = dN_xy[i][1]; 
                }
                
                GaussPoint2D gp_fake{xi, eta, 1.0};
                
                // Appel de la physique !
                compute_physics(3, dV, gp_fake);
            }
        } 
        else if (n == 4) {
            const auto& gauss_points = Quadrature::get_gauss_2x2();
            for (const auto& gp : gauss_points) {
                auto N_std = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
                auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
                auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);

                if (std::isfinite(detJ) && std::abs(detJ) > 1e-12) {
                    double dV = gp.weight * std::abs(detJ);
                    for(int i = 0; i < 4; ++i) { 
                        N_buffer(i) = N_std[i]; 
                        grad_N_buffer(0, i) = dN_xy[0][i]; 
                        grad_N_buffer(1, i) = dN_xy[1][i]; 
                    }
                    
                    // Appel de la physique !
                    compute_physics(4, dV, gp);
                }
            }
        }
    }
};