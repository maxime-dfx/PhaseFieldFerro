#pragma once
#include <array>
#include <vector>

class ShapeFunctions {
public:
    static std::array<double,4> get_shape_functions(double xi, double eta) {
        return {
            0.25 * (1.0 - xi) * (1.0 - eta),
            0.25 * (1.0 + xi) * (1.0 - eta),
            0.25 * (1.0 + xi) * (1.0 + eta),
            0.25 * (1.0 - xi) * (1.0 + eta)
        };
    }

    static std::array<std::array<double,4>,2> get_shape_function_gradients(double xi, double eta) {
        return {{
            { -0.25*(1.0 - eta),  0.25*(1.0 - eta),  0.25*(1.0 + eta), -0.25*(1.0 + eta) }, 
            { -0.25*(1.0 - xi), -0.25*(1.0 + xi),  0.25*(1.0 + xi),  0.25*(1.0 - xi) }    
        }};
    }

    static std::pair<std::array<std::array<double,4>,2>, double>
    compute_physical_derivatives(const std::vector<std::array<double,2>>& coords,
                                 const std::array<std::array<double,4>,2>& dN_xi_eta) {
        double dx_dxi = 0, dy_dxi = 0, dx_deta = 0, dy_deta = 0;
        for (int i = 0; i < 4; ++i) {
            dx_dxi += dN_xi_eta[0][i] * coords[i][0];
            dy_dxi += dN_xi_eta[0][i] * coords[i][1];
            dx_deta += dN_xi_eta[1][i] * coords[i][0];
            dy_deta += dN_xi_eta[1][i] * coords[i][1];
        }
        double detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

        std::array<std::array<double,4>,2> dN_xy;
        for (int i = 0; i < 4; ++i) {
            dN_xy[0][i] = ( dy_deta * dN_xi_eta[0][i] - dy_dxi * dN_xi_eta[1][i]) / detJ;
            dN_xy[1][i] = (-dx_deta * dN_xi_eta[0][i] + dx_dxi * dN_xi_eta[1][i]) / detJ;
        }
        return {dN_xy, detJ};
    }
};