#pragma once
#include <array>
#include <vector>
#include <utility>

class ShapeFunctions {
public:
    // --- 1. TRIANGLE (T3) ---
    static std::array<double, 3> get_shape_functions_tri(double xi, double eta) { 
        return { 1.0 - xi - eta, xi, eta }; 
    }
    static std::array<std::array<double, 3>, 2> get_shape_function_gradients_tri() { 
        return {{ {-1.0, 1.0, 0.0}, {-1.0, 0.0, 1.0} }}; 
    }

    // --- 2. QUADRANGLE (Q4) ---
    static std::array<double, 4> get_shape_functions(double xi, double eta) {
        return { 0.25*(1-xi)*(1-eta), 0.25*(1+xi)*(1-eta), 0.25*(1+xi)*(1+eta), 0.25*(1-xi)*(1+eta) };
    }
    static std::array<std::array<double, 4>, 2> get_shape_function_gradients(double xi, double eta) {
        return {{ { -0.25*(1-eta), 0.25*(1-eta), 0.25*(1+eta), -0.25*(1+eta) }, 
                  { -0.25*(1-xi), -0.25*(1+xi), 0.25*(1+xi), 0.25*(1-xi) } }};
    }

    // --- 3. TRIANGLE (T6) ---
    static std::array<double, 6> get_shape_functions_t6(double xi, double eta) {
        double L1 = 1.0 - xi - eta;
        return { L1*(2*L1-1), xi*(2*xi-1), eta*(2*eta-1), 4*L1*xi, 4*xi*eta, 4*L1*eta };
    }
    static std::array<std::array<double, 6>, 2> get_shape_function_gradients_t6(double xi, double eta) {
        double L1 = 1.0 - xi - eta;
        return {{ {1-4*L1, 4*xi-1, 0, 4*(L1-xi), 4*eta, -4*eta}, 
                  {1-4*L1, 0, 4*eta-1, -4*xi, 4*xi, 4*(L1-eta)} }};
    }

    // --- 4. QUADRANGLE (Q8) ---
    static std::array<double, 8> get_shape_functions_q8(double xi, double eta) {
        double x2 = xi*xi, e2 = eta*eta;
        return { 0.25*(1-xi)*(1-eta)*(-xi-eta-1), 0.25*(1+xi)*(1-eta)*(xi-eta-1), 0.25*(1+xi)*(1+eta)*(xi+eta-1), 0.25*(1-xi)*(1+eta)*(-xi+eta-1), 0.5*(1-x2)*(1-eta), 0.5*(1+xi)*(1-e2), 0.5*(1-x2)*(1+eta), 0.5*(1-xi)*(1-e2) };
    }
    static std::array<std::array<double, 8>, 2> get_shape_function_gradients_q8(double xi, double eta) {
        return {{
            { 0.25*(1-eta)*(2*xi+eta), 0.25*(1-eta)*(2*xi-eta), 0.25*(1+eta)*(2*xi+eta), 0.25*(1+eta)*(2*xi-eta), -xi*(1-eta), 0.5*(1-eta*eta), -xi*(1+eta), -0.5*(1-eta*eta) },
            { 0.25*(1-xi)*(2*eta+xi), 0.25*(1+xi)*(2*eta-xi), 0.25*(1+xi)*(2*eta+xi), 0.25*(1-xi)*(2*eta-xi), -0.5*(1-xi*xi), -eta*(1+xi), 0.5*(1-xi*xi), -eta*(1-xi) }
        }};
    }

    // --- AJOUTS : CALCULATEURS DE DÉRIVÉES PHYSIQUES ---

    static std::pair<std::array<std::array<double,2>, 3>, double>
    compute_physical_derivatives_tri(const std::array<std::array<double,2>, 8>& coords) {
        auto dN_xi_eta = get_shape_function_gradients_tri();
        double dx_dxi = 0, dy_dxi = 0, dx_deta = 0, dy_deta = 0;
        for (int i = 0; i < 3; ++i) {
            dx_dxi  += dN_xi_eta[0][i] * coords[i][0];
            dy_dxi  += dN_xi_eta[0][i] * coords[i][1];
            dx_deta += dN_xi_eta[1][i] * coords[i][0];
            dy_deta += dN_xi_eta[1][i] * coords[i][1];
        }
        double detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
        std::array<std::array<double,2>, 3> dN_xy;
        for (int i = 0; i < 3; ++i) {
            dN_xy[i][0] = ( dy_deta * dN_xi_eta[0][i] - dy_dxi * dN_xi_eta[1][i]) / detJ;
            dN_xy[i][1] = (-dx_deta * dN_xi_eta[0][i] + dx_dxi * dN_xi_eta[1][i]) / detJ;
        }
        return {dN_xy, detJ};
    }

    template<size_t N>
    static std::pair<std::array<std::array<double,N>,2>, double>
    compute_physical_derivatives(const std::array<std::array<double, 2>, 8>& coords, 
                                 const std::array<std::array<double,N>,2>& dN_xi_eta) {
        double dx_dxi = 0, dy_dxi = 0, dx_deta = 0, dy_deta = 0;
        for (size_t i = 0; i < N; ++i) {
            dx_dxi  += dN_xi_eta[0][i] * coords[i][0];
            dy_dxi  += dN_xi_eta[0][i] * coords[i][1];
            dx_deta += dN_xi_eta[1][i] * coords[i][0];
            dy_deta += dN_xi_eta[1][i] * coords[i][1];
        }
        double detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
        std::array<std::array<double,N>,2> dN_xy;
        for (size_t i = 0; i < N; ++i) {
            dN_xy[0][i] = ( dy_deta * dN_xi_eta[0][i] - dy_dxi * dN_xi_eta[1][i]) / detJ;
            dN_xy[1][i] = (-dx_deta * dN_xi_eta[0][i] + dx_dxi * dN_xi_eta[1][i]) / detJ;
        }
        return {dN_xy, detJ};
    }
};