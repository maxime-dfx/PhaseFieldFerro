#pragma once
#include <array>
#include <vector>
#include <utility>

class ShapeFunctions {
public:
    // ==========================================
    // 1. TRIANGLE LINÉAIRE (T3 / CST)
    // ==========================================
    static std::vector<double> get_shape_functions_tri(double xi, double eta) {
        return { 1.0 - xi - eta, xi, eta };
    }

    static std::array<std::vector<double>, 2> get_shape_function_gradients_tri() {
        return {{
            { -1.0, 1.0, 0.0 },   // dN/dxi
            { -1.0, 0.0, 1.0 }    // dN/deta
        }};
    }

    // ==========================================
    // 2. QUADRANGLE BILINÉAIRE (Q4)
    // ==========================================
    static std::array<double, 4> get_shape_functions(double xi, double eta) {
        return {
            0.25 * (1.0 - xi) * (1.0 - eta),
            0.25 * (1.0 + xi) * (1.0 - eta),
            0.25 * (1.0 + xi) * (1.0 + eta),
            0.25 * (1.0 - xi) * (1.0 + eta)
        };
    }

    static std::array<std::array<double, 4>, 2> get_shape_function_gradients(double xi, double eta) {
        return {{
            { -0.25*(1.0 - eta),  0.25*(1.0 - eta),  0.25*(1.0 + eta), -0.25*(1.0 + eta) },
            { -0.25*(1.0 - xi),  -0.25*(1.0 + xi),   0.25*(1.0 + xi),   0.25*(1.0 - xi) }
        }};
    }

    // ==========================================
    // 3. TRIANGLE QUADRATIQUE (T6)
    // ==========================================
    static std::vector<double> get_shape_functions_t6(double xi, double eta) {
        double L1 = 1.0 - xi - eta;
        double L2 = xi;
        double L3 = eta;
        return {
            L1 * (2.0 * L1 - 1.0), 
            L2 * (2.0 * L2 - 1.0), 
            L3 * (2.0 * L3 - 1.0),
            4.0 * L1 * L2, 
            4.0 * L2 * L3, 
            4.0 * L3 * L1
        };
    }

    static std::array<std::vector<double>, 2> get_shape_function_gradients_t6(double xi, double eta) {
        double L1 = 1.0 - xi - eta;
        return {{
            // dN/dxi
            {1.0 - 4.0 * L1, 4.0 * xi - 1.0, 0.0, 4.0 * (L1 - xi), 4.0 * eta, -4.0 * eta},
            // dN/deta
            {1.0 - 4.0 * L1, 0.0, 4.0 * eta - 1.0, -4.0 * xi, 4.0 * xi, 4.0 * (L1 - eta)}
        }};
    }

    // ==========================================
    // 4. QUADRANGLE QUADRATIQUE (Q8 - Serendipity)
    // ==========================================
    static std::vector<double> get_shape_functions_q8(double xi, double eta) {
        double xi2 = xi * xi;
        double eta2 = eta * eta;
        return {
            0.25 * (1.0 - xi) * (1.0 - eta) * (-xi - eta - 1.0),
            0.25 * (1.0 + xi) * (1.0 - eta) * (xi - eta - 1.0),
            0.25 * (1.0 + xi) * (1.0 + eta) * (xi + eta - 1.0),
            0.25 * (1.0 - xi) * (1.0 + eta) * (-xi + eta - 1.0),
            0.5 * (1.0 - xi2) * (1.0 - eta),
            0.5 * (1.0 + xi) * (1.0 - eta2),
            0.5 * (1.0 - xi2) * (1.0 + eta),
            0.5 * (1.0 - xi) * (1.0 - eta2)
        };
    }

    static std::array<std::vector<double>, 2> get_shape_function_gradients_q8(double xi, double eta) {
        double xi2 = xi * xi;
        double eta2 = eta * eta;
        return {{
            // dN/dxi
            {
                0.25 * (1.0 - eta) * (2.0 * xi + eta),
                0.25 * (1.0 - eta) * (2.0 * xi - eta),
                0.25 * (1.0 + eta) * (2.0 * xi + eta),
                0.25 * (1.0 + eta) * (2.0 * xi - eta),
                -xi * (1.0 - eta),
                0.5 * (1.0 - eta2),
                -xi * (1.0 + eta),
                -0.5 * (1.0 - eta2)
            },
            // dN/deta
            {
                0.25 * (1.0 - xi) * (2.0 * eta + xi),
                0.25 * (1.0 + xi) * (2.0 * eta - xi),
                0.25 * (1.0 + xi) * (2.0 * eta + xi),
                0.25 * (1.0 - xi) * (2.0 * eta - xi),
                -0.5 * (1.0 - xi2),
                -eta * (1.0 + xi),
                0.5 * (1.0 - xi2),
                -eta * (1.0 - xi)
            }
        }};
    }
};