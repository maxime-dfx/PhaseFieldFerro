#pragma once
#include <vector>
#include <stdexcept>
#include "FEM/Quadrature.h"
#include "FEM/ShapeFunctions.h"

class ElementT3 {
public:
    static constexpr int get_num_nodes() { return 3; }
    static constexpr double get_weight_multiplier() { return 0.5; }
    
    static const std::vector<GaussPoint2D>& get_gauss_points() {
        static const std::vector<GaussPoint2D> gp = {{1.0/3.0, 1.0/3.0, 1.0}};
        return gp;
    }
    static void compute_shape_functions(double xi, double eta, double* N_out) {
        auto N = ShapeFunctions::get_shape_functions_tri(xi, eta);
        for(int i = 0; i < 3; ++i) N_out[i] = N[i];
    }
    static void compute_shape_gradients(double /*xi*/, double /*eta*/, double* dN_xi_out, double* dN_eta_out) {
        auto dN = ShapeFunctions::get_shape_function_gradients_tri(); 
        for(int i = 0; i < 3; ++i) {
            dN_xi_out[i] = dN[0][i];
            dN_eta_out[i] = dN[1][i];
        }
    }
};

class ElementQ4 {
public:
    static constexpr int get_num_nodes() { return 4; }
    static constexpr double get_weight_multiplier() { return 1.0; }
    
    static const std::vector<GaussPoint2D>& get_gauss_points() {
        static const std::vector<GaussPoint2D> gp = Quadrature::get_gauss_2x2();
        return gp;
    }
    static void compute_shape_functions(double xi, double eta, double* N_out) {
        auto N = ShapeFunctions::get_shape_functions(xi, eta);
        for(int i = 0; i < 4; ++i) N_out[i] = N[i];
    }
    static void compute_shape_gradients(double xi, double eta, double* dN_xi_out, double* dN_eta_out) {
        auto dN = ShapeFunctions::get_shape_function_gradients(xi, eta);
        for(int i = 0; i < 4; ++i) {
            dN_xi_out[i] = dN[0][i];
            dN_eta_out[i] = dN[1][i];
        }
    }
};

class ElementT6 {
public:
    static constexpr int get_num_nodes() { return 6; }
    static constexpr double get_weight_multiplier() { return 0.5; }
    
    static const std::vector<GaussPoint2D>& get_gauss_points() {
        static const std::vector<GaussPoint2D> gp = Quadrature::get_gauss_3_points_tri();
        return gp;
    }
    static void compute_shape_functions(double xi, double eta, double* N_out) {
        auto N = ShapeFunctions::get_shape_functions_t6(xi, eta);
        for(int i = 0; i < 6; ++i) N_out[i] = N[i];
    }
    static void compute_shape_gradients(double xi, double eta, double* dN_xi_out, double* dN_eta_out) {
        auto dN = ShapeFunctions::get_shape_function_gradients_t6(xi, eta);
        for(int i = 0; i < 6; ++i) {
            dN_xi_out[i] = dN[0][i];
            dN_eta_out[i] = dN[1][i];
        }
    }
};

class ElementQ8 {
public:
    static constexpr int get_num_nodes() { return 8; }
    static constexpr double get_weight_multiplier() { return 1.0; }
    
    static const std::vector<GaussPoint2D>& get_gauss_points() {
        static const std::vector<GaussPoint2D> gp = Quadrature::get_gauss_3x3();
        return gp;
    }
    static void compute_shape_functions(double xi, double eta, double* N_out) {
        auto N = ShapeFunctions::get_shape_functions_q8(xi, eta);
        for(int i = 0; i < 8; ++i) N_out[i] = N[i];
    }
    static void compute_shape_gradients(double xi, double eta, double* dN_xi_out, double* dN_eta_out) {
        auto dN = ShapeFunctions::get_shape_function_gradients_q8(xi, eta);
        for(int i = 0; i < 8; ++i) {
            dN_xi_out[i] = dN[0][i];
            dN_eta_out[i] = dN[1][i];
        }
    }
};