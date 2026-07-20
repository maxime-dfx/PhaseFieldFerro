#pragma once
#include <vector>
#include <stdexcept>
#include "FEM/Quadrature.h"
#include "FEM/ShapeFunctions.h"

// ============================================================
// INTERFACE DE BASE
// ============================================================
class ReferenceElement {
public:
    virtual ~ReferenceElement() = default;
    
    virtual int get_num_nodes() const = 0;
    virtual double get_weight_multiplier() const { return 1.0; }
    virtual const std::vector<GaussPoint2D>& get_gauss_points() const = 0;
    
    // Remplissage direct des buffers locaux pour éviter les allocations (zéro overhead)
    virtual void compute_shape_functions(double xi, double eta, double* N_out) const = 0;
    virtual void compute_shape_gradients(double xi, double eta, double* dN_xi_out, double* dN_eta_out) const = 0;
};

// ============================================================
// IMPLÉMENTATIONS SPÉCIFIQUES
// ============================================================

class ElementT3 : public ReferenceElement {
public:
    int get_num_nodes() const override { return 3; }
    double get_weight_multiplier() const override { return 0.5; }
    
    const std::vector<GaussPoint2D>& get_gauss_points() const override {
        static const std::vector<GaussPoint2D> gp = {{1.0/3.0, 1.0/3.0, 1.0}};
        return gp;
    }

    void compute_shape_functions(double xi, double eta, double* N_out) const override {
        auto N = ShapeFunctions::get_shape_functions_tri(xi, eta);
        for(int i = 0; i < 3; ++i) N_out[i] = N[i];
    }

    void compute_shape_gradients(double /*xi*/, double /*eta*/, double* dN_xi_out, double* dN_eta_out) const override {
        auto dN = ShapeFunctions::get_shape_function_gradients_tri(); // Ne prend pas de paramètres
        for(int i = 0; i < 3; ++i) {
            dN_xi_out[i] = dN[0][i];
            dN_eta_out[i] = dN[1][i];
        }
    }
};

class ElementQ4 : public ReferenceElement {
public:
    int get_num_nodes() const override { return 4; }
    
    const std::vector<GaussPoint2D>& get_gauss_points() const override {
        static const std::vector<GaussPoint2D> gp = Quadrature::get_gauss_2x2();
        return gp;
    }

    void compute_shape_functions(double xi, double eta, double* N_out) const override {
        auto N = ShapeFunctions::get_shape_functions(xi, eta);
        for(int i = 0; i < 4; ++i) N_out[i] = N[i];
    }

    void compute_shape_gradients(double xi, double eta, double* dN_xi_out, double* dN_eta_out) const override {
        auto dN = ShapeFunctions::get_shape_function_gradients(xi, eta);
        for(int i = 0; i < 4; ++i) {
            dN_xi_out[i] = dN[0][i];
            dN_eta_out[i] = dN[1][i];
        }
    }
};

class ElementT6 : public ReferenceElement {
public:
    int get_num_nodes() const override { return 6; }
    double get_weight_multiplier() const override { return 0.5; }
    
    const std::vector<GaussPoint2D>& get_gauss_points() const override {
        static const std::vector<GaussPoint2D> gp = Quadrature::get_gauss_3_points_tri();
        return gp;
    }

    void compute_shape_functions(double xi, double eta, double* N_out) const override {
        auto N = ShapeFunctions::get_shape_functions_t6(xi, eta);
        for(int i = 0; i < 6; ++i) N_out[i] = N[i];
    }

    void compute_shape_gradients(double xi, double eta, double* dN_xi_out, double* dN_eta_out) const override {
        auto dN = ShapeFunctions::get_shape_function_gradients_t6(xi, eta);
        for(int i = 0; i < 6; ++i) {
            dN_xi_out[i] = dN[0][i];
            dN_eta_out[i] = dN[1][i];
        }
    }
};

class ElementQ8 : public ReferenceElement {
public:
    int get_num_nodes() const override { return 8; }
    
    const std::vector<GaussPoint2D>& get_gauss_points() const override {
        static const std::vector<GaussPoint2D> gp = Quadrature::get_gauss_3x3();
        return gp;
    }

    void compute_shape_functions(double xi, double eta, double* N_out) const override {
        auto N = ShapeFunctions::get_shape_functions_q8(xi, eta);
        for(int i = 0; i < 8; ++i) N_out[i] = N[i];
    }

    void compute_shape_gradients(double xi, double eta, double* dN_xi_out, double* dN_eta_out) const override {
        auto dN = ShapeFunctions::get_shape_function_gradients_q8(xi, eta);
        for(int i = 0; i < 8; ++i) {
            dN_xi_out[i] = dN[0][i];
            dN_eta_out[i] = dN[1][i];
        }
    }
};

// ============================================================
// FACTORY (Distributeur d'éléments statiques)
// ============================================================
class ElementFactory {
public:
    static const ReferenceElement& get_element(int n_nodes) {
        // Initialisation paresseuse statique (Thread-safe depuis C++11, 0 allocation après le 1er appel)
        static const ElementT3 t3;
        static const ElementQ4 q4;
        static const ElementT6 t6;
        static const ElementQ8 q8;

        switch(n_nodes) {
            case 3: return t3;
            case 4: return q4;
            case 6: return t6;
            case 8: return q8;
            default: throw std::runtime_error("Nombre de noeuds non supporte pour l'integration : " + std::to_string(n_nodes));
        }
    }
};