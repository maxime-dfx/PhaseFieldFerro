#pragma once
// FEM.h
// ------------------
// Fusion des concepts mathématiques FEM auparavant éparpillés dans
// Quadrature.h, ShapeFunctions.h, ReferenceElement.h et ElementIntegrator.h.
//
// Ordre logique du fichier :
//   1. GaussPoint2D / Quadrature
//   2. ShapeFunctions
//   3. Définitions des éléments (ElementT3, ElementQ4, ElementT6, ElementQ8)
//   4. ElementIntegrator

#include "Mesh/include/Mesh.h"
#include "Utils/include/Profiling.h"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#if defined(_MSC_VER)
    #define FORCE_INLINE inline __forceinline
#elif defined(__GNUC__) || defined(__clang__)
    #define FORCE_INLINE inline __attribute__((always_inline))
#else
    #define FORCE_INLINE inline
#endif

// ======================================================================
// 1. GaussPoint2D / Quadrature
// ======================================================================

struct GaussPoint2D { double xi, eta, weight; };

class Quadrature {
public:
    static std::vector<GaussPoint2D> get_gauss_2x2() {
        double pt = 1.0 / std::sqrt(3.0);
        return {{-pt, -pt, 1.0}, {pt, -pt, 1.0}, {pt, pt, 1.0}, {-pt, pt, 1.0}};
    }
    static std::vector<GaussPoint2D> get_gauss_3_points_tri() {
        return {{1.0/6.0, 1.0/6.0, 1.0/6.0}, {2.0/3.0, 1.0/6.0, 1.0/6.0}, {1.0/6.0, 2.0/3.0, 1.0/6.0}};
    }
    static std::vector<GaussPoint2D> get_gauss_3x3() {
        double pt = std::sqrt(3.0 / 5.0);
        double w_center = 8.0 / 9.0;
        double w_edge = 5.0 / 9.0;
        return {
            {-pt, -pt, w_edge*w_edge}, {0.0, -pt, w_center*w_edge}, {pt, -pt, w_edge*w_edge},
            {-pt, 0.0, w_edge*w_center}, {0.0, 0.0, w_center*w_center}, {pt, 0.0, w_edge*w_center},
            {-pt, pt, w_edge*w_edge}, {0.0, pt, w_center*w_edge}, {pt, pt, w_edge*w_edge}
        };
    }
};

// ======================================================================
// 2. ShapeFunctions
// ======================================================================

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
        return { 0.25*(1-xi)*(1-eta)*(-xi-eta-1), 0.25*(1+xi)*(1-eta)*(xi-eta-1), 0.25*(1+xi)*(1+eta)*(xi+eta-1),
                 0.25*(1-xi)*(1+eta)*(-xi+eta-1), 0.5*(1-x2)*(1-eta), 0.5*(1+xi)*(1-e2), 0.5*(1-x2)*(1+eta), 0.5*(1-xi)*(1-e2) };
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

// ======================================================================
// 3. Définitions des éléments de référence
// ======================================================================

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

// ======================================================================
// 4. ElementIntegrator
// ======================================================================

class ElementIntegrator {
private:
    // Implémentation templatisée qui sera inlinée à 100% par le compilateur
    template<typename RefElem, typename Func>
    FORCE_INLINE static void integrate_impl(const std::array<std::array<double, 2>, 8>& coords,
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
    template<typename Func>
    FORCE_INLINE static void integrate(const Element& elem,
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
