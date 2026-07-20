#pragma once
#include "Physics/MaterialModel.h"
#include "IO/ConfigTypes.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

class FerroelectricMaterial : public MaterialModel {
private:
    const MaterialConfig& material;

public:
    explicit FerroelectricMaterial(const MaterialConfig& config) 
        : material(config) {}

    double get_eta_k() const override { 
        return material.eta_k; 
    }
    
    double get_eps0() const override { 
        return material.eps0; 
    }

    // ==========================================
    // 0. CALCUL DES ÉNERGIES
    // ==========================================
    
    double U_energy(const Eigen::Matrix2d& grad_P) const override {
        double grad_P_norm_sq = grad_P(0,0)*grad_P(0,0) + grad_P(0,1)*grad_P(0,1) +
                                grad_P(1,0)*grad_P(1,0) + grad_P(1,1)*grad_P(1,1);
        return 0.5 * material.a0 * grad_P_norm_sq;
    }

    double W_energy(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain) const override {
        double p1 = P(0), p2 = P(1);
        double p1_2 = p1 * p1, p2_2 = p2 * p2;
        
        // Énergie électrostrictive (Couplage)
        double W_mec = -0.5 * (material.b1 * strain(0,0) * p1_2 +
                               material.b2 * strain(1,1) * p1_2 +
                               material.b1 * strain(1,1) * p2_2 +
                               material.b2 * strain(0,0) * p2_2 +
                               4.0 * material.b3 * strain(0,1) * p1 * p2);

        // Énergie élastique pure (Anisotrope)
        double c1 = material.c1, c2 = material.c2, c3 = material.c3;
        double W_elas = 0.5 * c1 * (strain(0,0)*strain(0,0) + strain(1,1)*strain(1,1))
                      + c2 * strain(0,0)*strain(1,1)
                      + 0.5 * c3 * (strain(0,1)*strain(0,1) + strain(1,0)*strain(1,0));

        return W_mec + W_elas;
    }

    double chi_energy(const Eigen::Vector2d& P) const override {
        double p1 = P(0), p2 = P(1);
        double p1_2 = p1 * p1, p1_4 = p1_2 * p1_2, p1_6 = p1_4 * p1_2, p1_8 = p1_4 * p1_4;
        double p2_2 = p2 * p2, p2_4 = p2_2 * p2_2, p2_6 = p2_4 * p2_2, p2_8 = p2_4 * p2_4;
        
        return material.alpha_1 * (p1_2 + p2_2) +
               material.alpha_11 * (p1_4 + p2_4) +
               material.alpha_12 * (p1_2 * p2_2) +
               material.alpha_111 * (p1_6 + p2_6) +
               material.alpha_112 * (p1_2 * p2_4 + p2_2 * p1_4) +
               material.alpha_1111 * (p1_8 + p2_8) +
               material.alpha_1112 * (p1_6 * p2_2 + p2_6 * p1_2) +
               material.alpha_1122 * (p1_4 * p2_4);
    }

    // ==========================================
    // 1. MÉCANIQUE
    // ==========================================
    
    Eigen::Matrix3d get_elastic_matrix() const override {
        Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
        C(0,0) = material.c1; 
        C(1,1) = material.c1; 
        C(0,1) = material.c2; 
        C(1,0) = material.c2;
        C(2,2) = material.c3 / 2.0;
        return C;
    }

    Eigen::Vector3d compute_sigma_0(const Eigen::Vector2d& P) const override {
        Eigen::Vector3d sigma_0;
        double p1_2 = P(0) * P(0);
        double p2_2 = P(1) * P(1);
        
        sigma_0(0) = -0.5 * (material.b1 * p1_2 + material.b2 * p2_2);
        sigma_0(1) = -0.5 * (material.b1 * p2_2 + material.b2 * p1_2);
        sigma_0(2) = -material.b3 * P(0) * P(1);
        
        return sigma_0;
    }

    // ==========================================
    // 2. ÉLECTROSTATIQUE
    // ==========================================
    
    double compute_effective_permittivity(double v, double eta_k_in, bool is_impermeable) const override {
        (void)eta_k_in; 
        double phase_factor = is_impermeable ? (v * v + material.eta_k) : 1.0;
        return material.eps0 * phase_factor;
    }

    Eigen::Vector2d compute_effective_polarization(const Eigen::Vector2d& P, double v, double eta_k_in, bool is_impermeable) const override {
        (void)eta_k_in;
        double phase_factor = is_impermeable ? (v * v + material.eta_k) : 1.0;
        return P * phase_factor;
    }

    // ==========================================
    // 3. FRACTURE
    // ==========================================
    
    double compute_H_drive(const Eigen::Matrix2d& grad_P, const Eigen::Vector2d& P,
                           const Eigen::Matrix2d& strain, const Eigen::Vector2d& E,
                           bool is_impermeable) const override {
        double U = U_energy(grad_P);
        double W_mec = W_energy(P, strain);
        double W_elec = 0.0;
        
        if (is_impermeable) {
            W_elec = - (P(0) * E(0) + P(1) * E(1)) - 0.5 * material.eps0 * (E(0)*E(0) + E(1)*E(1));
        }
        return std::max(0.0, U + W_mec + W_elec);
    }

    // ==========================================
    // 4. POLARISATION (Ginzburg-Landau)
    // ==========================================
    
    GinzburgLandauTerms compute_GL_terms(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain,
                                         const Eigen::Vector2d& E, double penalite_fracture,
                                         bool is_impermeable) const override {
        (void)is_impermeable;
        GinzburgLandauTerms GL;
        
        double a1 = material.alpha_1, a11 = material.alpha_11, a12 = material.alpha_12;
        double a111 = material.alpha_111, a112 = material.alpha_112;
        double a1111 = material.alpha_1111, a1112 = material.alpha_1112, a1122 = material.alpha_1122;
        
        double p1 = P(0), p2 = P(1);
        double p1_2 = p1 * p1, p1_3 = p1_2 * p1, p1_4 = p1_2 * p1_2, p1_5 = p1_4 * p1, p1_6 = p1_4 * p1_2, p1_7 = p1_6 * p1;
        double p2_2 = p2 * p2, p2_3 = p2_2 * p2, p2_4 = p2_2 * p2_2, p2_5 = p2_4 * p2, p2_6 = p2_4 * p2_2, p2_7 = p2_6 * p2;
        
        double dchi_dp1 = 2.0*a1*p1 + 4.0*a11*p1_3 + 2.0*a12*p1*p2_2 + 6.0*a111*p1_5 + 4.0*a112*p1_3*p2_2 + 2.0*a112*p1*p2_4 + 8.0*a1111*p1_7 + 6.0*a1112*p1_5*p2_2 + 2.0*a1112*p1*p2_6 + 4.0*a1122*p1_3*p2_4;
        double dchi_dp2 = 2.0*a1*p2 + 4.0*a11*p2_3 + 2.0*a12*p2*p1_2 + 6.0*a111*p2_5 + 4.0*a112*p2_3*p1_2 + 2.0*a112*p2*p1_4 + 8.0*a1111*p2_7 + 6.0*a1112*p2_5*p1_2 + 2.0*a1112*p2*p1_6 + 4.0*a1122*p2_3*p1_4;
        
        double eps_11 = strain(0,0), eps_22 = strain(1,1), eps_12 = 0.5 * (strain(0,1) + strain(1,0));
        
        double dW_dp1 = - (material.b1 * eps_11 + material.b2 * eps_22) * p1 - 2.0 * material.b3 * eps_12 * p2;
        double dW_dp2 = - (material.b1 * eps_22 + material.b2 * eps_11) * p2 - 2.0 * material.b3 * eps_12 * p1;
        
        GL.force_px = penalite_fracture * dW_dp1 + dchi_dp1 - E(0);
        GL.force_py = penalite_fracture * dW_dp2 + dchi_dp2 - E(1);
        
        double d2chi_dp12 = 2.0*a1 + 12.0*a11*p1_2 + 2.0*a12*p2_2 + 30.0*a111*p1_4 + 12.0*a112*p1_2*p2_2 + 2.0*a112*p2_4 + 56.0*a1111*p1_6 + 30.0*a1112*p1_4*p2_2 + 2.0*a1112*p2_6 + 12.0*a1122*p1_2*p2_4;
        double d2chi_dp22 = 2.0*a1 + 12.0*a11*p2_2 + 2.0*a12*p1_2 + 30.0*a111*p2_4 + 12.0*a112*p2_2*p1_2 + 2.0*a112*p1_4 + 56.0*a1111*p2_6 + 30.0*a1112*p2_4*p1_2 + 2.0*a1112*p1_6 + 12.0*a1122*p2_2*p1_4;
        double d2chi_dp1dp2 = 4.0*a12*p1*p2 + 8.0*a112*p1*p2*(p1_2 + p2_2) + 12.0*a1112*p1*p2*(p1_4 + p2_4) + 16.0*a1122*p1_2*p1*p2_2*p2;
        
        double d2W_dp12 = -material.b1 * eps_11 - material.b2 * eps_22;
        double d2W_dp22 = -material.b1 * eps_22 - material.b2 * eps_11;
        double d2W_dp1dp2 = -2.0 * material.b3 * eps_12;
        
        GL.J_11 = penalite_fracture * d2W_dp12 + d2chi_dp12;
        GL.J_22 = penalite_fracture * d2W_dp22 + d2chi_dp22;
        GL.J_12 = penalite_fracture * d2W_dp1dp2 + d2chi_dp1dp2;
        
        return GL;
    }
};