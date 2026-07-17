#ifndef PHYSICS_MATH_H
#define PHYSICS_MATH_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include "IO/Datafile.h"

// Structure pour regrouper les termes de l'équation de Ginzburg-Landau
struct GinzburgLandauTerms {
    double J_11;
    double J_12;
    double J_22;
    double force_px;
    double force_py;
};

class Math {
private:
    const Datafile& config;

public:
    // Variables publiques nécessaires pour la classe Diagnostics
    const double eta_k;
    const double eps0;

    explicit Math(const Datafile& config) 
        : config(config), eta_k(config.material.eta_k), eps0(config.material.eps0) {}

    // ==========================================
    // 0. CALCUL DES ÉNERGIES (Pour Diagnostics et Fracture)
    // ==========================================
    
    inline double U_energy(const Eigen::Matrix2d& grad_P) const {
        double grad_P_norm_sq = grad_P(0,0)*grad_P(0,0) + grad_P(0,1)*grad_P(0,1) + 
                                grad_P(1,0)*grad_P(1,0) + grad_P(1,1)*grad_P(1,1);
        return 0.5 * config.material.a0 * grad_P_norm_sq;
    }

    // W(p, eps) = -b1/2(eps11 p1^2 + eps22 p2^2) - b2/2(eps22 p1^2 + eps11 p2^2)
    //             - b3(eps21 + eps12) p1 p2 + c1/2(eps11^2+eps22^2) + c2 eps11 eps22
    //             + c3/2(eps12^2 + eps21^2)                              [Eq. (4), papier 2011]
    // strain(0,1) et strain(1,0) sont ici des composantes tensorielles egales
    // (deja symetrisees dans Mechanics.cpp), donc (eps21+eps12) = 2*strain(0,1).
    // => le prefacteur correct sur le terme croise est b3 * 2 * strain(0,1),
    //    et comme il est deja multiplie par -0.5 exterieur, il faut un facteur
    //    interne de 4.0 (et non 2.0) pour obtenir au final -b3*(eps21+eps12)*p1*p2.
    inline double W_energy(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain) const {
        double p1 = P(0), p2 = P(1);
        double p1_2 = p1 * p1, p2_2 = p2 * p2;
        
        // Énergie électrostrictive (Couplage)
        double W_mec = -0.5 * (config.material.b1 * strain(0,0) * p1_2 + 
                               config.material.b2 * strain(1,1) * p1_2 + 
                               config.material.b1 * strain(1,1) * p2_2 + 
                               config.material.b2 * strain(0,0) * p2_2 + 
                               4.0 * config.material.b3 * strain(0,1) * p1 * p2);   // <-- 2.0 -> 4.0

        // Énergie élastique pure (Anisotrope)
        double c1 = config.material.c1, c2 = config.material.c2, c3 = config.material.c3;
        double W_elas = 0.5 * c1 * (strain(0,0)*strain(0,0) + strain(1,1)*strain(1,1))
                      + c2 * strain(0,0)*strain(1,1)
                      + 0.5 * c3 * (strain(0,1)*strain(0,1) + strain(1,0)*strain(1,0));

        return W_mec + W_elas;
    }

    inline double chi_energy(const Eigen::Vector2d& P) const {
        double p1 = P(0), p2 = P(1);
        double p1_2 = p1 * p1, p1_4 = p1_2 * p1_2, p1_6 = p1_4 * p1_2, p1_8 = p1_4 * p1_4;
        double p2_2 = p2 * p2, p2_4 = p2_2 * p2_2, p2_6 = p2_4 * p2_2, p2_8 = p2_4 * p2_4;

        return config.material.alpha_1 * (p1_2 + p2_2) +
               config.material.alpha_11 * (p1_4 + p2_4) +
               config.material.alpha_12 * (p1_2 * p2_2) +
               config.material.alpha_111 * (p1_6 + p2_6) +
               config.material.alpha_112 * (p1_2 * p2_4 + p2_2 * p1_4) +
               config.material.alpha_1111 * (p1_8 + p2_8) +
               config.material.alpha_1112 * (p1_6 * p2_2 + p2_6 * p1_2) +
               config.material.alpha_1122 * (p1_4 * p2_4);
    }

    // ==========================================
    // 1. MÉCANIQUE : K_u * du = F_piezo + F_ext
    // ==========================================
    
    inline Eigen::Matrix3d get_elastic_matrix() const {
        Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
        C(0,0) = config.material.c1; 
        C(1,1) = config.material.c1; 
        C(0,1) = config.material.c2; 
        C(1,0) = config.material.c2;
        C(2,2) = config.material.c3 / 2.0;
        return C;
    }

    // sigma_0 = dW_mec/d(eps), exprime dans la base de Voigt utilisee par B
    // (B utilise le cisaillement "ingenieur" gamma_12 = 2*eps_12_tensoriel,
    //  donc la composante conjuguee sigma_0(2) doit valoir dW/dgamma_12
    //  = 0.5 * dW/deps_12_tensoriel).
    // Avec le W_energy corrige : dW_mec/deps_12_tensoriel = -2*b3*p1*p2
    //  => sigma_0(2) = 0.5 * (-2*b3*p1*p2) = -b3*p1*p2.
    inline Eigen::Vector3d compute_sigma_0(const Eigen::Vector2d& P) const {
        Eigen::Vector3d sigma_0;
        double p1_2 = P(0) * P(0);
        double p2_2 = P(1) * P(1);

        sigma_0(0) = -0.5 * (config.material.b1 * p1_2 + config.material.b2 * p2_2);
        sigma_0(1) = -0.5 * (config.material.b1 * p2_2 + config.material.b2 * p1_2);
        sigma_0(2) = -config.material.b3 * P(0) * P(1);   // <-- -0.5*b3 -> -1.0*b3

        return sigma_0;
    }

    // ==========================================
    // 2. ÉLECTROSTATIQUE : K_phi * d_phi = F_pol + F_w
    // ==========================================
    
    // (void)eta_k_in; permet de supprimer le warning du compilateur tout en gardant 
    // la compatibilité avec l'appel de la fonction dans l'assembleur.
    inline double compute_effective_permittivity(double v, double eta_k_in, bool is_impermeable) const {
        (void)eta_k_in; 
        double phase_factor = is_impermeable ? (v * v + eta_k) : 1.0;
        return eps0 * phase_factor;
    }

    inline Eigen::Vector2d compute_effective_polarization(const Eigen::Vector2d& P, double v, double eta_k_in, bool is_impermeable) const {
        (void)eta_k_in;
        double phase_factor = is_impermeable ? (v * v + eta_k) : 1.0;
        return P * phase_factor;
    }

    // ==========================================
    // 3. FRACTURE : K_v * d_v = F_drive
    // ==========================================
    
    inline double compute_H_drive(const Eigen::Matrix2d& grad_P, const Eigen::Vector2d& P, 
                                  const Eigen::Matrix2d& strain, const Eigen::Vector2d& E, 
                                  bool is_impermeable) const {
        double U = U_energy(grad_P);
        double W_mec = W_energy(P, strain);
        
        double W_elec = 0.0;
        // Seule la fissure imperméable pénalise l'énergie électrique totale
        if (is_impermeable) {
            W_elec = - (P(0) * E(0) + P(1) * E(1)) - 0.5 * eps0 * (E(0)*E(0) + E(1)*E(1));
        }

        return std::max(0.0, U + W_mec + W_elec);
    }

    // ==========================================
    // 4. POLARISATION : K_P * d_P = F_GL
    // ==========================================
    
    inline GinzburgLandauTerms compute_GL_terms(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain, 
                                                const Eigen::Vector2d& E, double penalite_fracture, 
                                                bool is_impermeable) const {
        (void)is_impermeable;
        GinzburgLandauTerms GL;
        
        double a1 = config.material.alpha_1, a11 = config.material.alpha_11, a12 = config.material.alpha_12;
        double a111 = config.material.alpha_111, a112 = config.material.alpha_112;
        double a1111 = config.material.alpha_1111, a1112 = config.material.alpha_1112, a1122 = config.material.alpha_1122;

        double p1 = P(0), p2 = P(1);
        double p1_2 = p1 * p1, p1_3 = p1_2 * p1, p1_4 = p1_2 * p1_2, p1_5 = p1_4 * p1, p1_6 = p1_4 * p1_2, p1_7 = p1_6 * p1;
        double p2_2 = p2 * p2, p2_3 = p2_2 * p2, p2_4 = p2_2 * p2_2, p2_5 = p2_4 * p2, p2_6 = p2_4 * p2_2, p2_7 = p2_6 * p2;

        double dchi_dp1 = 2.0*a1*p1 + 4.0*a11*p1_3 + 2.0*a12*p1*p2_2 + 6.0*a111*p1_5 + 4.0*a112*p1_3*p2_2 + 2.0*a112*p1*p2_4 + 8.0*a1111*p1_7 + 6.0*a1112*p1_5*p2_2 + 2.0*a1112*p1*p2_6 + 4.0*a1122*p1_3*p2_4;
        double dchi_dp2 = 2.0*a1*p2 + 4.0*a11*p2_3 + 2.0*a12*p2*p1_2 + 6.0*a111*p2_5 + 4.0*a112*p2_3*p1_2 + 2.0*a112*p2*p1_4 + 8.0*a1111*p2_7 + 6.0*a1112*p2_5*p1_2 + 2.0*a1112*p2*p1_6 + 4.0*a1122*p2_3*p1_4;

        double eps_11 = strain(0,0), eps_22 = strain(1,1), eps_12 = 0.5 * (strain(0,1) + strain(1,0));

        // dW/dp1, dW/dp2 : coefficient 2*b3 — cohérent avec W_energy() corrigée ci-dessus
        // (dW_mec/dp1 = -0.5*(2*b1*eps11*p1 + 4*b3*eps12*p2) = -(b1*eps11)*p1 - 2*b3*eps12*p2)
        double dW_dp1 = - (config.material.b1 * eps_11 + config.material.b2 * eps_22) * p1 - 2.0 * config.material.b3 * eps_12 * p2;
        double dW_dp2 = - (config.material.b1 * eps_22 + config.material.b2 * eps_11) * p2 - 2.0 * config.material.b3 * eps_12 * p1;

        GL.force_px = penalite_fracture * dW_dp1 + dchi_dp1 - E(0);
        GL.force_py = penalite_fracture * dW_dp2 + dchi_dp2 - E(1);

        double d2chi_dp12 = 2.0*a1 + 12.0*a11*p1_2 + 2.0*a12*p2_2 + 30.0*a111*p1_4 + 12.0*a112*p1_2*p2_2 + 2.0*a112*p2_4 + 56.0*a1111*p1_6 + 30.0*a1112*p1_4*p2_2 + 2.0*a1112*p2_6 + 12.0*a1122*p1_2*p2_4;
        double d2chi_dp22 = 2.0*a1 + 12.0*a11*p2_2 + 2.0*a12*p1_2 + 30.0*a111*p2_4 + 12.0*a112*p2_2*p1_2 + 2.0*a112*p1_4 + 56.0*a1111*p2_6 + 30.0*a1112*p2_4*p1_2 + 2.0*a1112*p1_6 + 12.0*a1122*p2_2*p1_4;
        double d2chi_dp1dp2 = 4.0*a12*p1*p2 + 8.0*a112*p1*p2*(p1_2 + p2_2) + 12.0*a1112*p1*p2*(p1_4 + p2_4) + 16.0*a1122*p1_2*p1*p2_2*p2;

        double d2W_dp12 = -config.material.b1 * eps_11 - config.material.b2 * eps_22;
        double d2W_dp22 = -config.material.b1 * eps_22 - config.material.b2 * eps_11;
        // coefficient 2*b3 — cohérent avec dW_dp1/dW_dp2 corrigés ci-dessus
        double d2W_dp1dp2 = -2.0 * config.material.b3 * eps_12;

        GL.J_11 = penalite_fracture * d2W_dp12 + d2chi_dp12;
        GL.J_22 = penalite_fracture * d2W_dp22 + d2chi_dp22;
        GL.J_12 = penalite_fracture * d2W_dp1dp2 + d2chi_dp1dp2;

        return GL;
    }
};

#endif // PHYSICS_MATH_H