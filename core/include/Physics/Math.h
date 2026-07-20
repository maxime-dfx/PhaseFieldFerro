#ifndef PHYSICS_MATH_H
#define PHYSICS_MATH_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include "Physics/MaterialModel.h"

class Math {
private:
    const MaterialModel& material;

public:
    // Variables publiques nécessaires pour la classe Diagnostics
    const double eta_k;
    const double eps0;

    explicit Math(const MaterialModel& material_in)
        : material(material_in), eta_k(material_in.get_eta_k()), eps0(material_in.get_eps0()) {}

    // ==========================================
    // 0. CALCUL DES ÉNERGIES (Pour Diagnostics et Fracture)
    // ==========================================
    
    inline double U_energy(const Eigen::Matrix2d& grad_P) const {
        return material.U_energy(grad_P);
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
        return material.W_energy(P, strain);
    }

    inline double chi_energy(const Eigen::Vector2d& P) const {
        return material.chi_energy(P);
    }

    // ==========================================
    // 1. MÉCANIQUE : K_u * du = F_piezo + F_ext
    // ==========================================
    
    inline Eigen::Matrix3d get_elastic_matrix() const {
        return material.get_elastic_matrix();
    }

    // sigma_0 = dW_mec/d(eps), exprime dans la base de Voigt utilisee par B
    // (B utilise le cisaillement "ingenieur" gamma_12 = 2*eps_12_tensoriel,
    //  donc la composante conjuguee sigma_0(2) doit valoir dW/dgamma_12
    //  = 0.5 * dW/deps_12_tensoriel).
    // Avec le W_energy corrige : dW_mec/deps_12_tensoriel = -2*b3*p1*p2
    //  => sigma_0(2) = 0.5 * (-2*b3*p1*p2) = -b3*p1*p2.
    inline Eigen::Vector3d compute_sigma_0(const Eigen::Vector2d& P) const {
        return material.compute_sigma_0(P);
    }

    // ==========================================
    // 2. ÉLECTROSTATIQUE : K_phi * d_phi = F_pol + F_w
    // ==========================================
    
    // (void)eta_k_in; permet de supprimer le warning du compilateur tout en gardant 
    // la compatibilité avec l'appel de la fonction dans l'assembleur.
    inline double compute_effective_permittivity(double v, double eta_k_in, bool is_impermeable) const {
        return material.compute_effective_permittivity(v, eta_k_in, is_impermeable);
    }

    inline Eigen::Vector2d compute_effective_polarization(const Eigen::Vector2d& P, double v, double eta_k_in, bool is_impermeable) const {
        return material.compute_effective_polarization(P, v, eta_k_in, is_impermeable);
    }

    // ==========================================
    // 3. FRACTURE : K_v * d_v = F_drive
    // ==========================================
    
    inline double compute_H_drive(const Eigen::Matrix2d& grad_P, const Eigen::Vector2d& P, 
                                  const Eigen::Matrix2d& strain, const Eigen::Vector2d& E, 
                                  bool is_impermeable) const {
        return material.compute_H_drive(grad_P, P, strain, E, is_impermeable);
    }

    // ==========================================
    // 4. POLARISATION : K_P * d_P = F_GL
    // ==========================================
    
    inline GinzburgLandauTerms compute_GL_terms(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain, 
                                                const Eigen::Vector2d& E, double penalite_fracture, 
                                                bool is_impermeable) const {
        return material.compute_GL_terms(P, strain, E, penalite_fracture, is_impermeable);
    }
};

#endif // PHYSICS_MATH_H