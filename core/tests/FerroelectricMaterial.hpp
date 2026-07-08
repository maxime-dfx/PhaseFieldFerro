#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>

using namespace Eigen;

class FerroelectricMaterial {
private:
    double Gc_, kappa_, eta_, eps0_, a0_;
    double mu_p_, mu_v_; // <-- AJOUT ICI

    double a1_, a11_, a12_, a111_, a112_, a1111_, a1112_, a1122_;
    Matrix3d C_;
    double q11_, q12_, q44_;

public:
    FerroelectricMaterial(double Gc, double kappa, double eta, double eps0, double a0, double mu_p, double mu_v)
        : Gc_(Gc), kappa_(kappa), eta_(eta), eps0_(eps0), a0_(a0), mu_p_(mu_p), mu_v_(mu_v) {
        
        a1_ = -0.0023; a11_ = -0.0029; a12_ = -0.0011; 
        C_ << 185, 111, 0, 111, 185, 0, 0, 0, 74; 
        q11_ = 0.05; q12_ = -0.02; q44_ = 0.04; 
    }

    // ========================================================================
    // 1. GESTION DE LA RUPTURE (Couplage Phase-Field)
    // ========================================================================

    // Fonction de dégradation : g(v) = v^2 + eta
    inline double degradation(double v) const {
        return (v * v) + eta_;
    }

    // Dérivée de la fonction de dégradation : g'(v) = 2v
    inline double d_degradation_dv(double v) const {
        return 2.0 * v;
    }

    // ========================================================================
    // 2. COMPORTEMENT ÉLECTRO-MÉCANIQUE
    // ========================================================================

    // Déformation spontanée eps^0 en fonction de la polarisation p (Notation de Voigt)
    // eps^0 = [q11*p1^2 + q12*p2^2,  q12*p1^2 + q11*p2^2,  q44*p1*p2]^T
    inline Vector3d spontaneous_strain(const Vector2d& p) const {
        Vector3d eps0;
        eps0(0) = q11_ * p(0)*p(0) + q12_ * p(1)*p(1);
        eps0(1) = q12_ * p(0)*p(0) + q11_ * p(1)*p(1);
        eps0(2) = q44_ * p(0) * p(1);
        return eps0;
    }

    // Matrice de dérivée d_eps^0 / dp (Taille 3x2)
    // Utile pour la force mécanique qui pousse les domaines à basculer
    inline Matrix<double, 3, 2> d_spontaneous_strain_dp(const Vector2d& p) const {
        Matrix<double, 3, 2> deps0_dp;
        deps0_dp(0, 0) = 2.0 * q11_ * p(0);  deps0_dp(0, 1) = 2.0 * q12_ * p(1);
        deps0_dp(1, 0) = 2.0 * q12_ * p(0);  deps0_dp(1, 1) = 2.0 * q11_ * p(1);
        deps0_dp(2, 0) = q44_ * p(1);        deps0_dp(2, 1) = q44_ * p(0);
        return deps0_dp;
    }

    // Calcul des contraintes effectives : sigma = C * (eps - eps^0)
    inline Vector3d stress(const Vector2d& p, const Vector3d& eps) const {
        return C_ * (eps - spontaneous_strain(p));
    }

    // Énergie élastique : W = 1/2 * (eps - eps^0)^T * C * (eps - eps^0)
    inline double elastic_energy(const Vector2d& p, const Vector3d& eps) const {
        Vector3d elastic_strain = eps - spontaneous_strain(p);
        return 0.5 * elastic_strain.dot(C_ * elastic_strain);
    }

    // ========================================================================
    // 3. COMPORTEMENT FERROÉLECTRIQUE (Landau-Devonshire)
    // ========================================================================

    // Dérivée du potentiel de Landau d(chi)/dp (Vecteur 2x1)
    // C'est la force chimique qui crée les puits de potentiel
    inline Vector2d landau_derivative(const Vector2d& p) const {
        double p1 = p(0); double p2 = p(1);
        double p1_2 = p1*p1; double p2_2 = p2*p2;
        
        Vector2d dchi;
        // Dérivée par rapport à p1 (Polynôme d'ordre 8 du papier d'Arias)
        dchi(0) = 2.0*a1_*p1 + 4.0*a11_*p1*p1_2 + 2.0*a12_*p1*p2_2 /* + Ajouter les termes d'ordre 6 et 8 ici */;
        
        // Dérivée par rapport à p2
        dchi(1) = 2.0*a1_*p2 + 4.0*a11_*p2*p2_2 + 2.0*a12_*p2*p1_2 /* + Ajouter les termes d'ordre 6 et 8 ici */;
        
        return dchi;
    }

    // ========================================================================
    // 4. LES FORCES MOTEURS DU SYSTÈME (Seconds membres des schémas)
    // ========================================================================

    // Force qui pousse la polarisation à changer sous l'effet des contraintes
    // dW/dp = - (d_eps0/dp)^T * sigma
    inline Vector2d mechanical_driving_force_on_p(const Vector2d& p, const Vector3d& eps) const {
        Vector3d sig = stress(p, eps);
        Matrix<double, 3, 2> deps0_dp = d_spontaneous_strain_dp(p);
        return -deps0_dp.transpose() * sig; // Renvoie un vecteur 2x1
    }

    // Force qui fait avancer la fissure : l'énergie accumulée
    // Force = U(grad_p) + W(p, eps)
    inline double crack_driving_force(const Vector2d& p, const Vector3d& eps, const Matrix2d& grad_p) const {
        double U_wall = 0.5 * a0_ * (grad_p(0,0)*grad_p(0,0) + grad_p(0,1)*grad_p(0,1) + 
                                     grad_p(1,0)*grad_p(1,0) + grad_p(1,1)*grad_p(1,1));
        double W_elast = elastic_energy(p, eps);
        return U_wall + W_elast;
    }

    // Getters pour les constantes d'intégration
    inline double Gc() const { return Gc_; }
    inline double kappa() const { return kappa_; }
    inline double eps0() const { return eps0_; }
    inline double a0() const { return a0_; }
    inline double mu_p() const { return mu_p_; }
    inline double mu_v() const { return mu_v_; }
    inline Matrix3d C() const { return C_; }
};