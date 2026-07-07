#include "Physics/Math.h"
#include <cmath>

Math::Math(const Datafile& config) {
    xi = config.material.xi;
    c0 = config.material.c0;
    mu_p = config.material.mu_p;
    mu_v = config.material.mu_v;
    a0 = config.material.a0;
    P0 = config.material.P0;
    t = config.material.t;
    eta_k = config.material.eta_k;
    eps0 = config.material.eps0;
    alpha1 = config.material.alpha_1;
    alpha11 = config.material.alpha_11;
    alpha12 = config.material.alpha_12;
    alpha111 = config.material.alpha_111;
    alpha112 = config.material.alpha_112;
    alpha1111 = config.material.alpha_1111;
    alpha1112 = config.material.alpha_1112;
    alpha1122 = config.material.alpha_1122;
    b1 = config.material.b1;
    b2 = config.material.b2;
    b3 = config.material.b3;
    c1 = config.material.c1;
    c2 = config.material.c2;
    c3 = config.material.c3;
}

void Math::DimensionLess() {
    // Vide pour le moment
}

double Math::U_energy(const Eigen::Matrix2d& Pij) const {
    return 0.5 * a0 * (Pij(0,0)*Pij(0,0) + Pij(1,1)*Pij(1,1) + Pij(0,1)*Pij(0,1) + Pij(1,0)*Pij(1,0));  
}

double Math::W_energy(const Eigen::Vector2d& Pi, const Eigen::Matrix2d& epsjk) const {
    double p1 = Pi(0);
    double p2 = Pi(1);
    double p1_2 = p1*p1;
    double p2_2 = p2*p2;

    return -0.5 * b1 *(epsjk(0,0)*p1_2 + epsjk(1,1)*p2_2)
           -0.5 * b2 *(epsjk(0,0)*p2_2 + epsjk(1,1)*p1_2)
           - b3 * (epsjk(1,0) + epsjk(0,1)) * p1 * p2
           + c1 * (epsjk(0,0)*epsjk(0,0) + epsjk(1,1)*epsjk(1,1))
           + c2 * (epsjk(0,0)*epsjk(1,1))
           + c3 * (epsjk(0,1)*epsjk(0,1) + epsjk(1,0)*epsjk(1,0));
}

double Math::dW_dp1(const Eigen::Vector2d& Pi, const Eigen::Matrix2d& epsjk) const {
    return - b1 * epsjk(0,0) * Pi(0) 
           - b2 * epsjk(1,1) * Pi(0) 
           - b3 * (epsjk(1,0) + epsjk(0,1)) * Pi(1);
}

double Math::dW_dp2(const Eigen::Vector2d& Pi, const Eigen::Matrix2d& epsjk) const {
    return - b1 * epsjk(1,1) * Pi(1) 
           - b2 * epsjk(0,0) * Pi(1) 
           - b3 * (epsjk(1,0) + epsjk(0,1)) * Pi(0);
}

double Math::chi_energy(const Eigen::Vector2d& Pi) const {
    double p1 = Pi(0); double p2 = Pi(1);    
    double p1_2 = p1*p1; double p1_4 = p1_2*p1_2; double p1_6 = p1_4*p1_2; double p1_8 = p1_4*p1_4;
    double p2_2 = p2*p2; double p2_4 = p2_2*p2_2; double p2_6 = p2_4*p2_2; double p2_8 = p2_4*p2_4;

    return alpha1 * (p1_2 + p2_2)
         + alpha11 * (p1_4 + p2_4)
         + alpha12 * (p1_2*p2_2)       
         + alpha111 * (p1_6 + p2_6)
         + alpha112 * (p1_2*p2_4 + p2_2*p1_4)
         + alpha1111 * (p1_8 + p2_8)
         + alpha1112 * (p1_6*p2_2 + p2_6*p1_2)
         + alpha1122 * (p1_4 * p2_4);
}

double Math::dchi_dp1(const Eigen::Vector2d& Pi) const {
    double p1 = Pi(0); double p2 = Pi(1);    
    double p1_2 = p1*p1; double p1_3 = p1_2*p1; double p1_4 = p1_2*p1_2; double p1_5 = p1_4*p1; double p1_7 = p1_4*p1_3;
    double p2_2 = p2*p2; double p2_4 = p2_2*p2_2; double p2_6 = p2_4*p2_2;

    return 2.0 * alpha1 * p1 
         + 4.0 * alpha11 * p1_3 
         + 2.0 * alpha12 * p1 * p2_2 
         + 6.0 * alpha111 * p1_5 
         + 4.0 * alpha112 * p1_3 * p2_2 
         + 2.0 * alpha112 * p1 * p2_4 
         + 8.0 * alpha1111 * p1_7 
         + 6.0 * alpha1112 * p1_5 * p2_2 
         + 2.0 * alpha1112 * p1 * p2_6 
         + 4.0 * alpha1122 * p1_3 * p2_4;
}

double Math::dchi_dp2(const Eigen::Vector2d& Pi) const {
    double p1 = Pi(0); double p2 = Pi(1);
    double p1_2 = p1*p1; double p1_4 = p1_2*p1_2; double p1_6 = p1_4*p1_2;
    double p2_2 = p2*p2; double p2_3 = p2_2*p2; double p2_4 = p2_2*p2_2; double p2_5 = p2_4*p2; double p2_7 = p2_4*p2_3;

    return 2.0 * alpha1 * p2 
         + 4.0 * alpha11 * p2_3 
         + 2.0 * alpha12 * p1_2 * p2
         + 6.0 * alpha111 * p2_5
         + 4.0 * alpha112 * p1_2 * p2_3
         + 2.0 * alpha112 * p1_4 * p2 
         + 8.0 * alpha1111 * p2_7
         + 6.0 * alpha1112 * p1_2 * p2_5 
         + 2.0 * alpha1112 * p1_6 * p2 
         + 4.0 * alpha1122 * p1_4 * p2_3;
}

double Math::h_enthalpy_density(double U, double W, double chi, const Eigen::Vector2d& E, const Eigen::Vector2d& Pi) const {
    return U + W + chi - eps0/2.0*E.squaredNorm() - E.dot(Pi);
}

Eigen::Matrix3d Math::get_elastic_matrix() const {
    Eigen::Matrix3d C;
    C << c1, c2, 0.0,
         c2, c1, 0.0,
         0.0, 0.0, c3 / 2.0;
    return C;
}

Eigen::Vector3d Math::compute_sigma_0(const Eigen::Vector2d& Pi) const {
    double p1 = Pi(0); double p2 = Pi(1);
    Eigen::Vector3d sigma_0;
    sigma_0(0) = - (b1 / 2.0) * p1 * p1 - (b2 / 2.0) * p2 * p2;
    sigma_0(1) = - (b2 / 2.0) * p1 * p1 - (b1 / 2.0) * p2 * p2;
    sigma_0(2) = - b3 * p1 * p2;
    return sigma_0;
}

double Math::compute_H_drive(const Eigen::Matrix2d& Pij, const Eigen::Vector2d& Pi, 
                             const Eigen::Matrix2d& eps, const Eigen::Vector2d& E, 
                             bool is_impermeable) const {
    double U = U_energy(Pij);
    double W = W_energy(Pi, eps);
    double H_drive = U + W; 

    if (is_impermeable) {
        H_drive += -0.5 * eps0 * E.squaredNorm() - E.dot(Pi);
    }
    return H_drive;
}

double Math::compute_polarization_force(int component, const Eigen::Vector2d& Pi, 
                                       const Eigen::Matrix2d& eps, double E_gp, 
                                       double penalite_fracture, bool is_impermeable) const {
    double dW   = (component == 0) ? dW_dp1(Pi, eps) : dW_dp2(Pi, eps);
    double dchi = (component == 0) ? dchi_dp1(Pi)   : dchi_dp2(Pi);

    if (!is_impermeable) {
        return - (penalite_fracture * dW + dchi - E_gp);
    } else {
        return - (penalite_fracture * (dW - E_gp) + dchi);
    }
}