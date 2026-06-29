#include "Physics/Math.h"
#include "Physics/Material.h"
#include "Physics/Electrostatics.h"
#include "Physics/Polarization.h"

Math::Math(const Material& mat_in) : mat(mat_in) {}

void Math::DimensionLess() {
    xi = mat.xi;
    c0 = mat.c0;
    mu_p = mat.mu_p;
    a0 = mat.a0;
    P0 = mat.P0;
    t = mat.t;
    eps0 = mat.eps0;
    alpha1 = mat.alpha_1;
    alpha11 = mat.alpha_11;
    alpha12 = mat.alpha_12;
    alpha111 = mat.alpha_111;
    alpha112 = mat.alpha_112;
    alpha1111 = mat.alpha_1111;
    alpha1112 = mat.alpha_1112;
    alpha1122 = mat.alpha_1122;
    b1 = mat.b1;
    b2 = mat.b2;
    b3 = mat.b3;
    c1 = mat.c1;
    c2 = mat.c2;
    c3 = mat.c3;

    // xi_adim = xi*sqrt(c0/a0)/P0;
    // pi_adim = Pi/P0;
    // t_adim = t*c0/(mu*P0);
    // eps0_adim = eps0*c0/(P0*P0);
    // phi_adim = phi*sqrt(a0*c0)/(P0*P0);
    // alpha1_adim = alpha1*P0*P0/c0;
    // alpha11_adim = alpha11*P0*P0*P0*P0/c0;
    // alpha12_adim = alpha12*P0*P0*P0*P0/c0;
    // alpha111_adim = alpha111*P0*P0*P0*P0*P0*P0/c0;
    // alpha112_adim = alpha112*P0*P0*P0*P0*P0*P0/c0;
    // alpha1111_adim = alpha1111*P0*P0*P0*P0*P0*P0*P0*P0/c0;  
    // alpha1112_adim = alpha1112*P0*P0*P0*P0*P0*P0*P0*P0/c0;
    // alpha1122_adim = alpha1122*P0*P0*P0*P0*P0*P0*P0*P0/c0;
    // b1_adim = b1*P0*P0/c0;
    // b2_adim = b2*P0*P0/c0;
    // b3_adim = b3*P0*P0/c0;
    // c1_adim = c1/c0;
    // c2_adim = c2/c0;
    // c3_adim = c3/c0;
}

double Math::U_energy(const Eigen::Matrix2d& Pij) const {
    return 0.5 * mat.a0 * (Pij(0,0)*Pij(0,0) + Pij(1,1)*Pij(1,1) + Pij(0,1)*Pij(0,1) + Pij(1,0)*Pij(1,0));  
}

double Math::W_energy(const Eigen::Vector2d& Pi, const Eigen::Matrix2d& epsjk) const {
    double p1 = Pi(0);
    double p2 = Pi(1);
    double p1_2 = p1*p1;
    double p2_2 = p2*p2;

    return -0.5 * mat.b1 *(epsjk(0,0)*p1_2 + epsjk(1,1)*p2_2)
            -0.5 * mat.b2 *(epsjk(0,0)*p2_2 + epsjk(1,1)*p1_2)
            - mat.b3 * (epsjk(1,0) + epsjk(0,1))
            + mat.c1 * (epsjk(0,0)*epsjk(0,0) + epsjk(1,1)*epsjk(1,1))
            + mat.c2 * (epsjk(0,0)*epsjk(1,1));
            + mat.c3 * (epsjk(0,1)*epsjk(0,1) + epsjk(1,0)*epsjk(1,0));
}

double Math::dW_dp1(const Eigen::Vector2d& Pi, const Eigen::Matrix2d& epsjk) const {
    double p1 = Pi(0);
    double p2 = Pi(1);
    return - mat.b1 * epsjk(0,0) * p1 
           - mat.b2 * epsjk(1,1) * p1 
           - mat.b3 * (epsjk(1,0) + epsjk(0,1)) * p2;
}

double Math::dW_dp2(const Eigen::Vector2d& Pi, const Eigen::Matrix2d& epsjk) const {
    double p1 = Pi(0);
    double p2 = Pi(1);
    return - mat.b1 * epsjk(1,1) * p2 
           - mat.b2 * epsjk(0,0) * p2 
           - mat.b3 * (epsjk(1,0) + epsjk(0,1)) * p1;
}


double Math::chi_energy(const Eigen::Vector2d& Pi) const {
    double p1 = Pi(0);
    double p2 = Pi(1);    
    double p1_2 = p1*p1; double p1_3 = p1_2*p1; double p1_4 = p1_3*p1; double p1_5 = p1_4*p1; double p1_6 = p1_5*p1; double p1_7 = p1_6*p1; double p1_8 = p1_7*p1;
    double p2_2 = p2*p2; double p2_3 = p2_2*p2; double p2_4 = p2_3*p2; double p2_5 = p2_4*p2; double p2_6 = p2_5*p2; double p2_7 = p2_6*p2; double p2_8 = p2_7*p2;

    return mat.alpha_1 * (p1_2 + p2_2)
            + mat.alpha_11 * (p1_4 + p2_4)
            + mat.alpha_12 * (p1_2*p2_2)       
            + mat.alpha_111 * (p1_6 + p2_6)
            + mat.alpha_112 * (p1_2*p2_4 + p2_2*p1_4)
            + mat.alpha_1111 * (p1_8 + p2_8)
            + mat.alpha_1112 * (p1_6*p2_2 + p2_6*p1_2)
            + mat.alpha_1122 * (p1_4 * p2_4);
}

double Math::dchi_dp1(const Eigen::Vector2d& Pi) const {
    double p1 = Pi(0);
    double p2 = Pi(1);    
    double p1_2 = p1*p1; double p1_3 = p1_2*p1; double p1_4 = p1_3*p1; double p1_5 = p1_4*p1; double p1_6 = p1_5*p1; double p1_7 = p1_6*p1;
    double p2_2 = p2*p2; double p2_3 = p2_2*p2; double p2_4 = p2_3*p2; double p2_5 = p2_4*p2; double p2_6 = p2_5*p2; double p2_7 = p2_6*p2;

    return 2.0 * mat.alpha_1 * p1 
         + 4.0 * mat.alpha_11 * p1_3 
         + 2.0 * mat.alpha_12 * p1 * p2_2 
         + 6.0 * mat.alpha_111 * p1_5 
         + 4.0 * mat.alpha_112 * p1_3 * p2_2 
         + 2.0 * mat.alpha_112 * p1 * p2_4 
         + 8.0 * mat.alpha_1111 * p1_7 
         + 6.0 * mat.alpha_1112 * p1_5 * p2_2 
         + 2.0 * mat.alpha_1112 * p1 * p2_6 
         + 4.0 * mat.alpha_1122 * p1_3 * p2_4;
}

double Math::dchi_dp2(const Eigen::Vector2d& Pi) const {
    double p1 = Pi(0);
    double p2 = Pi(1);
    double p1_2 = p1*p1; double p1_3 = p1_2*p1; double p1_4 = p1_3*p1; double p1_5 = p1_4*p1; double p1_6 = p1_5*p1; double p1_7 = p1_6*p1;
    double p2_2 = p2*p2; double p2_3 = p2_2*p2; double p2_4 = p2_3*p2; double p2_5 = p2_4*p2; double p2_6 = p2_5*p2; double p2_7 = p2_6*p2;

    return 2.0 * mat.alpha_1 * p2 
         + 4.0 * mat.alpha_11 * p2_3 
         + 2.0 * mat.alpha_12 * p1_2 * p2
         + 6.0 * mat.alpha_111 * p2_5
         + 4.0 * mat.alpha_112 * p1_2 * p2_3
         + 2.0 * mat.alpha_112 * p1_4 * p2 
         + 8.0 * mat.alpha_1111 * p2_7
         + 6.0 * mat.alpha_1112 * p1_2 * p2_5 
         + 2.0 * mat.alpha_1112 * p1_6 * p2 
         + 4.0 * mat.alpha_1122 * p1_4 * p2_3;
}

double Math::h_enthalpy_density(double U, double W, double chi, const Eigen::Vector2d& E, const Eigen::Vector2d& Pi) const {
    return U + W + chi - mat.eps0/2*E.squaredNorm() - E.dot(Pi);
}


