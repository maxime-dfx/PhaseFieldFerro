#pragma once
#include <Eigen/Dense>
#include "Physics/Material.h"
#include "Physics/Electrostatics.h"
#include "Physics/Polarization.h"

class Math {
    private :
    double p1, p2;
    double p1_2, p1_3, p1_5, p1_7;
    double p2_2, p2_4, p2_6;

    public :
    const Material& mat;

    double xi_adim, pi_adim, t_adim, eps0_adim, phi_adim, alpha1_adim, alpha11_adim, alpha12_adim, alpha111_adim, alpha112_adim, alpha1111_adim, alpha1112_adim, alpha1122_adim, b1_adim, b2_adim, b3_adim, c1_adim, c2_adim, c3_adim;  
    double xi, c0, mu_p, mu_v, a0, P0, t, eps0, phi, alpha1, alpha11, alpha12, alpha111, alpha112, alpha1111, alpha1112, alpha1122, b1, b2, b3, c1, c2, c3, eta_k;

    Math(const Material& mat);
    void DimensionLess();
    double U_energy(const Eigen::Matrix2d& Pij) const;
    double W_energy(const Eigen::Vector2d& Pi, const Eigen::Matrix2d& epsjk) const;
    double dW_dp1(const Eigen::Vector2d& Pi, const Eigen::Matrix2d& epsjk) const;
    double dW_dp2(const Eigen::Vector2d& Pi, const Eigen::Matrix2d& epsjk) const;
    double chi_energy(const Eigen::Vector2d& Pi) const;
    double dchi_dp1(const Eigen::Vector2d& Pi) const;
    double dchi_dp2(const Eigen::Vector2d& Pi) const;
    double h_enthalpy_density(double U_energy, double W_energy, double chi_energy, const Eigen::Vector2d& E, const Eigen::Vector2d& Pi) const;
    double H_enthalpy(double u, double p, double phi) const;

};