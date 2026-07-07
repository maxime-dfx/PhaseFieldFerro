#pragma once
#include <Eigen/Dense>
#include "IO/Datafile.h"

class Math {
public:
    // Variables du matériau stockées localement en toute sécurité
    double xi, c0, mu_p, mu_v, a0, P0, t, eta_k, eps0;
    double alpha1, alpha11, alpha12, alpha111, alpha112, alpha1111, alpha1112, alpha1122;
    double b1, b2, b3, c1, c2, c3;

    // Le constructeur prend désormais le Datafile
    Math(const Datafile& config);
    
    void DimensionLess();

    double U_energy(const Eigen::Matrix2d& Pij) const;
    double W_energy(const Eigen::Vector2d& Pi, const Eigen::Matrix2d& epsjk) const;
    double chi_energy(const Eigen::Vector2d& Pi) const;
    double dW_dp1(const Eigen::Vector2d& Pi, const Eigen::Matrix2d& epsjk) const;
    double dW_dp2(const Eigen::Vector2d& Pi, const Eigen::Matrix2d& epsjk) const;
    double dchi_dp1(const Eigen::Vector2d& Pi) const;
    double dchi_dp2(const Eigen::Vector2d& Pi) const;
    double h_enthalpy_density(double U, double W, double chi, const Eigen::Vector2d& E, const Eigen::Vector2d& Pi) const;

    Eigen::Matrix3d get_elastic_matrix() const;
    Eigen::Vector3d compute_sigma_0(const Eigen::Vector2d& Pi) const;

    double compute_H_drive(const Eigen::Matrix2d& Pij, const Eigen::Vector2d& Pi, 
                           const Eigen::Matrix2d& eps, const Eigen::Vector2d& E, 
                           bool is_impermeable) const;

    double compute_polarization_force(int component, const Eigen::Vector2d& Pi, 
                                      const Eigen::Matrix2d& eps, double E_gp, 
                                      double penalite_fracture, bool is_impermeable) const;
};