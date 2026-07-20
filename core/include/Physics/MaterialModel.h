#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

// On conserve la structure de retour pour Ginzburg-Landau (ou on la rend plus générique)
struct GinzburgLandauTerms {
    double J_11;
    double J_12;
    double J_22;
    double force_px;
    double force_py;
}; // Conservé de l'ancienne implémentation

class MaterialModel {
public:
    virtual ~MaterialModel() = default;

    // Constantes physiques de base accessibles à tous
    virtual double get_eta_k() const = 0;
    virtual double get_eps0() const = 0;

    // Calculs des énergies
    virtual double U_energy(const Eigen::Matrix2d& grad_P) const = 0;
    virtual double W_energy(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain) const = 0;
    virtual double chi_energy(const Eigen::Vector2d& P) const = 0;

    // Calculs mécaniques
    virtual Eigen::Matrix3d get_elastic_matrix() const = 0;
    virtual Eigen::Vector3d compute_sigma_0(const Eigen::Vector2d& P) const = 0;

    // Calculs électrostatiques et couplages
    virtual double compute_effective_permittivity(double v, double eta_k_in, bool is_impermeable) const = 0;
    virtual Eigen::Vector2d compute_effective_polarization(const Eigen::Vector2d& P, double v, double eta_k_in, bool is_impermeable) const = 0;

    // Calculs de forces motrices (Fracture & Polarisation)
    virtual double compute_H_drive(const Eigen::Matrix2d& grad_P, const Eigen::Vector2d& P, 
                                   const Eigen::Matrix2d& strain, const Eigen::Vector2d& E, 
                                   bool is_impermeable) const = 0;
                                   
    virtual GinzburgLandauTerms compute_GL_terms(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain,
                                                 const Eigen::Vector2d& E, double penalite_fracture,
                                                 bool is_impermeable) const = 0;
};