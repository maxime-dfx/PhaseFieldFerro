#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include "IO/include/ConfigTypes.h"

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

    // Type concret (cf. MaterialType) : utilisé pour les diagnostics/logs et
    // par MaterialManager pour retrouver quelle classe a été instanciée.
    // Chaque classe Models/*.h renvoie sa propre valeur fixe.
    virtual MaterialType get_type() const = 0;

    // Constantes physiques de base accessibles à tous
    virtual double get_eta_k() const = 0;
    virtual double get_eps0() const = 0;

    // Paramètres matériau désormais résolus par élément (via MaterialManager
    // + Element::ref_tag) plutôt que lus globalement sur Datafile::material.
    // Indispensable pour un piézocomposite : la résine (matrice passive) et
    // la céramique (inclusions actives) n'ont pas la même ténacité à la
    // fracture Gc, ni la même longueur de régularisation kappa, ni les mêmes
    // mobilités/coefficients de gradient de polarisation.
    virtual double get_Gc() const = 0;    // Ténacité à la fracture (module Fracture)
    virtual double get_kappa() const = 0; // Longueur de régularisation AT (module Fracture)
    virtual double get_mu_v() const = 0;  // Mobilité du champ de phase v (module Fracture)
    virtual double get_mu_p() const = 0;  // Mobilité de la polarisation (module Polarization)
    virtual double get_a0() const = 0;    // Coefficient du gradient de polarisation (module Polarization)

    // Indique si ce matériau porte une physique ferroélectrique active
    // (Landau, gradient, couplage électrostrictif, DOFs de polarisation
    // libres). Vrai par défaut pour ne pas casser les matériaux existants ;
    // les matériaux passifs (ex. matrice polymère diélectrique) doivent
    // le surcharger à false.
    virtual bool is_ferroelectric() const { return true; }

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
