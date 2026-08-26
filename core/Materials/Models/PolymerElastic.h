#pragma once
#include "Materials/Core/MaterialConcepts.h"
#include "Materials/Models/internal/ElasticMath.h"
#include "IO/include/ConfigTypes.h"

// Matériau passif pour la matrice polymère d'un piézocomposite : un isolant
// diélectrique linéaire, ni conducteur ni ferroélectrique.
//
//  - Permittivité : constante (eps0 de la config), valide pour l'équation de
//    Poisson (électrostatique). Contrairement à PureElasticMaterial (qui
//    renvoyait 0), on garde ici une vraie permittivité de condensateur.
//  - Mécanique : purement élastique (pas de couplage électrostrictif avec P).
//  - Polarisation : P doit rester bloquée à 0 dans ce matériau. Ce modèle
//    annule toutes les énergies/forces motrices qui en dépendent (Landau,
//    gradient, couplage), mais le VERROUILLAGE effectif des DOFs de
//    polarisation (P=0 imposé) est une responsabilité du DofMapper de
//    polarisation (cf. PolarizationDofMapper::get_bcs), qui interroge
//    is_ferroelectric() pour savoir quels noeuds bloquer.
class PolymerElasticMaterial final : public MaterialModel {
private:
    const MaterialConfig& material;
    ElasticMaterial m_elastic;

public:
    explicit PolymerElasticMaterial(const MaterialConfig& config)
        : material(config), m_elastic(config) {}

    MaterialType get_type() const override { return MaterialType::PolymerElastic; }

    bool is_ferroelectric() const override { return false; }
    double get_eta_k() const override { return material.eta_k; }
    double get_eps0() const override { return material.eps0; }
    double get_Gc() const override { return material.Gc; }
    double get_kappa() const override { return material.kappa; }
    double get_mu_v() const override { return material.mu_v; }
    double get_mu_p() const override { return 0.0; }
    double get_a0() const override { return 0.0; }


    double U_energy(const Eigen::Matrix2d& /*grad_P*/) const override { return 0.0; }
    double chi_energy(const Eigen::Vector2d& /*P*/) const override { return 0.0; }

    double W_energy(const Eigen::Vector2d& /*P*/, const Eigen::Matrix2d& strain) const override {
        return m_elastic.elastic_energy(strain);
    }

    // Mécanique : loi de Hooke standard, aucune contrainte spontanée liée à P.
    Eigen::Matrix3d get_elastic_matrix() const override {
        return m_elastic.elastic_matrix();
    }
    Eigen::Vector3d compute_sigma_0(const Eigen::Vector2d& /*P*/) const override {
        return Eigen::Vector3d::Zero();
    }

    double compute_effective_permittivity(double v, double /*eta_k_in*/, bool is_impermeable) const override {
        const double phase_factor = is_impermeable ? (v * v + material.eta_k) : 1.0;
        return material.eps0 * phase_factor;
    }

    Eigen::Vector2d compute_effective_polarization(const Eigen::Vector2d& /*P*/, double, double, bool) const override {
        return Eigen::Vector2d::Zero();
    }

    double compute_H_drive(const Eigen::Matrix2d& /*grad_P*/, const Eigen::Vector2d& /*P*/,
                        const Eigen::Matrix2d& strain, const Eigen::Vector2d& /*E*/, bool /*is_impermeable*/) const override {
        return m_elastic.elastic_energy(strain);
    }


    GinzburgLandauTerms compute_GL_terms(const Eigen::Vector2d&, const Eigen::Matrix2d&,
                                          const Eigen::Vector2d&, double, bool) const override {
        return GinzburgLandauTerms{0.0, 0.0, 0.0, 0.0, 0.0};
    }
};
