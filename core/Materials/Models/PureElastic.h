#pragma once
#include "Materials/Core/MaterialConcepts.h"
#include "Materials/Models/internal/ElasticMath.h"

class PureElasticMaterial final : public MaterialModel {
private:
    const MaterialConfig& material;
    ElasticMaterial m_elastic;

public:
    explicit PureElasticMaterial(const MaterialConfig& config) 
        : material(config), m_elastic(config) {}

    MaterialType get_type() const override { return MaterialType::PureElastic; }

    bool is_ferroelectric() const override { return false; }

    // Paramètres inactifs pour un matériau non-ferroélectrique
    double get_eta_k() const override { return 0.0; }
    double get_eps0() const override { return 0.0; } 
    double get_Gc() const override { return material.Gc; }
    double get_kappa() const override { return material.kappa; }
    double get_mu_v() const override { return material.mu_v; }
    double get_mu_p() const override { return 0.0; }
    double get_a0() const override { return 0.0; }

    // Énergies : Seule l'énergie mécanique existe
    double U_energy(const Eigen::Matrix2d& /*grad_P*/) const override { return 0.0; }
    double chi_energy(const Eigen::Vector2d& /*P*/) const override { return 0.0; }
    
    double W_energy(const Eigen::Vector2d& /*P*/, const Eigen::Matrix2d& strain) const override {
        return m_elastic.elastic_energy(strain); // On ignore P
    }

    // Mécanique
    Eigen::Matrix3d get_elastic_matrix() const override {
        return m_elastic.elastic_matrix();
    }
    Eigen::Vector3d compute_sigma_0(const Eigen::Vector2d& /*P*/) const override {
        return Eigen::Vector3d::Zero(); 
    }

    // Couplages annulés
    double compute_effective_permittivity(double, double, bool) const override { return 0.0; }
    Eigen::Vector2d compute_effective_polarization(const Eigen::Vector2d&, double, double, bool) const override {
        return Eigen::Vector2d::Zero();
    }
    double compute_H_drive(const Eigen::Matrix2d& /*grad_P*/, const Eigen::Vector2d& /*P*/,
                        const Eigen::Matrix2d& strain, const Eigen::Vector2d& /*E*/, bool /*is_impermeable*/) const override {
        return m_elastic.elastic_energy(strain);
    }

    // Ginzburg-Landau annulé
    GinzburgLandauTerms compute_GL_terms(const Eigen::Vector2d&, const Eigen::Matrix2d&, const Eigen::Vector2d&, double, bool) const override {
        return GinzburgLandauTerms{0.0, 0.0, 0.0, 0.0, 0.0};
    }
};
