#pragma once

#include "Materials/Core/MaterialConcepts.h"
#include "Materials/Models/internal/ElasticMath.h"
#include "Materials/Models/internal/LandauMath.h"
#include "Materials/Models/internal/Math/Math.h"
#include "IO/include/ConfigTypes.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

// SingleCrystalMaterial (matériau Ferroelectric)
// -----------------------------------------------------------------------
// Anciennement Models/FerroelectricMaterial.h, qui déléguait toute la
// physique à internal/FerroelectricCore.h (elle-même une classe complète
// implémentant MaterialModel, dupliquant donc en grande partie
// SingleCrystalMaterial). Fusionné ici : SingleCrystalMaterial implémente
// directement l'interface en s'appuyant sur les toolboxes mathématiques
// pures ElasticMath (couplage électrostrictif + élasticité) et LandauMath
// (énergie de Landau-Devonshire + dérivées), sans classe intermédiaire.
class SingleCrystalMaterial final : public MaterialModel {
private:
    const MaterialConfig& material;
    ElasticMaterial elastic_family;

    double gradient_energy(const Eigen::Matrix2d& grad_P) const {
        return 0.5 * material.a0 * grad_P.squaredNorm();
    }

public:
    explicit SingleCrystalMaterial(const MaterialConfig& config)
        : material(config), elastic_family(config) {}

    MaterialType get_type() const override { return MaterialType::Ferroelectric; }

    double get_eta_k() const override { return material.eta_k; }
    double get_eps0() const override { return material.eps0; }
    double get_Gc() const override { return material.Gc; }
    double get_kappa() const override { return material.kappa; }
    double get_mu_v() const override { return material.mu_v; }
    double get_mu_p() const override { return material.mu_p; }
    double get_a0() const override { return material.a0; }

    // DENSITÉS D'ÉNERGIE
    double U_energy(const Eigen::Matrix2d& grad_P) const override {
        return gradient_energy(grad_P);
    }
    double W_energy(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain) const override {
        return elastic_family.total_energy(P, strain);
    }
    double chi_energy(const Eigen::Vector2d& P) const override {
        return LandauMath::compute(P, material).energy;
    }
    double electrostatic_energy(const Eigen::Vector2d& P, const Eigen::Vector2d& E) const {
        return -(P.dot(E)) - 0.5 * material.eps0 * E.squaredNorm();
    }

    // MÉCANIQUE
    Eigen::Matrix3d get_elastic_matrix() const override {
        return elastic_family.elastic_matrix();
    }
    Eigen::Vector3d compute_sigma_0(const Eigen::Vector2d& P) const override {
        return elastic_family.sigma_0(P);
    }

    // ÉLECTROSTATIQUE & COUPLAGE (Branches polymorphiques)
    double compute_effective_permittivity(double v, double eta_k_in, bool is_impermeable) const override {
        (void)eta_k_in;
        const double phase_factor = is_impermeable ? (v * v + material.eta_k) : 1.0;
        return material.eps0 * phase_factor;
    }
    Eigen::Vector2d compute_effective_polarization(const Eigen::Vector2d& P, double v, double eta_k_in, bool is_impermeable) const override {
        (void)eta_k_in;
        const double phase_factor = is_impermeable ? (v * v + material.eta_k) : 1.0;
        return P * phase_factor;
    }

    // FRACTURE DRIVE
    double compute_H_drive(const Eigen::Matrix2d& grad_P, const Eigen::Vector2d& P,
                           const Eigen::Matrix2d& strain, const Eigen::Vector2d& E,
                           bool is_impermeable) const override {
        const double U = gradient_energy(grad_P);
        const double W = elastic_family.total_energy(P, strain);
        const double H = is_impermeable ? (U + W + electrostatic_energy(P, E))
                                        : (U + W);
        return H;
    }

    // GINZBURG-LANDAU
    GinzburgLandauTerms compute_GL_terms(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain,
                                         const Eigen::Vector2d& E, double penalite_fracture,
                                         bool is_impermeable) const override {
        GinzburgLandauTerms GL{};
        const LandauTerms ld = LandauMath::compute(P, material);

        const double eps_11 = strain(0, 0);
        const double eps_22 = strain(1, 1);
        const double eps_12 = MaterialMath::symmetrized_shear(strain);

        Eigen::Vector2d dW_dP;
        dW_dP(0) = - (material.b1 * eps_11 + material.b2 * eps_22) * P(0) - 2.0 * material.b3 * eps_12 * P(1);
        dW_dP(1) = - (material.b1 * eps_22 + material.b2 * eps_11) * P(1) - 2.0 * material.b3 * eps_12 * P(0);

        Eigen::Matrix2d d2W_dP2 = Eigen::Matrix2d::Zero();
        d2W_dP2(0, 0) = -material.b1 * eps_11 - material.b2 * eps_22;
        d2W_dP2(1, 1) = -material.b1 * eps_22 - material.b2 * eps_11;
        d2W_dP2(0, 1) = d2W_dP2(1, 0) = -2.0 * material.b3 * eps_12;

        // CORRECTION : chi_factor ne doit pas inclure la pénalité de fracture
        const double chi_factor = 1.0;
        const double E_factor = is_impermeable ? penalite_fracture : 1.0;

        GL.force_px = penalite_fracture * dW_dP(0) + chi_factor * ld.dchi_dP(0) - E_factor * E(0);
        GL.force_py = penalite_fracture * dW_dP(1) + chi_factor * ld.dchi_dP(1) - E_factor * E(1);

        GL.J_11 = penalite_fracture * d2W_dP2(0, 0) + chi_factor * ld.d2chi_dP2(0, 0);
        GL.J_22 = penalite_fracture * d2W_dP2(1, 1) + chi_factor * ld.d2chi_dP2(1, 1);
        GL.J_12 = penalite_fracture * d2W_dP2(0, 1) + chi_factor * ld.d2chi_dP2(0, 1);

        return GL;
    }
};
