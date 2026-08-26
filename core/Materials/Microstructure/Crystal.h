#pragma once

#include "Materials/Core/MaterialConcepts.h"
#include "Utils/include/Rotation2D.h"

// Crystal
// -----------------------------------------------------------------------
// Gère 1 grain d'une microstructure polycristalline : son germe (seed
// point, issu de VoronoiTessellation), son orientation cristalline (angle
// theta, en radians, entre l'axe x1 du cristal local et l'axe x global) et
// la matrice de rotation qui en découle.
//
// Fusionne les anciens Microstructure/Grain.h (pure structure de données)
// et Models/OrientedMaterial.h (décorateur de MaterialModel appliquant la
// rotation aux champs physiques) : les deux ne faisaient sens que l'un par
// rapport à l'autre (l'orientation d'un grain n'a d'utilité que pour
// tourner un matériau), ils vivent donc désormais dans le même fichier.
struct Crystal {
    int id = -1;
    double seed_x = 0.0;
    double seed_y = 0.0;
    double theta = 0.0; // orientation cristalline, radians
};

// CrystalMaterial
// -----------------------------------------------------------------------
// Décore un MaterialModel de base (typiquement un SingleCrystalMaterial
// monocristallin) avec la rotation d'un Crystal donné : les champs entrant
// (P, strain, E, grad_P) sont tournés vers le repère local du cristal
// avant d'être passés au matériau de base, puis les résultats (forces,
// contraintes, tenseur de rigidité) sont tournés en retour vers le repère
// global. Anciennement OrientedMaterial.
class CrystalMaterial final : public MaterialModel {
private:
    const MaterialModel& base_material;
    Rotation2D rotation;

public:
    CrystalMaterial(const MaterialModel& material, double orientation_theta)
        : base_material(material), rotation(orientation_theta) {}

    CrystalMaterial(const MaterialModel& material, const Crystal& crystal)
        : base_material(material), rotation(crystal.theta) {}

    MaterialType get_type() const override { return base_material.get_type(); }

    double get_eta_k() const override { return base_material.get_eta_k(); }
    double get_eps0() const override { return base_material.get_eps0(); }

    double get_Gc() const override { return base_material.get_Gc(); }
    double get_kappa() const override { return base_material.get_kappa(); }
    double get_mu_v() const override { return base_material.get_mu_v(); }
    double get_mu_p() const override { return base_material.get_mu_p(); }
    double get_a0() const override { return base_material.get_a0(); }

    bool is_ferroelectric() const override { return base_material.is_ferroelectric(); }

    double U_energy(const Eigen::Matrix2d& grad_P) const override {
        return base_material.U_energy(rotation.to_local(grad_P));
    }

    double W_energy(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain) const override {
        return base_material.W_energy(rotation.to_local(P), rotation.to_local(strain));
    }

    double chi_energy(const Eigen::Vector2d& P) const override {
        return base_material.chi_energy(rotation.to_local(P));
    }

    Eigen::Matrix3d get_elastic_matrix() const override {
        return rotation.to_global(base_material.get_elastic_matrix());
    }

    Eigen::Vector3d compute_sigma_0(const Eigen::Vector2d& P) const override {
        return rotation.to_global(base_material.compute_sigma_0(rotation.to_local(P)));
    }

    double compute_effective_permittivity(double v, double eta_k_in, bool is_impermeable) const override {
        return base_material.compute_effective_permittivity(v, eta_k_in, is_impermeable);
    }

    Eigen::Vector2d compute_effective_polarization(const Eigen::Vector2d& P, double v, double eta_k_in, bool is_impermeable) const override {
        return rotation.to_global(base_material.compute_effective_polarization(rotation.to_local(P), v, eta_k_in, is_impermeable));
    }

    double compute_H_drive(const Eigen::Matrix2d& grad_P, const Eigen::Vector2d& P,
                           const Eigen::Matrix2d& strain, const Eigen::Vector2d& E,
                           bool is_impermeable) const override {
        return base_material.compute_H_drive(rotation.to_local(grad_P), rotation.to_local(P), rotation.to_local(strain), rotation.to_local(E), is_impermeable);
    }

    GinzburgLandauTerms compute_GL_terms(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain,
                                         const Eigen::Vector2d& E, double penalite_fracture,
                                         bool is_impermeable) const override {
        GinzburgLandauTerms terms = base_material.compute_GL_terms(rotation.to_local(P), rotation.to_local(strain), rotation.to_local(E), penalite_fracture, is_impermeable);
        Eigen::Vector2d force = rotation.to_global(Eigen::Vector2d(terms.force_px, terms.force_py));
        Eigen::Matrix2d J_local;
        J_local << terms.J_11, terms.J_12,
                   terms.J_12, terms.J_22;
        Eigen::Matrix2d J_global = rotation.to_global(J_local);
        terms.force_px = force.x();
        terms.force_py = force.y();
        terms.J_11 = J_global(0, 0);
        terms.J_12 = J_global(0, 1);
        terms.J_22 = J_global(1, 1);
        return terms;
    }
};

// Alias conservé pour compatibilité de nommage avec le code/tests existants
// qui référencent le décorateur d'orientation sous son ancien nom.
using OrientedMaterial = CrystalMaterial;
