#pragma once
#include "Materials/Core/MaterialConcepts.h"
#include "Materials/Models/internal/ElasticMath.h"
#include "Materials/Models/internal/LandauMath.h"
#include "IO/include/ConfigTypes.h"

// Matériau du JOINT DE GRAIN, à l'intérieur d'une phase ferroélectrique
// polycristalline : une couche fine, structurellement désordonnée, avec ses
// PROPRES constantes (Gc, kappa, mu_v, module élastique, permittivité,
// et EVENTUELLEMENT ses propres constantes de Landau a0/b1/b2/b3/alpha_*),
// entièrement indépendantes de celles du grain massif.
//
// Ce matériau n'est JAMAIS associé à un Physical Tag Gmsh (Element::ref_tag) :
// il est résolu par élément via la détection de joint de grain
// (Polycrystal::is_grain_boundary_element) et référencé par
// CrystalConfig::grain_boundary_material_id. Il doit donc être déclaré comme
// n'importe quelle entrée [[materials]], avec type = "GrainBoundary", mais
// son `id` n'est utilisé que via ce champ crystal, jamais matché sur un
// ref_tag de maillage.
//
// DEUX COMPORTEMENTS PHYSIQUES, pilotés par `is_ferroelectric` dans le TOML
// de CETTE entrée (cf. Datafile::parse_materials_list, qui - CONTRAIREMENT
// aux autres types - respecte la valeur explicite de l'utilisateur pour
// GrainBoundary au lieu de la déduire de `type`) :
//
//  - is_ferroelectric = false (défaut historique) : couche désordonnée sans
//    polarisation propre. P=0 est verrouillé en amont par
//    PolarizationDofMapper (cf. is_ferroelectric() ci-dessous), et tous les
//    termes de Landau/couplage électromécanique sont nuls ici - seules Gc,
//    kappa, l'élasticité linéaire et la permittivité diélectrique du GB
//    interviennent (comme PolymerElasticMaterial).
//
//  - is_ferroelectric = true : le GB devient un second matériau
//    ferroélectrique à part entière (même physique que SingleCrystalMaterial
//    dans Ferroelectric.h : gradient d'énergie U, énergie de Landau chi,
//    couplage électrostrictif W, force motrice de fracture H_drive avec
//    terme électrostatique, polarisation effective), mais avec SES PROPRES
//    constantes a0/mu_p/b1/b2/b3/alpha_* lues depuis SA propre entrée
//    [[materials]] (héritées de [material] si non surchargées, cf.
//    parse_materials_list). Utile pour représenter un joint "ferroélectrique
//    affaibli" plutôt qu'un joint totalement désordonné.
class GrainBoundaryMaterial final : public MaterialModel {
private:
    const MaterialConfig& material;
    ElasticMaterial m_elastic;

    double gradient_energy(const Eigen::Matrix2d& grad_P) const {
        return 0.5 * material.a0 * grad_P.squaredNorm();
    }

    double electrostatic_energy(const Eigen::Vector2d& P, const Eigen::Vector2d& E) const {
        return -(P.dot(E)) - 0.5 * material.eps0 * E.squaredNorm();
    }

public:
    explicit GrainBoundaryMaterial(const MaterialConfig& config)
        : material(config), m_elastic(config) {}

    MaterialType get_type() const override { return MaterialType::GrainBoundary; }

    // Piloté par le TOML (cf. commentaire de classe ci-dessus) : ce n'est
    // plus figé à `false`. C'est cette valeur que consulte
    // PolarizationDofMapper::element_is_ferroelectric() pour décider si les
    // nœuds du joint de grain doivent avoir P verrouillé à 0 ou non.
    bool is_ferroelectric() const override { return material.is_ferroelectric; }

    double get_eta_k() const override { return material.eta_k; }
    double get_eps0() const override { return material.eps0; }
    double get_Gc() const override { return material.Gc; }
    double get_kappa() const override { return material.kappa; }
    double get_mu_v() const override { return material.mu_v; }
    // mu_p/a0 ne pilotent que la dynamique de P (mass_coeff = mu_p/dt,
    // stiff = a0*penalite dans Polarization.h) : sans intérêt (et
    // dangereux, mass_coeff nul -> système singulier) si is_ferroelectric
    // est false, cas où ce nœud est de toute façon verrouillé en amont et
    // ces coefficients ne sont jamais consultés par l'assemblage.
    double get_mu_p() const override { return material.is_ferroelectric ? material.mu_p : 0.0; }
    double get_a0()   const override { return material.is_ferroelectric ? material.a0   : 0.0; }

    // DENSITÉS D'ÉNERGIE
    double U_energy(const Eigen::Matrix2d& grad_P) const override {
        return material.is_ferroelectric ? gradient_energy(grad_P) : 0.0;
    }
    double chi_energy(const Eigen::Vector2d& P) const override {
        return material.is_ferroelectric ? LandauMath::compute(P, material).energy : 0.0;
    }
    double W_energy(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain) const override {
        // Le couplage électrostrictif (b1/b2/b3) n'a de sens que si P peut
        // réellement varier ; sinon on retombe sur l'énergie élastique pure
        // (comportement historique, identique à PolymerElasticMaterial).
        return material.is_ferroelectric ? m_elastic.total_energy(P, strain)
                                          : m_elastic.elastic_energy(strain);
    }

    // MÉCANIQUE
    Eigen::Matrix3d get_elastic_matrix() const override {
        return m_elastic.elastic_matrix();
    }
    Eigen::Vector3d compute_sigma_0(const Eigen::Vector2d& P) const override {
        return material.is_ferroelectric ? m_elastic.sigma_0(P) : Eigen::Vector3d::Zero();
    }

    // ÉLECTROSTATIQUE & COUPLAGE
    double compute_effective_permittivity(double v, double /*eta_k_in*/, bool is_impermeable) const override {
        const double phase_factor = is_impermeable ? (v * v + material.eta_k) : 1.0;
        return material.eps0 * phase_factor;
    }
    Eigen::Vector2d compute_effective_polarization(const Eigen::Vector2d& P, double v, double /*eta_k_in*/, bool is_impermeable) const override {
        if (!material.is_ferroelectric) return Eigen::Vector2d::Zero();
        const double phase_factor = is_impermeable ? (v * v + material.eta_k) : 1.0;
        return P * phase_factor;
    }

    // FRACTURE DRIVE
    double compute_H_drive(const Eigen::Matrix2d& grad_P, const Eigen::Vector2d& P,
                        const Eigen::Matrix2d& strain, const Eigen::Vector2d& E, bool is_impermeable) const override {
        if (!material.is_ferroelectric) {
            // La force motrice de rupture est l'énergie élastique stockée.
            return m_elastic.elastic_energy(strain);
        }
        const double U = gradient_energy(grad_P);
        const double W = m_elastic.total_energy(P, strain);
        return is_impermeable ? (U + W + electrostatic_energy(P, E)) : (U + W);
    }

    // GINZBURG-LANDAU
    GinzburgLandauTerms compute_GL_terms(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain,
                                          const Eigen::Vector2d& E, double penalite_fracture,
                                          bool is_impermeable) const override {
        if (!material.is_ferroelectric) {
            return GinzburgLandauTerms{0.0, 0.0, 0.0, 0.0, 0.0};
        }

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
