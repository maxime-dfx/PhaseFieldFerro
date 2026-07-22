#pragma once
#include "Physics/MaterialModel.h"
#include "IO/ConfigTypes.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <tracy/Tracy.hpp>

// =============================================================================
// FerroelectricMaterial
// -----------------------------------------------------------------------------
// Implementation of the coupled ferroelectric / brittle-fracture constitutive
// model of Abdollahi & Arias, Acta Materialia 59 (2011) 4733-4746.
//
// Conventions used throughout this file:
//   - P = (p1, p2)                : normalized polarization components
//   - strain(i,j)                 : full (symmetrized) tensorial strain, i.e.
//                                    strain(0,1) == strain(1,0) == eps_12
//   - E = (E1, E2)                 : electric field, E = -grad(phi)
//   - v                            : fracture phase-field (1 = intact, 0 = broken)
//   - degradation / penalite       : (v^2 + eta_k), the jump-set function
//   - is_impermeable               : selects between Eq. (11) [permeable] and
//                                    Eq. (13) [impermeable] of the paper
//
// IMPORTANT PHYSICS NOTE ON THE TWO CRACK MODELS (Section 2.3 of the paper):
//
//   Permeable   (Eq. 11): h = (v^2+eta_k)[U(grad p) + W(p,eps) + chi(p)]
//                             - eps0/2 |E|^2 - E.p
//     => the Landau-Devonshire energy chi(p) IS degraded by the crack,
//        the electrostatic terms are NOT.
//
//   Impermeable (Eq. 13): h = (v^2+eta_k)[U(grad p) + W(p,eps) - eps0/2|E|^2 - E.p]
//                             + chi(p)
//     => the electrostatic terms ARE degraded by the crack,
//        chi(p) is NOT.
//
// Both compute_GL_terms() (the driving force for the polarization gradient
// flow, Eq. 15) and compute_H_drive() (the fracture driving force feeding
// Eq. 16) must branch on is_impermeable accordingly. Prior versions of this
// file silently dropped this branch, which made both crack models identical
// and wrong. This version restores the distinction explicitly.
// =============================================================================
class FerroelectricMaterial : public MaterialModel {
private:
    const MaterialConfig& material;

public:
    explicit FerroelectricMaterial(const MaterialConfig& config)
        : material(config) {}

    double get_eta_k() const override {
        return material.eta_k;
    }

    double get_eps0() const override {
        return material.eps0;
    }

    // =========================================================================
    // 0. ENERGY DENSITIES (Eqs. 1, 3, 4, 5)
    // =========================================================================

    // U(grad p) = a0/2 * (p1,1^2 + p1,2^2 + p2,1^2 + p2,2^2)   [Eq. (3)]
    double U_energy(const Eigen::Matrix2d& grad_P) const override {
        double grad_P_norm_sq = grad_P(0,0)*grad_P(0,0) + grad_P(0,1)*grad_P(0,1) +
                                 grad_P(1,0)*grad_P(1,0) + grad_P(1,1)*grad_P(1,1);
        return 0.5 * material.a0 * grad_P_norm_sq;
    }

    // W(p, eps) = W_mec(p, eps) + W_elas(eps)    [Eq. (4)]
    // W_mec = -b1/2(eps11 p1^2 + eps22 p2^2) - b2/2(eps22 p1^2 + eps11 p2^2)
    //         - b3(eps21+eps12) p1 p2
    // W_elas = c1/2(eps11^2+eps22^2) + c2 eps11 eps22 + c3/2(eps12^2+eps21^2)
    double W_energy(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain) const override {
        return W_mec_energy(P, strain) + W_elas_energy(strain);
    }

    // Electrostrictive coupling part of W only (useful for diagnostics/logging).
    double W_mec_energy(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain) const {
        double p1 = P(0), p2 = P(1);
        double p1_2 = p1 * p1, p2_2 = p2 * p2;

        return -0.5 * (material.b1 * strain(0,0) * p1_2 +
                       material.b2 * strain(1,1) * p1_2 +
                       material.b1 * strain(1,1) * p2_2 +
                       material.b2 * strain(0,0) * p2_2 +
                       4.0 * material.b3 * strain(0,1) * p1 * p2);
    }

    // Pure anisotropic elastic part of W only (useful for diagnostics/logging).
    double W_elas_energy(const Eigen::Matrix2d& strain) const {
        double c1 = material.c1, c2 = material.c2, c3 = material.c3;
        return 0.5 * c1 * (strain(0,0)*strain(0,0) + strain(1,1)*strain(1,1))
             + c2 * strain(0,0)*strain(1,1)
             + 0.5 * c3 * (strain(0,1)*strain(0,1) + strain(1,0)*strain(1,0));
    }

    // chi(p) : Landau-Devonshire phase-separation energy    [Eq. (5)]
    double chi_energy(const Eigen::Vector2d& P) const override {
        double p1 = P(0), p2 = P(1);
        double p1_2 = p1 * p1, p1_4 = p1_2 * p1_2, p1_6 = p1_4 * p1_2, p1_8 = p1_4 * p1_4;
        double p2_2 = p2 * p2, p2_4 = p2_2 * p2_2, p2_6 = p2_4 * p2_2, p2_8 = p2_4 * p2_4;

        return material.alpha_1    * (p1_2 + p2_2) +
               material.alpha_11   * (p1_4 + p2_4) +
               material.alpha_12   * (p1_2 * p2_2) +
               material.alpha_111  * (p1_6 + p2_6) +
               material.alpha_112  * (p1_2 * p2_4 + p2_2 * p1_4) +
               material.alpha_1111 * (p1_8 + p2_8) +
               material.alpha_1112 * (p1_6 * p2_2 + p2_6 * p1_2) +
               material.alpha_1122 * (p1_4 * p2_4);
    }

    // Electrostatic contribution -eps0/2|E|^2 - E.p appearing in h (Eq. 2).
    // Exposed here (rather than inlined) so it can be logged/diagnosed and
    // reused identically by compute_H_drive().
    double electrostatic_energy(const Eigen::Vector2d& P, const Eigen::Vector2d& E) const {
        return -(P(0) * E(0) + P(1) * E(1)) - 0.5 * material.eps0 * (E(0)*E(0) + E(1)*E(1));
    }

    // Full electromechanical enthalpy density h at a Gauss point, accounting
    // for the crack jump-set function according to the permeable (Eq. 11) or
    // impermeable (Eq. 13) model. Handy single entry point for diagnostics.
    double enthalpy_density(const Eigen::Matrix2d& grad_P, const Eigen::Vector2d& P,
                             const Eigen::Matrix2d& strain, const Eigen::Vector2d& E,
                             double v, bool is_impermeable) const {
        double degradation = v * v + material.eta_k;
        double U = U_energy(grad_P);
        double W = W_energy(P, strain);
        double chi = chi_energy(P);
        double W_elec = electrostatic_energy(P, E);

        if (is_impermeable) {
            // h = (v^2+eta_k)[U + W + W_elec] + chi        [Eq. (13)]
            return degradation * (U + W + W_elec) + chi;
        }
        // h = (v^2+eta_k)[U + W + chi] + W_elec             [Eq. (11)]
        return degradation * (U + W + chi) + W_elec;
    }

    // =========================================================================
    // 1. MECHANICS : K_u * du = F_piezo + F_ext
    // =========================================================================

    // Elastic stiffness matrix C in Voigt notation, engineering shear
    // convention (gamma_12 = 2*eps_12), consistent with the B-matrix used
    // by MechanicsAssembler.
    Eigen::Matrix3d get_elastic_matrix() const override {
        Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
        C(0,0) = material.c1;
        C(1,1) = material.c1;
        C(0,1) = material.c2;
        C(1,0) = material.c2;
        C(2,2) = material.c3 / 2.0;
        return C;
    }

    // sigma_0 = dW_mec/d(eps), expressed in the Voigt basis used by B.
    // B uses the engineering shear gamma_12 = 2*eps_12_tensorial, so the
    // conjugate component sigma_0(2) must equal dW/dgamma_12
    // = 0.5 * dW/deps_12_tensorial.
    // With W_mec as implemented: dW_mec/deps_12_tensorial = -2*b3*p1*p2
    //  => sigma_0(2) = 0.5 * (-2*b3*p1*p2) = -b3*p1*p2.
    //
    // NOTE (assembly-side bug, not fixed here): MechanicsAssembler currently
    // adds this term with a '+=' to F_local instead of '-=' . The correct
    // weak form (Eq. 17) is  K*u = F_traction - integral(B^T sigma_0), so the
    // assembler must be corrected to subtract this contribution.
    Eigen::Vector3d compute_sigma_0(const Eigen::Vector2d& P) const override {
        Eigen::Vector3d sigma_0;
        double p1_2 = P(0) * P(0);
        double p2_2 = P(1) * P(1);

        sigma_0(0) = -0.5 * (material.b1 * p1_2 + material.b2 * p2_2);
        sigma_0(1) = -0.5 * (material.b1 * p2_2 + material.b2 * p1_2);
        sigma_0(2) = -material.b3 * P(0) * P(1);

        return sigma_0;
    }

    // =========================================================================
    // 2. ELECTROSTATICS : K_phi * d_phi = F_pol + F_w
    // =========================================================================

    // Effective permittivity multiplying grad(phi) in K_phi.
    // Impermeable crack (Eq. 13): the -eps0/2|E|^2 - E.p terms ARE degraded
    //   => eps_eff = eps0 * (v^2 + eta_k)
    // Permeable crack (Eq. 11): these terms are untouched by the crack
    //   => eps_eff = eps0
    double compute_effective_permittivity(double v, double eta_k_in, bool is_impermeable) const override {
        (void)eta_k_in;
        double phase_factor = is_impermeable ? (v * v + material.eta_k) : 1.0;
        return material.eps0 * phase_factor;
    }

    // Effective polarization entering the electrostatic RHS (-E.p term),
    // degraded identically to the permittivity above.
    Eigen::Vector2d compute_effective_polarization(const Eigen::Vector2d& P, double v, double eta_k_in, bool is_impermeable) const override {
        (void)eta_k_in;
        double phase_factor = is_impermeable ? (v * v + material.eta_k) : 1.0;
        return P * phase_factor;
    }

    // =========================================================================
    // 3. FRACTURE : K_v * d_v = F_drive       (Eq. 16)
    // =========================================================================

    // Driving force H = max(0, dh/d(v^2)|_relevant) feeding the Allen-Cahn
    // fracture evolution. Only the energy contributions that are actually
    // multiplied by the jump-set function (v^2+eta_k) in h contribute to
    // this driving force - see the branch below, which mirrors Eqs. (11)/(13).
    //
    // Permeable   (Eq. 11): degraded part = U + W + chi   -> H = U + W + chi
    // Impermeable (Eq. 13): degraded part = U + W + W_elec -> H = U + W + W_elec
    //
    // Previous implementation always used U + W + W_elec regardless of
    // is_impermeable, which is only correct for the impermeable case and
    // silently dropped chi(p) from the permeable case (the regime used in
    // the majority of the paper's simulations, Section 3.2.1).
    double compute_H_drive(const Eigen::Matrix2d& grad_P, const Eigen::Vector2d& P,
                            const Eigen::Matrix2d& strain, const Eigen::Vector2d& E,
                            bool is_impermeable) const override {
        double U = U_energy(grad_P);
        double W = W_energy(P, strain);

        double H;
        if (is_impermeable) {
            double W_elec = electrostatic_energy(P, E);
            H = U + W + W_elec;
        } else {
            double chi = chi_energy(P);
            H = U + W + chi;
        }
        return std::max(0.0, H);
    }

    // =========================================================================
    // 4. POLARIZATION (Ginzburg-Landau gradient flow)      (Eq. 15)
    // =========================================================================

    // Effective single-field weak-form terms driving the polarization
    // evolution:
    //   mu_p * dP/dt = -dh/dP = -[ degradation-weighted dW/dP
    //                              + (branch-dependent) dchi/dP
    //                              - (branch-dependent) E ]
    //
    // Permeable   (Eq. 11): chi IS degraded, E is NOT
    //   => dh/dp_i = penalite*dW/dp_i + penalite*dchi/dp_i - E_i
    // Impermeable (Eq. 13): E IS degraded, chi is NOT
    //   => dh/dp_i = penalite*dW/dp_i + dchi/dp_i - penalite*E_i
    //
    // Previous implementation always used
    //   dh/dp_i = penalite*dW/dp_i + dchi/dp_i - E_i
    // which is neither the permeable nor the impermeable form (it silently
    // discarded the is_impermeable flag), and made both crack models
    // indistinguishable in the polarization dynamics.
    GinzburgLandauTerms compute_GL_terms(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain,
                                          const Eigen::Vector2d& E, double penalite_fracture,
                                          bool is_impermeable) const override {
        GinzburgLandauTerms GL;

        // --- dchi/dP and d2chi/dP2 ---
        double a1 = material.alpha_1, a11 = material.alpha_11, a12 = material.alpha_12;
        double a111 = material.alpha_111, a112 = material.alpha_112;
        double a1111 = material.alpha_1111, a1112 = material.alpha_1112, a1122 = material.alpha_1122;

        double p1 = P(0), p2 = P(1);
        double p1_2 = p1 * p1, p1_3 = p1_2 * p1, p1_4 = p1_2 * p1_2, p1_5 = p1_4 * p1, p1_6 = p1_4 * p1_2, p1_7 = p1_6 * p1;
        double p2_2 = p2 * p2, p2_3 = p2_2 * p2, p2_4 = p2_2 * p2_2, p2_5 = p2_4 * p2, p2_6 = p2_4 * p2_2, p2_7 = p2_6 * p2;

        double dchi_dp1 = 2.0*a1*p1 + 4.0*a11*p1_3 + 2.0*a12*p1*p2_2 + 6.0*a111*p1_5 + 4.0*a112*p1_3*p2_2 + 2.0*a112*p1*p2_4 + 8.0*a1111*p1_7 + 6.0*a1112*p1_5*p2_2 + 2.0*a1112*p1*p2_6 + 4.0*a1122*p1_3*p2_4;
        double dchi_dp2 = 2.0*a1*p2 + 4.0*a11*p2_3 + 2.0*a12*p2*p1_2 + 6.0*a111*p2_5 + 4.0*a112*p2_3*p1_2 + 2.0*a112*p2*p1_4 + 8.0*a1111*p2_7 + 6.0*a1112*p2_5*p1_2 + 2.0*a1112*p2*p1_6 + 4.0*a1122*p2_3*p1_4;

        double d2chi_dp12    = 2.0*a1 + 12.0*a11*p1_2 + 2.0*a12*p2_2 + 30.0*a111*p1_4 + 12.0*a112*p1_2*p2_2 + 2.0*a112*p2_4 + 56.0*a1111*p1_6 + 30.0*a1112*p1_4*p2_2 + 2.0*a1112*p2_6 + 12.0*a1122*p1_2*p2_4;
        double d2chi_dp22    = 2.0*a1 + 12.0*a11*p2_2 + 2.0*a12*p1_2 + 30.0*a111*p2_4 + 12.0*a112*p2_2*p1_2 + 2.0*a112*p1_4 + 56.0*a1111*p2_6 + 30.0*a1112*p2_4*p1_2 + 2.0*a1112*p1_6 + 12.0*a1122*p2_2*p1_4;
        double d2chi_dp1dp2  = 4.0*a12*p1*p2 + 8.0*a112*p1*p2*(p1_2 + p2_2) + 12.0*a1112*p1*p2*(p1_4 + p2_4) + 16.0*a1122*p1_2*p1*p2_2*p2;

        // --- dW/dP and d2W/dP2 ---
        double eps_11 = strain(0,0), eps_22 = strain(1,1), eps_12 = 0.5 * (strain(0,1) + strain(1,0));

        double dW_dp1 = - (material.b1 * eps_11 + material.b2 * eps_22) * p1 - 2.0 * material.b3 * eps_12 * p2;
        double dW_dp2 = - (material.b1 * eps_22 + material.b2 * eps_11) * p2 - 2.0 * material.b3 * eps_12 * p1;

        double d2W_dp12   = -material.b1 * eps_11 - material.b2 * eps_22;
        double d2W_dp22   = -material.b1 * eps_22 - material.b2 * eps_11;
        double d2W_dp1dp2 = -2.0 * material.b3 * eps_12;

        // --- branch-dependent factors, cf. Eqs. (11)/(13) ---
        // Permeable:   chi degraded (factor = penalite), E not degraded (factor = 1)
        // Impermeable: chi not degraded (factor = 1),      E degraded (factor = penalite)
        double chi_factor = is_impermeable ? 1.0 : penalite_fracture;
        double E_factor    = is_impermeable ? penalite_fracture : 1.0;

        GL.force_px = penalite_fracture * dW_dp1 + chi_factor * dchi_dp1 - E_factor * E(0);
        GL.force_py = penalite_fracture * dW_dp2 + chi_factor * dchi_dp2 - E_factor * E(1);

        GL.J_11 = penalite_fracture * d2W_dp12   + chi_factor * d2chi_dp12;
        GL.J_22 = penalite_fracture * d2W_dp22   + chi_factor * d2chi_dp22;
        GL.J_12 = penalite_fracture * d2W_dp1dp2 + chi_factor * d2chi_dp1dp2;

        return GL;
    }

    // =========================================================================
    // 5. DIAGNOSTICS / LOGGING HELPERS
    // =========================================================================
    // These do not participate in any assembly; they exist purely so that
    // callers (Diagnostics, post-processing, unit tests) can log individual
    // physical quantities without duplicating the constitutive formulas
    // above or risking a permeable/impermeable inconsistency.

    // Bundle of every energy density term at a Gauss point, useful for
    // energy-balance plots such as Fig. 6 of the paper (surface energy vs.
    // load step) or for sanity-checking single-phase vs. multi-phase runs.
    struct EnergyBreakdown {
        double U;          // domain-wall energy
        double W_mec;       // electrostrictive coupling energy
        double W_elas;      // pure elastic energy
        double chi;         // Landau-Devonshire phase-separation energy
        double W_elec;      // electrostatic energy (-eps0/2|E|^2 - E.p)
        double h;           // full enthalpy density (branch-dependent)
    };

    EnergyBreakdown compute_energy_breakdown(const Eigen::Matrix2d& grad_P, const Eigen::Vector2d& P,
                                              const Eigen::Matrix2d& strain, const Eigen::Vector2d& E,
                                              double v, bool is_impermeable) const {
        EnergyBreakdown eb;
        eb.U       = U_energy(grad_P);
        eb.W_mec   = W_mec_energy(P, strain);
        eb.W_elas  = W_elas_energy(strain);
        eb.chi     = chi_energy(P);
        eb.W_elec  = electrostatic_energy(P, E);
        eb.h       = enthalpy_density(grad_P, P, strain, E, v, is_impermeable);
        return eb;
    }

    // Magnitude of the polarization vector, |P| = sqrt(p1^2+p2^2). Useful to
    // log detwinning/switching progress (cf. Figs. 13-15 of the paper).
    double polarization_magnitude(const Eigen::Vector2d& P) const {
        return P.norm();
    }

    // Angle (radians, in [-pi, pi]) of the polarization vector relative to
    // the x1 axis. Useful to track 90-degree ferroelastic domain switching.
    double polarization_angle(const Eigen::Vector2d& P) const {
        return std::atan2(P(1), P(0));
    }

    // Degradation (jump-set) factor (v^2 + eta_k) at a given fracture field
    // value, exposed for logging/consistency checks against the assemblers.
    double degradation_factor(double v) const {
        return v * v + material.eta_k;
    }

    // Full Cauchy stress sigma = degradation * (C:eps + sigma_0), in Voigt
    // notation (sigma_11, sigma_22, sigma_12). Useful to reproduce stress
    // cross-sections such as Fig. 8 of the paper.
    Eigen::Vector3d compute_stress(const Eigen::Matrix2d& strain, const Eigen::Vector2d& P, double v) const {
        Eigen::Vector3d eps_voigt;
        eps_voigt(0) = strain(0,0);
        eps_voigt(1) = strain(1,1);
        eps_voigt(2) = strain(0,1) + strain(1,0); // engineering shear gamma_12

        Eigen::Matrix3d C = get_elastic_matrix();
        Eigen::Vector3d sigma_0 = compute_sigma_0(P);

        return degradation_factor(v) * (C * eps_voigt + sigma_0);
    }

    // Full electric displacement D = -dh/dE = degradation-branch-dependent
    // (eps0*E + P). Useful to reproduce Figs. 11-12 of the paper.
    Eigen::Vector2d compute_electric_displacement(const Eigen::Vector2d& E, const Eigen::Vector2d& P,
                                                   double v, bool is_impermeable) const {
        double factor = is_impermeable ? degradation_factor(v) : 1.0;
        return factor * (material.eps0 * E + P);
    }
};