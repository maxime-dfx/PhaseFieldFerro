#include "Tests/include/TestFramework.h"
#include "Tests/include/TestFixtures.h"
#include "Materials/Models/Ferroelectric.h"
#include <Eigen/Core>
#include <cmath>

// =========================================================================
// Verification des parametres (Arias 2011) et couplages physiques.
// Migre depuis l'ancien main.cpp pour ne pas polluer la production.
// =========================================================================

TEST_CASE("Validation Parametres : Coefficients du papier (Arias 2011)") {
    MaterialConfig mat = TestFixtures::make_paper_material_config();
    
    // Parametres de Fracture et Phase-Field
    CHECK_NEAR(mat.Gc, 4.0, 1e-5);
    CHECK_NEAR(mat.kappa, 2.0, 1e-5);
    CHECK_NEAR(mat.eta_k, 1e-6, 1e-5);
    CHECK_NEAR(mat.a0, 0.1, 1e-5);
    CHECK_NEAR(mat.mu_p, 1.0, 1e-5);
    CHECK_NEAR(mat.mu_v, 15.0, 1e-5);

    // Coefficients de Landau
    CHECK_NEAR(mat.c1, 185.0, 1e-5);
    CHECK_NEAR(mat.c2, 111.0, 1e-5);
    CHECK_NEAR(mat.c3, 74.0, 1e-5);
    CHECK_NEAR(mat.b1, 1.4282, 1e-5);
    CHECK_NEAR(mat.b2, -0.185, 1e-5);
    CHECK_NEAR(mat.b3, 0.8066, 1e-5);
    CHECK_NEAR(mat.alpha_1, -0.0023, 1e-5);
    CHECK_NEAR(mat.alpha_11, -0.0029, 1e-5);
    CHECK_NEAR(mat.alpha_12, -0.0011, 1e-5);
    CHECK_NEAR(mat.alpha_111, 0.003, 1e-5);
    CHECK_NEAR(mat.alpha_1111, 0.001, 1e-5);
    CHECK_NEAR(mat.alpha_1122, 1.24, 1e-5);
    CHECK_NEAR(mat.eps0, 0.131, 1e-5);
}

TEST_CASE("Validation Multiphysique : Forces et contraintes (Differences Finies)") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Vector2d P(0.8, 0.3);
    Eigen::Matrix2d strain;
    strain << 0.01, -0.005,
             -0.005, 0.02;
    Eigen::Vector2d E(0.001, -0.002);
    
    double v = 0.5;
    double penalite = (v * v) + material.get_eta_k();
    bool is_impermeable = true;
    double delta = 1e-6;
    double tol = 1e-4;

    // Fonction Enthalpie Totale
    auto compute_H_tot = [&](const Eigen::Vector2d& p) {
        double W = material.W_energy(p, strain);
        double chi = material.chi_energy(p);
        double W_elec = -p.dot(E);
        if (is_impermeable) {
            return penalite * (W + W_elec) + chi;
        }
        return penalite * (W + chi) + W_elec;
    };

    // 1. Forces Ginzburg-Landau
    double dH_dp1_FD = (compute_H_tot(P + Eigen::Vector2d(delta, 0)) - compute_H_tot(P - Eigen::Vector2d(delta, 0))) / (2.0 * delta);
    double dH_dp2_FD = (compute_H_tot(P + Eigen::Vector2d(0, delta)) - compute_H_tot(P - Eigen::Vector2d(0, delta))) / (2.0 * delta);
    GinzburgLandauTerms GL = material.compute_GL_terms(P, strain, E, penalite, is_impermeable);
    
    CHECK_NEAR(dH_dp1_FD, GL.force_px, tol);
    CHECK_NEAR(dH_dp2_FD, GL.force_py, tol);

    // 2. Jacobien Ginzburg-Landau
    auto GL_p1_plus  = material.compute_GL_terms(P + Eigen::Vector2d(delta, 0), strain, E, penalite, is_impermeable);
    auto GL_p1_minus = material.compute_GL_terms(P - Eigen::Vector2d(delta, 0), strain, E, penalite, is_impermeable);
    double J11_FD = (GL_p1_plus.force_px - GL_p1_minus.force_px) / (2.0 * delta);
    double J12_FD = (GL_p1_plus.force_py - GL_p1_minus.force_py) / (2.0 * delta);
    
    auto GL_p2_plus  = material.compute_GL_terms(P + Eigen::Vector2d(0, delta), strain, E, penalite, is_impermeable);
    auto GL_p2_minus = material.compute_GL_terms(P - Eigen::Vector2d(0, delta), strain, E, penalite, is_impermeable);
    double J22_FD = (GL_p2_plus.force_py - GL_p2_minus.force_py) / (2.0 * delta);

    CHECK_NEAR(J11_FD, GL.J_11, tol);
    CHECK_NEAR(J22_FD, GL.J_22, tol);
    CHECK_NEAR(J12_FD, GL.J_12, tol);

    // 3. Contrainte Elastique
    Eigen::Matrix2d strain_eps11_plus = strain; strain_eps11_plus(0,0) += delta;
    Eigen::Matrix2d strain_eps11_minus = strain; strain_eps11_minus(0,0) -= delta;
    double dW_deps11_FD = (material.W_energy(P, strain_eps11_plus) - material.W_energy(P, strain_eps11_minus)) / (2.0 * delta);
    
    Eigen::Vector3d sigma_0 = material.compute_sigma_0(P);
    Eigen::Matrix3d C = material.get_elastic_matrix();
    double sigma_11_total_analytic = sigma_0(0) + C(0,0) * strain(0,0) + C(0,1) * strain(1,1);
    
    CHECK_NEAR(dW_deps11_FD, sigma_11_total_analytic, tol);
}