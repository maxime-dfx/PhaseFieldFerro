#include "Tests/include/TestFramework.h"
#include "Tests/include/TestFixtures.h"
#include "Materials/Models/Ferroelectric.h"
#include <algorithm>
#include <cmath>

// =========================================================================
// Niveau 1 : dynamique locale 0D (relaxation de polarisation par descente
// de gradient adaptative, sans maillage) + patch test mecanique.
// Ported depuis Validation.cpp (voir note dans
// test_validation_level0_derivatives.cpp sur l'origine de ce portage).
// =========================================================================

namespace {

// Descente de gradient adaptative (line search a pas decroissant) sur
// l'enthalpie W + chi - P.E, utilisee pour relaxer P vers un minimum
// local a strain/E fixes.
Eigen::Vector2d integrate_polarization_point(const SingleCrystalMaterial& material,
                                              Eigen::Vector2d P, const Eigen::Matrix2d& strain,
                                              const Eigen::Vector2d& E, int n_steps) {
    auto enthalpy = [&](const Eigen::Vector2d& p_val) {
        return material.W_energy(p_val, strain) + material.chi_energy(p_val) - p_val(0) * E(0) - p_val(1) * E(1);
    };

    double alpha = 0.05;

    for (int step = 0; step < n_steps; ++step) {
        auto GL = material.compute_GL_terms(P, strain, E, 1.0, false);
        Eigen::Vector2d grad(GL.force_px, GL.force_py);

        if (grad.norm() < 1e-8) break;

        double current_h = enthalpy(P);
        bool step_accepted = false;

        for (int ls = 0; ls < 20; ++ls) {
            Eigen::Vector2d step_vec = alpha * grad;

            // Securite critique : on borne la taille du pas pour eviter
            // l'explosion de P^8 (termes de Landau d'ordre eleve).
            if (step_vec.norm() > 0.1) step_vec = step_vec.normalized() * 0.1;

            Eigen::Vector2d P_new = P - step_vec;
            double new_h = enthalpy(P_new);

            if (new_h < current_h) {
                P = P_new;
                step_accepted = true;
                break;
            }
            alpha *= 0.5;
        }

        if (!step_accepted) break;
        alpha *= 1.2;
        if (alpha > 1.0) alpha = 1.0;
    }
    return P;
}

} // namespace

TEST_CASE("Validation Niv1 : dynamique 0D - relaxation vers un equilibre stable") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Vector2d P0(0.01, 0.001);
    Eigen::Vector2d P_final = integrate_polarization_point(material, P0, Eigen::Matrix2d::Zero(), Eigen::Vector2d::Zero(), 20000);

    auto GL_final = material.compute_GL_terms(P_final, Eigen::Matrix2d::Zero(), Eigen::Vector2d::Zero(), 1.0, false);
    double force_norm = std::sqrt(GL_final.force_px * GL_final.force_px + GL_final.force_py * GL_final.force_py);
    double max_stiffness = std::max({std::abs(GL_final.J_11), std::abs(GL_final.J_22), 1e-6});

    CHECK_NEAR(force_norm / max_stiffness, 0.0, 1e-4);
}

TEST_CASE("Validation Niv1 : dynamique 0D - stabilite locale du puits de potentiel") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Vector2d P_min = integrate_polarization_point(material, Eigen::Vector2d(0.01, 0.0), Eigen::Matrix2d::Zero(), Eigen::Vector2d::Zero(), 20000);

    Eigen::Vector2d P_perturbed = P_min;
    P_perturbed(0) *= 0.8;
    P_perturbed(1) += 0.05;

    // L'algorithme doit retrouver (quasi) exactement le meme point.
    Eigen::Vector2d P_recovered = integrate_polarization_point(material, P_perturbed, Eigen::Matrix2d::Zero(), Eigen::Vector2d::Zero(), 20000);
    double drift = (P_recovered - P_min).norm() / std::max(1e-6, P_min.norm());

    CHECK_NEAR(drift, 0.0, 1e-4);
}

TEST_CASE("Validation Niv1 : dynamique 0D - basculement ferroelastique sous traction") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Vector2d P = integrate_polarization_point(material, Eigen::Vector2d(0.01, 0.001), Eigen::Matrix2d::Zero(), Eigen::Vector2d::Zero(), 20000);
    bool switched = false;

    // Rampe de traction selon Y pour forcer le basculement a 90 degres.
    for (int i = 1; i <= 200; ++i) {
        Eigen::Matrix2d strain; strain << 0.0, 0.0, 0.0, i * 0.005;
        P = integrate_polarization_point(material, P, strain, Eigen::Vector2d::Zero(), 2000);

        if (std::abs(P(1)) > 1.5 * std::abs(P(0))) {
            switched = true;
            break;
        }
    }

    CHECK_TRUE(switched);
}

TEST_CASE("Validation Niv1 : mecanique - patch test W_elas avec matrice C") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    double a = 0.01, b = -0.005;
    Eigen::Matrix3d C = material.get_elastic_matrix();
    Eigen::Vector3d eps_voigt(a, b, 0.0);
    double W_analytic = 0.5 * eps_voigt.transpose() * C * eps_voigt;

    Eigen::Matrix2d strain; strain << a, 0.0, 0.0, b;
    double W_numeric = material.W_energy(Eigen::Vector2d::Zero(), strain);

    CHECK_NEAR(W_analytic, W_numeric, 1e-8);
}
