#include "Tests/include/TestFramework.h"
#include "Tests/include/TestFixtures.h"
#include "Materials/Models/Ferroelectric.h"
#include <algorithm>
#include <cmath>

// =========================================================================
// Niveau 2 : couplage de la fracture (champ de phase v), toujours en 0D
// (aucun maillage requis). Ported depuis Validation.cpp (voir note dans
// test_validation_level0_derivatives.cpp).
// =========================================================================

namespace {
constexpr double kFdEpsilon = 1e-6;
constexpr double kDefaultTol = 1e-4;
} // namespace

TEST_CASE("Validation Niv2 : fracture - contrainte liberee (traction-free) pour v=0") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Vector2d P(0.5, 0.3);
    Eigen::Matrix2d strain; strain << 0.01, 0.005, 0.005, -0.008;
    Eigen::Vector2d E = Eigen::Vector2d::Zero();

    auto GL_intact = material.compute_GL_terms(P, strain, E, 1.0, false);
    auto GL_broken = material.compute_GL_terms(P, strain, E, material.get_eta_k(), false);
    auto GL_chi_only = material.compute_GL_terms(P, Eigen::Matrix2d::Zero(), E, 0.0, false);

    double J11_mec_intact = GL_intact.J_11 - GL_chi_only.J_11;
    double J11_mec_broken = GL_broken.J_11 - GL_chi_only.J_11;

    double ratio = std::abs(J11_mec_broken) / std::max(1e-12, std::abs(J11_mec_intact));

    CHECK_NEAR(ratio, material.get_eta_k(), 1e-10);
}

TEST_CASE("Validation Niv2 : fracture - champ D (condition permeable vs impermeable)") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Vector2d P(0.0, 0.5);
    Eigen::Vector2d E(0.0, 0.01);

    Eigen::Vector2d D_perm = material.compute_effective_permittivity(0.0, material.get_eta_k(), false) * E
                            + material.compute_effective_polarization(P, 0.0, material.get_eta_k(), false);

    Eigen::Vector2d D_imperm = material.compute_effective_permittivity(0.0, material.get_eta_k(), true) * E
                              + material.compute_effective_polarization(P, 0.0, material.get_eta_k(), true);

    double expected_perm = (material.get_eps0() * E + P).norm();
    double expected_imperm = material.get_eta_k() * expected_perm;

    double err_perm = std::abs(D_perm.norm() - expected_perm);
    double err_imperm = std::abs(D_imperm.norm() - expected_imperm);

    CHECK_TRUE(err_perm < 1e-10);
    CHECK_TRUE(err_imperm < 1e-10);
}

TEST_CASE("Validation Niv2 : fracture - force thermodynamique dh/dv (H_drive)") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Matrix2d grad_P; grad_P << 0.1, -0.05, 0.02, 0.08;
    Eigen::Vector2d P(0.4, 0.4);
    Eigen::Matrix2d strain; strain << 0.01, 0.0, 0.0, -0.01;
    Eigen::Vector2d E(0.01, 0.01);
    double v_test = 0.5;

    auto enthalpy = [&](double v_val, bool imperm) {
        double penalty = v_val * v_val + material.get_eta_k();
        double U = material.U_energy(grad_P);
        double W = material.W_energy(P, strain);
        double chi = material.chi_energy(P);
        double W_elec_base = -P.dot(E);
        double W_elec_imperm = W_elec_base - 0.5 * material.get_eps0() * E.squaredNorm();

        return imperm ? penalty * (U + W + W_elec_imperm) + chi
                      : penalty * (U + W) + chi + W_elec_base - 0.5 * material.get_eps0() * E.squaredNorm();
    };

    double dh_dv_perm_num = (enthalpy(v_test + kFdEpsilon, false) - enthalpy(v_test - kFdEpsilon, false)) / (2.0 * kFdEpsilon);
    double dh_dv_imperm_num = (enthalpy(v_test + kFdEpsilon, true) - enthalpy(v_test - kFdEpsilon, true)) / (2.0 * kFdEpsilon);

    double dh_dv_perm_analy = 2.0 * v_test * material.compute_H_drive(grad_P, P, strain, E, false);
    double dh_dv_imperm_analy = 2.0 * v_test * material.compute_H_drive(grad_P, P, strain, E, true);

    double err_perm = std::abs(dh_dv_perm_num - dh_dv_perm_analy) / std::max(1e-12, std::abs(dh_dv_perm_num));
    double err_imperm = std::abs(dh_dv_imperm_num - dh_dv_imperm_analy) / std::max(1e-12, std::abs(dh_dv_imperm_num));

    CHECK_TRUE(err_perm < kDefaultTol);
    CHECK_TRUE(err_imperm < kDefaultTol);
}
