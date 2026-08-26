#include "Tests/include/TestFramework.h"
#include "Tests/include/TestFixtures.h"
#include "Materials/Models/Ferroelectric.h"
#include <algorithm>
#include <cmath>
#include <functional>

// =========================================================================
// Niveau 0 : Verification mathematique des derivees analytiques (forces,
// Jacobien, contrainte) par comparaison a des differences finies.
// Ported depuis Validation.cpp (le fichier d'origine etait entierement
// commente et referencait une classe "Math" obsolete ; la logique est
// reprise ici contre l'API publique actuelle, SingleCrystalMaterial /
// MaterialModel, dont les methodes correspondent 1:1 a l'ancien "math.*").
// =========================================================================

namespace {

constexpr double kFdEpsilon = 1e-6;
constexpr double kFdEpsilon2nd = 1e-4;
constexpr double kDefaultTol = 1e-4;

double finite_diff_dP(const std::function<double(const Eigen::Vector2d&)>& f,
                       const Eigen::Vector2d& P, int component) {
    Eigen::Vector2d P_plus = P, P_minus = P;
    P_plus(component)  += kFdEpsilon;
    P_minus(component) -= kFdEpsilon;
    return (f(P_plus) - f(P_minus)) / (2.0 * kFdEpsilon);
}

double finite_diff_d2P(const std::function<double(const Eigen::Vector2d&)>& f,
                        const Eigen::Vector2d& P) {
    const double h = kFdEpsilon2nd;
    Eigen::Vector2d Pxy_pp = P, Pxy_pm = P, Pxy_mp = P, Pxy_mm = P;
    Pxy_pp(0) += h; Pxy_pp(1) += h;
    Pxy_pm(0) += h; Pxy_pm(1) -= h;
    Pxy_mp(0) -= h; Pxy_mp(1) += h;
    Pxy_mm(0) -= h; Pxy_mm(1) -= h;

    return (f(Pxy_pp) - f(Pxy_pm) - f(Pxy_mp) + f(Pxy_mm)) / (4.0 * h * h);
}

} // namespace

TEST_CASE("Validation Niv0 : dchi_dp1 (Landau-Devonshire, strain=0)") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Vector2d P(0.6, 0.3);
    auto f = [&](const Eigen::Vector2d& p) { return material.chi_energy(p); };
    double numeric = finite_diff_dP(f, P, 0);

    Eigen::Matrix2d zero_strain = Eigen::Matrix2d::Zero();
    Eigen::Vector2d zero_E = Eigen::Vector2d::Zero();
    auto GL = material.compute_GL_terms(P, zero_strain, zero_E, 1.0, false);

    CHECK_NEAR(GL.force_px, numeric, kDefaultTol);
}

TEST_CASE("Validation Niv0 : dchi_dp2 (Landau-Devonshire, strain=0)") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Vector2d P(0.6, 0.3);
    auto f = [&](const Eigen::Vector2d& p) { return material.chi_energy(p); };
    double numeric = finite_diff_dP(f, P, 1);

    Eigen::Matrix2d zero_strain = Eigen::Matrix2d::Zero();
    Eigen::Vector2d zero_E = Eigen::Vector2d::Zero();
    auto GL = material.compute_GL_terms(P, zero_strain, zero_E, 1.0, false);

    CHECK_NEAR(GL.force_py, numeric, kDefaultTol);
}

TEST_CASE("Validation Niv0 : d2chi_dp1dp2 (Jacobien croise, strain=0)") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Vector2d P(0.5, 0.4);
    auto f = [&](const Eigen::Vector2d& p) { return material.chi_energy(p); };
    double numeric = finite_diff_d2P(f, P);

    Eigen::Matrix2d zero_strain = Eigen::Matrix2d::Zero();
    Eigen::Vector2d zero_E = Eigen::Vector2d::Zero();
    auto GL = material.compute_GL_terms(P, zero_strain, zero_E, 1.0, false);

    CHECK_NEAR(GL.J_12, numeric, kDefaultTol);
}

TEST_CASE("Validation Niv0 : dW_dp1 (couplage electrostrictif P <-> Eps)") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Vector2d P(0.5, 0.3);
    Eigen::Matrix2d strain; strain << 0.01, 0.005, 0.005, -0.008;

    auto f = [&](const Eigen::Vector2d& p) { return material.W_energy(p, strain); };
    double numeric = finite_diff_dP(f, P, 0);

    Eigen::Vector2d zero_E = Eigen::Vector2d::Zero();
    auto GL = material.compute_GL_terms(P, strain, zero_E, 1.0, false);
    auto GL_chi_only = material.compute_GL_terms(P, Eigen::Matrix2d::Zero(), zero_E, 1.0, false);
    double analytic_dW_dp1 = GL.force_px - GL_chi_only.force_px;

    CHECK_NEAR(analytic_dW_dp1, numeric, kDefaultTol);
}

TEST_CASE("Validation Niv0 : dW_dp2 (couplage electrostrictif P <-> Eps)") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Vector2d P(0.5, 0.3);
    Eigen::Matrix2d strain; strain << 0.01, 0.005, 0.005, -0.008;

    auto f = [&](const Eigen::Vector2d& p) { return material.W_energy(p, strain); };
    double numeric = finite_diff_dP(f, P, 1);

    Eigen::Vector2d zero_E = Eigen::Vector2d::Zero();
    auto GL = material.compute_GL_terms(P, strain, zero_E, 1.0, false);
    auto GL_chi_only = material.compute_GL_terms(P, Eigen::Matrix2d::Zero(), zero_E, 1.0, false);
    double analytic_dW_dp2 = GL.force_py - GL_chi_only.force_py;

    CHECK_NEAR(analytic_dW_dp2, numeric, kDefaultTol);
}

TEST_CASE("Validation Niv0 : d2W_dp1dp2 (Jacobien croise electrostrictif)") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Vector2d P(0.4, 0.4);
    Eigen::Matrix2d strain; strain << 0.01, 0.006, 0.006, -0.004;

    auto f = [&](const Eigen::Vector2d& p) { return material.W_energy(p, strain); };
    double numeric = finite_diff_d2P(f, P);

    Eigen::Vector2d zero_E = Eigen::Vector2d::Zero();
    auto GL = material.compute_GL_terms(P, strain, zero_E, 1.0, false);
    auto GL_chi_only = material.compute_GL_terms(P, Eigen::Matrix2d::Zero(), zero_E, 1.0, false);
    double analytic = GL.J_12 - GL_chi_only.J_12;

    CHECK_NEAR(analytic, numeric, kDefaultTol);
}

TEST_CASE("Validation Niv0 : sigma_0 conjugue a W_mec (test thermodynamique)") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Vector2d P(0.5, 0.3);
    Eigen::Vector2d P_zero(0.0, 0.0);
    double eps11 = 0.02, eps22 = -0.015, eps12 = 0.008;

    auto W_mec_only = [&](double e11, double e22, double e12) {
        Eigen::Matrix2d strain; strain << e11, e12, e12, e22;
        return material.W_energy(P, strain) - material.W_energy(P_zero, strain);
    };

    double dW_deps11_num = (W_mec_only(eps11 + kFdEpsilon, eps22, eps12) - W_mec_only(eps11 - kFdEpsilon, eps22, eps12)) / (2.0 * kFdEpsilon);
    double dW_deps22_num = (W_mec_only(eps11, eps22 + kFdEpsilon, eps12) - W_mec_only(eps11, eps22 - kFdEpsilon, eps12)) / (2.0 * kFdEpsilon);
    double dW_deps12_num = (W_mec_only(eps11, eps22, eps12 + kFdEpsilon) - W_mec_only(eps11, eps22, eps12 - kFdEpsilon)) / (2.0 * kFdEpsilon);

    Eigen::Vector3d sigma_0 = material.compute_sigma_0(P);
    double max_err = std::max({std::abs(sigma_0(0) - dW_deps11_num), std::abs(sigma_0(1) - dW_deps22_num), std::abs(sigma_0(2) - 0.5 * dW_deps12_num)});
    double scale = std::max({std::abs(dW_deps11_num), std::abs(dW_deps22_num), std::abs(dW_deps12_num), 1e-12});

    CHECK_NEAR(max_err / scale, 0.0, kDefaultTol);
}

TEST_CASE("Validation Niv0 : H_drive strictement positif (critere d'irreversibilite)") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Matrix2d grad_P; grad_P << 0.1, -0.05, 0.02, 0.08;
    Eigen::Vector2d P(0.7, 0.2);
    Eigen::Matrix2d strain; strain << 0.01, 0.002, 0.002, -0.005;
    Eigen::Vector2d E(0.01, -0.02);

    double H_perm = material.compute_H_drive(grad_P, P, strain, E, false);
    double H_imperm = material.compute_H_drive(grad_P, P, strain, E, true);

    CHECK_TRUE(H_perm >= 0.0);
    CHECK_TRUE(H_imperm >= 0.0);
}

TEST_CASE("Validation Niv0 : H_drive nul a l'equilibre homogene") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(cfg);

    Eigen::Vector2d P(1.0, 0.0);
    double H = material.compute_H_drive(Eigen::Matrix2d::Zero(), P, Eigen::Matrix2d::Zero(), Eigen::Vector2d::Zero(), false);

    CHECK_NEAR(H, 0.0, 1e-10);
}
