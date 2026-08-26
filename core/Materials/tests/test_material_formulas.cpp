#include "Tests/include/TestFramework.h"
#include "Tests/include/TestFixtures.h"
#include "Materials/Models/Ferroelectric.h"
#include <cmath>

// =========================================================================
// Tests des formules de base de SingleCrystalMaterial
// -------------------------------------------------------------------------
// Verifications minimales de sante physique/numerique (positivite de U,
// absence de NaN dans les termes de Ginzburg-Landau). Ported depuis
// TestSuite::run_section_A_material_formulas, en passant par l'API
// publique (SingleCrystalMaterial) plutot que la classe interne
// FerroelectricMaterial.
// =========================================================================

TEST_CASE("SingleCrystalMaterial : U_energy est strictement positive") {
    MaterialConfig config = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial mat(config);

    Eigen::Matrix2d grad_P;
    grad_P << 0.1, 0.0,
              0.0, 0.1;

    CHECK_TRUE(mat.U_energy(grad_P) > 0.0);
}

TEST_CASE("SingleCrystalMaterial : GL_terms sont finis (pas de NaN)") {
    MaterialConfig config = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial mat(config);

    Eigen::Vector2d P(0.1, 0.05);
    Eigen::Matrix2d strain;
    strain << 0.01, 0.005,
              0.005, -0.01;
    Eigen::Vector2d E(1.0, 0.0);

    auto gl = mat.compute_GL_terms(P, strain, E, 1.0, false);
    CHECK_TRUE(std::isfinite(gl.J_11));
    CHECK_TRUE(std::isfinite(gl.force_px));
}
