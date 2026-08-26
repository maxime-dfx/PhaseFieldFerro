#include "Tests/include/TestFramework.h"
#include "Mesh/include/FEM.h"
#include <vector>
#include <utility>

// =========================================================================
// Tests FEM (fonctions de forme + quadrature) - proprietes independantes
// de tout maillage : partition de l'unite, poids de quadrature.
// =========================================================================

TEST_CASE("ShapeFunctions Q4 : partition de l'unite en plusieurs points") {
    const std::vector<std::pair<double,double>> pts = {{0.0,0.0}, {0.5,-0.3}, {-1.0,1.0}, {0.7,0.7}};
    for (const auto& [xi, eta] : pts) {
        auto N = ShapeFunctions::get_shape_functions(xi, eta);
        double sum = N[0] + N[1] + N[2] + N[3];
        CHECK_NEAR(sum, 1.0, 1e-12);
    }
}

TEST_CASE("ShapeFunctions Q4 : valeur nodale = delta de Kronecker aux sommets de reference") {
    // Sommet 0 = (-1,-1) : seule N0 doit valoir 1, les autres 0.
    auto N = ShapeFunctions::get_shape_functions(-1.0, -1.0);
    CHECK_NEAR(N[0], 1.0, 1e-12);
    CHECK_NEAR(N[1], 0.0, 1e-12);
    CHECK_NEAR(N[2], 0.0, 1e-12);
    CHECK_NEAR(N[3], 0.0, 1e-12);
}

TEST_CASE("ShapeFunctions T3 : partition de l'unite") {
    auto N = ShapeFunctions::get_shape_functions_tri(0.3, 0.4);
    CHECK_NEAR(N[0] + N[1] + N[2], 1.0, 1e-12);
}

TEST_CASE("Quadrature : la somme des poids Gauss 2x2 vaut l'aire du carre de reference (4)") {
    double sum_w = 0.0;
    for (const auto& gp : Quadrature::get_gauss_2x2()) sum_w += gp.weight;
    CHECK_NEAR(sum_w, 4.0, 1e-12);
}

TEST_CASE("Quadrature : la somme des poids Gauss 3x3 vaut l'aire du carre de reference (4)") {
    double sum_w = 0.0;
    for (const auto& gp : Quadrature::get_gauss_3x3()) sum_w += gp.weight;
    CHECK_NEAR(sum_w, 4.0, 1e-12);
}

TEST_CASE("Quadrature : la somme des poids Gauss triangle vaut l'aire du triangle de reference (0.5)") {
    double sum_w = 0.0;
    for (const auto& gp : Quadrature::get_gauss_3_points_tri()) sum_w += gp.weight;
    CHECK_NEAR(sum_w, 0.5, 1e-12);
}
