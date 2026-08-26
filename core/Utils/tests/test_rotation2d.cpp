#include "Tests/include/TestFramework.h"
#include "Utils/include/Rotation2D.h"
#include <Eigen/Core>

// =========================================================================
// Tests Rotation2D : to_local/to_global doivent etre inverses l'une de
// l'autre (rotation orthogonale), quelle que soit la representation
// (vecteur, tenseur 2x2, notation de Voigt, matrice de rigidite 3x3).
// =========================================================================

TEST_CASE("Rotation2D : to_global(to_local(v)) == v pour un vecteur") {
    Rotation2D rot(0.37);
    Eigen::Vector2d v(1.3, -0.7);
    Eigen::Vector2d round_trip = rot.to_global(rot.to_local(v));
    CHECK_NEAR(round_trip.x(), v.x(), 1e-12);
    CHECK_NEAR(round_trip.y(), v.y(), 1e-12);
}

TEST_CASE("Rotation2D : theta=0 est l'identite") {
    Rotation2D rot(0.0);
    Eigen::Vector2d v(0.42, -1.1);
    Eigen::Vector2d out = rot.to_local(v);
    CHECK_NEAR(out.x(), v.x(), 1e-12);
    CHECK_NEAR(out.y(), v.y(), 1e-12);
}

TEST_CASE("Rotation2D : to_global(to_local(C)) == C pour une matrice de rigidite 3x3") {
    Rotation2D rot(1.1);
    Eigen::Matrix3d C;
    C << 200.0,  80.0, 0.0,
          80.0, 200.0, 0.0,
           0.0,   0.0, 60.0;
    Eigen::Matrix3d round_trip = rot.to_global(rot.to_local(C));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            CHECK_NEAR(round_trip(i, j), C(i, j), 1e-9);
}
