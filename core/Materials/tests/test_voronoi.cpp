#include "Tests/include/TestFramework.h"
#include "Materials/Microstructure/Voronoi.h"

// =========================================================================
// Tests VoronoiTessellation : reproductibilite du tirage aleatoire.
// Ported depuis TestSuite::run_section_D_voronoi.
// =========================================================================

TEST_CASE("VoronoiTessellation : meme seed => meme microstructure") {
    VoronoiTessellation vt1(10.0, 10.0, 5, 4242);
    VoronoiTessellation vt2(10.0, 10.0, 5, 4242);

    CHECK_TRUE(vt1.grains()[0].seed_x == vt2.grains()[0].seed_x);
    CHECK_TRUE(vt1.grains()[0].theta == vt2.grains()[0].theta);
}

TEST_CASE("VoronoiTessellation : seed differente => microstructure differente") {
    VoronoiTessellation vt1(10.0, 10.0, 5, 4242);
    VoronoiTessellation vt_diff(10.0, 10.0, 5, 9999);

    CHECK_TRUE(vt1.grains()[0].seed_x != vt_diff.grains()[0].seed_x);
}
