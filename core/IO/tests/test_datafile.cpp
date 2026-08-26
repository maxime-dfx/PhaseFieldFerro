#include "Tests/include/TestFramework.h"
#include "IO/include/Datafile.h"
#include <fstream>
#include <cstdio>

// =========================================================================
// Tests Datafile (parsing TOML -> ConfigTypes).
// =========================================================================

TEST_CASE("Datafile : parse correctement une configuration minimale") {
    std::string path = "test_io_datafile_minimal.toml";
    std::ofstream out(path);
    out << "[simulation]\noutput_dir = \".\"\ntotal_time = 1\ndt = 0.1\n"
        << "[mesh]\nL_x=2.0\nL_y=3.0\nn_x=4\nn_y=5\nelement_type=\"Q4\"\n"
        << "[material]\n[crystal]\nnum_grains=1\n";
    out.close();

    Datafile config(path);
    std::remove(path.c_str());

    CHECK_NEAR(config.mesh.Lx, 2.0, 1e-12);
    CHECK_NEAR(config.mesh.Ly, 3.0, 1e-12);
    CHECK_TRUE(config.mesh.nx == 4);
    CHECK_TRUE(config.mesh.ny == 5);
    CHECK_TRUE(config.crystal.num_grains == 1);
}