#include "Tests/include/TestFramework.h"
#include "Mesh/include/Mesh.h"
#include "Mesh/include/MeshGenerators.h"
#include "IO/include/ConfigTypes.h"

// =========================================================================
// Tests Mesh / MeshManager (generation structuree Q4).
// =========================================================================

TEST_CASE("MeshGenerators::generate_structured : nombre de noeuds/elements Q4") {
    MeshConfig cfg;
    cfg.calcul_mesh = true;
    cfg.Lx = 2.0; cfg.Ly = 1.0;
    cfg.nx = 4; cfg.ny = 2;
    cfg.dx = cfg.Lx / cfg.nx; cfg.dy = cfg.Ly / cfg.ny;
    cfg.element_type = ElementType::QUAD4;

    Mesh mesh = MeshGenerators::generate_structured(cfg);

    CHECK_TRUE(mesh.get_num_nodes() == (cfg.nx + 1) * (cfg.ny + 1));
    CHECK_TRUE(mesh.get_num_elements() == cfg.nx * cfg.ny);
    CHECK_NEAR(mesh.get_Lx(), 2.0, 1e-12);
    CHECK_NEAR(mesh.get_Ly(), 1.0, 1e-12);
}

TEST_CASE("Mesh : chaque element Q4 a 4 noeuds valides") {
    MeshConfig cfg;
    cfg.calcul_mesh = true;
    cfg.Lx = 1.0; cfg.Ly = 1.0;
    cfg.nx = 3; cfg.ny = 3;
    cfg.dx = cfg.Lx / cfg.nx; cfg.dy = cfg.Ly / cfg.ny;
    cfg.element_type = ElementType::QUAD4;

    Mesh mesh = MeshGenerators::generate_structured(cfg);

    for (int e = 0; e < mesh.get_num_elements(); ++e) {
        const Element& elem = mesh.get_elements()[e];
        CHECK_TRUE(elem.get_num_nodes() == 4);
        for (int n = 0; n < elem.get_num_nodes(); ++n) {
            CHECK_TRUE(mesh.get_node_index(e, n) >= 0 && mesh.get_node_index(e, n) < mesh.get_num_nodes());
        }
    }
}
