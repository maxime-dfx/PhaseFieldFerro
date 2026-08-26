#include "Tests/include/TestFramework.h"
#include "Tests/include/TestFixtures.h"
#include "Materials/Microstructure/Polycrystal.h"
#include "Mesh/include/Mesh.h"
#include "Mesh/include/MeshGenerators.h"

// =========================================================================
// Tests Polycrystal (mode monocristal)
// Ported depuis TestSuite::run_section_C_materials_manager.
// NOTE : l'ancienne classe unifiee MaterialsManager (materiau + grains) a
// ete scindee en MaterialRegistry (materiaux) et Polycrystal
// (microstructure). Ces tests portent donc desormais sur Polycrystal
// directement, qui est le seul des deux a exposer num_grains()/
// grain_id_for_element().
// =========================================================================

namespace {
    // Fonction utilitaire pour générer un vrai maillage 20x20 pour les tests
    Mesh make_test_mesh() {
        MeshConfig mesh_cfg;
        mesh_cfg.calcul_mesh = true;
        mesh_cfg.Lx = 1.0; 
        mesh_cfg.Ly = 1.0;
        mesh_cfg.nx = 20; 
        mesh_cfg.ny = 20;
        mesh_cfg.dx = mesh_cfg.Lx / mesh_cfg.nx; 
        mesh_cfg.dy = mesh_cfg.Ly / mesh_cfg.ny;
        mesh_cfg.element_type = ElementType::QUAD4;
        return MeshGenerators::generate_structured(mesh_cfg);
    }
}

TEST_CASE("Polycrystal : mode monocristal (num_grains=1) ne cree aucun grain explicite") {
    CrystalConfig crystal_cfg; // num_grains = 1 par defaut -> monocristal
    Mesh mesh = make_test_mesh();
    Polycrystal poly(crystal_cfg, mesh);

    CHECK_TRUE(poly.num_grains() == 0);
}

TEST_CASE("Polycrystal : grain_id_for_element renvoie 0 en monocristal") {
    CrystalConfig crystal_cfg;
    Mesh mesh = make_test_mesh();

    Polycrystal poly(crystal_cfg, mesh);

    CHECK_TRUE(poly.grain_id_for_element(100) == 0);
}

TEST_CASE("Polycrystal : en monocristal, aucun element n'est marque joint de grain") {
    CrystalConfig crystal_cfg; // num_grains = 1 -> pas de joint possible
    Mesh mesh = make_test_mesh();
    Polycrystal poly(crystal_cfg, mesh);

    for (int e = 0; e < mesh.get_num_elements(); ++e) {
        CHECK_TRUE(!poly.is_grain_boundary_element(e));
    }
}

TEST_CASE("Polycrystal : en polycristal, certains elements sont marques joint de grain") {

    CrystalConfig crystal_cfg;
    crystal_cfg.num_grains = 6;
    crystal_cfg.use_random_seed = false;
    crystal_cfg.seed = 42;
    Mesh mesh = make_test_mesh();
    Polycrystal poly(crystal_cfg, mesh);

    int gb_count = 0;
    for (int e = 0; e < mesh.get_num_elements(); ++e) {
        if (poly.is_grain_boundary_element(e)) ++gb_count;
    }
    CHECK_TRUE(gb_count > 0);

    for (int e = 0; e < mesh.get_num_elements(); ++e) {
        if (!poly.is_grain_boundary_element(e)) continue;
        const Element& elem = mesh.get_elements()[static_cast<std::size_t>(e)];
        const int g = poly.grain_id_for_element(e);
        bool found_foreign_neighbor = false;
        for (int i = 0; i < elem.get_num_nodes() && !found_foreign_neighbor; ++i) {
            const int idx = mesh.get_node_index(e, i);
            for (int e2 = 0; e2 < mesh.get_num_elements(); ++e2) {
                if (e2 == e) continue;
                const Element& elem2 = mesh.get_elements()[static_cast<std::size_t>(e2)];
                for (int i2 = 0; i2 < elem2.get_num_nodes(); ++i2) {
                    if (mesh.get_node_index(e2, i2) == idx && poly.grain_id_for_element(e2) != g) {
                        found_foreign_neighbor = true;
                        break;
                    }
                }
                if (found_foreign_neighbor) break;
            }
        }
        CHECK_TRUE(found_foreign_neighbor);
        break; 
    }
}