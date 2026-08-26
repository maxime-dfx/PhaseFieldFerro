#include "Tests/include/TestFramework.h"
#include "Mesh/include/Mesh.h"
#include "Mesh/include/MeshGenerators.h"
#include "Physics/Core/BoundaryManager.h"
#include "IO/include/ConfigTypes.h"
#include <memory>

// =========================================================================
// Tests BoundaryManager (regles de bord programmatiques, sans passer par le
// parsing TOML de Datafile - cf IO/tests pour ce dernier).
// =========================================================================

TEST_CASE("BoundaryManager : une regle Dirichlet constante sur un bord peuple bien get_bcs(phi)") {
    MeshConfig cfg;
    cfg.calcul_mesh = true;
    cfg.Lx = 1.0; cfg.Ly = 1.0;
    cfg.nx = 4; cfg.ny = 4;
    cfg.dx = cfg.Lx / cfg.nx; cfg.dy = cfg.Ly / cfg.ny;
    cfg.element_type = ElementType::QUAD4;
    Mesh mesh = MeshGenerators::generate_structured(cfg);

    std::vector<BoundaryRuleConfig> empty_rules;
    BoundaryManager bc(mesh, empty_rules);

    auto left_edge = std::make_shared<EdgeShape>("left", cfg.Lx, cfg.Ly);
    bc.add_rule_constant("phi", left_edge, BCType::DIRICHLET, 1.0, 0.0);
    bc.initialize_all_boundaries();

    // (ny+1) noeuds sur le bord gauche d'un maillage structure nx x ny.
    CHECK_TRUE(bc.get_bcs("phi").size() == static_cast<std::size_t>(mesh.get_num_nodes()));

    int num_dirichlet = 0;
    for (const auto& node_bc : bc.get_bcs("phi")) {
        if (node_bc.type == BCType::DIRICHLET) ++num_dirichlet;
    }
    CHECK_TRUE(num_dirichlet == cfg.ny + 1);
}

TEST_CASE("BoundaryManager : EdgeShape/RectShape/CircleShape - contains() coherent") {
    EdgeShape left("left", 2.0, 2.0);
    CHECK_TRUE(left.contains(0.0, 1.0));
    CHECK_TRUE(!left.contains(2.0, 1.0));

    RectShape rect(0.0, 1.0, 0.0, 1.0);
    CHECK_TRUE(rect.contains(0.5, 0.5));
    CHECK_TRUE(!rect.contains(1.5, 0.5));

    CircleShape circle(0.0, 0.0, 1.0);
    CHECK_TRUE(circle.contains(0.5, 0.5));
    CHECK_TRUE(!circle.contains(1.0, 1.0));
}
