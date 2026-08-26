#include "Tests/include/TestFramework.h"
#include "Core/include/SimulationManager.h"
#include "IO/include/Datafile.h"
#include "Mesh/include/Mesh.h"
#include "Mesh/include/MeshGenerators.h"
#include "IO/include/ResultsExporter.h"
#include "Materials/Core/MaterialManager.h"
#include "Materials/Models/Ferroelectric.h"
#include <fstream>
#include <cstdio>
#include <memory>

// =========================================================================
// Tests SimulationManager : instanciation de bout en bout (equivalent de
// ce que fait Core/src/main.cpp), sans executer la boucle temporelle
// complete (couverte plus tard par les tests globaux dans Tests/, avec
// donnees experimentales/reference).
// =========================================================================

TEST_CASE("SimulationManager : instanciation + initialize_physics() ne crashe pas") {
    std::string path = "test_simulation_manager_minimal.toml";
    std::ofstream out(path);
    out << "[simulation]\noutput_dir = \".\"\ntotal_time = 1\ndt = 0.1\n"
        << "[mesh]\nLx=1.0\nLy=1.0\nnx=2\nny=2\nelement_type=\"Q4\"\n"
        << "[material]\n[crystal]\nnum_grains=1\n";
    out.close();

    Datafile config(path);
    std::remove(path.c_str());

    Mesh mesh = MeshGenerators::generate_structured(config.mesh);
    ResultsExporter exporter(mesh);

    MaterialManager materials;
    materials.register_material(0, std::make_unique<SingleCrystalMaterial>(config.material));

    SimulationManager sim(config, mesh, exporter, materials, config.simulation.output_dir);
    sim.initialize_physics();

    std::remove((config.simulation.output_dir + "/energies.csv").c_str());
    CHECK_TRUE(true);
}