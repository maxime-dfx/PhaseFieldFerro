#include "Core/include/SimulationApp.h"
#include "Mesh/include/MeshGenerators.h"
#include "Materials/Core/MaterialManager.h"
#include "Utils/include/Logger.h"
#include "Utils/include/RunId.h"
#include "Utils/include/EnsureDir.h"
#include "Utils/include/Profiling.h"

SimulationApp::SimulationApp(int argc, char* argv[]) {
    m_config_file = (argc > 1) ? argv[1] : "../config.toml";
}

void SimulationApp::setup() {
    PROFILE_ZONE_NC("SimulationApp::setup", PROFILE_COLOR_SEQUENTIAL);
    
    // 1. Fichier de Configuration & Logs
    m_config = std::make_unique<Datafile>(m_config_file);
    Logger::set_level(m_config->simulation.debug_enabled ? LogLevel::DEBUG : LogLevel::INFO);
    Logger::info("Configuration chargée : " + m_config_file);

    // 2. Dossier de sortie
    m_run_output_dir = m_config->simulation.output_dir + "/run_" + RunId::make();
    EnsureDir::prepareOutputDirectory(m_run_output_dir);
    Logger::info("Sortie de ce run : " + m_run_output_dir);

    // 3. Maillage
    if (m_config->mesh.calcul_mesh) {
        m_mesh = std::make_unique<Mesh>(MeshGenerators::generate_structured(m_config->mesh));
    } else {
        m_mesh = std::make_unique<Mesh>(MeshGenerators::load_from_gmsh(m_config->mesh.get_mesh_file, m_config->mesh.Lx, m_config->mesh.Ly));
    }
    m_mesh->apply_rcm();
    m_mesh->compute_coloring();

    // 4. Matériaux
    m_materials_registry = std::make_unique<MaterialManager>();
    for (const auto& mat_cfg : m_config->materials) {
        m_materials_registry->create_and_register(mat_cfg.id, mat_cfg);
    }

    // 5. Exporter et Manager principal
    m_exporter = std::make_unique<ResultsExporter>(*m_mesh);
    
    m_sim_manager = std::make_unique<SimulationManager>(
        *m_config, *m_mesh, *m_exporter, *m_materials_registry, m_run_output_dir
    );
    
    m_sim_manager->initialize_physics();
}

void SimulationApp::run() {
    if (m_sim_manager) {
        m_sim_manager->run();
    }
}