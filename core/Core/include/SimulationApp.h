#pragma once

#include <string>
#include <memory>

#include "IO/include/Datafile.h"
#include "Mesh/include/Mesh.h"
#include "Materials/Core/MaterialManager.h"
#include "IO/include/ResultsExporter.h"
#include "Core/include/SimulationManager.h"

class SimulationApp {
private:
    std::string m_config_file;
    std::string m_run_output_dir;

    // L'ordre de déclaration est important : ils seront détruits dans l'ordre inverse
    std::unique_ptr<Datafile> m_config;
    std::unique_ptr<Mesh> m_mesh;
    std::unique_ptr<MaterialManager> m_materials_registry;
    std::unique_ptr<ResultsExporter> m_exporter;
    std::unique_ptr<SimulationManager> m_sim_manager;

public:
    SimulationApp(int argc, char* argv[]);
    void setup();
    void run();
};