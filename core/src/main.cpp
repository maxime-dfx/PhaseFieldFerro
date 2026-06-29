#include "Simulation.h"
#include "IO/Datafile.h"
#include "Core/Mesh.h"
#include "IO/ResultsExporter.h"
#include "Utils/Logger.h"

int main() {
    // 1. Initialiser la configuration
    Datafile config("../input/config.toml");
    Logger::debug("Configuration loaded successfully.");
    // 2. Initialiser le maillage

    Mesh mesh(config);
    Logger::debug("Mesh initialized successfully.");

    // 3. Déléguer tout le reste à la classe Simulation !
    ResultsExporter exporter(mesh);
    Logger::debug("Results exporter initialized successfully.");
    Simulation sim(config, mesh, exporter);
    Logger::debug("Simulation initialized successfully.");

    // 4. Lancer les calculs

    sim.run();
    Logger::debug("Simulation completed successfully.");
    
    return 0;
}