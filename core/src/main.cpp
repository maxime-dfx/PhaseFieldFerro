#include "Simulation.h"
#include "IO/Datafile.h"
#include "Core/Mesh.h"
#include "IO/ResultsExporter.h"
#include "Utils/Logger.h"
#include <string>

int main(int argc, char** argv) {
    // 0. Detection de l'option --restart en ligne de commande
    bool do_restart = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--restart") {
            do_restart = true;
        }
    }

    // 1. Initialiser la configuration
    Datafile config("../input/config.toml");
    Logger::debug("Configuration loaded successfully.", config.debug_enabled());

    // 2. Initialiser le maillage
    Mesh mesh(config);
    Logger::debug("Mesh initialized successfully.", config.debug_enabled());

    // 3. Déléguer tout le reste à la classe Simulation !
    ResultsExporter exporter(mesh);
    Logger::debug("Results exporter initialized successfully.", config.debug_enabled());

    Simulation sim(config, mesh, exporter);
    Logger::debug("Simulation initialized successfully.", config.debug_enabled());

    if (do_restart) {
        std::string checkpoint_path = config.getOutputDir() + "/checkpoint.bin";
        sim.set_restart_path(checkpoint_path);
        Logger::info("Mode restart active : reprise depuis " + checkpoint_path);
    }

    // --- LIGNES MANQUANTES RAJOUTÉES ---
    sim.initializeMesh();
    sim.initializePhysics();
    // -----------------------------------

    // 4. Lancer les calculs 
    sim.run();
    Logger::debug("Simulation completed successfully.", config.debug_enabled());

    return 0;
}