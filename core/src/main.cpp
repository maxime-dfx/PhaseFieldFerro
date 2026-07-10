#include "Simulation.h"
#include "IO/Datafile.h"
#include "Core/Mesh.h"
#include "Core/MeshGenerator.h" 
#include "Core/MeshGeneratorGmsh.h"
#include "IO/ResultsExporter.h"
#include "Utils/Logger.h"
#include <string>

int main(int argc, char** argv) {
    try {
        (void)argc; 
        (void)argv;
        // 1. Initialiser la configuration (Parsing TOML)
        Datafile config("../input/config.toml");
        Logger::debug("Configuration loaded successfully.", config.simulation.debug_enabled);

        // 2. Initialiser le maillage (Factory Pattern)
        Mesh mesh = [&]() -> Mesh {
            if (config.mesh.calcul_mesh) {
                Logger::info("Mesh generated successfully.");
                return MeshGenerator::generate_structured_mesh(config.mesh);
            } else {
                Logger::info("Mesh loaded from file: " + config.mesh.get_mesh_file);
                return MeshGeneratorGmsh::load_from_msh(config.mesh.get_mesh_file, config.mesh.Lx, config.mesh.Ly);
            }
        }(); // <-- L'appel immédiat de la lambda est ici

        // 'mesh' est maintenant parfaitement initialisé et a la bonne portée !
        ResultsExporter exporter(mesh);
        Logger::debug("Results exporter initialized successfully.", config.simulation.debug_enabled);

        // 4. Initialiser l'Orchestrateur (Simulation)
        Simulation sim(config, mesh, exporter);
        Logger::debug("Simulation initialized successfully.", config.simulation.debug_enabled);

        // 5. Initialisations pré-run
        sim.initialize_physics();

        // 6. Lancer la boucle temporelle
        sim.run();
        Logger::debug("Simulation completed successfully.", config.simulation.debug_enabled);
        
    } catch (const std::exception& e) {
        // Gestion propre des erreurs fatales (ex: fichier TOML introuvable)
        std::cerr << "[ERREUR FATALE] " << e.what() << std::endl;
        return 1;
    }

    return 0;
}