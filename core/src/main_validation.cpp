#include "Validation.h"
#include "IO/Datafile.h"
#include "Physics/Math.h"
#include "Core/Mesh.h"
#include "Core/MeshGenerator.h"
#include "Core/BoundaryManager.h"
#include "Physics/Polarization.h"
#include "Physics/Fracture.h"
#include "Physics/Mechanics.h"       // INCLUSION AJOUTÉE
#include "Physics/Electrostatics.h"
#include <iostream>
#include <cstdlib>
#include <vector>

int main(int argc, char** argv) {
    // // Vérification de la présence du fichier de configuration en argument
    // if (argc < 2) {
    //     std::cerr << "Usage: " << argv[0] << " <config.toml>\n";
    //     std::cerr << "(Necessaire pour charger les parametres materiau utilises par Math.h)\n";
    //     return EXIT_FAILURE;
    // }

    // try {
    //     // 1. Chargement des paramètres et des outils mathématiques
    //     Datafile config(argv[1]);
    //     Math math(config);

    //     // 2. Création du maillage spatial
    //     Mesh mesh = MeshGenerator::generate_structured_mesh(config.mesh);

    //     // 3. Initialisation des conditions aux limites
    //     BoundaryManager boundary_manager(mesh, config);
    //     boundary_manager.initialize_all_boundaries();
    //     boundary_manager.update_time(0.0);

    //     // 4. Instanciation de TOUTES les physiques (Nécessaire pour le couplage)
    //     Polarization pol(config, mesh, boundary_manager);
    //     Fracture frac(config, mesh);
    //     Mechanics mec(config, mesh, boundary_manager); // INSTANCIATION AJOUTÉE
    //     Electrostatics elec(config, mesh, boundary_manager);

    //     // 5. Initialisation du moteur de validation
    //     Validation val(math);

    //     // 6. Exécution des tests découplés (Niveaux 0, 1 et 2)
    //     std::vector<ValidationResult> results = val.run_all();

    //     // 7. Exécution des tests spatiaux couplés (Niveau 3)
    //     // Les signatures correspondent maintenant à celles de Simulation.cpp
    //     results.push_back(val.test_domain_wall_profile(config, mesh, pol, frac, mec, elec));
    //     results.push_back(val.test_phase_field_crack_profile(config, mesh, frac, pol, mec, elec));
    //     results.push_back(val.test_electrostatics_patch(config, mesh, elec, pol, frac, mec));

    //     // 8. Affichage du bilan dans la console
    //     val.print_report(results);

    //     // 9. Le programme retourne 0 (succès) uniquement si tous les tests passent
    //     return Validation::all_passed(results) ? EXIT_SUCCESS : EXIT_FAILURE;

    // } catch (const std::exception& e) {
    //     // Capture propre des erreurs (ex: fichier TOML introuvable, maillage invalide)
    //     std::cerr << "[ERREUR FATALE] " << e.what() << "\n";
    //     return EXIT_FAILURE;
    // }
}