#include "Core/include/SimulationApp.h"
#include "Tests/include/TestFramework.h"
#include "Utils/include/PetscGuard.h"
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    // Initialisation RAII (S'occupe de Initialize/Finalize automatiquement)
    PetscGuard petsc_guard(argc, argv);

    // Mode Test Unitaire
    if (argc > 1 && std::string(argv[1]) == "--test") {
        return Testing::run_all() ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    // Mode Production
    try {
        SimulationApp app(argc, argv);
        app.setup();
        app.run();
    } catch (const std::exception& e) {
        std::cerr << "[ERREUR FATALE] " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}