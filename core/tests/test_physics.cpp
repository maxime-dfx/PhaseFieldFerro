#include <iostream>
#include <cassert>
#include "Utils/Logger.h"

void test_laplacian_zero() {
    Logger::info("Test: Laplacien d'un champ constant...");
    double laplacien = 10.0; 
    if (laplacien != 0.0) {
        Logger::error("Echec du test : Laplacien attendu 0.0, obtenu " + std::to_string(laplacien));
        exit(1); 
    }
    Logger::success("Test reussi !");
}

int main() {
    test_laplacian_zero();
    return 0;
}