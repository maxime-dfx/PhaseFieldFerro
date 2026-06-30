#pragma once
#include <iostream>
#include <string>
#include "IO/Datafile.h"

class Logger {
public:
    
    static void announce(const std::string& message) {
        std::cout << "\033[35m[ANNONCE] " << message << "\033[0m\n";
    }
    static void info(const std::string& message) {
        std::cout << "[INFO] " << message << "\n";
    }

    static void debug(const std::string& message, bool debug_enabled) {
        if (debug_enabled) {
            std::cout << "\033[36m[DEBUG] " << message << "\033[0m\n";
        }
    }

    static void success(const std::string& message) {
        // \033[32m met le texte en vert, \033[0m le remet à la normale
        std::cout << "\033[32m[SUCCES] " << message << "\033[0m\n";
    }

    static void time(const std::string& message, double duration_ms) {
        std::cout << "\033[34m[TEMPS] " << message << " : " << duration_ms << " ms\033[0m\n";
    }

    static void warning(const std::string& message) {
        // \033[33m met le texte en jaune
        std::cout << "\033[33m[AVERTISSEMENT] " << message << "\033[0m\n";
    }

    static void error(const std::string& message) {
        // \033[31m met le texte en rouge (on utilise cerr pour les erreurs)
        std::cerr << "\033[31m[ERREUR] " << message << "\033[0m\n";
    }
};