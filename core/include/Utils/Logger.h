#pragma once
#include <iostream>
#include <string>
#include <mutex>

// 1. Définition des niveaux d'importance
enum class LogLevel {
    DEBUG = 0,
    INFO = 1,
    WARNING = 2,
    ERROR = 3,
    NONE = 4
};

class Logger {
private:
    // 2. Variables globales (inline permet de les définir directement dans le .h en C++17)
    inline static LogLevel current_level = LogLevel::INFO;
    inline static std::mutex log_mutex; // Empêche les conflits si OpenMP écrit en parallèle

public:
    // Permet de configurer le niveau de verbosité une seule fois dans le main
    static void set_level(LogLevel level) {
        current_level = level;
    }

    static void announce(const std::string& message) {
        if (current_level > LogLevel::INFO) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "\033[35m[ANNONCE] " << message << "\033[0m\n";
    }

    static void info(const std::string& message) {
        if (current_level > LogLevel::INFO) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "[INFO] " << message << "\n";
    }

    // NOUVELLE VERSION : S'appuie sur le niveau global (plus besoin de passer le booléen)
    static void debug(const std::string& message) {
        if (current_level > LogLevel::DEBUG) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "\033[36m[DEBUG] " << message << "\033[0m\n";
    }

    // ANCIENNE VERSION (Conservée pour ne pas casser ton code existant)
    static void debug(const std::string& message, bool debug_enabled) {
        if (debug_enabled || current_level <= LogLevel::DEBUG) {
            std::lock_guard<std::mutex> lock(log_mutex);
            std::cout << "\033[36m[DEBUG] " << message << "\033[0m\n";
        }
    }

    static void success(const std::string& message) {
        if (current_level > LogLevel::INFO) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "\033[32m[SUCCES] " << message << "\033[0m\n";
    }

    static void time(const std::string& message, double duration_ms) {
        if (current_level > LogLevel::INFO) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "\033[34m[TEMPS] " << message << " : " << duration_ms << " ms\033[0m\n";
    }

    static void warning(const std::string& message) {
        if (current_level > LogLevel::WARNING) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "\033[33m[AVERTISSEMENT] " << message << "\033[0m\n";
    }

    static void error(const std::string& message) {
        if (current_level > LogLevel::ERROR) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cerr << "\033[31m[ERREUR] " << message << "\033[0m\n";
    }
};