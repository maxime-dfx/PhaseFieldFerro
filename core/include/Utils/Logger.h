#pragma once
#include <iostream>
#include <string>
#include <mutex>

// Support conditionnel pour OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

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
    inline static LogLevel current_level = LogLevel::INFO;
    inline static std::mutex log_mutex;

    // <-- NOUVEAU : Filtre pour empêcher les conflits OpenMP
    static inline bool should_ignore() {
#ifdef _OPENMP
        // Si on est dans une boucle parallèle, seul le thread 0 loggue
        if (omp_in_parallel() && omp_get_thread_num() != 0) {
            return true;
        }
#endif
        return false;
    }

public:
    static void set_level(LogLevel level) {
        current_level = level;
    }

    static void announce(const std::string& message) {
        if (current_level > LogLevel::INFO || should_ignore()) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "\033[35m[ANNONCE] " << message << "\033[0m\n";
    }

    static void info(const std::string& message) {
        if (current_level > LogLevel::INFO || should_ignore()) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "[INFO] " << message << "\n";
    }

    static void debug(const std::string& message) {
        if (current_level > LogLevel::DEBUG || should_ignore()) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "\033[36m[DEBUG] " << message << "\033[0m\n";
    }

    static void debug(const std::string& message, bool debug_enabled) {
        if ((debug_enabled || current_level <= LogLevel::DEBUG) && !should_ignore()) {
            std::lock_guard<std::mutex> lock(log_mutex);
            std::cout << "\033[36m[DEBUG] " << message << "\033[0m\n";
        }
    }

    static void success(const std::string& message) {
        if (current_level > LogLevel::INFO || should_ignore()) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "\033[32m[SUCCES] " << message << "\033[0m\n";
    }

    static void time(const std::string& message, double duration_ms) {
        if (current_level > LogLevel::INFO || should_ignore()) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "\033[34m[TEMPS] " << message << " : " << duration_ms << " ms\033[0m\n";
    }

    static void warning(const std::string& message) {
        if (current_level > LogLevel::WARNING || should_ignore()) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "\033[33m[AVERTISSEMENT] " << message << "\033[0m\n";
    }

    static void error(const std::string& message) {
        if (current_level > LogLevel::ERROR || should_ignore()) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cerr << "\033[31m[ERREUR] " << message << "\033[0m\n";
    }
};