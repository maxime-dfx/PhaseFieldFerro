#pragma once

#include <iostream>
#include <string>
#include <mutex>
#include <utility>

// Support conditionnel pour OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

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

    // --- TEMPLATES VARIADIQUES C++17 ---

    template<typename... Args>
    static void info(Args&&... args) {
        if (current_level > LogLevel::INFO || should_ignore()) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "[INFO] ";
        // Fold expression : applique "std::cout <<" à chaque argument
        (std::cout << ... << std::forward<Args>(args));
        std::cout << "\n";
    }

    template<typename... Args>
    static void debug(Args&&... args) {
        if (current_level > LogLevel::DEBUG || should_ignore()) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "\033[36m[DEBUG] ";
        (std::cout << ... << std::forward<Args>(args));
        std::cout << "\033[0m\n";
    }

    // Surcharge optionnelle de debug (si forcée par un paramètre externe)
    template<typename... Args>
    static void debug_force(bool force_enabled, Args&&... args) {
        if ((force_enabled || current_level <= LogLevel::DEBUG) && !should_ignore()) {
            std::lock_guard<std::mutex> lock(log_mutex);
            std::cout << "\033[36m[DEBUG] ";
            (std::cout << ... << std::forward<Args>(args));
            std::cout << "\033[0m\n";
        }
    }

    template<typename... Args>
    static void warning(Args&&... args) {
        if (current_level > LogLevel::WARNING || should_ignore()) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "\033[33m[AVERTISSEMENT] ";
        (std::cout << ... << std::forward<Args>(args));
        std::cout << "\033[0m\n";
    }

    template<typename... Args>
    static void error(Args&&... args) {
        if (current_level > LogLevel::ERROR || should_ignore()) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cerr << "\033[31m[ERREUR] ";
        (std::cerr << ... << std::forward<Args>(args));
        std::cerr << "\033[0m\n";
    }

    static void time(const std::string& message, double duration_ms) {
        if (current_level > LogLevel::INFO || should_ignore()) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "\033[34m[TEMPS] " << message << " : " << duration_ms << " ms\033[0m\n";
    }
};