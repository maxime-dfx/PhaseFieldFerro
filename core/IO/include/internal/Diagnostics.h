#pragma once
#include <string>
#include <vector>
#include "Physics/include/Core/SystemMetrics.h" // On inclut la structure au lieu de la redéfinir !

struct EnergyRecord {
    double load_step = 0.0;
    double time = 0.0;
    SystemMetrics metrics;
};

class Diagnostics {
public:
    Diagnostics() = default; // Constructeur par défaut

    void record(double load_step, double time, const SystemMetrics& metrics);
    void write_csv(const std::string& path) const;
    void append_csv(const std::string& path);

    const std::vector<EnergyRecord>& get_history() const { return history; }

private:
    std::vector<EnergyRecord> history;
    bool append_header_written = false;
};