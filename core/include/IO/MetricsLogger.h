#pragma once
#include <fstream>
#include <string>

class MetricsLogger {
public:
    MetricsLogger(const std::string& output_dir);
    ~MetricsLogger();

    // On lui passe directement les valeurs, il ne calcule plus rien !
    void record_step(int step, double time, double load_w, double surface_energy, double elastic_energy);

private:
    std::ofstream m_file;
};