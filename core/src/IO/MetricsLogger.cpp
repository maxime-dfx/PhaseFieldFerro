#include "IO/MetricsLogger.h"
#include <iostream>
#include <iomanip>

MetricsLogger::MetricsLogger(const std::string& output_dir) {
    std::string filename = output_dir + "/simulation_metrics.dat";
    m_file.open(filename);
    if (m_file.is_open()) {
        m_file << "Step,Time,Load_W,Surface_Energy,Elastic_Energy\n";
    }
}

MetricsLogger::~MetricsLogger() {
    if (m_file.is_open()) m_file.close();
}

void MetricsLogger::record_step(int step, double time, double load_w, double surface_energy, double elastic_energy) {
    if (m_file.is_open()) {
        m_file << step << "," 
               << std::scientific << std::setprecision(6) << time << ","
               << load_w << ","
               << surface_energy << ","
               << elastic_energy << "\n";
        m_file.flush(); 
    }
}