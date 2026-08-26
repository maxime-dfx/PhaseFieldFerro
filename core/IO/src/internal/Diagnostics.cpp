#include "IO/include/internal/Diagnostics.h"
#include "Utils/include/Profiling.h"
#include <fstream>

void Diagnostics::record(double load_step, double time, const SystemMetrics& metrics) {
    PROFILE_ZONE_NC("Diagnostics::record", PROFILE_COLOR_SEQUENTIAL);
    EnergyRecord rec;
    rec.load_step = load_step;
    rec.time = time;
    rec.metrics = metrics;
    history.push_back(rec);
}

void Diagnostics::write_csv(const std::string& path) const {
    PROFILE_ZONE_NC("Diagnostics::write_csv", PROFILE_COLOR_SEQUENTIAL);
    std::ofstream file(path);
    file << "LoadStep,Time,U,W,Chi,Elec,Bulk,Surface,Total,Vmin,Vmax,SolveFailed\n";
    for (const auto& r : history) {
        file << r.load_step << "," << r.time << "," 
             << r.metrics.U_total << "," << r.metrics.W_total << ","
             << r.metrics.chi_total << "," << r.metrics.electric_total << "," 
             << r.metrics.bulk_enthalpy << "," << r.metrics.surface_energy << "," 
             << r.metrics.total_energy << "," << r.metrics.v_min << "," 
             << r.metrics.v_max << "," << (r.metrics.solve_failed ? 1 : 0) << "\n";
    }
}

void Diagnostics::append_csv(const std::string& path) {
    PROFILE_ZONE_NC("Diagnostics::append_csv", PROFILE_COLOR_SEQUENTIAL);
    std::ofstream file(path, std::ios::app);
    if (!append_header_written) {
        file << "LoadStep,Time,U,W,Chi,Elec,Bulk,Surface,Total,Vmin,Vmax,SolveFailed\n";
        append_header_written = true;
    }
    const auto& r = history.back();
    file << r.load_step << "," << r.time << "," 
         << r.metrics.U_total << "," << r.metrics.W_total << ","
         << r.metrics.chi_total << "," << r.metrics.electric_total << "," 
         << r.metrics.bulk_enthalpy << "," << r.metrics.surface_energy << "," 
         << r.metrics.total_energy << "," << r.metrics.v_min << "," 
         << r.metrics.v_max << "," << (r.metrics.solve_failed ? 1 : 0) << "\n";
}