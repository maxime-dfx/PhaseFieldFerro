#pragma once
#include "IO/ResultsExporter.h"
#include "IO/Diagnostics.h"
#include "IO/Datafile.h"
#include <string>

// Forward declarations pour éviter les inclusions circulaires
class Polarization;
class Mechanics;
class Fracture;
class Electrostatics;
class MaterialModel;

class IOManager {
private:
    ResultsExporter& exporter;
    Diagnostics& diagnostics;
    std::string energy_csv_path;
    const Datafile& config;

public:
    IOManager(ResultsExporter& exp, Diagnostics& diag, const std::string& csv_path, const Datafile& conf);
    
    void extract_and_save_results(double time, int step, const std::string& initial_time_str,
                                  Polarization& polarization, Mechanics& mechanics, 
                                  Fracture& fracture, Electrostatics& electrostatics, 
                                  const MaterialModel& material);
};