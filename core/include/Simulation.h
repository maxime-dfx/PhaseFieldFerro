#pragma once

#include <string>
#include "IO/Datafile.h"
#include "Mesh/Mesh.h"
#include "IO/ResultsExporter.h"
#include "IO/Diagnostics.h"
#include "Utils/Chrono.h"
#include "BC/BoundaryManager.h"

// Modules physiques
#include "Physics/Math.h"
#include "Physics/Mechanics.h"
#include "Physics/Electrostatics.h"
#include "Physics/Fracture.h"
#include "Physics/Polarization.h"

class Simulation {
private:
    const Datafile& config;
    const Mesh& mesh;
    ResultsExporter& exporter;
    
    // Utilitaires et gestionnaires
    BoundaryManager boundary_manager;   
    Math math;
    Diagnostics diagnostics;   
    Chrono chrono;
    std::string energy_csv_path; 

    // Modules Physiques
    Mechanics mechanics;
    Electrostatics electrostatics;
    Fracture fracture;
    Polarization polarization;

public:
    Simulation(const Datafile& config, const Mesh& mesh, ResultsExporter& exporter);
    
    void initialize_mesh();
    void initialize_physics();
    void run();

private:
    // Sous-routines de la boucle temporelle (SRP)
    int compute_one_step_physics(double time, double dt);
    
    void save_previous_states();
    void restore_previous_states();
    void update_physics_history();
    void extract_and_save_results(double time, int step, const std::string& initial_time_str);
};