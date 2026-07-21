#pragma once
#include <string>
#include "IO/Datafile.h"
#include "Mesh/Mesh.h"
#include "IO/ResultsExporter.h"
#include "IO/Diagnostics.h"
#include "Utils/Chrono.h"
#include "BC/BoundaryManager.h"
#include "Physics/MaterialModel.h"

// Modules physiques
#include "Physics/Mechanics.h"
#include "Physics/Electrostatics.h"
#include "Physics/Fracture.h"
#include "Physics/Polarization.h"

// Nouveaux managers
#include "Simulation/TimeManager.h"
#include "IO/IOManager.h"

class Simulation {
private:
    // /!\ L'ORDRE DE DÉCLARATION DOIT CORRESPONDRE AU CONSTRUCTEUR POUR ÉVITER -Wreorder
    const Datafile& config;
    const Mesh& mesh;
    ResultsExporter& exporter;
    BoundaryManager boundary_manager;
    const MaterialModel& material;
    Diagnostics diagnostics;
    Chrono chrono;

    // Modules Physiques (explicites pour la boucle de Picard)
    Mechanics mechanics;
    Electrostatics electrostatics;
    Fracture fracture;
    Polarization polarization;

    // Managers délégués
    TimeManager time_manager;
    IOManager io_manager;

public:
    Simulation(const Datafile& config, const Mesh& mesh, ResultsExporter& exporter, const MaterialModel& material);
    
    void initialize_mesh();
    void initialize_physics();
    void run();

private:
    // Sous-routines de la boucle temporelle
    int compute_one_step_physics(double time, double dt);
    
    void save_previous_states();
    void restore_previous_states();
    void update_physics_history();
};