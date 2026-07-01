#pragma once
#include "IO/Datafile.h"
#include "Core/Mesh.h"
#include "IO/ResultsExporter.h"
#include "Utils/Chrono.h"
#include "Utils/Logger.h"
#include "Utils/ProgressBar.h"
#include "Physics/Material.h"
#include "Physics/Math.h"
#include "Physics/Mechanics.h"
#include "Physics/Electrostatics.h"
#include "Physics/Fracture.h"
#include "Physics/Polarization.h"
#include "IO/Diagnostics.h"

class Simulation {
private:
    const Datafile& config;
    const Mesh& mesh;
    ResultsExporter& exporter;
    Material material;
    Math math;
    BoundaryManager boundary_manager;   

    Mechanics mechanics;
    Electrostatics electrostatics;
    Fracture fracture;
    Polarization polarization;
    Chrono chrono;
    Diagnostics diagnostics;   
    std::string energy_csv_path; 

public:
    Simulation(const Datafile& config, const Mesh& mesh, ResultsExporter& exporter);
    void initializeMesh();
    void initializePhysics();
    int ComputeOneStepPhysics(double time, double dt);
    void run();
};