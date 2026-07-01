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
<<<<<<< HEAD
#include "IO/Diagnostics.h"
=======
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da

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
<<<<<<< HEAD
    Diagnostics diagnostics;   
    std::string energy_csv_path; 
=======
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da

public:
    Simulation(const Datafile& config, const Mesh& mesh, ResultsExporter& exporter);
    void initializeMesh();
    void initializePhysics();
<<<<<<< HEAD
    int ComputeOneStepPhysics(double time, double dt);
=======
    void ComputeOneStepPhysics(double time);
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
    void run();
};