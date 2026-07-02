#pragma once
#include "IO/Datafile.h"
#include "Core/Mesh.h"
#include "IO/ResultsExporter.h"
#include "IO/Diagnostics.h"
#include "Utils/Chrono.h"
#include "Utils/Logger.h"
#include "Utils/ProgressBar.h"
#include "Physics/Material.h"
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
    Material material;
    Math math;
    BoundaryManager boundary_manager;   

    Mechanics mechanics;
    Electrostatics electrostatics;
    Fracture fracture;
    Polarization polarization;
    Diagnostics diagnostics;
    Chrono chrono;

    std::string diagnostics_csv_path;
    std::string restart_path_;

public:
    Simulation(const Datafile& config, const Mesh& mesh, ResultsExporter& exporter);
    void initializeMesh();
    void initializePhysics();
    void ComputeOneStepPhysics(double time);
    void set_restart_path(const std::string& path) { restart_path_ = path; }
    void save_checkpoint(const std::string& path, double time, int step) const;
    void load_checkpoint(const std::string& path, double& time, int& step);
    void run();
};