#include "Simulation.h"
#include "Utils/Logger.h"

// Constructeur : respect strict de l'ordre de déclaration du fichier .h
Simulation::Simulation(const Datafile& config_in, const Mesh& mesh_in, ResultsExporter& exporter_in)
    : config(config_in),
      mesh(mesh_in),
      exporter(exporter_in),
      material(config_in),
      math(material),
      boundary_manager(mesh_in, config_in),
      mechanics(config_in, mesh_in, boundary_manager), 
      electrostatics(config_in, mesh_in),
      fracture(config_in, mesh_in),
      polarization(config_in, mesh_in, boundary_manager), 
      chrono()
{
}

void Simulation::initializeMesh() {
    if (config.get_chrono_mesh()) { 
        chrono.start(); 
        Logger::info("Starting mesh generation..."); 
    }
    
    exporter.exportToMesh(config.getOutputDir() + "/" + config.getMeshCreateFile() + ".mesh");
    
    if (config.get_chrono_mesh()) { 
        chrono.stop(); 
        Logger::time("Mesh generation completed in ", chrono.elapsed_ms()); 
    }
}

void Simulation::initializePhysics() {
    Logger::info("Starting multi-physics initialization...");

    if (config.get_chrono__Polarization()) { chrono.start(); Logger::info("Polarization init..."); chrono.stop(); }
    if (config.get_chrono_Mechanics()) { chrono.start(); Logger::info("Mechanics init..."); chrono.stop(); }
    if (config.get_chrono_Electrostatics()) { chrono.start(); Logger::info("Electrostatics init..."); chrono.stop(); }
    if (config.get_chrono_Fracture()) { chrono.start(); Logger::info("Fracture init..."); chrono.stop(); }

    Logger::info("Starting multi-physics VTK export...");
    
    exporter.exportMultiPhysicsVTK(
        config.getOutputDir() + "/initial_state/multiphysics_initial.vtk", 
        {"v"}, 
        {&fracture.get_v()},
        {"P", "U", "E"}, 
        {&polarization.get_Px(), &mechanics.get_ux(), &electrostatics.get_Ex()}, 
        {&polarization.get_Py(), &mechanics.get_uy(), &electrostatics.get_Ey()}
    );
    
    Logger::info("Initialization phase complete.");
}

void Simulation::ComputeOneStepPhysics(double time) {
    boundary_manager.update_time(time);
    // 1. Sauvegarde de l'état n (itération m=0)
    // Cela permet aux physiques d'avoir accès à u_{n}, p_{n}, v_{n}
    if (config.enable_polarization()) polarization.save_previous_iteration();
    if (config.enable_mecanics()) mechanics.save_previous_iteration();
    if (config.enable_electrostatics()) electrostatics.save_previous_iteration();
    if (config.enable_fracture()) fracture.save_previous_iteration();

    int m = 0;
    double err_p = 1.0;
    double err_v = 1.0;
    const double tol_ferro = config.get_tol_ferro();  // À remplacer par config.get_tol_ferro()
    const double tol_vfield = config.get_tol_vfield(); // À remplacer par config.get_tol_vfield()
    const int MAX_ITER = 50;

    // 2. Boucle repeat ... until (L'algorithme de ton papier)
    do {
        m++;

        // Étape 6: Compute p^m
        if (config.enable_polarization()) { 
            polarization.update_polarization_component(time, fracture, mechanics, electrostatics, math, 0);
            polarization.update_polarization_component(time, fracture, mechanics, electrostatics, math, 1);
        }

        // Étape 7: Compute u^m
        if (config.enable_mecanics()) {
            mechanics.update_u(time, polarization, fracture, math);
        }

        // Étape 8: Compute phi^m
        if (config.enable_electrostatics()) {
            electrostatics.update_Ex(time, polarization, mechanics, fracture);
            electrostatics.update_Ey(time, polarization, mechanics, fracture);
        }

        // Étape 9: Compute v^m
        if (config.enable_fracture()) {
            fracture.update_v(time, polarization, mechanics, electrostatics);
        }

        // Étape 10: Vérification de la convergence
        err_p = 0.0;
        err_v = 0.0;
        
        if (config.enable_polarization()) {
            err_p = polarization.calculate_error(); // || p^m - p^{m-1} ||
            polarization.save_previous_iteration(); // p^{m-1} devient p^m pour le tour suivant
        }
        
        if (config.enable_fracture()) {
            err_v = fracture.calculate_error();     // || v^m - v^{m-1} ||
            fracture.save_previous_iteration();     // v^{m-1} devient v^m
        }

    } while ((err_p > tol_ferro || err_v > tol_vfield) && m < MAX_ITER);

    if (m >= MAX_ITER) {
        Logger::warning("Non-convergence au temps " + std::to_string(time) + 
                        " (err_p=" + std::to_string(err_p) + ", err_v=" + std::to_string(err_v) + ")");
    }
}

void Simulation::run() {
    int total_steps = static_cast<int>(config.get_total_time() / config.get_dt());
    ProgressBar progressBar(total_steps, "[SIMULATION]");
    
    if (config.get_chrono_run()) { chrono.start(); }
    
    Logger::info("Starting simulation...");
    std::string initial_time = chrono.get_datetime_string();
    double time = 0.0;
    
    for (int step = 1; step <= total_steps; ++step) {
        time += config.get_dt();
        
        ComputeOneStepPhysics(time);

        if (step % config.get_save_frequency() == 0) {
            std::string filename = config.getOutputDir() + "/VTK_" + initial_time + 
                                   "/multiphysics_results" + std::to_string(step) + ".vtk";
            
            exporter.exportMultiPhysicsVTK(
                filename, 
                {"v"}, 
                {&fracture.get_v()},
                {"P", "U", "E"}, 
                {&polarization.get_Px(), &mechanics.get_ux(), &electrostatics.get_Ex()}, 
                {&polarization.get_Py(), &mechanics.get_uy(), &electrostatics.get_Ey()}
            );
        }
        progressBar.update(step, 0.0); 
    }
    
    progressBar.finish();
    
    if (config.get_chrono_run()) { 
        chrono.stop(); 
        Logger::time("Simulation completed in ", chrono.elapsed_ms()); 
    }
    Logger::info("Simulation completed.");
}