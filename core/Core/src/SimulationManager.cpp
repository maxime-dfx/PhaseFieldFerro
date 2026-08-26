#include "../include/SimulationManager.h"
#include "Materials/Core/MaterialConcepts.h"
#include "Utils/include/Logger.h"
#include "Utils/include/ProgressBar.h"
#include "Utils/include/Profiling.h"
#include <algorithm> 
#include <stdexcept> 

// =========================================================================
// Constructeur : Orchestration et injection des dépendances
// =========================================================================
SimulationManager::SimulationManager(const Datafile& config_in, const Mesh& mesh_in, ResultsExporter& exporter_in,
                                      const MaterialManager& materials_manager_in, const std::string& run_output_dir_in) :
    config(config_in),
    mesh(mesh_in),
    exporter(exporter_in),
    boundary_manager(mesh_in, config_in.boundary_rules),
    materials_manager(materials_manager_in),
    chrono(),
    run_output_dir(run_output_dir_in),
    m_polycrystal(config_in.crystal, mesh_in),
    physics_manager(config_in, mesh_in, boundary_manager, materials_manager_in, m_polycrystal),
    time_manager(config_in),
    io_manager(exporter_in, run_output_dir, config_in)
{
    boundary_manager.initialize_all_boundaries();
    m_polycrystal.attach_materials(materials_manager, mesh, config.crystal);
}

void SimulationManager::initialize_mesh() {
    PROFILE_ZONE_NC("initialize_mesh", PROFILE_COLOR_SEQUENTIAL);
    if (config.chrono.mesh) { 
        chrono.start(); 
        Logger::info("Starting mesh generation..."); 
    }
    
    exporter.exportToMesh(run_output_dir + "/" + config.simulation.mesh_create_file + ".mesh");
    
    if (config.chrono.mesh) { 
        chrono.stop(); 
        Logger::time("Mesh generation completed in ", chrono.elapsed_ms()); 
    }
}

void SimulationManager::initialize_physics() {
    PROFILE_ZONE_NC("initialize_physics", PROFILE_COLOR_SEQUENTIAL);
    Logger::info("Starting multi-physics initialization...");

    if (m_polycrystal.num_grains() > 1) {
        Eigen::VectorXd grain_ids(mesh.get_num_elements());
        for (int e = 0; e < mesh.get_num_elements(); ++e) {
            grain_ids(e) = static_cast<double>(m_polycrystal.grain_id_for_element(e));
        }
        exporter.exportCellScalarVTK(run_output_dir + "/grains.vtk", "grain_id", grain_ids);
        Logger::info("Microstructure polycristalline exportée : ", run_output_dir, "/grains.vtk");
    }

    io_manager.extract_and_save_results(0.0, 0, physics_manager, m_polycrystal);
    Logger::info("Initialization phase complete.");
}

// =========================================================================
// BOUCLE TEMPORELLE - pas de temps FIXE (Algorithm 1, Abdollahi & Arias 2011)
// n = 100 increments de charge, Delta t^n = 3e-2, sans rollback ni
// adaptation de dt. La robustesse pres du saut instable de propagation est
// geree par la boucle de Picard interne (tolerance / max_iter dans
// PhysicsManager::compute_one_step_physics), pas par une reduction du pas de charge.
// =========================================================================
void SimulationManager::run() {
    PROFILE_ZONE_NC("SimulationManager::run", PROFILE_COLOR_SEQUENTIAL);
    if (config.chrono.run) chrono.start();
    Logger::info("Starting simulation... (sortie : ", run_output_dir, ")");
    
    ProgressBar progressBar(100, "[SIMULATION]");
    
    while (!time_manager.is_finished()) {
        PROFILE_FRAME_MARK_N("TimeStep");
        PROFILE_ZONE_NC("SimulationManager::run::step", PROFILE_COLOR_SEQUENTIAL);
        double next_time = time_manager.get_time() + time_manager.get_dt();
        Logger::debug("==== [SimulationManager] debut etape n=", time_manager.get_step(),
                       " t=", time_manager.get_time(),
                       " -> next_time=", next_time,
                       " dt_load=", time_manager.get_dt(), " ====");
        int iters = physics_manager.compute_one_step_physics(next_time, time_manager.get_dt());
        Logger::debug("==== [SimulationManager] fin etape n=", time_manager.get_step(),
                       " : compute_one_step_physics a retourne ", iters,
                       " (>=0 : converge en 'iters' iterations Picard ; -1 : max_iter atteint sans convergence ; -2 : echec de resolution lineaire) ====");

        if (iters < 0) {
            Logger::warning("Non-convergence a t=", next_time,
                             " : etat non stabilise conserve, la charge avance quand meme.");
        }

        physics_manager.update_physics_history();
        time_manager.advance_step();

        PROFILE_ZONE_NC("SimulationManager::run::save_results_call", PROFILE_COLOR_SEQUENTIAL);
        io_manager.extract_and_save_results(
            time_manager.get_time(), 
            time_manager.get_step(), 
            physics_manager,
            m_polycrystal
        );

        // Mise à jour de l'interface
        int progress = static_cast<int>((time_manager.get_time() / config.simulation.total_time) * 100.0);
        Logger::debug("[SimulationManager] progress = time(", time_manager.get_time(),
                       ") / total_time(", config.simulation.total_time,
                       ") * 100 = ", progress, "%");
        progressBar.update(std::min(progress, 100), 0.0); 
    }
    
    progressBar.finish();
    io_manager.finalize();
    
    if (config.chrono.run) { 
        chrono.stop(); 
        Logger::time("SimulationManager completed in ", chrono.elapsed_ms()); 
    }
}
