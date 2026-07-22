#include "Simulation/Simulation.h"
#include "Physics/MaterialModel.h"
#include "Utils/Logger.h"
#include "Utils/ProgressBar.h"
#include <algorithm> 
#include <stdexcept> 
#include <tracy/Tracy.hpp>

// =========================================================================
// Constructeur : Orchestration et injection des dépendances
// =========================================================================
Simulation::Simulation(const Datafile& config_in, const Mesh& mesh_in, ResultsExporter& exporter_in, const MaterialModel& material_in)
    : config(config_in),
      mesh(mesh_in),
      exporter(exporter_in),
      boundary_manager(mesh_in, config_in.boundary_rules),
      material(material_in), 
      diagnostics(config_in, mesh_in),
      chrono(),
      // Initialisation des modules physiques (Couplage explicite conservé)
      mechanics(config_in, mesh_in, boundary_manager),
      electrostatics(config_in, mesh_in, boundary_manager),
      fracture(config_in, mesh_in),
      polarization(config_in, mesh_in, boundary_manager),
      // Initialisation des nouveaux managers délégués
      time_manager(config_in),
      io_manager(exporter_in, diagnostics, config_in.simulation.output_dir + "/energies.csv", config_in)
{
    boundary_manager.initialize_all_boundaries();
}

void Simulation::initialize_mesh() {
    if (config.chrono.mesh) { 
        chrono.start(); 
        Logger::info("Starting mesh generation..."); 
    }
    
    exporter.exportToMesh(config.simulation.output_dir + "/" + config.simulation.mesh_create_file + ".mesh");
    
    if (config.chrono.mesh) { 
        chrono.stop(); 
        Logger::time("Mesh generation completed in ", chrono.elapsed_ms()); 
    }
}

void Simulation::initialize_physics() {
    Logger::info("Starting multi-physics initialization...");

    // Enregistrement de l'état initial (t=0) délégué au IOManager
    io_manager.extract_and_save_results(0.0, 0, "initial", polarization, mechanics, fracture, electrostatics, material);
    
    Logger::info("Initialization phase complete.");
}

// =========================================================================
// RÉSOLUTION NON-LINÉAIRE (BOUCLE DE PICARD) - Algorithm 1 (Abdollahi & Arias)
// =========================================================================
int Simulation::compute_one_step_physics(double time, double dt) {
    ZoneScoped;
    boundary_manager.update_time(time);
    
    // 1. Sauvegarde de l'état n (itération m=0) pour le calcul d'erreur intra-pas
    if (config.physics_toggle.polarization) polarization.save_previous_iteration();
    if (config.physics_toggle.mechanics) mechanics.save_previous_iteration();
    if (config.physics_toggle.electrostatics) electrostatics.save_previous_iteration();

    fracture.update_precrack_geometry(time);

    int m = 0;
    double err_p = 1.0, err_v = 1.0;
    const double tol_ferro = config.simulation.tol_ferro;  
    const double tol_vfield = config.simulation.tol_vfield; 
    const int MAX_ITER = config.simulation.max_iter; 
    
    // Temps de relaxation pseudo-temporel (section 3.1 du papier, eq 15 et 16)
    const double dt_relax = config.simulation.dt_relax;
    
    // 2. Boucle Repeat-Until (Algorithme couplé itératif staggered)
    do {
        m++;

        // Ligne 6 : P_m utilise P_m-1, Phi_m-1, V_m-1
        if (config.physics_toggle.polarization) { 
            polarization.update_P(time, dt_relax, fracture, mechanics, electrostatics, material); 
        }
        
        // Ligne 7 : u_m utilise P_m et V_m-1 (Instantané)
        if (config.physics_toggle.mechanics) {
            mechanics.update_u(time, polarization, fracture, material);
        }
        
        // Ligne 8 : Phi_m utilise P_m et V_m-1 (Instantané)
        if (config.physics_toggle.electrostatics) {
            electrostatics.update_phi(time, polarization, fracture, material); 
        }
        
        // Ligne 9 : V_m utilise P_m, u_m, Phi_m et V_m-1
        if (config.physics_toggle.fracture) {
            fracture.update_v(dt_relax, polarization, mechanics, electrostatics, material);
        }

        // 3. Vérification de la convergence (Ligne 10)
        err_p = 0.0; err_v = 0.0;
        
        if (config.physics_toggle.polarization) {
            err_p = polarization.calculate_error();
            polarization.save_previous_iteration(); 
        }
        if (config.physics_toggle.fracture) {
            err_v = fracture.calculate_error(); 
            fracture.save_previous_iteration(); 
        }
        
        Logger::debug("[Convergence] t=" + std::to_string(time) + " m=" + std::to_string(m) + 
                      " err_p=" + std::to_string(err_p) + " err_v=" + std::to_string(err_v), config.simulation.debug_enabled);
                      
    } while ((err_p > tol_ferro || err_v > tol_vfield) && m < MAX_ITER);

    if (m >= MAX_ITER) {
        Logger::debug("Non-convergence au temps " + std::to_string(time), config.simulation.debug_enabled);
        return -1; // Échec
    }

    return m; // Succès
}

// =========================================================================
// BOUCLE TEMPORELLE - pas de temps FIXE (Algorithm 1, Abdollahi & Arias 2011)
// n = 100 increments de charge, Delta t^n = 3e-2, sans rollback ni
// adaptation de dt. La robustesse pres du saut instable de propagation est
// geree par la boucle de Picard interne (tolerance / max_iter dans
// compute_one_step_physics), pas par une reduction du pas de charge.
// =========================================================================
void Simulation::run() {
    if (config.chrono.run) chrono.start();
    Logger::info("Starting simulation...");
    std::string initial_time_str = chrono.get_datetime_string();
    
    ProgressBar progressBar(100, "[SIMULATION]");
    
    while (!time_manager.is_finished()) {
        double next_time = time_manager.get_time() + time_manager.get_dt();
        int iters = compute_one_step_physics(next_time, time_manager.get_dt());

        if (iters < 0) {
            // Pas de rollback : comme dans le papier, on garde le meilleur
            // etat obtenu apres MAX_ITER et on avance quand meme. On logue
            // pour pouvoir remonter la tolerance / max_iter si ca arrive trop.
            Logger::warning("Non-convergence a t=" + std::to_string(next_time) +
                             " : etat non stabilise conserve, la charge avance quand meme.");
        }

        // On valide le pas : les etats "_current" deviennent les etats "_n"
        update_physics_history();
        time_manager.advance_step();

        // Extraction des données déléguée au IOManager
        io_manager.extract_and_save_results(
            time_manager.get_time(), 
            time_manager.get_step(), 
            initial_time_str,
            polarization, mechanics, fracture, electrostatics, material
        );

        // Mise à jour de l'interface
        int progress = static_cast<int>((time_manager.get_time() / config.simulation.total_time) * 100.0);
        progressBar.update(std::min(progress, 100), 0.0); 
    }
    
    progressBar.finish();
    
    // Forcer la dernière écriture du CSV si nécessaire
    diagnostics.write_csv(config.simulation.output_dir + "/energies_final.csv");
    
    if (config.chrono.run) { 
        chrono.stop(); 
        Logger::time("Simulation completed in ", chrono.elapsed_ms()); 
    }
}

// =========================================================================
// UTILITAIRES PRIVÉS DE GESTION D'ÉTAT
// =========================================================================

void Simulation::update_physics_history() {
    // Cette fonction valide t_n, indispensable pour l'irréversibilité v_n
    if (config.physics_toggle.mechanics) mechanics.update_history();
    if (config.physics_toggle.fracture) fracture.update_history();
    if (config.physics_toggle.polarization) polarization.update_history();
    if (config.physics_toggle.electrostatics) electrostatics.update_history();
}