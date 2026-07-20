#include "Simulation.h"
#include "Physics/MaterialModel.h"
#include "Utils/Logger.h"
#include "Utils/ProgressBar.h"
#include <algorithm> 
#include <stdexcept> 

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
// BOUCLE TEMPORELLE ADAPTATIVE (ORCHESTRÉE PAR LE TIMEMANAGER)
// =========================================================================
void Simulation::run() {
    if (config.chrono.run) chrono.start();
    Logger::info("Starting simulation...");
    std::string initial_time_str = chrono.get_datetime_string();
    
    ProgressBar progressBar(100, "[SIMULATION]");
    
    // Remplacement du while(time < config.simulation.total_time) par le manager
    while (!time_manager.is_finished()) {
        bool step_accepted = false;
        
        while (!step_accepted) {
            // 1. Sauvegarde pour éventuel rollback en cas de divergence
            save_previous_states();

            // 2. Tentative de résolution
            double next_time = time_manager.get_time() + time_manager.get_dt();
            int iters = compute_one_step_physics(next_time, time_manager.get_dt()); 
            
            // 3. Analyse du résultat
            if (iters > 0) {
                step_accepted = true;
                
                // On valide le pas de temps, les états "_current" deviennent les états "_n"
                update_physics_history();
                time_manager.advance_step();
                time_manager.adapt_dt_after_success(); // Si implémentée

            } else {
                // ÉCHEC : Rollback strict et réduction du pas de temps géré par le manager
                restore_previous_states();
                time_manager.adapt_dt_after_failure(); 
                // Note : adapt_dt_after_failure() lèvera une exception si dt < dt_min
            }
        } 

        // 4. Extraction des données déléguée au IOManager
        io_manager.extract_and_save_results(
            time_manager.get_time(), 
            time_manager.get_step(), 
            initial_time_str,
            polarization, mechanics, fracture, electrostatics, material
        );

        // 5. Mise à jour de l'interface
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

void Simulation::save_previous_states() {
    if (config.physics_toggle.mechanics) mechanics.save_previous_state();
    if (config.physics_toggle.fracture) fracture.save_previous_state();
    if (config.physics_toggle.polarization) polarization.save_previous_state();
    if (config.physics_toggle.electrostatics) electrostatics.save_previous_state();
}

void Simulation::restore_previous_states() {
    if (config.physics_toggle.mechanics) mechanics.restore_previous_state();
    if (config.physics_toggle.fracture) fracture.restore_previous_state();
    if (config.physics_toggle.polarization) polarization.restore_previous_state();
    if (config.physics_toggle.electrostatics) electrostatics.restore_previous_state();
}

void Simulation::update_physics_history() {
    // Cette fonction valide t_n, indispensable pour l'irréversibilité v_n
    if (config.physics_toggle.mechanics) mechanics.update_history();
    if (config.physics_toggle.fracture) fracture.update_history();
    if (config.physics_toggle.polarization) polarization.update_history();
    if (config.physics_toggle.electrostatics) electrostatics.update_history();
}