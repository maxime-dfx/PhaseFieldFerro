#include "Simulation.h"
#include "Utils/Logger.h"
#include "Utils/ProgressBar.h"
#include <algorithm> 
#include <stdexcept> 

// Constructeur : Orchestration et injection des dépendances
Simulation::Simulation(const Datafile& config_in, const Mesh& mesh_in, ResultsExporter& exporter_in)
    : config(config_in),
      mesh(mesh_in),
      exporter(exporter_in),
      boundary_manager(mesh_in, config_in),
      math(config_in), 
      diagnostics(config_in, mesh_in),
      chrono(),
      mechanics(config_in, mesh_in, boundary_manager),
      electrostatics(config_in, mesh_in, boundary_manager),
      fracture(config_in, mesh_in),
      polarization(config_in, mesh_in, boundary_manager){
    boundary_manager.initialize_all_boundaries();
    energy_csv_path = config.simulation.output_dir + "/energies.csv";
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

    // Enregistrement de l'état initial (t=0) dans le CSV d'énergies
    diagnostics.record(0.0, 0.0, polarization, mechanics, fracture, electrostatics, math);
    diagnostics.append_csv(energy_csv_path);
    
    Logger::info("Initialization phase complete.");
}

// =========================================================================
// RÉSOLUTION NON-LINÉAIRE (BOUCLE DE PICARD) - Algorithm 1 (Abdollahi & Arias)
// =========================================================================
int Simulation::compute_one_step_physics(double time, double dt) {
    boundary_manager.update_time(time);
    
    // 1. Sauvegarde de l'état n (itération m=0) pour le calcul d'erreur
    if (config.physics_toggle.polarization) polarization.save_previous_iteration();
    if (config.physics_toggle.mechanics) mechanics.save_previous_iteration();
    if (config.physics_toggle.electrostatics) electrostatics.save_previous_iteration();
    if (config.physics_toggle.fracture) fracture.save_previous_iteration();

    int m = 0;
    double err_p = 1.0, err_v = 1.0;
    const double tol_ferro = config.simulation.tol_ferro;  
    const double tol_vfield = config.simulation.tol_vfield; 
    const int MAX_ITER = 50;

    // 2. Boucle Repeat-Until (Algorithme couplé itératif staggered)
    do {
        m++;

        // Ligne 6 de l'Algorithme 1 : P_m utilise P_m-1, Phi_m-1, V_m-1
        if (config.physics_toggle.polarization) { 
            // Appel au nouveau solveur monolithique
            polarization.update_P(time, fracture, mechanics, electrostatics, math); // CORRIGÉ
        }
        
        // Ligne 7 de l'Algorithme 1 : u_m utilise P_m et V_m-1
        if (config.physics_toggle.mechanics) {
            mechanics.update_u(time, polarization, fracture, math);
        }
        
        // Ligne 8 de l'Algorithme 1 : Phi_m utilise P_m et V_m-1
        if (config.physics_toggle.electrostatics) {
            electrostatics.update_phi(time, polarization, fracture, math); // CORRIGÉ
        }
        
        // Ligne 9 de l'Algorithme 1 : V_m utilise P_m, u_m, Phi_m et V_m-1
        if (config.physics_toggle.fracture) {
            fracture.update_v(dt, polarization, mechanics, electrostatics, math);
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
// BOUCLE TEMPORELLE ADAPTATIVE (TIME STEPPER)
// =========================================================================
void Simulation::run() {
    if (config.chrono.run) chrono.start();
    Logger::info("Starting simulation...");
    std::string initial_time = chrono.get_datetime_string();
    
    double time = 0.0;
    double dt = config.simulation.dt;
    const double dt_min = 1e-8;
    const double dt_max = config.simulation.dt * 5.0;
    const int MAX_ITERS = 50; 
    
    int step = 0;
    ProgressBar progressBar(100, "[SIMULATION]");
    
    while (time < config.simulation.total_time) {
        bool step_accepted = false;
        
        while (!step_accepted) {
            // 1. Sauvegarde pour éventuel rollback
            save_previous_states();

            // 2. Tentative de résolution
            int iters = compute_one_step_physics(time + dt, dt); 
            
            // 3. Analyse du résultat
            if (iters > 0 && iters <= MAX_ITERS) {
                step_accepted = true;
                time += dt;
                step++;
                
                // On valide le pas de temps, les états "_current" deviennent les états "_n"
                update_physics_history();

            } else {
                // ÉCHEC : Rollback strict et réduction du pas de temps
                Logger::warning("Non-convergence a t=" + std::to_string(time + dt) + ". Reduction de dt...");
                
                restore_previous_states();
                dt *= 0.5;
                
                if (dt < dt_min) {
                    throw std::runtime_error("Erreur fatale : dt est devenu trop petit (< dt_min) ! Rupture numerique.");
                }
            }
        } 

        // 4. Extraction des données
        extract_and_save_results(time, step, initial_time);

        // 5. Mise à jour de l'interface
        int progress = static_cast<int>((time / config.simulation.total_time) * 100.0);
        progressBar.update(std::min(progress, 100), 0.0); 
    }
    
    progressBar.finish();
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
    if (config.physics_toggle.electrostatics) electrostatics.save_previous_state(); // Ajout essentiel
}

void Simulation::restore_previous_states() {
    if (config.physics_toggle.mechanics) mechanics.restore_previous_state();
    if (config.physics_toggle.fracture) fracture.restore_previous_state();
    if (config.physics_toggle.polarization) polarization.restore_previous_state();
    if (config.physics_toggle.electrostatics) electrostatics.restore_previous_state(); // Ajout essentiel
}

void Simulation::update_physics_history() {
    // Cette fonction valide t_n, indispensable pour l'irréversibilité v_n (Algorithme 1, Ligne 11)
    if (config.physics_toggle.mechanics) mechanics.update_history();
    if (config.physics_toggle.fracture) fracture.update_history();
    if (config.physics_toggle.polarization) polarization.update_history();
    if (config.physics_toggle.electrostatics) electrostatics.update_history(); // Ajout essentiel
}

void Simulation::extract_and_save_results(double time, int step, const std::string& initial_time_str) {
    diagnostics.record(time, time, polarization, mechanics, fracture, electrostatics, math);
    diagnostics.append_csv(energy_csv_path);
    diagnostics.compute_nodal_energies(polarization, mechanics, fracture, electrostatics, math);

    if (step % config.simulation.save_frequency == 0) {
        std::string filename = config.simulation.output_dir + "/VTK_" + initial_time_str + 
                               "/multiphysics_results" + std::to_string(step / config.simulation.save_frequency) + ".vtk";            
        
        exporter.exportMultiPhysicsVTK(
            filename, 
            {"v", "phi", "Energy_Gradient", "Energy_Elastic", "Energy_Landau", "Energy_Electric", "Energy_Surface"}, 
            {&fracture.get_v(), &electrostatics.get_phi(),
             &diagnostics.get_U_nodal(), &diagnostics.get_W_nodal(), 
             &diagnostics.get_chi_nodal(), &diagnostics.get_elec_nodal(), 
             &diagnostics.get_surf_nodal()},
            {"P", "U", "E"}, 
            {&polarization.get_Px(), &mechanics.get_ux(), &electrostatics.get_Ex()}, 
            {&polarization.get_Py(), &mechanics.get_uy(), &electrostatics.get_Ey()}
        );
    }
}