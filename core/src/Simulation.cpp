#include "Simulation.h"
#include "Utils/Logger.h"
<<<<<<< HEAD
#include <algorithm> 
#include <stdexcept> 
=======
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da

// Constructeur : respect strict de l'ordre de déclaration du fichier .h
Simulation::Simulation(const Datafile& config_in, const Mesh& mesh_in, ResultsExporter& exporter_in)
    : config(config_in),
      mesh(mesh_in),
      exporter(exporter_in),
      material(config_in),
      math(material),
      boundary_manager(mesh_in, config_in),
      mechanics(config_in, mesh_in, boundary_manager),
      electrostatics(config_in, mesh_in, boundary_manager),
      fracture(config_in, mesh_in),
      polarization(config_in, mesh_in, boundary_manager),
<<<<<<< HEAD
      diagnostics(config_in, mesh_in),
      chrono()
{
    boundary_manager.initialize_all_boundaries();
    energy_csv_path = config.getOutputDir() + "/energies.csv";
=======
      chrono()
{
    boundary_manager.initialize_all_boundaries();
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
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
<<<<<<< HEAD

    // Enregistrement de l'état initial (t=0) dans le CSV d'énergies, pour
    // avoir un point de départ cohérent avec les graphes du papier (Fig. 6
    // part de zéro à l'état intact).
    diagnostics.record(0.0, 0.0, polarization, mechanics, fracture, electrostatics, math);
    diagnostics.append_csv(energy_csv_path);
=======
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
    
    Logger::info("Initialization phase complete.");
}

<<<<<<< HEAD
// =========================================================================
// RÉSOLUTION D'UN PAS DE TEMPS (Retourne le nombre d'itérations, ou -1 si échec)
// =========================================================================
int Simulation::ComputeOneStepPhysics(double time, double dt) {
    boundary_manager.update_time(time);
    
=======
void Simulation::ComputeOneStepPhysics(double time) {
    boundary_manager.update_time(time);
    double dt = config.get_dt();
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
    // 1. Sauvegarde de l'état n (itération m=0)
    // Cela permet aux physiques d'avoir accès à u_{n}, p_{n}, v_{n}
    if (config.enable_polarization()) polarization.save_previous_iteration();
    if (config.enable_mecanics()) mechanics.save_previous_iteration();
    if (config.enable_electrostatics()) electrostatics.save_previous_iteration();
    if (config.enable_fracture()) fracture.save_previous_iteration();

    // Fige P_n / v_n (etat de debut de pas de temps physique), utilises
<<<<<<< HEAD
    // uniquement dans les termes de masse implicites.
=======
    // uniquement dans les termes de masse implicites. Doit etre appele une
    // seule fois ici, JAMAIS a l'interieur de la boucle do...while ci-dessous
    // (contrairement a save_previous_iteration(), qui elle est rappelee a
    // chaque sous-iteration pour le calcul de l'erreur de Picard).
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
    if (config.enable_polarization()) polarization.freeze_time_step();
    if (config.enable_fracture()) fracture.freeze_time_step();

    int m = 0;
    double err_p = 1.0;
    double err_v = 1.0;
<<<<<<< HEAD
    const double tol_ferro = config.get_tol_ferro();  
    const double tol_vfield = config.get_tol_vfield(); 
=======
    const double tol_ferro = config.get_tol_ferro();  // À remplacer par config.get_tol_ferro()
    const double tol_vfield = config.get_tol_vfield(); // À remplacer par config.get_tol_vfield()
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
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
            electrostatics.update_electric_potential(dt, polarization, fracture, math);
        }

        // Étape 9: Compute v^m
        if (config.enable_fracture()) {
            fracture.update_v(dt, polarization, mechanics, electrostatics, math);
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
<<<<<<< HEAD
        
        Logger::debug("[DEBUG][Convergence] t=" + std::to_string(time) + " m=" + std::to_string(m) + 
                      " err_p=" + std::to_string(err_p) + " err_v=" + std::to_string(err_v), config.debug_enabled());
                      
    } while ((err_p > tol_ferro || err_v > tol_vfield) && m < MAX_ITER);

    // Si on n'a pas convergé dans le temps imparti, on retourne -1
    if (m >= MAX_ITER) {
        Logger::debug("Non-convergence au temps " + std::to_string(time) + 
                      " (err_p=" + std::to_string(err_p) + ", err_v=" + std::to_string(err_v) + ")", config.debug_enabled());
        return -1;
    }

    // Succès
    return m;
}

// =========================================================================
// BOUCLE PRINCIPALE AVEC PAS DE TEMPS ADAPTATIF
// =========================================================================
void Simulation::run() {
=======
        Logger::debug("[DEBUG][Convergence] t=" + std::to_string(time) + " m=" + std::to_string(m) + 
                      " err_p=" + std::to_string(err_p) + " err_v=" + std::to_string(err_v), config.debug_enabled());
    } while ((err_p > tol_ferro || err_v > tol_vfield) && m < MAX_ITER);

    if (m >= MAX_ITER) {
        Logger::warning("Non-convergence au temps " + std::to_string(time) + 
                        " (err_p=" + std::to_string(err_p) + ", err_v=" + std::to_string(err_v) + ")");
    }
}

void Simulation::run() {
    int total_steps = static_cast<int>(config.get_total_time() / config.get_dt());
    ProgressBar progressBar(total_steps, "[SIMULATION]");
    
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
    if (config.get_chrono_run()) { chrono.start(); }
    
    Logger::info("Starting simulation...");
    std::string initial_time = chrono.get_datetime_string();
<<<<<<< HEAD
    
    double time = 0.0;
    double dt = config.get_dt();
    const double dt_min = 1e-8;
    const double dt_max = config.get_dt() * 5.0;
    
    int step = 0;
    int vtk_counter = 1;
    
    ProgressBar progressBar(100, "[SIMULATION]");
    
    while (time < config.get_total_time()) {
        bool step_accepted = false;
        
        while (!step_accepted) {
            // 1. SAUVEGARDE de l'état (Mémoire N-1 pour rollback si échec)
            if (config.enable_mecanics()) mechanics.save_previous_state();
            if (config.enable_fracture()) fracture.save_previous_state();
            if (config.enable_polarization()) polarization.save_previous_state();

            // 2. RÉSOLUTION du pas de temps
            int iters = ComputeOneStepPhysics(time + dt, dt); 
            
            const int MAX_ITERS = 50; 
            
            // 3. VÉRIFICATION DE LA CONVERGENCE
            if (iters > 0 && iters <= MAX_ITERS) {
                // ---> SUCCÈS : Le pas est accepté
                step_accepted = true;
                time += dt;
                step++;
                
                // On valide les nouvelles valeurs (v_n = v_current, etc.)
                if (config.enable_mecanics()) mechanics.update_history();
                if (config.enable_fracture()) fracture.update_history();
                if (config.enable_polarization()) polarization.update_history();

                // Accélération si convergence très rapide
                if (iters <= 4) {
                    dt = std::min(dt * 1.2, dt_max);
                }
            } else {
                // ---> ÉCHEC : Divergence ou convergence trop lente
                Logger::warning("Non-convergence a t = " + std::to_string(time + dt) + " (Iters: " + std::to_string(iters) + "). Reduction de dt...");
                
                // Restauration de l'état précédent
                if (config.enable_mecanics()) mechanics.restore_previous_state();
                if (config.enable_fracture()) fracture.restore_previous_state();
                if (config.enable_polarization()) polarization.restore_previous_state();

                // Ralentissement
                dt *= 0.5;
                
                if (dt < dt_min) {
                    throw std::runtime_error("Erreur fatale : dt est devenu trop petit (< dt_min) ! Rupture numerique.");
                }
            }
        } // Fin du while (!step_accepted)

        // 4. EXTRACTIONS ET SAUVEGARDES (Uniquement si le pas a été accepté)
        diagnostics.record(time, time, polarization, mechanics, fracture, electrostatics, math);
        diagnostics.append_csv(energy_csv_path);
        diagnostics.compute_nodal_energies(polarization, mechanics, fracture, electrostatics, math);

        if (step % config.get_save_frequency() == 0) {
            std::string filename = config.getOutputDir() + "/VTK_" + initial_time + 
                                   "/multiphysics_results" + std::to_string(vtk_counter++) + ".vtk";            
            exporter.exportMultiPhysicsVTK(
                filename, 
                // Liste de TOUS les champs scalaires
                {"v", "phi", "Energy_Gradient", "Energy_Elastic", "Energy_Landau", "Energy_Electric", "Energy_Surface"}, 
                {&fracture.get_v(), &electrostatics.get_phi(),
                &diagnostics.get_U_nodal(), &diagnostics.get_W_nodal(), 
                &diagnostics.get_chi_nodal(), &diagnostics.get_elec_nodal(), 
                &diagnostics.get_surf_nodal()},
                // Liste de TOUS les champs vectoriels
=======
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
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
                {"P", "U", "E"}, 
                {&polarization.get_Px(), &mechanics.get_ux(), &electrostatics.get_Ex()}, 
                {&polarization.get_Py(), &mechanics.get_uy(), &electrostatics.get_Ey()}
            );
        }
<<<<<<< HEAD
        
        int progress = static_cast<int>((time / config.get_total_time()) * 100.0);
        progressBar.update(std::min(progress, 100), 0.0); 
    }
    
    progressBar.finish();

    // Réécriture complète du CSV en fin de simulation
    diagnostics.write_csv(config.getOutputDir() + "/energies_final.csv");
=======
        progressBar.update(step, 0.0); 
    }
    
    progressBar.finish();
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
    
    if (config.get_chrono_run()) { 
        chrono.stop(); 
        Logger::time("Simulation completed in ", chrono.elapsed_ms()); 
    }
    Logger::info("Simulation completed.");
}