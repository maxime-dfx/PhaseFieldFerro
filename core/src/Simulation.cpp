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
      electrostatics(config_in, mesh_in, boundary_manager),
      fracture(config_in, mesh_in),
      polarization(config_in, mesh_in, boundary_manager),
      diagnostics(config_in, mesh_in),
      chrono(),
      diagnostics_csv_path(config_in.getOutputDir() + "/diagnostics.csv")
{
    boundary_manager.initialize_all_boundaries();
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

    // Enregistre l'etat initial (pas de charge 0, t=0) dans le suivi des
    // grandeurs globales, avant que la boucle temporelle ne demarre.
    diagnostics.record(0.0, 0.0, polarization, mechanics, fracture, electrostatics, math);
    diagnostics.append_csv(diagnostics_csv_path);
    
    Logger::info("Initialization phase complete.");
}

void Simulation::ComputeOneStepPhysics(double time) {
    boundary_manager.update_time(time);
    double dt_relax = config.get_dt_relax();
    // 1. Sauvegarde de l'état n (itération m=0)
    // Cela permet aux physiques d'avoir accès à u_{n}, p_{n}, v_{n}
    if (config.enable_polarization()) polarization.save_previous_iteration();
    if (config.enable_mecanics()) mechanics.save_previous_iteration();
    if (config.enable_electrostatics()) electrostatics.save_previous_iteration();
    if (config.enable_fracture()) fracture.save_previous_iteration();

    // Fige P_n / v_n (etat de debut de pas de temps physique), utilises
    // uniquement dans les termes de masse implicites. Doit etre appele une
    // seule fois ici, JAMAIS a l'interieur de la boucle do...while ci-dessous
    // (contrairement a save_previous_iteration(), qui elle est rappelee a
    // chaque sous-iteration pour le calcul de l'erreur de Picard).
    if (config.enable_polarization()) polarization.freeze_time_step();
    if (config.enable_fracture()) fracture.freeze_time_step();

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
        // Un seul appel : Px et Py sont résolus ensemble (même matrice A,
        // factorisation Cholesky partagée, résolution multi-RHS) au lieu de
        // deux appels séparés qui refactorisaient deux fois la même matrice.
        if (config.enable_polarization()) { 
            polarization.update_polarization(time, fracture, mechanics, electrostatics, math);
        }

        // Étape 7: Compute u^m
        if (config.enable_mecanics()) {
            mechanics.update_u(time, polarization, fracture, math);
        }

        // Étape 8: Compute phi^m
        if (config.enable_electrostatics()) {
            electrostatics.update_electric_potential(time, polarization, fracture, math);
        }

        // Étape 9: Compute v^m
        if (config.enable_fracture()) {
            fracture.update_v(dt_relax, polarization, mechanics, electrostatics, math);
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
        Logger::debug("[DEBUG][Convergence] t=" + std::to_string(time) + " m=" + std::to_string(m) + 
                      " err_p=" + std::to_string(err_p) + " err_v=" + std::to_string(err_v), config.debug_enabled());
    } while ((err_p > tol_ferro || err_v > tol_vfield) && m < MAX_ITER);

    if (m >= MAX_ITER) {
        Logger::warning("Non-convergence au temps " + std::to_string(time) + 
                        " (err_p=" + std::to_string(err_p) + ", err_v=" + std::to_string(err_v) + ")");
    }
}

void Simulation::save_checkpoint(const std::string& path, double time, int step) const {
    std::ofstream out(path, std::ios::binary);
    out.write(reinterpret_cast<const char*>(&time), sizeof(double));
    out.write(reinterpret_cast<const char*>(&step), sizeof(int));

    auto write_vec = [&](const Eigen::VectorXd& v) {
        int n = static_cast<int>(v.size());
        out.write(reinterpret_cast<const char*>(&n), sizeof(int));
        out.write(reinterpret_cast<const char*>(v.data()), n * sizeof(double));
    };

    write_vec(polarization.get_Px());   write_vec(polarization.get_Py());
    write_vec(polarization.get_Px_n()); write_vec(polarization.get_Py_n());  
    write_vec(mechanics.get_ux());      write_vec(mechanics.get_uy());
    write_vec(electrostatics.get_phi());
    write_vec(fracture.get_v());        write_vec(fracture.get_v_n());       
}

void Simulation::load_checkpoint(const std::string& path, double& time, int& step) {
    std::ifstream in(path, std::ios::binary);
    in.read(reinterpret_cast<char*>(&time), sizeof(double));
    in.read(reinterpret_cast<char*>(&step), sizeof(int));

    auto read_vec = [&](Eigen::VectorXd& v) {
        int n; in.read(reinterpret_cast<char*>(&n), sizeof(int));
        v.resize(n);
        in.read(reinterpret_cast<char*>(v.data()), n * sizeof(double));
    };

    Eigen::VectorXd Px, Py, Px_n, Py_n, ux, uy, phi, v, v_n;
    read_vec(Px); read_vec(Py); read_vec(Px_n); read_vec(Py_n);
    read_vec(ux); read_vec(uy); read_vec(phi); read_vec(v); read_vec(v_n);

    polarization.set_state(Px, Py, Px_n, Py_n);
    mechanics.set_state(ux, uy);
    electrostatics.set_state(phi);
    fracture.set_state(v, v_n);
}


void Simulation::run() {
    int total_steps = static_cast<int>(config.get_total_time() / config.get_dt());
 
    double time = 0.0;
    int start_step = 0;
 
    // --- Restart eventuel ---
    if (!restart_path_.empty()) {
        load_checkpoint(restart_path_, time, start_step);
    }
 
    ProgressBar progressBar(total_steps, "[SIMULATION]");
 
    if (config.get_chrono_run()) { chrono.start(); }
 
    Logger::info("Starting simulation...");
    std::string initial_time = chrono.get_datetime_string();
 
    for (int step = start_step + 1; step <= total_steps; ++step) {
        time += config.get_dt();
 
        ComputeOneStepPhysics(time);
 
        diagnostics.record(static_cast<double>(step), time, polarization, mechanics, fracture, electrostatics, math);
        diagnostics.append_csv(diagnostics_csv_path);
 
        if (step % config.get_save_frequency() == 0) {
            std::string filename = config.getOutputDir() + "/VTK_" + initial_time +
                                   "/multiphysics_results" + std::to_string(step) + ".vtk";
 
            diagnostics.compute_nodal_energies(polarization, mechanics, fracture, electrostatics, math);
 
            exporter.exportMultiPhysicsVTK(
                filename,
                {"v", "U_energy", "W_energy", "chi_energy", "elec_energy", "surf_energy"},
                {&fracture.get_v(),
                 &diagnostics.get_U_nodal(), &diagnostics.get_W_nodal(), &diagnostics.get_chi_nodal(),
                 &diagnostics.get_elec_nodal(), &diagnostics.get_surf_nodal()},
                {"P", "U", "E"},
                {&polarization.get_Px(), &mechanics.get_ux(), &electrostatics.get_Ex()},
                {&polarization.get_Py(), &mechanics.get_uy(), &electrostatics.get_Ey()}
            );
        }
 
        // --- Checkpoint periodique : toutes les 20 pas, ecrase le meme fichier
        // (evite de saturer le disque). Adapte la frequence a ton besoin. ---
        if (step % 20 == 0 || step == total_steps) {
            save_checkpoint(config.getOutputDir() + "/checkpoint.bin", time, step);
        }
 
        progressBar.update(step, 0.0);
    }
 
    progressBar.finish();
 
    diagnostics.write_csv(diagnostics_csv_path);
 
    if (config.get_chrono_run()) {
        chrono.stop();
        Logger::time("Simulation completed in ", chrono.elapsed_ms());
    }
    Logger::info("Simulation completed.");
}
