#include "Physics/include/Core/PhysicsManager.h"
#include "Physics/include/Core/PhysicsConcepts.h"
#include "Utils/include/Profiling.h"

#include "Physics/include/Modules/Mechanics.h"
#include "Physics/include/Modules/Electrostatics.h"
#include "Physics/include/Modules/Fracture.h"
#include "Physics/include/Modules/Polarization.h"

PhysicsManager::PhysicsManager(const Datafile& config_in, const Mesh& mesh_in, BoundaryManager& boundary_manager_in,
                               const MaterialManager& materials_manager_in, const Polycrystal& polycrystal_in)
    : config(config_in),
      mesh(mesh_in),
      materials_manager(materials_manager_in),
      polycrystal(polycrystal_in),
      boundary_manager(boundary_manager_in),
      m_metrics_calculator(config_in, mesh_in, materials_manager_in)
{
    m_current_state.materials_manager = &materials_manager;
    m_current_state.polycrystal = &polycrystal;
    m_current_state.config = &config;
    m_energy_post.register_fields(m_current_state);
    
    // --- Utilise la deduction de types (C++17) ---
    auto register_module = [&](const std::string& name, const auto& mod_config, auto eq, auto mapper) {
        SolverConfig cfg;
        cfg.solver_type = mod_config.solver_type;
        cfg.preconditioner_type = mod_config.preconditioner_type;
        cfg.tolerance = mod_config.solver_tol;
        cfg.max_iterations = mod_config.solver_max_iter;

        // Deduction automatique des types pour forcer l'instanciation du template
        using EqType = decltype(eq);
        using MapType = decltype(mapper);

        auto mod = std::make_unique<GenericPhysicsModule<EqType, MapType>>(
            name, cfg, mesh_in, std::move(eq), std::move(mapper)
        );

        mod->get_mapper().register_fields(m_current_state);
        
        auto* raw_ptr = mod.get();
        m_active_modules.push_back(std::move(mod));
        return raw_ptr;
    };

    // --- INSTANCIATIONS DIRECTES SUR LA PILE ---
    if (config.physics_toggle.polarization) {
        register_module("Polarization", config.polarization,
                        PolarizationEquation(), // Instancie par valeur !
                        PolarizationDofMapper(mesh_in, boundary_manager, polycrystal, config_in, materials_manager_in));
    }
    
    if (config.physics_toggle.mechanics) {
        auto* mod = register_module("Mechanics", config.mechanics,
                                    MechanicsEquation(),
                                    MechanicsDofMapper(mesh_in, boundary_manager));
        mod->add_post_processor(std::make_unique<StressPostProcessor>(), m_current_state);
    }
    
    if (config.physics_toggle.electrostatics) {
        auto* mod = register_module("Electrostatics", config.electrostatics,
                                    ElectrostaticsEquation(),
                                    ElectrostaticsDofMapper(mesh_in, boundary_manager));
        mod->add_post_processor(std::make_unique<ElectrostaticsPostProcessor>(), m_current_state);
    }
    
    if (config.physics_toggle.fracture) {
        register_module("Fracture", config.fracture,
                        FractureEquation(),
                        FractureDofMapper(mesh_in, boundary_manager, config));
    }
}

void PhysicsManager::update_state_pointers(double time, double dt, double dt_load) {
    m_current_state.time = time;
    m_current_state.dt = dt;
    m_current_state.dt_load = dt_load;
}

int PhysicsManager::compute_one_step_physics(double time, double dt) {
    PROFILE_ZONE_NC("compute_one_step_physics", PROFILE_COLOR_SEQUENTIAL);
    
    boundary_manager.update_time(time);
    update_state_pointers(time, config.simulation.dt_relax, dt);
    
    for (auto& module : m_active_modules) module->save_previous_iteration();
    
    int m = 0;
    bool any_solve_failed = false;
    bool all_converged = false;
    
    do {
        PROFILE_ZONE_NC("Picard_Iteration", PROFILE_COLOR_SEQUENTIAL);
        m++;
        any_solve_failed = false;
        all_converged = true;
        
        Logger::debug("---- [PhysicsManager] Picard m=", m,
                       " (t=", time, ", dt_relax=", config.simulation.dt_relax,
                       ") : debut resolution des ", m_active_modules.size(), " modules ----");
        
        for (auto& module : m_active_modules) {
            module->compute_step(m_current_state);
            bool failed = module->last_solve_failed();
            if (failed) any_solve_failed = true;
            Logger::debug("  [PhysicsManager] m=", m, " module=", module->get_name(),
                           " solve_failed=", (failed ? "OUI" : "non"));
        }
        
        for (const auto& kv : m_current_state.fields) {
            if (!kv.second) continue;
            const Eigen::VectorXd& f = *kv.second;
            if (f.size() == 0) continue;
            Logger::debug("  [PhysicsManager] m=", m, " champ=", kv.first,
                           " min=", f.minCoeff(), " max=", f.maxCoeff(),
                           " norm=", f.norm());
        }
        
        for (auto& module : m_active_modules) {
            double err = module->calculate_error();
            bool this_module_converged = true;

            if (module->get_name() == "Polarization") {
                this_module_converged = (err <= config.simulation.tol_ferro);
                if (!this_module_converged) all_converged = false;
            } else if (module->get_name() == "Fracture") {
                this_module_converged = (err <= config.simulation.tol_vfield);
                if (!this_module_converged) all_converged = false;
            }
            Logger::debug("  [PhysicsManager] m=", m, " module=", module->get_name(),
                           " error=", err,
                           " converged=", (this_module_converged ? "oui" : "NON"));

            module->save_previous_iteration();
        }
        
        Logger::debug("---- [PhysicsManager] Picard m=", m,
                       " : all_converged=", (all_converged ? "oui" : "non"),
                       " (min_iter=", config.simulation.min_iter,
                       ", max_iter=", config.simulation.max_iter, ") ----");
        
    } while ((!all_converged || m < config.simulation.min_iter) && m < config.simulation.max_iter);
    
    m_last_step_had_solve_failure = any_solve_failed;
    Logger::debug("==== [PhysicsManager] fin compute_one_step_physics : m final=", m,
                   " any_solve_failed=", (any_solve_failed ? "oui" : "non"),
                   " all_converged=", (all_converged ? "oui" : "non"), " ====");

    if (m >= config.simulation.max_iter) return -1;
    if (any_solve_failed) return -2;
    
    return m;
}

void PhysicsManager::update_physics_history() {
    for (auto& module : m_active_modules) module->update_history();
}
