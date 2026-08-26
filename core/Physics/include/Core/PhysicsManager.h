#pragma once
#include <vector>
#include <memory>
#include "IO/include/Datafile.h"
#include "Mesh/include/Mesh.h"
#include "Physics/Core/BoundaryManager.h"
#include "Materials/Core/MaterialManager.h"
#include "Materials/Microstructure/Polycrystal.h"
#include "Physics/include/Core/PhysicsConcepts.h"
#include "Physics/include/Core/SystemMetrics.h"
#include "Physics/include/Core/EnergyPostProcessor.h"

class PhysicsManager {
private:
    const Datafile& config;
    const Mesh& mesh;
    const MaterialManager& materials_manager;
    const Polycrystal& polycrystal;
    BoundaryManager& boundary_manager;

    // LE REGISTRE POLYMORPHIQUE : Fini les instances en dur !
    std::vector<std::unique_ptr<IPhysicsModule>> m_active_modules;
    PhysicsState m_current_state;

    SystemMetricsCalculator m_metrics_calculator;
    bool m_last_step_had_solve_failure = false;

    mutable EnergyPostProcessor m_energy_post;

    void update_state_pointers(double time, double dt, double dt_load);

public:
    PhysicsManager(const Datafile& config, const Mesh& mesh, BoundaryManager& boundary_manager,
                   const MaterialManager& materials_manager, const Polycrystal& polycrystal);

    int compute_one_step_physics(double time, double dt);
    void update_physics_history();
    SystemMetrics compute_system_metrics() const {
        m_energy_post.compute(m_current_state, mesh);
        return m_metrics_calculator.compute(m_current_state, m_last_step_had_solve_failure);
    }
    bool last_step_had_solve_failure() const { return m_last_step_had_solve_failure; }

    // Expose l'etat global (le Blackboard) pour ResultsExporter
    const PhysicsState& get_current_state() const { return m_current_state; }

    // Expose le maillage pour la couche d'export/visualisation (IOManager)
    const Mesh& get_mesh() const { return mesh; }
};
