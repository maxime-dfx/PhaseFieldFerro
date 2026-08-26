#pragma once
//
// SystemMetrics.h
// -----------------
// Diagnostics : le POD SystemMetrics (valeurs scalaires exportees vers
// energies.csv) et le SystemMetricsCalculator qui les calcule a partir du
// PhysicsState courant. Fusionnes ici car SystemMetricsCalculator n'a
// aucune raison d'exister sans SystemMetrics.

#include "Physics/include/Core/PhysicsConcepts.h"
#include "Materials/Core/MaterialManager.h"
#include "IO/include/Datafile.h"

struct SystemMetrics {
    double U_total = 0.0;
    double W_total = 0.0;
    double chi_total = 0.0;
    double electric_total = 0.0;
    double bulk_enthalpy = 0.0;
    double surface_energy = 0.0;
    double total_energy = 0.0;
    double v_min = 1.0;
    double v_max = 1.0;
    bool solve_failed = false;
};

class SystemMetricsCalculator {
private:
    const Datafile& config;
    const Mesh& mesh;
    const MaterialManager& materials_manager;

public:
    SystemMetricsCalculator(const Datafile& config, const Mesh& mesh, const MaterialManager& materials_manager)
        : config(config), mesh(mesh), materials_manager(materials_manager) {}

    SystemMetrics compute(const PhysicsState& state, bool solve_failed) const;
};
