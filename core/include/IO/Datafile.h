#pragma once

#include <string>
#include <vector>
#include <iostream>
#include "Utils/toml.hpp"
#include "Utils/Logger.h"
#include "IO/ConfigTypes.h"

// ------------------------------------------------------------
//  CLASSE DATAFILE REFACTORISÉE
// ------------------------------------------------------------
class Datafile {
public:
    SimulationConfig simulation;
    MeshConfig mesh;
    ChronoConfig chrono;
    PhysicsToggle physics_toggle;
    PhysicsInitConfig polarization;
    PhysicsInitConfig mechanics;
    PhysicsInitConfig electrostatics;
    FractureConfig fracture;          
    MaterialConfig material;
    std::vector<BoundaryRuleConfig> boundary_rules;

    Datafile(const std::string& filename);

private:
    void parse_simulation(const toml::table& config);
    void parse_mesh(const toml::table& config);
    void parse_chrono(const toml::table& config);
    void parse_physics(const toml::table& config);
    void parse_material(const toml::table& config);
    void parse_fracture(const toml::table& config); 
    void parse_boundaries(const toml::table& config);

    template <typename T>
    T get_val(const toml::table& config, const std::string& section, const std::string& key, T default_value);
};