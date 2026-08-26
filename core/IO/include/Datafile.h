#pragma once

#include <string>
#include <vector>
#include <iostream>
#include "librairies/toml.hpp"
#include "Utils/include/Logger.h"
#include "IO/include/ConfigTypes.h"

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
    std::vector<MaterialConfig> materials;
    CrystalConfig crystal;
    std::vector<BoundaryRuleConfig> boundary_rules;

    Datafile(const std::string& filename);

private:
    void parse_simulation(const toml::table& config);
    void parse_mesh(const toml::table& config);
    void parse_chrono(const toml::table& config);
    void parse_physics(const toml::table& config);
    void parse_material(const toml::table& config);
    void parse_materials_list(const toml::table& config);
    void parse_crystal(const toml::table& config);
    void parse_fracture(const toml::table& config); 
    void parse_boundaries(const toml::table& config);

    template <typename T>
    T get_val(const toml::table& config, const std::string& section, const std::string& key, T default_value);

    template <typename T>
    T get_val_from_table(const toml::table& tbl, const std::string& key, T default_value);
};