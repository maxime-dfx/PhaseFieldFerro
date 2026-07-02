/*
 * ============================================================
 *  Datafile.h  —  Fichier de configuration pour le code de simulation
 *  Basé sur l'article : Acta Materialia 2011, Abdollahi & Arias
 *  
 *  Paramètres normalisés conformément au Tableau 1 et Section 3.1
 * ============================================================
 */

#pragma once

#include <string>
#include <vector>
#include <iostream>
#include "Utils/toml.hpp"
#include "Utils/Logger.h"

// ------------------------------------------------------------
//  ENUMÉRATIONS
// ------------------------------------------------------------

enum class CrackBCType {
    PERMEABLE,      
    IMPERMEABLE,    
    CONDUCTIVE      
};

enum class ElementType {
    TRIANGLE,
    QUAD
};

enum class PolarizationInitializationType {
    UNIFORM,
    RANDOM
};

enum class MechanicsInitializationType {
    UNIFORM,
    RANDOM
};

enum class ElectrostaticsInitializationType {
    UNIFORM,
    RANDOM
};

// ------------------------------------------------------------
//  STRUCTURE DE CONDITION AUX LIMITES
// ------------------------------------------------------------


struct BoundaryRuleConfig {
    std::string field;       
    std::string bc_type;     
    
    std::string shape;       
    std::string edge_name;   
    double cx = 0.0, cy = 0.0, radius = 0.0;    
    double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;  
    double px = 0.0, py = 0.0;
    
    std::string profile;     
    double val = 0.0;       
    
    // Profil "ramp"
    double val_start = 0.0, val_end = 0.0, t_end = 1.0;
    
    // Profil "sinusoidal"
    double amplitude = 0.0, frequency = 0.0, offset = 0.0;
    
    // Profil "gaussian"
    double max_val = 0.0, y_center = 0.0, width = 0.0;
    
    // Profil "step" (fonction sigmoïde)
    double P0 = 0.0, x0 = 0.0, epsilon = 1e-3;
    double k = 0.0, omega = 0.0;
};

// ------------------------------------------------------------
//  FONCTION UTILITAIRE TOML
// ------------------------------------------------------------

template <typename T>
T get_val(const toml::table& config, const std::string& section, 
          const std::string& key, T default_value) {
    if (auto node = config[section][key]) {
        return node.value<T>().value_or(default_value);
    }
    return default_value;
}

// ------------------------------------------------------------
//  CLASSE DATAFILE
// ------------------------------------------------------------

class Datafile {
private:
    toml::table m_config;

public:
    // ============================================================
    //  CONSTRUCTEUR
    // ============================================================
    Datafile(const std::string& filename) {
        try {
            m_config = toml::parse_file(filename);
        } catch (const toml::parse_error& err) {
            // Logger::error("Erreur de parsing du fichier TOML : " + std::string(err.description()));
            throw;
        }
    }

    // ============================================================
    //  SORTIE ET RÉSULTATS
    // ============================================================
    std::string getOutputDir() const {
        return get_val<std::string>(m_config, "output", "OutputDir", "../results");
    }
    std::string getMeshCreateFile() const {
        return get_val<std::string>(m_config, "output", "MeshCreateFile", "mesh");
    }

    // ============================================================
    //  DEBUG
    // ============================================================
    bool debug_enabled() const {
        return get_val<bool>(m_config, "debug", "enabled", false);
    }

    // ============================================================
    //  CHRONOMÉTRAGE
    // ============================================================

    bool get_chrono_mesh() const {
        return get_val<bool>(m_config, "chrono", "chrono_mesh", false);
    }
    bool get_chrono__Polarization() const {
        return get_val<bool>(m_config, "chrono", "chrono_Polarization", false);
    }
    bool get_chrono_Mechanics() const {
        return get_val<bool>(m_config, "chrono", "chrono_Mechanics", false);
    }
    bool get_chrono_Electrostatics() const {
        return get_val<bool>(m_config, "chrono", "chrono_Electrostatics", false);
    }
    bool get_chrono_Fracture() const {
        return get_val<bool>(m_config, "chrono", "chrono_Fracture", false);
    }
    bool get_chrono_run() const {
        return get_val<bool>(m_config, "chrono", "chrono_run", false);
    }

    // ============================================================
    //  PHYSIQUE ACTIVÉE
    // ============================================================
    bool enable_polarization() const {
        return get_val<bool>(m_config, "physics", "enable_polarization", true);
    }
    bool enable_electrostatics() const {
        return get_val<bool>(m_config, "physics", "enable_electrostatics", true);
    }
    bool enable_mecanics() const {
        return get_val<bool>(m_config, "physics", "enable_mecanics", true);
    }
    bool enable_fracture() const {
        return get_val<bool>(m_config, "physics", "enable_fracture", true);
    }

    // ============================================================
    //  SIMULATION
    // ============================================================
    int get_total_time() const {
        return get_val<int>(m_config, "simulation", "total_time", 100);
    }
    double get_dt() const {
        return get_val<double>(m_config, "simulation", "dt", 0.1);
    }
    double get_dt_relax() const { 
        return get_val<double>(m_config, "simulation", "dt_relax", 0.1);
    }
    int get_save_frequency() const {
        return get_val<int>(m_config, "simulation", "save_frequency", 10);
    }
    double get_tol_ferro() const {
        return get_val<double>(m_config, "simulation", "tolerance_ferro", 1e-3);
    }
    double get_tol_vfield() const {
        return get_val<double>(m_config, "simulation", "tolerance_vfield", 1e-6);
    }

    // ============================================================
    //  MAILLAGE
    // ============================================================
    bool calcul_mesh() const {
        return get_val<bool>(m_config, "mesh", "calcul_mesh", true);
    }
    double get_L_x() const {
        return get_val<double>(m_config, "mesh", "L_x", 200.0);   
    }
    double get_L_y() const {
        return get_val<double>(m_config, "mesh", "L_y", 200.0);   
    }
    int get_n_x() const {
        return get_val<int>(m_config, "mesh", "n_x", 200);
    }
    int get_n_y() const {
        return get_val<int>(m_config, "mesh", "n_y", 200);
    }
    double get_dx() const {
        return get_L_x() / get_n_x();
    }
    double get_dy() const {
        return get_L_y() / get_n_y();
    }    
    ElementType get_element_type() const {
        std::string t = get_val<std::string>(m_config, "mesh", "element_type", "QUAD");
        if (t == "TRIANGLE") return ElementType::TRIANGLE;
        return ElementType::QUAD;
    }

    // ============================================================
    //  POLARISATION
    // ============================================================

    PolarizationInitializationType get_initial_polarization() const {
        std::string t = get_val<std::string>(m_config, "polarization", "initial_polarization", "UNIFORM");
        if (t == "UNIFORM") return PolarizationInitializationType::UNIFORM;
        if (t == "RANDOM") return PolarizationInitializationType::RANDOM;
        // Logger::error("Type d'initialisation de polarisation non reconnu.");
        return PolarizationInitializationType::UNIFORM;
    }
    double get_Px_0() const {
        return get_val<double>(m_config, "polarization", "Px_0", 0.0);
    }
    double get_Py_0() const {
        return get_val<double>(m_config, "polarization", "Py_0", 0.0);
    }

    // ============================================================
    //  MECHANICS
    // ============================================================

    MechanicsInitializationType get_initial_mechanics() const {
        std::string t = get_val<std::string>(m_config, "mechanics", "initial_mechanics", "UNIFORM");
        if (t == "UNIFORM") return MechanicsInitializationType::UNIFORM;
        if (t == "RANDOM") return MechanicsInitializationType::RANDOM;
        // Logger::error("Type d'initialisation de mechanics non reconnu.");
        return MechanicsInitializationType::UNIFORM;
    }
    double get_ux_0() const {
        return get_val<double>(m_config, "mechanics", "ux_0", 0.0);
    }
    double get_uy_0() const {
        return get_val<double>(m_config, "mechanics", "uy_0", 0.0);
    }

    // ============================================================
    //  ELECTROSTATICS
    // ============================================================
    ElectrostaticsInitializationType get_initial_electrostatics() const {
        std::string t = get_val<std::string>(m_config, "electrostatics", "initial_electrostatics", "UNIFORM");
        if (t == "UNIFORM") return ElectrostaticsInitializationType::UNIFORM;
        if (t == "RANDOM") return ElectrostaticsInitializationType::RANDOM;
        // Logger::error("Type d'initialisation de electrostatics non reconnu.");
        return ElectrostaticsInitializationType::UNIFORM;
    }
    double get_Ex_0() const {
        return get_val<double>(m_config, "electrostatics", "Ex_0", 0.0);
    }
    double get_Ey_0() const {
        return get_val<double>(m_config, "electrostatics", "Ey_0", 0.0);
    }

    // ============================================================
    //  FRACTURE 
    // ============================================================

    CrackBCType get_fracture_mode() const {
        std::string t = get_val<std::string>(m_config, "fracture", "mode", "PERMEABLE");
        if (t == "PERMEABLE") return CrackBCType::PERMEABLE;
        if (t == "IMPERMEABLE") return CrackBCType::IMPERMEABLE;
        // Logger::error("Type de condition aux limites de fracture non reconnu.");
        return CrackBCType::PERMEABLE;
    }

    // ============================================================
    //  MATÉRIAU
    // ============================================================
    double get_a0() const {
        return get_val<double>(m_config, "material", "a0", 1.0);
    }
    double get_b1() const {
        return get_val<double>(m_config, "material", "b1", 1.0);
    }
    double get_b2() const {
        return get_val<double>(m_config, "material", "b2", 1.0);
    }
    double get_b3() const {
        return get_val<double>(m_config, "material", "b3", 1.0);
    }
    double get_c1() const {
        return get_val<double>(m_config, "material", "c1", 1.0);
    }
    double get_c2() const {
        return get_val<double>(m_config, "material", "c2", 1.0);
    }
    double get_c3() const {
        return get_val<double>(m_config, "material", "c3", 1.0);
    }
    double get_eps0() const {
        return get_val<double>(m_config, "material", "eps0", 8.854e-12);
    }
    double get_t() const {
        return get_val<double>(m_config, "material", "t", 300.0);
    }
    double get_P0() const {
        return get_val<double>(m_config, "material", "P0", 1.0);
    }
    double get_mu_p() const {
        return get_val<double>(m_config, "material", "mu_p", 1.0);
    }
    double get_mu_v() const {
        return get_val<double>(m_config, "material", "mu_v", 1.0);
    }
    double get_xi() const {
        return get_val<double>(m_config, "material", "xi", 1.0);
    }
    double get_c0() const {
        return get_val<double>(m_config, "material", "c0", 1.0);
    }
    double get_alpha_1() const {
        return get_val<double>(m_config, "material", "alpha_1", 1.0);
    }
    double get_alpha_11() const {
        return get_val<double>(m_config, "material", "alpha_11", 1.0);
    }
    double get_alpha_111() const {
        return get_val<double>(m_config, "material", "alpha_111", 1.0);
    }
    double get_alpha_1111() const {
        return get_val<double>(m_config, "material", "alpha_1111", 1.0);
    }
    double get_alpha_12() const {
        return get_val<double>(m_config, "material", "alpha_12", 1.0);
    }
    double get_alpha_112() const {
        return get_val<double>(m_config, "material", "alpha_112", 1.0);
    }
    double get_alpha_1112() const {
        return get_val<double>(m_config, "material", "alpha_1112", 1.0);
    }
    double get_alpha_1122() const {
        return get_val<double>(m_config, "material", "alpha_1122", 1.0);
    }
    double get_eta_k() const {
        return get_val<double>(m_config, "material", "eta_k", 1.0);
    }
    double get_Gc() const {
        return get_val<double>(m_config, "material", "Gc", 1.0);
    }
    double get_kappa() const {
        return get_val<double>(m_config, "material", "kappa", 1.0);
    }

    // ============================================================
    //  CONDITIONS AUX LIMITES
    // ============================================================
    
    std::vector<BoundaryRuleConfig> get_boundary_rules() const;

    bool get_fix_top_y_ramp() const {
        return false;
    }
};