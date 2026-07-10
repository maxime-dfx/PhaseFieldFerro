#include "IO/Datafile.h"
#include "Utils/Logger.h" 
#include <stdexcept>

// Implémentation du template utilitaire (placé ici car utilisé uniquement en interne par cette classe)
template <typename T>
T Datafile::get_val(const toml::table& config, const std::string& section, const std::string& key, T default_value) {
    if (auto node = config[section][key]) {
        return node.value<T>().value_or(default_value);
    }
    return default_value;
}

// Constructeur : orchestration du parsing
Datafile::Datafile(const std::string& filename) {
    try {
        toml::table raw_config = toml::parse_file(filename);
        
        parse_simulation(raw_config);
        parse_mesh(raw_config);
        parse_chrono(raw_config);
        parse_physics(raw_config);
        parse_material(raw_config);
        parse_boundaries(raw_config);

    } catch (const toml::parse_error& err) {
        throw std::runtime_error("Erreur de parsing TOML (" + filename + ") : " + std::string(err.description()));
    }
}

// ------------------------------------------------------------
//  MÉTHODES DE PARSING DÉDIÉES
// ------------------------------------------------------------

void Datafile::parse_simulation(const toml::table& config) {
    simulation.output_dir = get_val<std::string>(config, "output", "OutputDir", "../results");
    simulation.mesh_create_file = get_val<std::string>(config, "output", "MeshCreateFile", "mesh");
    simulation.debug_enabled = get_val<bool>(config, "debug", "enabled", false);
    
    simulation.total_time = get_val<int>(config, "simulation", "total_time", 100);
    simulation.dt = get_val<double>(config, "simulation", "dt", 0.1);
    simulation.save_frequency = get_val<int>(config, "simulation", "save_frequency", 10);
    simulation.tol_ferro = get_val<double>(config, "simulation", "tolerance_ferro", 1e-3);
    simulation.tol_vfield = get_val<double>(config, "simulation", "tolerance_vfield", 1e-6);
}

void Datafile::parse_mesh(const toml::table& config) {
    mesh.calcul_mesh = get_val<bool>(config, "mesh", "calcul_mesh", true);
    mesh.get_mesh_file = get_val<std::string>(config, "mesh", "get_mesh_file", "../input/rectangle.msh");
    mesh.Lx = get_val<double>(config, "mesh", "L_x", 200.0);
    mesh.Ly = get_val<double>(config, "mesh", "L_y", 200.0);
    mesh.nx = get_val<int>(config, "mesh", "n_x", 200);
    mesh.ny = get_val<int>(config, "mesh", "n_y", 200);
    
    // Pré-calcul des deltas
    mesh.dx = mesh.Lx / (mesh.nx - 1); 
    mesh.dy = mesh.Ly / (mesh.ny - 1);
    
    std::string type = get_val<std::string>(config, "mesh", "element_type", "QUAD");
    mesh.element_type = (type == "TRIANGLE") ? ElementType::TRIANGLE : ElementType::QUAD;
}

void Datafile::parse_chrono(const toml::table& config) {
    chrono.mesh = get_val<bool>(config, "chrono", "chrono_mesh", false);
    chrono.polarization = get_val<bool>(config, "chrono", "chrono_Polarization", false);
    chrono.mechanics = get_val<bool>(config, "chrono", "chrono_Mechanics", false);
    chrono.electrostatics = get_val<bool>(config, "chrono", "chrono_Electrostatics", false);
    chrono.fracture = get_val<bool>(config, "chrono", "chrono_Fracture", false);
    chrono.run = get_val<bool>(config, "chrono", "chrono_run", false);
}

void Datafile::parse_physics(const toml::table& config) {
    // 1. Activation des modules
    physics_toggle.polarization = get_val<bool>(config, "physics", "enable_polarization", true);
    physics_toggle.electrostatics = get_val<bool>(config, "physics", "enable_electrostatics", true);
    physics_toggle.mechanics = get_val<bool>(config, "physics", "enable_mecanics", true);
    physics_toggle.fracture = get_val<bool>(config, "physics", "enable_fracture", true);

    // 2. Fonction lambda pour factoriser le code d'initialisation (DRY)
    auto parse_init = [&](const std::string& section, const std::string& key_type, const std::string& key_x, const std::string& key_y) {
        PhysicsInitConfig c;
        std::string t = get_val<std::string>(config, section, key_type, "UNIFORM");
        c.type = (t == "RANDOM") ? InitializationType::RANDOM : InitializationType::UNIFORM;
        c.val_x_0 = get_val<double>(config, section, key_x, 0.0);
        c.val_y_0 = get_val<double>(config, section, key_y, 0.0);
        return c;
    };

    // 3. Récupération des initialisations
    polarization = parse_init("polarization", "initial_polarization", "Px_0", "Py_0");
    mechanics = parse_init("mechanics", "initial_mechanics", "ux_0", "uy_0");
    electrostatics = parse_init("electrostatics", "initial_electrostatics", "Ex_0", "Ey_0");

    // 4. Mode spécifique à la fracture
    std::string frac_mode = get_val<std::string>(config, "fracture", "mode", "PERMEABLE");
    if (frac_mode == "IMPERMEABLE") fracture_mode = CrackBCType::IMPERMEABLE;
    else if (frac_mode == "CONDUCTIVE") fracture_mode = CrackBCType::CONDUCTIVE;
    else fracture_mode = CrackBCType::PERMEABLE;
}

void Datafile::parse_material(const toml::table& config) {
    material.a0 = get_val<double>(config, "material", "a0", 1.0);
    material.b1 = get_val<double>(config, "material", "b1", 1.0);
    material.b2 = get_val<double>(config, "material", "b2", 1.0);
    material.b3 = get_val<double>(config, "material", "b3", 1.0);
    material.c1 = get_val<double>(config, "material", "c1", 1.0);
    material.c2 = get_val<double>(config, "material", "c2", 1.0);
    material.c3 = get_val<double>(config, "material", "c3", 1.0);
    material.eps0 = get_val<double>(config, "material", "eps0", 8.854e-12);
    material.t = get_val<double>(config, "material", "t", 300.0);
    material.P0 = get_val<double>(config, "material", "P0", 1.0);
    material.mu_p = get_val<double>(config, "material", "mu_p", 1.0);
    material.mu_v = get_val<double>(config, "material", "mu_v", 1.0);
    material.xi = get_val<double>(config, "material", "xi", 1.0);
    material.c0 = get_val<double>(config, "material", "c0", 1.0);
    
    // Paramètres Alpha
    material.alpha_1 = get_val<double>(config, "material", "alpha_1", 1.0);
    material.alpha_11 = get_val<double>(config, "material", "alpha_11", 1.0);
    material.alpha_111 = get_val<double>(config, "material", "alpha_111", 1.0);
    material.alpha_1111 = get_val<double>(config, "material", "alpha_1111", 1.0);
    material.alpha_12 = get_val<double>(config, "material", "alpha_12", 1.0);
    material.alpha_112 = get_val<double>(config, "material", "alpha_112", 1.0);
    material.alpha_1112 = get_val<double>(config, "material", "alpha_1112", 1.0);
    material.alpha_1122 = get_val<double>(config, "material", "alpha_1122", 1.0);
    
    // Autres paramètres fracture
    material.eta_k = get_val<double>(config, "material", "eta_k", 1.0);
    material.Gc = get_val<double>(config, "material", "Gc", 1.0);
    material.kappa = get_val<double>(config, "material", "kappa", 1.0);
}

void Datafile::parse_boundaries(const toml::table& config) {
    boundary_rules.clear();
    if (auto arr = config["boundary_rule"].as_array()) {
        // std::cout << "[DEBUG TOML] J'ai trouve " << arr->size() << " blocs [[boundary_rule]] dans le fichier." << std::endl;
        for (auto& node : *arr) {
            if (auto tbl = node.as_table()) {
                BoundaryRuleConfig r;
                r.field     = (*tbl)["field"].value_or<std::string>("");
                r.bc_type   = (*tbl)["bc_type"].value_or<std::string>("DIRICHLET");
                
                r.shape     = (*tbl)["shape"].value_or<std::string>("edge");
                r.edge_name = (*tbl)["edge_name"].value_or<std::string>("");
                r.cx        = (*tbl)["cx"].value_or<double>(0.0);
                r.cy        = (*tbl)["cy"].value_or<double>(0.0);
                r.radius    = (*tbl)["radius"].value_or<double>(0.0);
                r.px        = (*tbl)["px"].value_or<double>(0.0);
                r.py        = (*tbl)["py"].value_or<double>(0.0);
                r.xmin      = (*tbl)["xmin"].value_or<double>(0.0);
                r.xmax      = (*tbl)["xmax"].value_or<double>(0.0);
                r.ymin      = (*tbl)["ymin"].value_or<double>(0.0);
                r.ymax      = (*tbl)["ymax"].value_or<double>(0.0);

                r.profile   = (*tbl)["profile"].value_or<std::string>("constant");
                r.val       = (*tbl)["val"].value_or<double>(0.0);
                r.val_start = (*tbl)["val_start"].value_or<double>(0.0);
                r.val_end   = (*tbl)["val_end"].value_or<double>(0.0);
                r.t_end     = (*tbl)["t_end"].value_or<double>(1.0);
                r.amplitude = (*tbl)["amplitude"].value_or<double>(0.0);
                r.frequency = (*tbl)["frequency"].value_or<double>(0.0);
                r.offset    = (*tbl)["offset"].value_or<double>(0.0);
                r.P0        = (*tbl)["P0"].value_or<double>(0.0);
                r.x0        = (*tbl)["x0"].value_or<double>(0.0);
                r.epsilon   = (*tbl)["epsilon"].value_or<double>(1e-3);
                r.max_val   = (*tbl)["max_val"].value_or<double>(0.0);
                r.y_center  = (*tbl)["y_center"].value_or<double>(0.0);
                r.width     = (*tbl)["width"].value_or<double>(0.0);
                r.k         = (*tbl)["k"].value_or<double>(0.0);
                r.omega     = (*tbl)["omega"].value_or<double>(0.0);
                
                boundary_rules.push_back(r);
            }
        }
    }   else {
        // std::cout << "[DEBUG TOML] ERREUR : Impossible de trouver [[boundary_rule]] dans le fichier TOML !" << std::endl;
    }
}