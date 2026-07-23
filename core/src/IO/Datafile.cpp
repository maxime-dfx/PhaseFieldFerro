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
        parse_fracture(raw_config);       
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
    simulation.dt_relax = get_val<double>(config, "simulation", "dt_relax", 0.1);
    simulation.save_frequency = get_val<int>(config, "simulation", "save_frequency", 10);
    simulation.tol_ferro = get_val<double>(config, "simulation", "tolerance_ferro", 1e-3);
    simulation.tol_vfield = get_val<double>(config, "simulation", "tolerance_vfield", 1e-6);
    simulation.min_iter = get_val<int>(config, "simulation", "min_iter", 20);
    simulation.max_iter = get_val<int>(config, "simulation", "max_iter", 50);
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
    mesh.element_type = (type == "TRIANGLE") ? ElementType::TRIANGLE3 : ElementType::QUAD4;
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
    physics_toggle.polarization = get_val<bool>(config, "physics", "enable_polarization", true);
    physics_toggle.electrostatics = get_val<bool>(config, "physics", "enable_electrostatics", true);
    physics_toggle.mechanics = get_val<bool>(config, "physics", "enable_mecanics", true);
    physics_toggle.fracture = get_val<bool>(config, "physics", "enable_fracture", true);

    auto parse_init = [&](const std::string& section, const std::string& key_type, const std::string& key_x, const std::string& key_y) {
        PhysicsInitConfig c;
        std::string t = get_val<std::string>(config, section, key_type, "UNIFORM");
        c.type = (t == "RANDOM") ? InitializationType::RANDOM : InitializationType::UNIFORM;
        c.val_x_0 = get_val<double>(config, section, key_x, 0.0);
        c.val_y_0 = get_val<double>(config, section, key_y, 0.0);
        return c;
    };

    polarization = parse_init("polarization", "initial_polarization", "Px_0", "Py_0");
    mechanics = parse_init("mechanics", "initial_mechanics", "ux_0", "uy_0");
    electrostatics = parse_init("electrostatics", "initial_electrostatics", "Ex_0", "Ey_0");

}

void Datafile::parse_fracture(const toml::table& config) {
    // --- Mode de fissure (perméable / imperméable / conductive) ---
    std::string frac_mode = get_val<std::string>(config, "fracture", "mode", "PERMEABLE");
    if (frac_mode == "IMPERMEABLE") fracture.mode = CrackBCType::IMPERMEABLE;
    else if (frac_mode == "CONDUCTIVE") fracture.mode = CrackBCType::CONDUCTIVE;
    else fracture.mode = CrackBCType::PERMEABLE;

    // --- Pré-fissure : activation ---
    fracture.enable_precrack = get_val<bool>(config, "fracture", "enable_precrack", false);
    if (auto precrack_tbl = config["fracture"]["precrack"].as_table()) {
        fracture.enable_precrack = (*precrack_tbl)["enable"].value_or(bool(fracture.enable_precrack));
    }

    if (!fracture.enable_precrack) {
        return;
    }

    auto precrack_tbl = config["fracture"]["precrack"].as_table();
    if (!precrack_tbl) {
        Logger::error("[Datafile] enable_precrack=true mais section [fracture.precrack] introuvable. Pre-fissure desactivee.");
        fracture.enable_precrack = false;
        return;
    }

    // --- Forme de la pré-fissure ---
    std::string shape_str = (*precrack_tbl)["shape"].value_or<std::string>("segment");
    if (shape_str == "segment") fracture.precrack.shape = PrecrackShape::SEGMENT;
    else if (shape_str == "rect") fracture.precrack.shape = PrecrackShape::RECT;
    else {
        Logger::error("[Datafile] forme de precrack inconnue: '" + shape_str + "'. Utilisation de NONE.");
        fracture.precrack.shape = PrecrackShape::NONE;
    }

    // --- Paramètres géométriques SEGMENT ---
    fracture.precrack.x0         = (*precrack_tbl)["x0"].value_or<double>(0.0);
    fracture.precrack.y0         = (*precrack_tbl)["y0"].value_or<double>(mesh.Ly / 2.0);
    fracture.precrack.length     = (*precrack_tbl)["length"].value_or<double>(5.0);
    fracture.precrack.half_width = (*precrack_tbl)["half_width"].value_or<double>(0.5);

    // --- Paramètres géométriques RECT ---
    fracture.precrack.xmin = (*precrack_tbl)["xmin"].value_or<double>(0.0);
    fracture.precrack.xmax = (*precrack_tbl)["xmax"].value_or<double>(5.0);
    fracture.precrack.ymin = (*precrack_tbl)["ymin"].value_or<double>(mesh.Ly / 2.0 - 2.0);
    fracture.precrack.ymax = (*precrack_tbl)["ymax"].value_or<double>(mesh.Ly / 2.0 + 2.0);

    // --- Profil du champ v ---
    fracture.precrack.smooth           = (*precrack_tbl)["smooth"].value_or<bool>(true);
    fracture.precrack.smoothing_length = (*precrack_tbl)["smoothing_length"].value_or(double(material.kappa));

    // --- Croissance temporelle prescrite de la longueur (test de controle, v impose non resolu) ---
    // Si absents, growth_length_start/end valent length par defaut (=> pas de croissance, comportement statique).
    fracture.precrack.growth_enable = (*precrack_tbl)["growth_enable"].value_or<bool>(false);

    double default_length = fracture.precrack.length;
    double default_t_end  = simulation.total_time;

    fracture.precrack.growth_enable = (*precrack_tbl)["growth_enable"].value_or<bool>(false);

    fracture.precrack.growth_length_start = (*precrack_tbl)["growth_length_start"].value_or<double>(fracture.precrack.length + 0.0);
    fracture.precrack.growth_length_end   = (*precrack_tbl)["growth_length_end"].value_or<double>(fracture.precrack.length + 0.0);
    fracture.precrack.growth_t_end        = (*precrack_tbl)["growth_t_end"].value_or<double>(simulation.total_time + 0.0);

    if (fracture.precrack.growth_t_end <= 0.0) {
        Logger::error("[Datafile] precrack.growth_t_end doit etre strictement positif. Valeur recue: " +
                       std::to_string(fracture.precrack.growth_t_end) + ". Utilisation de simulation.total_time.");
        fracture.precrack.growth_t_end = simulation.total_time;
    }

    // --- Garde-fou : la demi-largeur ne doit pas être plus petite que la taille de maille ---
    double h_approx = std::min(mesh.dx, mesh.dy);
    if (fracture.precrack.half_width < h_approx) {
        Logger::error("[Datafile] precrack.half_width (" + std::to_string(fracture.precrack.half_width) +
                       ") plus petit que la taille de maille (" + std::to_string(h_approx) +
                       "). La pre-fissure risque de ne pas etre resolue par la discretisation.");
    }
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
                r.t_start  = (*tbl)["t_start"].value_or<double>(0.0);
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