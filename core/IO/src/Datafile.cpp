#include "IO/include/Datafile.h"
#include "Utils/include/Logger.h"
#include <stdexcept>
#include <algorithm>

template <typename T>
T Datafile::get_val(const toml::table& config, const std::string& section, const std::string& key, T default_value) {
    if (auto node = config[section][key]) {
        return node.value<T>().value_or(default_value);
    }
    return default_value;
} 

template <typename T>
T Datafile::get_val_from_table(const toml::table& tbl, const std::string& key, T default_value) {
    if (auto node = tbl[key]) {
        return node.value<T>().value_or(default_value);
    }
    return default_value;
}

Datafile::Datafile(const std::string& filename) {
    try {
        toml::table raw_config = toml::parse_file(filename);
        
        parse_simulation(raw_config);
        parse_mesh(raw_config);
        parse_chrono(raw_config);
        parse_physics(raw_config);
        parse_material(raw_config);
        parse_materials_list(raw_config);
        parse_crystal(raw_config);
        parse_fracture(raw_config);
        parse_boundaries(raw_config);
    } catch (const toml::parse_error& err) {
        throw std::runtime_error("Erreur de parsing TOML (" + filename + ") : " + std::string(err.description()));
    }
}

void Datafile::parse_simulation(const toml::table& config) {
    simulation.output_dir = get_val<std::string>(config, "output", "OutputDir", "../results");
    simulation.mesh_create_file = get_val<std::string>(config, "output", "MeshCreateFile", "mesh");
    simulation.debug_enabled = get_val<bool>(config, "debug", "enabled", false);
    
    simulation.total_time = get_val<double>(config, "simulation", "total_time", 3.0);
    simulation.dt = get_val<double>(config, "simulation", "dt", 3e-2);
    simulation.dt_relax = get_val<double>(config, "simulation", "dt_relax", 0.1);
    simulation.save_frequency = get_val<int>(config, "simulation", "save_frequency", 10);
    simulation.tol_ferro = get_val<double>(config, "simulation", "tolerance_ferro", 1e-3);
    simulation.tol_vfield = get_val<double>(config, "simulation", "tolerance_vfield", 1e-3);
    simulation.min_iter = get_val<int>(config, "simulation", "min_iter", 20);
    simulation.max_iter = get_val<int>(config, "simulation", "max_iter", 50);
    simulation.omega_polarization = get_val<double>(config, "simulation", "omega_polarization", 1.0);
}

void Datafile::parse_mesh(const toml::table& config) {
    mesh.calcul_mesh = get_val<bool>(config, "mesh", "calcul_mesh", true);
    mesh.get_mesh_file = get_val<std::string>(config, "mesh", "get_mesh_file", "../input/rectangle.msh");
    mesh.Lx = get_val<double>(config, "mesh", "L_x", 200.0);
    mesh.Ly = get_val<double>(config, "mesh", "L_y", 200.0);
    mesh.nx = get_val<int>(config, "mesh", "n_x", 200);
    mesh.ny = get_val<int>(config, "mesh", "n_y", 200);
    mesh.dx = mesh.Lx / mesh.nx;
    mesh.dy = mesh.Ly / mesh.ny;
    std::string type = get_val<std::string>(config, "mesh", "element_type", "TRIANGLE");
    mesh.element_type = (type == "QUAD") ? ElementType::QUAD4 : ElementType::TRIANGLE3;
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

    auto parse_init = [&](const std::string& section, const std::string& key_type, const std::string& key_x, const std::string& key_y, const std::string& key_norm = "") {
        PhysicsInitConfig c;
        std::string t = get_val<std::string>(config, section, key_type, "UNIFORM");
        c.type = (t == "RANDOM") ? InitializationType::RANDOM : InitializationType::UNIFORM;
        c.val_x_0 = get_val<double>(config, section, key_x, 0.0);
        c.val_y_0 = get_val<double>(config, section, key_y, 0.0);
        
        if (!key_norm.empty()) {
            c.norm_P0 = get_val<double>(config, section, key_norm, 0.0);
        } else {
            c.norm_P0 = 0.0;
        }
        
        std::string def_solver = (section == "mechanics") ? "DIRECT" : "CG";
        c.solver_type         = get_val<std::string>(config, section, "solver_type", def_solver);
        c.preconditioner_type = get_val<std::string>(config, section, "preconditioner_type", "DIAGONAL"); 
        c.solver_tol          = get_val<double>(config, section, "solver_tol", 1e-6);
        c.solver_max_iter     = get_val<int>(config, section, "solver_max_iter", (section == "mechanics" ? 2000 : 1000));
        return c;
    };

    polarization = parse_init("polarization", "initial_polarization", "Px_0", "Py_0", "norm_P0");
    mechanics = parse_init("mechanics", "initial_mechanics", "ux_0", "uy_0");
    electrostatics = parse_init("electrostatics", "initial_electrostatics", "Ex_0", "Ey_0");
}

void Datafile::parse_fracture(const toml::table& config) {
    fracture.solver_type     = get_val<std::string>(config, "fracture", "solver_type", "CG");
    fracture.solver_tol      = get_val<double>(config, "fracture", "solver_tol", 1e-6);
    fracture.solver_max_iter = get_val<int>(config, "fracture", "solver_max_iter", 1000);
    fracture.alpha_irreversibility = get_val<double>(config, "fracture", "alpha_irreversibility", 2e-2);

    std::string frac_mode = get_val<std::string>(config, "fracture", "mode", "PERMEABLE");
    if (frac_mode == "IMPERMEABLE") fracture.mode = CrackBCType::IMPERMEABLE;
    else if (frac_mode == "CONDUCTIVE") fracture.mode = CrackBCType::CONDUCTIVE;
    else fracture.mode = CrackBCType::PERMEABLE;

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

    std::string shape_str = (*precrack_tbl)["shape"].value_or<std::string>("segment");
    if (shape_str == "segment") fracture.precrack.shape = PrecrackShape::SEGMENT;
    else if (shape_str == "rect") fracture.precrack.shape = PrecrackShape::RECT;
    else {
        Logger::error("[Datafile] forme de precrack inconnue: '", shape_str, "'. Utilisation de NONE.");
        fracture.precrack.shape = PrecrackShape::NONE;
    }

    fracture.precrack.x0         = (*precrack_tbl)["x0"].value_or<double>(0.0);
    fracture.precrack.y0         = (*precrack_tbl)["y0"].value_or<double>(mesh.Ly / 2.0);
    fracture.precrack.length     = (*precrack_tbl)["length"].value_or<double>(5.0);
    fracture.precrack.half_width = (*precrack_tbl)["half_width"].value_or<double>(0.5);

    fracture.precrack.xmin = (*precrack_tbl)["xmin"].value_or<double>(0.0);
    fracture.precrack.xmax = (*precrack_tbl)["xmax"].value_or<double>(5.0);
    fracture.precrack.ymin = (*precrack_tbl)["ymin"].value_or<double>(mesh.Ly / 2.0 - 2.0);
    fracture.precrack.ymax = (*precrack_tbl)["ymax"].value_or<double>(mesh.Ly / 2.0 + 2.0);

    fracture.precrack.smooth           = (*precrack_tbl)["smooth"].value_or<bool>(true);
    fracture.precrack.smoothing_length = (*precrack_tbl)["smoothing_length"].value_or(double(material.kappa));

    fracture.precrack.growth_enable = (*precrack_tbl)["growth_enable"].value_or<bool>(false);
    fracture.precrack.growth_length_start = (*precrack_tbl)["growth_length_start"].value_or<double>(fracture.precrack.length + 0.0);
    fracture.precrack.growth_length_end   = (*precrack_tbl)["growth_length_end"].value_or<double>(fracture.precrack.length + 0.0);
    fracture.precrack.growth_t_end        = (*precrack_tbl)["growth_t_end"].value_or<double>(simulation.total_time + 0.0);

    if (fracture.precrack.growth_t_end <= 0.0) {
        Logger::error("[Datafile] precrack.growth_t_end doit etre strictement positif. Valeur recue: ",
                       fracture.precrack.growth_t_end, ". Utilisation de simulation.total_time.");
        fracture.precrack.growth_t_end = simulation.total_time;
    }

    double h_approx = std::min(mesh.dx, mesh.dy);
    if (fracture.precrack.half_width < h_approx) {
        Logger::error("[Datafile] precrack.half_width (", fracture.precrack.half_width,
                       ") plus petit que la taille de maille (", h_approx,
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
    
    material.alpha_1 = get_val<double>(config, "material", "alpha_1", 1.0);
    material.alpha_11 = get_val<double>(config, "material", "alpha_11", 1.0);
    material.alpha_111 = get_val<double>(config, "material", "alpha_111", 1.0);
    material.alpha_1111 = get_val<double>(config, "material", "alpha_1111", 1.0);
    material.alpha_12 = get_val<double>(config, "material", "alpha_12", 1.0);
    material.alpha_112 = get_val<double>(config, "material", "alpha_112", 1.0);
    material.alpha_1112 = get_val<double>(config, "material", "alpha_1112", 1.0);
    material.alpha_1122 = get_val<double>(config, "material", "alpha_1122", 1.0);
    
    material.eta_k = get_val<double>(config, "material", "eta_k", 1.0);
    material.Gc = get_val<double>(config, "material", "Gc", 1.0);
    material.kappa = get_val<double>(config, "material", "kappa", 1.0);

    // [material] historique -> id 0, actif (Ferroelectric) par défaut.
    material.id = 0;
    material.is_ferroelectric = true;
    material.type = MaterialType::Ferroelectric;
}

// Piézocomposites : bloc TOML optionnel [[materials]] (array-of-tables), un
// par Physical Tag Gmsh, ex. :
//
//   [[materials]]
//   id = 1
//   is_ferroelectric = false   # matrice polymère isolante
//   eps0 = 3.5e-11
//   c1 = 8.0e9
//   c2 = 3.0e9
//   c3 = 2.5e9
//
//   [[materials]]
//   id = 2
//   is_ferroelectric = true    # inclusions PZT/BaTiO3
//   a0 = ...
//   alpha_1 = ...
//   ...
//
// Chaque entrée démarre des valeurs par défaut de [material] (déjà parsé
// ci-dessus) puis les surcharge avec les clés présentes dans son propre
// sous-tableau, ce qui permet de ne spécifier que ce qui diffère de la
// matrice/du matériau de base. Si [[materials]] est absent, on retombe sur
// un unique matériau {material, id=0} pour rester compatible avec les .toml
// mono-matériau existants.
void Datafile::parse_materials_list(const toml::table& config) {
    materials.clear();

    const toml::array* arr = config["materials"].as_array();
    if (arr == nullptr || arr->empty()) {
        materials.push_back(material);
        return;
    }

    for (const auto& node : *arr) {
        const toml::table* tbl = node.as_table();
        if (tbl == nullptr) {
            Logger::error("[Datafile] Entrée [[materials]] invalide (pas une table), ignorée.");
            continue;
        }

        MaterialConfig m = this->material; // instance propre : chaque [[materials]] est 100% autonome (plus d'héritage implicite de [material])

        m.id = get_val_from_table<int>(*tbl, "id", m.id);
        // On lit d'abord la valeur explicite (ou heritee de [material]) AVANT
        // de resoudre "type", car le defaut de "type" en depend ci-dessous.
        const bool is_ferro_explicit = get_val_from_table<bool>(*tbl, "is_ferroelectric", m.is_ferroelectric);
        m.is_ferroelectric = is_ferro_explicit;

        // "type" prend le pas sur "is_ferroelectric" s'il est présent (nouveau
        // format explicite : Ferroelectric / PureElastic / PolymerElastic /
        // GrainBoundary). Sinon, dérivé du booléen legacy pour compatibilité
        // ascendante avec les .toml existants.
        const std::string default_type_str = m.is_ferroelectric ? "Ferroelectric" : "PolymerElastic";
        const std::string type_str = get_val_from_table<std::string>(*tbl, "type", default_type_str);
        try {
            m.type = material_type_from_string(type_str);
        } catch (const std::invalid_argument& err) {
            Logger::error("[Datafile] [[materials]] id=" + std::to_string(m.id) + " : " + err.what());
            throw;
        }
        // Garde is_ferroelectric cohérent avec type explicite pour
        // Ferroelectric/PureElastic/PolymerElastic (un seul choix possible,
        // is_ferroelectric == true ssi type == Ferroelectric). GrainBoundary
        // est un CAS A PART : c'est le seul type dont le comportement
        // ferroelectrique (P dynamique, energie de Landau, couplage
        // electromecanique) est piloté explicitement par l'utilisateur via
        // is_ferroelectric dans le TOML (cf. GrainBoundaryMaterial), et non
        // déduit du type - on ne l'écrase donc PAS ici, contrairement aux
        // autres types.
        if (m.type != MaterialType::GrainBoundary) {
            m.is_ferroelectric = (m.type == MaterialType::Ferroelectric);
        }

        m.a0 = get_val_from_table<double>(*tbl, "a0", m.a0);
        m.b1 = get_val_from_table<double>(*tbl, "b1", m.b1);
        m.b2 = get_val_from_table<double>(*tbl, "b2", m.b2);
        m.b3 = get_val_from_table<double>(*tbl, "b3", m.b3);
        m.c1 = get_val_from_table<double>(*tbl, "c1", m.c1);
        m.c2 = get_val_from_table<double>(*tbl, "c2", m.c2);
        m.c3 = get_val_from_table<double>(*tbl, "c3", m.c3);
        m.eps0 = get_val_from_table<double>(*tbl, "eps0", m.eps0);
        m.t = get_val_from_table<double>(*tbl, "t", m.t);
        m.P0 = get_val_from_table<double>(*tbl, "P0", m.P0);
        m.mu_p = get_val_from_table<double>(*tbl, "mu_p", m.mu_p);
        m.mu_v = get_val_from_table<double>(*tbl, "mu_v", m.mu_v);
        m.xi = get_val_from_table<double>(*tbl, "xi", m.xi);
        m.c0 = get_val_from_table<double>(*tbl, "c0", m.c0);

        m.alpha_1 = get_val_from_table<double>(*tbl, "alpha_1", m.alpha_1);
        m.alpha_11 = get_val_from_table<double>(*tbl, "alpha_11", m.alpha_11);
        m.alpha_111 = get_val_from_table<double>(*tbl, "alpha_111", m.alpha_111);
        m.alpha_1111 = get_val_from_table<double>(*tbl, "alpha_1111", m.alpha_1111);
        m.alpha_12 = get_val_from_table<double>(*tbl, "alpha_12", m.alpha_12);
        m.alpha_112 = get_val_from_table<double>(*tbl, "alpha_112", m.alpha_112);
        m.alpha_1112 = get_val_from_table<double>(*tbl, "alpha_1112", m.alpha_1112);
        m.alpha_1122 = get_val_from_table<double>(*tbl, "alpha_1122", m.alpha_1122);

        m.eta_k = get_val_from_table<double>(*tbl, "eta_k", m.eta_k);
        m.Gc = get_val_from_table<double>(*tbl, "Gc", m.Gc);
        m.kappa = get_val_from_table<double>(*tbl, "kappa", m.kappa);

        materials.push_back(m);
    }

    if (materials.empty()) {
        Logger::error("[Datafile] [[materials]] présent mais vide/invalide, retour au matériau [material] unique (id=0).");
        materials.push_back(material);
    }
}

void Datafile::parse_crystal(const toml::table& config) {
    crystal.num_grains = get_val<int>(config, "crystal", "num_grains", 1);
    crystal.seed = static_cast<unsigned int>(get_val<int>(config, "crystal", "seed", 42));
    crystal.use_random_seed = get_val<bool>(config, "crystal", "use_random_seed", false);
    // id (MaterialConfig::id) d'une entrée [[materials]] de type GrainBoundary,
    // à appliquer aux éléments situés à un joint de grain. -1 = désactivé.
    crystal.grain_boundary_material_id = get_val<int>(config, "crystal", "grain_boundary_material_id", -1);
    // Dilatation de la bande de joint de grain, cf. commentaire dans
    // ConfigTypes.h. Doit rester >= 1 (1 = comportement legacy).
    crystal.grain_boundary_band_rings = std::max(1, get_val<int>(config, "crystal", "grain_boundary_band_rings", 1));
    // Dilatation INDEPENDANTE de la bande de verrouillage P=0 (polarisation),
    // decouplee de grain_boundary_band_rings ci-dessus. Doit rester >= 1
    // (1 = comportement legacy, identique a la bande stricte).
    crystal.polarization_lock_band_rings = std::max(1, get_val<int>(config, "crystal", "polarization_lock_band_rings", 1));

    // Fail-fast : mieux vaut une erreur claire au chargement qu'un
    // std::runtime_error opaque de MaterialRegistry::get_material() au
    // premier pas de temps.
    if (crystal.grain_boundary_material_id != -1) {
        const auto it = std::find_if(materials.begin(), materials.end(), [&](const MaterialConfig& m) {
            return m.id == crystal.grain_boundary_material_id;
        });
        if (it == materials.end()) {
            throw std::runtime_error(
                "[Datafile] crystal.grain_boundary_material_id=" + std::to_string(crystal.grain_boundary_material_id) +
                " ne correspond à aucune entrée [[materials]].");
        }
    }
}

void Datafile::parse_boundaries(const toml::table& config) {
    boundary_rules.clear();
    if (auto arr = config["boundary_rule"].as_array()) {
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
                r.t_start   = (*tbl)["t_start"].value_or<double>(0.0);
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
    }
}