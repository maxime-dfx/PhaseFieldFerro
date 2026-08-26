#pragma once

#include <string>
#include "Utils/include/Types.h"

struct BoundaryRuleConfig {
    std::string field, bc_type, shape, edge_name, profile;
    double cx = 0.0, cy = 0.0, radius = 0.0;
    double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;
    double px = 0.0, py = 0.0, val = 0.0;
    double val_start = 0.0, val_end = 0.0, t_start = 0.0, t_end = 1.0;
    double amplitude = 0.0, frequency = 0.0, offset = 0.0;
    double max_val = 0.0, y_center = 0.0, width = 0.0;
    double P0 = 0.0, x0 = 0.0, epsilon = 1e-3, k = 0.0, omega = 0.0;
};

struct SimulationConfig {
    std::string output_dir;
    std::string mesh_create_file;
    // CORRECTIF : total_time doit etre un double, pas un int.
    // C'est un temps pseudo-total potentiellement fractionnaire (ex : 100 pas
    // de charge * Delta t_n = 3e-2 => total_time = 3.0, cf. Section 3.1 et
    // l'axe des Fig. 6/14/15 du papier). Avec un int, un TOML contenant
    // "total_time = 3.0" pouvait echouer au parsing strict (le noeud TOML est
    // de type float, pas integer) et retomber silencieusement sur la valeur
    // par defaut, faisant tourner la simulation ~100x trop longtemps.
    double total_time;
    double dt;
    double dt_relax;
    int save_frequency;
    double tol_ferro;
    double tol_vfield;
    int min_iter;
    int max_iter;
    bool debug_enabled;
    double omega_polarization;
};

struct MeshConfig {
    bool calcul_mesh;
    std::string get_mesh_file;
    double Lx, Ly;
    int nx, ny;
    double dx, dy;
    ElementType element_type;
};

struct ChronoConfig {
    bool mesh, polarization, mechanics, electrostatics, fracture, run;
};

struct PhysicsToggle {
    bool polarization, electrostatics, mechanics, fracture;
};

struct PhysicsInitConfig {
    InitializationType type;
    double val_x_0; 
    double val_y_0; 
    double norm_P0; 

    // --- Configuration Solveur ---
    std::string solver_type;
    std::string preconditioner_type;
    double solver_tol;
    int solver_max_iter;
};

// Type de matériau : détermine QUELLE classe MaterialModel MaterialManager
// instancie pour une entrée [[materials]] donnée. Un fichier Models/*.h par
// valeur (cf. Materials/include/Models/).
//  - Ferroelectric   -> SingleCrystalMaterial   (céramique active PZT/BaTiO3, Landau+gradient+couplage)
//  - PureElastic     -> PureElasticMaterial     (solide élastique pur, ni P ni permittivité)
//  - PolymerElastic  -> PolymerElasticMaterial  (matrice polymère : élastique + diélectrique linéaire, P=0)
//  - GrainBoundary   -> GrainBoundaryMaterial   (joint de grain : élastique + diélectrique linéaire + Gc/kappa/mu_v propres, P=0)
enum class MaterialType {
    Ferroelectric,
    PureElastic,
    PolymerElastic,
    GrainBoundary
};

// Conversion depuis la chaîne TOML "type = ...". Lève std::invalid_argument
// si la valeur est inconnue (fail-fast plutôt qu'un fallback silencieux).
MaterialType material_type_from_string(const std::string& s);
std::string material_type_to_string(MaterialType type);

struct MaterialConfig {
    double a0 = 0.1, b1 = 1.4282, b2 = -0.185, b3 = 0.8066, c1 = 185.0, c2 = 111.0, c3 = 74.0,
        eps0 = 0.131, t = 300.0, P0 = 1.0, mu_p = 1.0, mu_v = 15.0;
    double xi = 1.0, c0 = 1.0, alpha_1 = -0.0023, alpha_11 = -0.0029, alpha_111 = 0.003, alpha_1111 = 0.001;
    double alpha_12 = -0.0011, alpha_112 = -0.00068, alpha_1112 = 0.0093, alpha_1122 = 1.24;
    double eta_k = 1e-6, Gc = 4.0, kappa = 2.0;

    // --- Support multi-matériaux (piézocomposites + joints de grains) ---
    // id : Physical Tag Gmsh auquel ce matériau est associé (cf. Element::ref_tag
    // rempli par GmshReader), SAUF pour un matériau de type GrainBoundary, qui
    // n'est jamais associé à un ref_tag géométrique (cf. CrystalConfig::
    // grain_boundary_material_id) : il est résolu par élément à l'intérieur
    // d'une phase polycristalline via la détection de joint de grain.
    // Vaut 0 par défaut (comportement legacy mono-matériau).
    int id = 0;

    // type : pilote MaterialManager::create(). Si absent du TOML, dérivé de
    // is_ferroelectric pour compatibilité ascendante (cf. Datafile::parse_*).
    MaterialType type = MaterialType::Ferroelectric;

    // is_ferroelectric : conservé pour compatibilité (ancien format TOML
    // booléen) et car MaterialModel::is_ferroelectric() reste la source de
    // vérité consultée par PolarizationDofMapper pour verrouiller P=0. Les
    // classes concrètes (Models/*.h) renvoient leur propre valeur fixe et
    // n'ont plus besoin de lire ce champ.
    bool is_ferroelectric = true;
};

// -------------------------------------------------------------
// ATTENTION A L'ORDRE : PrecrackConfig AVANT FractureConfig !
// -------------------------------------------------------------
struct PrecrackConfig {
    PrecrackShape shape = PrecrackShape::NONE;
    double x0 = 0.0;
    double y0 = 50.0;
    double length = 5.0;
    double half_width = 0.5;
    double xmin = 0.0, xmax = 5.0, ymin = 48.0, ymax = 52.0;
    bool smooth = true;
    double smoothing_length = 1.0;
    bool growth_enable = false;
    double growth_length_start = 5.0;
    double growth_length_end   = 5.0;
    double growth_t_end        = 3.0;
};

struct CrystalConfig {
    int num_grains = 1;
    unsigned int seed = 42;
    bool use_random_seed = false;

    // id (MaterialConfig::id) de l'entrée [[materials]] de type GrainBoundary
    // à utiliser pour les éléments situés à un joint de grain (cf.
    // Polycrystal::is_grain_boundary_element). -1 = désactivé : les
    // joints de grains ne reçoivent alors aucun traitement matériau dédié,
    // ils utilisent le matériau de volume orienté comme n'importe quel
    // élément du grain (comportement legacy).
    int grain_boundary_material_id = -1;

    // Nombre d'anneaux de dilatation appliqués à la bande d'éléments
    // detectés comme "joint de grain" (cf. Polycrystal::
    // compute_grain_boundary_elements). 1 = comportement legacy (bande
    // d'un seul élément de large, autour des noeuds ambigus). A choisir
    // de sorte que largeur_bande ≈ N * h (taille d'élément) couvre au
    // moins la longueur de régularisation kappa (=l0) du champ de
    // fracture, sans quoi le champ diffus v moyenne le contraste de Gc
    // au lieu de le "sentir" nettement.
    int grain_boundary_band_rings = 1;

    // Nombre d'anneaux de dilatation appliqués à la bande d'éléments où
    // P=0 est verrouillé (cf. Polycrystal::is_polarization_lock_element,
    // consommé par PolarizationDofMapper). Découplé de
    // grain_boundary_band_rings : cette largeur pilote une zone physique
    // différente (l'étendue réelle de la couche non-ferroélectrique /
    // désordonnée), sans lien avec la largeur choisie pour que le champ
    // de fracture sente le contraste de Gc. 1 = comportement legacy
    // (bande d'un seul élément de large, identique à la détection stricte).
    int polarization_lock_band_rings = 1;
};

struct FractureConfig {
    CrackBCType mode = CrackBCType::PERMEABLE;
    bool enable_precrack = false;
    PrecrackConfig precrack; // Maintenant le compilateur sait ce que c'est !

    // --- Configuration Solveur ---
    std::string solver_type;
    std::string preconditioner_type;
    double solver_tol;
    int solver_max_iter;


    double alpha_irreversibility = 2e-2;
};