#pragma once

#include <map>
#include <memory>
#include <random>
#include <vector>

#include "Materials/Core/MaterialConcepts.h"
#include "Materials/Core/MaterialManager.h"
#include "Materials/Microstructure/Crystal.h"
#include "Materials/Microstructure/Voronoi.h"
#include "IO/include/ConfigTypes.h"
#include "Mesh/include/Mesh.h"

// Polycrystal
// -----------------------------------------------------------------------
// Assigne les Crystals au maillage : remplace PolycrystalBuilder, et
// absorbe également l'ancien OrientedMaterialCache. Les deux anciennes
// classes répondaient à la même question ("quel matériau/grain voit cet
// élément ?"), l'une pour la géométrie (grain_id, joint de grain), l'autre
// pour le matériau tourné qui en résultait : elles vivent maintenant
// ensemble dans une seule classe à deux étapes :
//
//  1) Construction (géométrie) : tessellation de Voronoi, assignation
//     element -> grain_id, détection des éléments de joint de grain.
//     Équivalent strict de l'ancien PolycrystalBuilder.
//
//  2) attach_materials() (physique) : une fois un MaterialManager
//     disponible, construit une instance CrystalMaterial par couple
//     (ref_tag, grain_id) réellement présent dans le maillage, puis expose
//     une table element -> MaterialModel oriente en O(1). Équivalent
//     strict de l'ancien OrientedMaterialCache, appelé en une passe
//     séquentielle (aucune synchronisation nécessaire lors des boucles
//     d'assemblage OpenMP, qui ne font que lire la table).
class Polycrystal final {
private:
    std::vector<Crystal> grains;
    std::vector<int> element_to_grain;
    // Un élément est "joint de grain" s'il possède au moins un nœud partagé
    // avec un élément voisin appartenant à un grain différent. Dimensionné à
    // num_elements, entièrement false si num_grains<=1 (monocristal : pas de
    // joint) ou si l'appelant n'a pas de matériau GrainBoundary configuré
    // (le flag reste calculé mais peut alors être ignoré, cf.
    // attach_materials()).
    std::vector<bool> element_is_gb_strict;
    // Version dilatee de element_is_gb_strict (BFS sur grain_boundary_band_rings
    // anneaux). Utilisee UNIQUEMENT pour affaiblir la tenacite Gc du champ de
    // fracture (cf. Fracture.h::compute_element_matrices), afin que le champ
    // diffus "sente" le contraste materiau sur une largeur de bande comparable
    // a la longueur de regularisation. La mecanique, elle, doit continuer a
    // utiliser la bande stricte via is_grain_boundary_element(), sans quoi
    // H_drive est verrouille sur une zone trop large et la fissure ne peut
    // plus se propager.
    std::vector<bool> element_is_gb_wide;
    // Version dilatee independante (BFS sur polarization_lock_band_rings
    // anneaux), utilisee UNIQUEMENT par PolarizationDofMapper pour verrouiller
    // P=0. Decouplee de element_is_gb_wide : la largeur de la zone
    // "depolarisee" physiquement (joint de grain "epais") n'a aucune raison
    // de coincider avec la largeur choisie pour que le champ de fracture
    // sente le contraste de Gc - les deux effets sont pilotes independamment
    // via deux parametres TOML distincts (cf. CrystalConfig).
    std::vector<bool> element_is_gb_polarization;
    unsigned int seed = 0;

    std::vector<std::unique_ptr<MaterialModel>> m_owned;
    std::vector<const MaterialModel*> m_material_for_element;

    // Un nœud est "ambigu" si les éléments incidents ne portent pas tous le
    // même grain_id ; un élément est alors marqué joint de grain s'il touche
    // au moins un nœud ambigu. Coût O(num_elements * num_nodes_per_elem),
    // exécuté une seule fois à la construction.
    void compute_grain_boundary_elements(const Mesh& mesh) {
        const int num_elements = static_cast<int>(element_to_grain.size());
        const int num_nodes = mesh.get_num_nodes();
        std::vector<int> node_grain(static_cast<std::size_t>(num_nodes), -1);
        std::vector<bool> node_is_ambiguous(static_cast<std::size_t>(num_nodes), false);

        for (int e = 0; e < num_elements; ++e) {
            const Element& elem = mesh.get_elements()[static_cast<std::size_t>(e)];
            const int g = element_to_grain[static_cast<std::size_t>(e)];
            for (int i = 0; i < elem.get_num_nodes(); ++i) {
                const int idx = mesh.get_node_index(e, i);
                if (node_grain[static_cast<std::size_t>(idx)] == -1) {
                    node_grain[static_cast<std::size_t>(idx)] = g;
                } else if (node_grain[static_cast<std::size_t>(idx)] != g) {
                    node_is_ambiguous[static_cast<std::size_t>(idx)] = true;
                }
            }
        }

        for (int e = 0; e < num_elements; ++e) {
            const Element& elem = mesh.get_elements()[static_cast<std::size_t>(e)];
            for (int i = 0; i < elem.get_num_nodes(); ++i) {
                const int idx = mesh.get_node_index(e, i);
                if (node_is_ambiguous[static_cast<std::size_t>(idx)]) {
                    element_is_gb_strict[static_cast<std::size_t>(e)] = true;
                    break;
                }
            }
        }
    }

    // Dilate le masque booleen `base` de (rings - 1) anneaux supplementaires
    // par propagation BFS via les noeuds partages entre elements, et renvoie
    // le resultat dans un nouveau vecteur (base n'est jamais modifie).
    // rings<=1 renvoie une copie inchangee de `base`.
    //
    // Fonction generique reutilisee pour construire independamment
    // element_is_gb_wide (bande fracture) et element_is_gb_polarization
    // (bande depolarisation), chacune avec son propre nombre d'anneaux.
    //
    // Utile car la largeur de bande necessaire pour que le champ de
    // fracture diffus "sente" reellement le contraste materiau (Gc, kappa,
    // elasticite) doit etre du meme ordre que la longueur de regularisation
    // kappa (=l0) du modele AT2, sans quoi le champ v moyenne le contraste
    // au lieu de le capturer nettement - et cette largeur n'a aucune raison
    // physique de coincider avec la largeur de la zone reellement depolarisee.
    static std::vector<bool> dilate_band(const Mesh& mesh, const std::vector<bool>& base, int rings) {
        std::vector<bool> result = base;
        if (rings <= 1) return result;

        const int num_elements = static_cast<int>(base.size());
        const int num_nodes = mesh.get_num_nodes();

        // node -> elements incidents (reconstruit a chaque appel : le cout
        // O(num_elements * num_nodes_par_elem) reste negligeable devant le
        // reste de la construction, et evite de mutualiser un etat partage
        // entre les deux bandes).
        std::vector<std::vector<int>> node_to_elements(static_cast<std::size_t>(num_nodes));
        for (int e = 0; e < num_elements; ++e) {
            const Element& elem = mesh.get_elements()[static_cast<std::size_t>(e)];
            for (int i = 0; i < elem.get_num_nodes(); ++i) {
                const int idx = mesh.get_node_index(e, i);
                node_to_elements[static_cast<std::size_t>(idx)].push_back(e);
            }
        }

        for (int ring = 1; ring < rings; ++ring) {
            std::vector<bool> next = result;
            for (int e = 0; e < num_elements; ++e) {
                if (!result[static_cast<std::size_t>(e)]) continue;
                const Element& elem = mesh.get_elements()[static_cast<std::size_t>(e)];
                for (int i = 0; i < elem.get_num_nodes(); ++i) {
                    const int idx = mesh.get_node_index(e, i);
                    for (int neighbor : node_to_elements[static_cast<std::size_t>(idx)]) {
                        next[static_cast<std::size_t>(neighbor)] = true;
                    }
                }
            }
            result.swap(next);
        }
        return result;
    }

public:
    Polycrystal() = default;

    unsigned int get_seed() const { return seed; }
    
    explicit Polycrystal(const CrystalConfig& crystal_config, const Mesh& mesh) {
        seed = crystal_config.seed;
        if (crystal_config.use_random_seed) {
            std::random_device rd;
            seed = rd();
        }

        const int num_elements = mesh.get_num_elements();
        element_to_grain.assign(static_cast<std::size_t>(num_elements), 0);
        element_is_gb_strict.assign(static_cast<std::size_t>(num_elements), false);
        element_is_gb_wide.assign(static_cast<std::size_t>(num_elements), false);
        element_is_gb_polarization.assign(static_cast<std::size_t>(num_elements), false);

        if (crystal_config.num_grains <= 1) {
            return;
        }

        VoronoiTessellation tessellation(mesh.get_Lx(), mesh.get_Ly(), crystal_config.num_grains, seed);
        grains = tessellation.grains();

        for (int e = 0; e < num_elements; ++e) {
            const Element& elem = mesh.get_elements()[static_cast<std::size_t>(e)];
            const auto coords = mesh.get_element_coords(e);
            double cx = 0.0;
            double cy = 0.0;
            for (int n = 0; n < elem.get_num_nodes(); ++n) {
                cx += coords[static_cast<std::size_t>(n)][0];
                cy += coords[static_cast<std::size_t>(n)][1];
            }
            cx /= elem.get_num_nodes();
            cy /= elem.get_num_nodes();
            element_to_grain[static_cast<std::size_t>(e)] = tessellation.nearest_grain(cx, cy);
        }

        compute_grain_boundary_elements(mesh);
        element_is_gb_wide = dilate_band(mesh, element_is_gb_strict, crystal_config.grain_boundary_band_rings);
        element_is_gb_polarization = dilate_band(mesh, element_is_gb_strict, crystal_config.polarization_lock_band_rings);
    }

    // Bande stricte (1 element de large) : utilisee par la mecanique et la
    // polarisation pour verrouiller P=0 sur le joint de grain lui-meme.
    bool is_grain_boundary_element(int elem_idx) const {
        return element_is_gb_strict.empty() ? false : element_is_gb_strict.at(static_cast<std::size_t>(elem_idx));
    }

    // Bande dilatee (grain_boundary_band_rings anneaux) : utilisee
    // uniquement par le champ de fracture pour affaiblir Gc sur une largeur
    // comparable a la longueur de regularisation, sans toucher a la
    // mecanique/polarisation.
    bool is_fracture_band_element(int elem_idx) const {
        return element_is_gb_wide.empty() ? false : element_is_gb_wide.at(static_cast<std::size_t>(elem_idx));
    }

    // Bande dilatee independante (polarization_lock_band_rings anneaux) :
    // utilisee uniquement par PolarizationDofMapper pour determiner la
    // largeur de la zone ou P=0 est verrouille (joint de grain "epais"),
    // sans aucun lien avec la largeur choisie pour l'affaiblissement de Gc.
    bool is_polarization_lock_element(int elem_idx) const {
        return element_is_gb_polarization.empty() ? false : element_is_gb_polarization.at(static_cast<std::size_t>(elem_idx));
    }

    int num_grains() const {
        return static_cast<int>(grains.size());
    }

    int num_elements() const {
        return static_cast<int>(element_to_grain.size());
    }

    int grain_id_for_element(int elem_idx) const {
        return element_to_grain.empty() ? 0 : element_to_grain.at(static_cast<std::size_t>(elem_idx));
    }

    double grain_angle(int grain_id) const {
        return grains.at(static_cast<std::size_t>(grain_id)).theta;
    }

    bool empty() const {
        return element_to_grain.empty();
    }

    // Construit, une seule fois et de maniere sequentielle, une instance
    // CrystalMaterial par couple (ref_tag, grain_id) reellement present
    // dans le maillage. A appeler une fois le MaterialManager peuplé
    // (registre des materiaux par ref_tag).
    //
    // crystal_config : uniquement consulté pour grain_boundary_material_id
    // (le reste, num_grains/seed, a déjà servi à construire la geometrie
    // ci-dessus).
    void attach_materials(const MaterialManager& manager, const Mesh& mesh,
                           const CrystalConfig& crystal_config = CrystalConfig{}) {
        // ATTENTION : empty() ne suffit pas ici. element_to_grain est
        // toujours dimensionne a num_elements (rempli de 0) meme quand
        // num_grains <= 1 (cf. constructeur ci-dessus), seul `grains` reste
        // vide dans ce cas. Se fier a empty() ferait passer le monocristal
        // comme "oriente" et grain_angle(0) crasherait sur `grains.at(0)`
        // vide.
        if (num_grains() <= 1) {
            return;
        }

        const int num_elements = mesh.get_num_elements();
        m_material_for_element.assign(static_cast<std::size_t>(num_elements), nullptr);

        std::map<std::pair<int, int>, MaterialModel*> cache_by_tag_and_grain;
        const int gb_id = crystal_config.grain_boundary_material_id;
        const MaterialModel* gb_material = (gb_id != -1) ? &manager.get_material(gb_id) : nullptr;

        for (int e = 0; e < num_elements; ++e) {
            // Joint de grain : matériau dédié (constantes indépendantes du
            // grain massif), SANS décoration d'orientation cristalline -
            // c'est une couche désordonnée, l'angle d'easy-axis n'a pas de
            // sens. Si aucun matériau GrainBoundary n'est configuré
            // (gb_material == nullptr), on retombe sur le comportement
            // legacy : l'élément reste traité comme un élément de volume
            // ordinaire de son grain.
            if (gb_material != nullptr && is_grain_boundary_element(e)) {
                m_material_for_element[static_cast<std::size_t>(e)] = gb_material;
                continue;
            }

            const int ref_tag = mesh.get_elements()[static_cast<std::size_t>(e)].ref_tag;
            const int grain_id = grain_id_for_element(e);
            const auto key = std::make_pair(ref_tag, grain_id);

            auto it = cache_by_tag_and_grain.find(key);
            if (it == cache_by_tag_and_grain.end()) {
                const MaterialModel& base = manager.get_material(ref_tag);
                const double theta = grain_angle(grain_id);
                m_owned.push_back(std::make_unique<CrystalMaterial>(base, theta));
                it = cache_by_tag_and_grain.emplace(key, m_owned.back().get()).first;
            }

            m_material_for_element[static_cast<std::size_t>(e)] = it->second;
        }
    }

    bool has_orientation() const {
        return !m_material_for_element.empty();
    }

    const MaterialModel& get_material(int elem_idx, int fallback_material_id, const MaterialManager& manager) const {
        if (!m_material_for_element.empty()) {
            const MaterialModel* oriented = m_material_for_element[static_cast<std::size_t>(elem_idx)];
            if (oriented != nullptr) {
                return *oriented;
            }
        }
        return manager.get_material(fallback_material_id);
    }
};
