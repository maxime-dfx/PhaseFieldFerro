#include "IO/include/IOManager.h"
#include "Physics/include/Core/PhysicsManager.h"
#include "Materials/Core/MaterialManager.h"
#include "Materials/Microstructure/Polycrystal.h"
#include "Utils/include/Profiling.h"
#include <iomanip>
#include <sstream>

IOManager::IOManager(ResultsExporter& exp, const std::string& run_output_dir_in, const Datafile& conf)
    : exporter(exp), diagnostics(), run_output_dir(run_output_dir_in), config(conf) {}

void IOManager::finalize() const {
    diagnostics.write_csv(run_output_dir + "/energies_final.csv");
} 

void IOManager::extract_and_save_results(double t, int step, const PhysicsManager& pm, const Polycrystal& poly) {
    PROFILE_ZONE_NC("IOManager::extract_and_save_results", PROFILE_COLOR_SEQUENTIAL);

    // 1. Diagnostics (Délègue le calcul au PhysicsManager)
    {
        PROFILE_ZONE_NC("Diagnostics_Record_And_Export", PROFILE_COLOR_SEQUENTIAL);
        SystemMetrics metrics = pm.compute_system_metrics();
        diagnostics.record(step, t, metrics);
        diagnostics.append_csv(run_output_dir + "/energies.csv");
    }

    // 2. Export VTK
    if (step % config.simulation.save_frequency == 0) {
        PROFILE_ZONE_NC("VTK_Export", PROFILE_COLOR_SEQUENTIAL);
        std::ostringstream idx;
        idx << std::setw(5) << std::setfill('0') << (step / config.simulation.save_frequency);
        std::string filename = run_output_dir + "/VTK/multiphysics_results" + idx.str() + ".vtk";

        // --- NOUVELLE ARCHITECTURE : Extraction via le Blackboard (PhysicsState) ---
        const PhysicsState& state = pm.get_current_state();

        // Déduction de la taille nodale via les champs actifs pour sécuriser les exports
        int n_nodes = 0;
        if (state.get_field("v")) n_nodes = state.get_field("v")->size();
        else if (state.get_field("phi")) n_nodes = state.get_field("phi")->size();
        else if (state.get_field("ux")) n_nodes = state.get_field("ux")->size();
        if (n_nodes == 0) n_nodes = 1000; // Fallback minimal de sécurité
        
        static Eigen::VectorXd dummy_zeros;
        if (dummy_zeros.size() != n_nodes) {
            dummy_zeros = Eigen::VectorXd::Zero(n_nodes);
        }

        // Fonction de récupération sécurisée (renvoie des zéros si la physique est inactive)
        auto get_field = [&](const std::string& name) -> const Eigen::VectorXd* {
            const Eigen::VectorXd* ptr = state.get_field(name);
            return (ptr && ptr->size() > 0) ? ptr : &dummy_zeros;
        };

        // --- Extraction des champs scalaires ---
        const Eigen::VectorXd* v_field          = get_field("v");
        const Eigen::VectorXd* phi_field        = get_field("phi");
        const Eigen::VectorXd* sigma_xx_field   = get_field("sigma_xx");
        const Eigen::VectorXd* sigma_yy_field   = get_field("sigma_yy");
        const Eigen::VectorXd* sigma_xy_field   = get_field("sigma_xy");
        const Eigen::VectorXd* von_mises_field  = get_field("von_mises");
        const Eigen::VectorXd* energy_grad_field    = get_field("energy_gradient");
        const Eigen::VectorXd* energy_elastic_field = get_field("energy_elastic");
        const Eigen::VectorXd* energy_landau_field  = get_field("energy_landau");
        const Eigen::VectorXd* energy_electric_field = get_field("energy_electric");
        const Eigen::VectorXd* energy_surface_field  = get_field("energy_surface");
        const Eigen::VectorXd* energy_bulk_field     = get_field("energy_bulk");

        std::vector<std::string> scalar_names = {
            "v", "phi", "sigma_xx", "sigma_yy", "sigma_xy", "von_mises",
            "energy_gradient", "energy_elastic", "energy_landau", "energy_electric",
            "energy_surface", "energy_bulk"
        };
        std::vector<const Eigen::VectorXd*> scalar_fields = {
            v_field, phi_field, sigma_xx_field, sigma_yy_field, sigma_xy_field, von_mises_field,
            energy_grad_field, energy_elastic_field, energy_landau_field, energy_electric_field,
            energy_surface_field, energy_bulk_field
        };

        // --- Extraction des champs vectoriels ---
        // Note : On extrait "Px" et "Py". S'ils nécessitent un moyennage nodal (Px_nodal),
        // il faudra l'implémenter dans un PolarizationPostProcessor qui expose "Px_nodal".
        const Eigen::VectorXd* Px_field = get_field("Px");
        const Eigen::VectorXd* Py_field = get_field("Py");
        const Eigen::VectorXd* ux_field = get_field("ux");
        const Eigen::VectorXd* uy_field = get_field("uy");
        const Eigen::VectorXd* Ex_field = get_field("Ex");
        const Eigen::VectorXd* Ey_field = get_field("Ey");
        const Eigen::VectorXd* Dx_field = get_field("Dx");
        const Eigen::VectorXd* Dy_field = get_field("Dy");

        std::vector<std::string> vector_names = {"P", "U", "E", "D"};
        std::vector<const Eigen::VectorXd*> vector_x = {Px_field, ux_field, Ex_field, Dx_field};
        std::vector<const Eigen::VectorXd*> vector_y = {Py_field, uy_field, Ey_field, Dy_field};

        // --- Extraction du maillage (Grains) ---
        if (!grain_id_field_ready && poly.num_grains() > 1) {
            int n_elem = poly.num_elements();
            grain_id_field.setZero(n_elem);
            for (int e = 0; e < n_elem; ++e) {
                grain_id_field(e) = static_cast<double>(poly.grain_id_for_element(e));
            }
            grain_id_field_ready = true;
        }

        if (grain_id_field_ready) {
            scalar_names.push_back("grain_id");
            scalar_fields.push_back(&grain_id_field);
        }

        // --- Extraction du champ de materiau (pour distinguer visuellement
        //     le joint de grain geometrique, si affecte, des grains normaux) ---
        {
            const Mesh& mesh = pm.get_mesh();
            const int n_elem = mesh.get_num_elements();
            const auto& elements = mesh.get_elements();

            material_id_field.setZero(n_elem);
            for (int e = 0; e < n_elem; ++e) {
                if (config.crystal.grain_boundary_material_id != -1 && poly.is_grain_boundary_element(e)) {
                    material_id_field(e) = static_cast<double>(config.crystal.grain_boundary_material_id);
                } else {
                    material_id_field(e) = static_cast<double>(elements[e].ref_tag);
                }
            }
            material_id_field_ready = true;
        }

        if (material_id_field_ready) {
            scalar_names.push_back("material_id");
            scalar_fields.push_back(&material_id_field);
        }

        // --- Extraction des bandes decouplees (debug ParaView) : bande de
        //     joint de grain "sentie" par le champ de fracture (Gc/mu_v/kappa
        //     affaiblis) et bande de verrouillage de polarisation (P=0),
        //     potentiellement de largeurs differentes (grain_boundary_band_rings
        //     vs polarization_lock_band_rings) et donc invisibles l'une par
        //     rapport a l'autre si l'on ne regarde que "material_id".
        {
            const int n_elem = pm.get_mesh().get_num_elements();

            frac_band_field.setZero(n_elem);
            polar_band_field.setZero(n_elem);
            for (int e = 0; e < n_elem; ++e) {
                frac_band_field(e) = poly.is_fracture_band_element(e) ? 1.0 : 0.0;
                polar_band_field(e) = poly.is_polarization_lock_element(e) ? 1.0 : 0.0;
            }
            frac_band_ready = true;
            polar_band_ready = true;
        }

        if (frac_band_ready) {
            scalar_names.push_back("fracture_band");
            scalar_fields.push_back(&frac_band_field);
        }
        if (polar_band_ready) {
            scalar_names.push_back("polarization_band");
            scalar_fields.push_back(&polar_band_field);
        }

        // Export VTK global
        exporter.exportMultiPhysicsVTK(filename, scalar_names, scalar_fields, vector_names, vector_x, vector_y);
    }
}