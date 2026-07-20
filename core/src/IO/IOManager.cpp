#include "IO/IOManager.h"
#include "Physics/Mechanics.h"
#include "Physics/Electrostatics.h"
#include "Physics/Fracture.h"
#include "Physics/Polarization.h"
#include "Physics/MaterialModel.h"

IOManager::IOManager(ResultsExporter& exp, Diagnostics& diag, const std::string& csv_path, const Datafile& conf)
    : exporter(exp), diagnostics(diag), energy_csv_path(csv_path), config(conf) 
{}

void IOManager::extract_and_save_results(double time, int step, const std::string& initial_time_str, 
                                         Polarization& polarization, Mechanics& mechanics, 
                                         Fracture& fracture, Electrostatics& electrostatics, 
                                         const MaterialModel& material) {
                                         
    // 1. Diagnostics et énergies
    diagnostics.record(time, time, polarization, mechanics, fracture, electrostatics, material);
    diagnostics.append_csv(energy_csv_path);
    diagnostics.compute_nodal_energies(polarization, mechanics, fracture, electrostatics, material);

    // 2. Projection nodale (Mécanique)
    if (config.physics_toggle.mechanics) {
        mechanics.compute_stress_field(polarization, fracture, material);
    }

    // 3. Export VTK
    if (step % config.simulation.save_frequency == 0) {
        std::string filename = config.simulation.output_dir + "/VTK_" + initial_time_str + 
                               "/multiphysics_results" + std::to_string(step / config.simulation.save_frequency) + ".vtk";            
        
        // Construction explicite des vecteurs pour satisfaire le typage
        std::vector<std::string> scalar_names = {"v", "phi", "Energy_Gradient", "Energy_Elastic", "Energy_Landau", "Energy_Electric", "Energy_Surface", "sigma_xx", "sigma_yy", "sigma_xy"};
        std::vector<const Eigen::VectorXd*> scalar_fields = {
            &fracture.get_v(), &electrostatics.get_phi(),
            &diagnostics.get_U_nodal(), &diagnostics.get_W_nodal(), 
            &diagnostics.get_chi_nodal(), &diagnostics.get_elec_nodal(), 
            &diagnostics.get_surf_nodal(),
            &mechanics.get_sigma_xx(), &mechanics.get_sigma_yy(), &mechanics.get_sigma_xy()
        };
        
        std::vector<std::string> vector_names = {"P", "U", "E"};
        std::vector<const Eigen::VectorXd*> vector_x = {&polarization.get_Px(), &mechanics.get_ux(), &electrostatics.get_Ex()};
        std::vector<const Eigen::VectorXd*> vector_y = {&polarization.get_Py(), &mechanics.get_uy(), &electrostatics.get_Ey()};

        exporter.exportMultiPhysicsVTK(
            filename, 
            scalar_names, 
            scalar_fields, 
            vector_names, 
            vector_x, 
            vector_y
        );
    }
}