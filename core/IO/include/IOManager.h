#pragma once
#include "IO/include/ResultsExporter.h"
#include "IO/include/internal/Diagnostics.h"
#include "IO/include/Datafile.h"
#include <Eigen/Dense>
#include <string>

class PhysicsManager;
class MaterialsManager;
class PolyCristal;

class IOManager {
private:
    ResultsExporter& exporter;
    Diagnostics diagnostics;
    std::string run_output_dir;
    const Datafile& config;

    Eigen::VectorXd grain_id_field;
    bool grain_id_field_ready = false;

    Eigen::VectorXd material_id_field;
    bool material_id_field_ready = false;

    // Champs de debug ParaView : visualisation des bandes decouplees de
    // joint de grain (fracture) et de verrouillage de polarisation.
    Eigen::VectorXd frac_band_field;
    bool frac_band_ready = false;

    Eigen::VectorXd polar_band_field;
    bool polar_band_ready = false;

public:
    IOManager(ResultsExporter& exp, const std::string& run_output_dir, const Datafile& conf);

    // Le PhysicsManager est maintenant passé en CONST (Garantie de non-modification par l'I/O)
    void extract_and_save_results(double t, int step, const PhysicsManager& pm, const Polycrystal& poly);

    // Ecriture finale (fin de run) de l'historique complet des energies,
    // dans run_output_dir/energies_final.csv.
    void finalize() const;

    const std::string& get_run_output_dir() const { return run_output_dir; }
};
