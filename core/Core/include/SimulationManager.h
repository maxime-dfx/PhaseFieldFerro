#pragma once
#include <string>
#include "IO/include/Datafile.h"
#include "Mesh/include/Mesh.h"
#include "IO/include/ResultsExporter.h"
#include "Utils/include/Chrono.h"
#include "Physics/Core/BoundaryManager.h"
#include "Materials/Core/MaterialManager.h"
#include "Materials/Microstructure/Polycrystal.h"

// Manager physique : seul point de contact avec les 4 modules physiques
// (Mechanics, Electrostatics, Polarization, Fracture) et la boucle de Picard.
#include "Physics/include/Core/PhysicsManager.h"

// Nouveaux managers
#include "Core/include/TimeManager.h"
#include "IO/include/IOManager.h"

class SimulationManager {
private:
    // /!\ L'ORDRE DE DÉCLARATION DOIT CORRESPONDRE AU CONSTRUCTEUR POUR ÉVITER -Wreorder
    const Datafile& config;
    const Mesh& mesh;
    ResultsExporter& exporter;
    BoundaryManager boundary_manager;
    const MaterialManager& materials_manager;
    Chrono chrono;

    // Dossier de sortie UNIQUE A CE RUN (cf. Utils/RunId.h, resolu par
    // l'appelant avant construction). Tous les exports (mesh, energies,
    // VTK) passent par lui : deux runs partageant le meme config.toml
    // n'ecrivent donc jamais dans les memes fichiers.
    std::string run_output_dir;

    // Microstructure polycristalline (Voronoi + affectation element -> grain),
    // construite une fois ici et partagee par PhysicsManager (via Polarization)
    // et par IOManager (export du champ grain_id).
    Polycrystal m_polycrystal;

    // Seul point de contact avec les 4 modules physiques et la boucle de Picard
    PhysicsManager physics_manager;

    // Managers délégués
    TimeManager time_manager;
    IOManager io_manager;

public:
    SimulationManager(const Datafile& config, const Mesh& mesh, ResultsExporter& exporter,
                       const MaterialManager& registry, const std::string& run_output_dir);
    void initialize_mesh();
    void initialize_physics();
    void run();

    const std::string& get_run_output_dir() const { return run_output_dir; }

    // NOTE : la boucle de Picard (compute_one_step_physics), la validation
    // de l'historique (update_physics_history) ainsi que save_previous_states/
    // restore_previous_states (declarees mais jamais appelees ni definies
    // dans l'ancienne version) sont desormais entierement encapsulees dans
    // PhysicsManager. SimulationManager ne fait plus que l'orchestrer.
};
