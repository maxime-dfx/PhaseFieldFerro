#pragma once
//
// PhysicsConcepts.h
// ------------------
// Point d'entree UNIQUE pour toutes les interfaces (concepts) du module
// Physics : ISolver, IElementEquation, IPhysicsDofMapper, IPostProcessor,
// IPhysicsModule, ainsi que la structure transverse PhysicsState, et le
// template GenericPhysicsModule qui les assemble.
//
// Un nouveau module physique n'a besoin d'inclure QUE ce header (voir
// Modules/Mechanics.h, Modules/Electrostatics.h, Modules/Fracture.h,
// Modules/Polarization.h pour des exemples).
//
// NOTE DE LAYERING : ElementContext et FlattenedBC vivent dans
// Assembly/SystemAssembler.h (couche plus bas niveau, n'a besoin que de
// Mesh.h). ISolver est defini ICI, AVANT l'inclusion de SolverFactory.h,
// pour que GenericPhysicsModule puisse a la fois definir ISolver et
// utiliser SolverFactory sans dependance circulaire reelle.

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include <cstring>

#include "Mesh/include/Mesh.h"
#include "Mesh/include/IDofMap.h"
#include "Physics/include/Assembly/SystemAssembler.h"
#include "Utils/include/Logger.h"
#include "Utils/include/Profiling.h"

class MaterialManager;
class Polycrystal;
class Datafile;

// =====================================================================
//  PhysicsState : le "blackboard" partage par tous les modules physiques.
// =====================================================================

// Mettre à jour le blackboard :
struct PhysicsState {
    double time = 0.0;
    double dt = 0.0;
    double dt_load = 0.0;
    std::unordered_map<std::string, const Eigen::VectorXd*> fields;

    const MaterialManager* materials_manager = nullptr;
    const Polycrystal* polycrystal = nullptr;
    const IDofMap* p_dof_map = nullptr;
    const Datafile* config = nullptr;

    const Eigen::VectorXd* get_field(const std::string& name) const {
        auto it = fields.find(name);
        return (it != fields.end()) ? it->second : nullptr;
    }
};

// =====================================================================
//  IElementEquation (ElementContext vient de Assembly/SystemAssembler.h)
// =====================================================================
class IElementEquation {
public:
    virtual ~IElementEquation() = default;
    virtual void compute_element_matrices(
        const ElementContext& ctx,
        const PhysicsState& state,
        Eigen::Ref<Eigen::MatrixXd> K_local,
        Eigen::Ref<Eigen::VectorXd> F_local,
        std::vector<int>& global_dofs) const = 0;
};

// =====================================================================
//  IPhysicsDofMapper (FlattenedBC vient de Assembly/SystemAssembler.h)
// =====================================================================
class IPhysicsDofMapper {
public:
    virtual ~IPhysicsDofMapper() = default;
    virtual int get_system_size() const = 0;
    virtual Eigen::VectorXd build_guess_vector() const = 0;
    virtual void map_solution_to_states(const Eigen::VectorXd& solution) = 0;
    virtual std::vector<FlattenedBC> get_bcs() const = 0;
    virtual void register_fields(PhysicsState& state) = 0;
    virtual void apply_initial_conditions() = 0;
    virtual void save_previous_iteration() = 0;
    virtual void save_previous_state() = 0;
    virtual void restore_previous_state() = 0;
    virtual void update_history() = 0;
    virtual double calculate_error() const = 0;
};

// =====================================================================
//  IPostProcessor
// =====================================================================
class IPostProcessor {
public:
    virtual ~IPostProcessor() = default;
    virtual void compute(const PhysicsState& state, const Mesh& mesh) = 0;
    // Optionnel : expose les champs internes du post-processeur dans le
    // blackboard (PhysicsState.fields) pour que l'export VTK (IOManager)
    // puisse les recuperer par nom, exactement comme les champs primaires
    // (Px, ux, phi...) le sont via IPhysicsDofMapper::register_fields.
    // Appele une seule fois, juste apres construction (cf. PhysicsManager).
    virtual void register_fields(PhysicsState& state) { (void)state; }
};

// =====================================================================
//  IPhysicsModule
// =====================================================================
class IPhysicsModule {
public:
    virtual ~IPhysicsModule() = default;
    virtual void compute_step(const PhysicsState& state) = 0;
    virtual double calculate_error() const = 0;
    virtual bool last_solve_failed() const = 0;
    virtual void save_previous_iteration() = 0;
    virtual void save_previous_state() = 0;
    virtual void restore_previous_state() = 0;
    virtual void update_history() = 0;

    // Permet au PhysicsManager d'identifier le module.
    virtual std::string get_name() const = 0;
};

// =====================================================================
//  ISolver
//  (defini ICI, AVANT l'include de SolverFactory.h : voir note de layering
//  en tete de fichier)
// =====================================================================
class ISolver {
public:
    virtual ~ISolver() = default;
    virtual void analyze_pattern(Eigen::SparseMatrix<double>& mat) = 0;
    virtual void factorize(Eigen::SparseMatrix<double>& mat) = 0;
    virtual Eigen::VectorXd solve(const Eigen::VectorXd& rhs, const Eigen::VectorXd& guess) = 0;
    virtual bool has_failed() const = 0;
    virtual int iterations() const = 0;
    virtual double error() const = 0;
};

#include "Physics/include/Solvers/SolverFactory.h"

// =====================================================================
//  GenericPhysicsModule<EquationType, MapperType>
//  Orchestration generique Equation + DofMapper + Solver + PostProcessors,
//  transforme en template pour l'inlining et l'effacement de type.
// =====================================================================
template <typename EquationType, typename MapperType>
class GenericPhysicsModule : public IPhysicsModule {
private:
    std::string m_name;
    const Mesh& m_mesh;
    std::unique_ptr<ISolver> m_solver;

    // Stockage PAR VALEUR du type exact (evite les allocations dynamiques et l'indirection)
    EquationType m_equation;
    MapperType m_mapper;

    std::vector<std::unique_ptr<IPostProcessor>> m_post_processors;

    Eigen::SparseMatrix<double> m_K_global;
    Eigen::VectorXd m_F_global;
    std::vector<long> m_csr_mapping;
    bool m_is_matrix_allocated = false;
    bool m_last_solve_failed = false;

public:
    GenericPhysicsModule(const std::string& name, const SolverConfig& solver_cfg, const Mesh& mesh,
                         EquationType eq, MapperType mapper)
        : m_name(name), m_mesh(mesh), m_equation(std::move(eq)), m_mapper(std::move(mapper))
    {
        // La Factory n'a plus besoin que de l'EquationType en tant que Tag
        m_solver = SolverFactory::create<EquationType>(solver_cfg);

        m_mapper.apply_initial_conditions();
    }

    void add_post_processor(std::unique_ptr<IPostProcessor> pp, PhysicsState& state) {
        pp->register_fields(state);
        m_post_processors.push_back(std::move(pp));
    }

    void compute_step(const PhysicsState& state) override {
        // Utilisation des variables statiques constexpr de l'equation !
        PROFILE_ZONE_NC(EquationType::ModuleName, EquationType::Color);

        int system_size = m_mapper.get_system_size();

        if (!m_is_matrix_allocated) {
            m_K_global.resize(system_size, system_size);
            m_F_global.resize(system_size);
        }

        std::vector<FlattenedBC> active_bcs = m_mapper.get_bcs();

        // Appel fortement optimise : SystemAssembler sait exactement ce qu'est m_equation
        SystemAssembler::assemble(
            m_mesh.get_num_elements(), system_size, !m_is_matrix_allocated,
            m_K_global, m_F_global, m_csr_mapping,
            m_equation, active_bcs, state,
            m_mesh
        );

        if (!m_is_matrix_allocated) {
            m_solver->analyze_pattern(m_K_global);
        }

        m_is_matrix_allocated = true;

        m_solver->factorize(m_K_global);

        if (m_solver->has_failed()) {
            Logger::error("[", m_name, "] Echec de la factorisation.");
            m_last_solve_failed = true;
            return;
        }

        Eigen::VectorXd guess = m_mapper.build_guess_vector();
        Eigen::VectorXd sol = m_solver->solve(m_F_global, guess);

        m_last_solve_failed = (m_solver->has_failed() || !sol.allFinite());

        if (!m_last_solve_failed) {
            Logger::debug("[", m_name, "] solve OK : ||sol||=", sol.norm(),
                           " min=", sol.minCoeff(), " max=", sol.maxCoeff(),
                           " ||F||=", m_F_global.norm());
            m_mapper.map_solution_to_states(sol);
        } else {
            Logger::error("[", m_name, "] Echec de la resolution.");
        }

        {
            // Sous-zone de post-processing avec traits statiques
            PROFILE_ZONE_NC(EquationType::ZonePostProcess, EquationType::Color);
            for (auto& pp : m_post_processors) {
                pp->compute(state, m_mesh);
            }
        }
    }

    bool last_solve_failed() const override { return m_last_solve_failed; }
    void save_previous_iteration() override { m_mapper.save_previous_iteration(); }
    void save_previous_state() override { m_mapper.save_previous_state(); }
    void restore_previous_state() override { m_mapper.restore_previous_state(); }
    void update_history() override { m_mapper.update_history(); }
    double calculate_error() const override { return m_mapper.calculate_error(); }
    std::string get_name() const override { return m_name; }

    MapperType& get_mapper() { return m_mapper; }
};
