# Architecture and Internal Mechanics

This section details the internal mechanics of PhaseFieldFerro. The engine is designed with a strictly modular architecture to separate time orchestration, equation solving, and data management.

## 1. Separation of Concerns

The main execution flow is divided into two management layers:
*   **`SimulationManager`**: Acts as the high-level orchestrator. It controls the global fixed-step time loop via `TimeManager` (following the Abdollahi & Arias algorithm), manages outputs through `IOManager`, and delegates physical computations to the `PhysicsManager`.
*   **`PhysicsManager`**: Exclusively orchestrates the Picard iteration loop. At each time step, it iterates through all active physical modules (Polarization, Mechanics, Electrostatics, Fracture) until global convergence criteria are met.

## 2. The Blackboard: `PhysicsState`

To prevent tight coupling between physical modules, data exchange is handled via a centralized "Blackboard" data structure named `PhysicsState`. 

This structure contains:
*   Current time and time steps (`time`, `dt`, `dt_load`).
*   Pointers to global managers (`MaterialManager`, `Polycrystal`, `Datafile`).
*   A dictionary of physical fields (`std::unordered_map<std::string, const Eigen::VectorXd*> fields`). 

For instance, the Mechanics module can read the polarization field by querying `state.get_field("Px_prev")` without requiring any direct dependency on the Polarization module.

## 3. Anatomy of a Physics Module

Each physical phenomenon is encapsulated within a `GenericPhysicsModule<EquationType, MapperType>`. This template-based approach allows the compiler to aggressively inline critical code, particularly during matrix assembly.

A physical module aggregates four main concepts (interfaces):

*   **`IElementEquation`**: Defines the pure mathematical formulation. The `compute_element_matrices` method takes an element's coordinates and state, and computes the local stiffness matrix (`K_local`) and right-hand side (`F_local`) via Gauss integration.
*   **`IPhysicsDofMapper`**: Manages the topology. It links geometric nodes to global degrees of freedom (DOFs), builds flattened boundary conditions (`FlattenedBC`), and stores the field's history (e.g., `current`, `prev_iter`, `n`).
*   **`ISolver`**: Solves the global linear system. Instantiated via the `SolverFactory`, it can utilize native iterative solvers (CG, BiCGSTAB), PARDISO, or PETSc/Hypre.
*   **`IPostProcessor`** (Optional): Computes derived fields after a step is solved. For example, `StressPostProcessor` calculates stresses (`sigma_xx`, `von_mises`) based on displacements.

## 4. Optimized Assembly (`SystemAssembler`)

To ensure high performance, the assembly of the global system `K * U = F` is delegated to the `SystemAssembler`. 

It operates in two phases:
1.  **First Assembly (Triplets):** The matrix is constructed using a `std::vector<Eigen::Triplet>`. Concurrently, a cache (`csr_mapping`) is generated to store the exact memory address of each local interaction within the global CSR matrix.
2.  **CSR Bypass (Subsequent Iterations):** For subsequent Picard iterations, the assembler uses the `csr_mapping` to write directly into the underlying arrays of the `Eigen::SparseMatrix` via `valuePtr()`. This avoids the substantial cost of reallocation and significantly accelerates the computation. The assembly is parallelized without data races using a mesh graph coloring strategy (`color_groups`).