# PhaseFieldFerro

PhaseFieldFerro is a multiphysics finite element simulation engine implementing a coupled ferroelectric/fracture phase-field model (Abdollahi & Arias, 2011).
It is written in C++17 and designed for high-performance computing with OpenMP parallelization.

## Supported Physics Modules

The codebase strongly couples four physical phenomena, solved iteratively:
- **Polarization (TDGL)**: Time-Dependent Ginzburg-Landau dynamics.
- **Mechanics**: Linear elasticity and electrostrictive coupling.
- **Electrostatics**: Poisson's equation for the electric potential.
- **Fracture (phase-field)**: Phase-field damage model driven by thermodynamic forces.

## Internal Architecture

The codebase is highly modular and separated into several key components:

- **`Core/` (Orchestration)**: `SimulationApp` loads the configuration and initializes the environment. `SimulationManager` orchestrates the fixed-step time loop (`TimeManager`).
- **`Physics/` (Equation Solving)**: `PhysicsManager` handles the Picard iteration loop to converge the 4 physics modules at each time step. Matrix assembly is highly optimized by `SystemAssembler` (with CSR format support).
- **`Materials/` (Constitutive Laws)**: Supports multiple models: Active ceramic (`Ferroelectric`), Passive matrix (`PolymerElastic`), Pure elasticity (`PureElastic`), and Grain boundaries (`GrainBoundary`).
- **`Microstructure/`**: Manages the creation of polycrystals via random Voronoi tessellation and assigns crystalline orientations (Euler angles) to the mesh elements.
- **`Mesh/` (Finite Elements)**: Manages reference finite elements (T3, Q4, T6, Q8), Gauss integration (`FEM.h`), and external mesh loading via Gmsh (`MeshGenerators.h`).
- **`IO/` (Input/Output)**: Parses TOML configuration files (`Datafile.h`) and exports regular results in `.vtk` (spatial visualization) and `.csv` (system metrics and energies) formats into a unique folder for each run.

## Dependencies

This project expects the following dependencies installed under a root folder pointed to by the `DEPS_ROOT` variable (default is `$HOME/dependances_pour_serveur`):
- **Eigen** (Matrix linear algebra).
- **Gmsh SDK** (Mesh generation and reading).
- **SuiteSparse / CHOLMOD** (Default direct solver).
- **mimalloc** (Optimized memory allocator).
- **PETSc / Hypre** (Optional, enabled via `--petsc` for advanced iterative solvers like BoomerAMG).
- **Intel MKL PARDISO** (Optional, enabled via `--pardiso` for an optimized multithreaded direct solver).
- **Tracy** (Optional, enabled via `--tracy` for real-time profiling).

## Configuration and Build

Ensure you have CMake (>= 3.15) and a C++17 compiler.

```bash
# Initial configuration
scripts/setup/configure.sh [--pardiso] [--petsc] [--tracy]

# Build in Release mode
cmake --build "$HOME/Min/build<corresponding_suffix>"