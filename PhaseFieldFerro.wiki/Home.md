# Getting Started

This guide outlines the installation, compilation, and execution of the PhaseFieldFerro simulation engine.

## 1. Project Overview

PhaseFieldFerro is a multiphysics finite element simulation engine developed in C++17. It implements a coupled ferroelectric/fracture phase-field model based on the reference work by Abdollahi & Arias (2011).

The engine is designed for high-performance computing, utilizing OpenMP parallelization, and iteratively solves four strongly coupled physical modules:
*   **Polarization (TDGL)**: Time-Dependent Ginzburg-Landau dynamics for the time evolution of ferroelectric domains.
*   **Mechanics**: Linear elasticity and electrostrictive coupling.
*   **Electrostatics**: Resolution of the electric potential via Poisson's equation.
*   **Fracture**: Phase-field damage model driven by thermodynamic forces.

## 2. Installation and Build

### Prerequisites
*   A compiler supporting the **C++17** standard.
*   **CMake** (version 3.15 or higher).

### Dependencies Management
The project expects to find its dependencies in a directory defined by the `DEPS_ROOT` environment variable (defaulting to `$HOME/dependances_pour_serveur`). 
The required dependencies are:
*   **Eigen** (Matrix linear algebra).
*   **Gmsh SDK** (Mesh generation and reading).
*   **SuiteSparse / CHOLMOD** (Default direct solver).
*   **mimalloc** (Optimized memory allocator).

*(Optional flags allow enabling **PETSc** via `--petsc`, **Intel MKL PARDISO** via `--pardiso`, or **Tracy Profiler** via `--tracy` if they are installed on the host system)*.

### Building the Project
A utility script is provided to generate the build directory and configure CMake. At the root of the project, execute:

```bash
# Standard configuration
./scripts/setup/configure.sh

# Configuration with advanced solvers or profiling enabled
./scripts/setup/configure.sh --pardiso --petsc --tracy