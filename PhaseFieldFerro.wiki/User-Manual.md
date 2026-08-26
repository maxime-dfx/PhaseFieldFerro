The entire simulation is controlled via a TOML configuration file parsed by the `Datafile` class. This document details the available parameters.

## 1. Global Simulation Parameters

The `[simulation]` section governs the time loop, convergence criteria, and output handling:

*   `output_dir`: Path to the directory where results will be saved.
*   `total_time`: Total simulation time (supports fractional values like 3.0).
*   `dt`: Macroscopic load step.
*   `dt_relax`: Relaxation step for the Picard iteration loop.
*   `max_iter` / `min_iter`: Maximum and minimum bounds for Picard iterations.
*   `tolerance_ferro` / `tolerance_vfield`: Convergence tolerances for the polarization and fracture fields, respectively.
*   `save_frequency`: Number of steps between each `.vtk` export.

## 2. Mesh Configuration

The `[mesh]` section dictates domain discretization:

*   `calcul_mesh`: Boolean. If `true`, generates an internal structured mesh. If `false`, loads an external `.msh` file.
*   `element_type`: `"QUAD"` (Q4 elements) or `"TRIANGLE"` (T3 elements).
*   `L_x`, `L_y`: Domain dimensions.
*   `n_x`, `n_y`: Number of elements along each axis (for internal generation).
*   `get_mesh_file`: Path to the external Gmsh file if `calcul_mesh = false`.

## 3. Materials and Microstructure

### Microstructure
The `[crystal]` section defines the polycrystalline topology using random Voronoi tessellation:
*   `num_grains`: Number of grains (1 = single crystal).
*   `seed`: Seed for the random number generator.
*   `grain_boundary_material_id`: Links the physical grain boundaries to a specific material ID defined in `[[materials]]`.
*   `grain_boundary_band_rings`: Thickness (in elements) of the fracture regularization band.
*   `polarization_lock_band_rings`: Thickness (in elements) of the depolarized zone.

### Material Definitions
Materials are defined as an array of tables using `[[materials]]`. They are mapped to the mesh via their `id` (corresponding to Gmsh physical tags). 

The `type` attribute defines the constitutive law:
*   `"Ferroelectric"`: Active piezoelectric ceramic with Landau dynamics.
*   `"PolymerElastic"`: Passive dielectric polymer (P is locked to 0).
*   `"PureElastic"`: Pure mechanical material without electrical properties.
*   `"GrainBoundary"`: Explicit grain boundary material.

Constants (`Gc`, `kappa`, `eps0`, `c1`, `c2`, `c3`, `a0`, etc.) can be overridden locally in each `[[materials]]` block.

## 4. Boundary Conditions (BCs)

Boundary conditions are strictly data-driven and defined using the `[[boundary_rule]]` array of tables. 

*   `field`: The target physical field (`"phi"`, `"ux"`, `"uy"`, `"Px"`, `"Py"`, `"v"`).
*   `bc_type`: `"DIRICHLET"` (prescribed value) or `"NEUMANN"` (flux).
*   `shape`: The geometric zone. Options include `"edge"`, `"point"`, `"circle"`, or `"rect"`.
    *   If `shape = "edge"`, specify `edge_name` (`"left"`, `"right"`, `"top"`, `"bottom"`).
*   `profile`: The temporal or spatial evolution function. Available profiles include:
    *   `"constant"`: Fixed value defined by `val`.
    *   `"time_ramp"`: Linear interpolation between `val_start` and `val_end` over the interval `[t_start, t_end]`.
    *   `"time_sine"`: Sinusoidal evolution using `amplitude`, `frequency`, and `offset`.
    *   `"spatial_tanh"`, `"spatial_parabola"`, `"traveling_wave"`: Advanced spatial profiles.