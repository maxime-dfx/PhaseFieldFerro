# Inputs, Outputs and Diagnostics

PhaseFieldFerro features a robust Input/Output management system to ensure data from different simulation runs are safely stored and easily visualizable.

## 1. Run Organization and Safety

To prevent accidental data overwriting when running multiple simulations simultaneously with the same `config.toml`, the `SimulationApp` utilizes a `RunId` generator. 
*   A unique output folder is created for each execution (e.g., `<output_dir>/run_<RunId>`).
*   All subsequent exports (meshes, energies, VTK files) are strictly directed into this isolated directory.

## 2. VTK Exports (`ResultsExporter`)

Results are exported in the `.vtk` format for spatial visualization in software like ParaView. The exports are handled by the `IOManager` based on the `save_frequency` parameter.

The codebase supports exporting:
*   **Point Data (Nodal fields):** Continuous physical fields such as electric potential (`phi`), displacements (`ux`, `uy`), and derived fields like stress components (`sigma_xx`, `von_mises`).
*   **Cell Data (Element fields):** Fields that are constant per element, which is critical for visualizing the microstructure. This includes the `grain_id` and the `material_id`.
*   **Decoupled Debug Bands:** Specific debug fields (`fracture_band`, `polarization_band`) are exported to visualize the exact elements affected by grain boundary modifications or polarization locking.

## 3. Diagnostics and Energy Tracking

The `Diagnostics` module continuously tracks the thermodynamic state of the system. 

*   At each load step, the `PhysicsManager` computes system metrics including the total elastic energy, electric energy, Landau energy, surface energy (fracture), and bulk enthalpy.
*   These metrics are appended in real-time to a CSV file (`energies.csv`) in the run directory. 
*   This file is invaluable for plotting load-displacement curves and monitoring the global energy dissipation during crack propagation.