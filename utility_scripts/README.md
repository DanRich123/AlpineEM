# Utility Scripts

Pre- and post-processing tools for the main FDTD solver. See [`examples/`](../examples) for demonstration workflows.

## Standard Workflow Scripts

Core utility scripts for everyday simulation setup, execution, and field analysis.

| Script | Description | Input / Dependencies | Output |
|---|---|---|---|
| `post_processor.py` | Extracts S-parameters, realized gain, and far-field radiation patterns from solver runs. | `.dat` simulation files | `.csv` performance metrics |
| `post_processor_kmax.py` | Post-processing variant customized for `kmax` simulation sweeps. | `.dat` files from `kmax` | `.csv` sweep data |
| `python_fields_viewer.py` | Generates 2-D field distribution plots and video animations. | 2-D field files (`.bin`) | Videos and plot images |
| `paraview_fields_viewer.py` | Exports 3-D volumetric field slices and timeframes for 3-D rendering. | 3-D field files (`.bin`) | ParaView VTK files |
| `make_optional_geom_bulk.py` | Constructs custom 3-D grid geometries in NumPy when native drawing is tedious. | NumPy | `optional_geom_bulk.bin` |
| `make_optional_geom_and_conform_bulk.py` | Generates fine-resolution sub-pixel conformal geometries ($\epsilon$, $\sigma$). | `conformal_builder.py` | `optional_geom_bulk.bin`, `materials_id_opfile.npy` |
| `conformal_builder.py` | Engine class (`ConformalGeometry`) for sub-cell material averaging and slice plotting. | `numpy`, `matplotlib` | Object class dependency (no Drude of permeability support yet) |


## Specialized & Non-Routine Scripts

Advanced scripts designed for specialized workflows (CFD-plasma coupling, parallel parameter sweeps, Touchstone conversions, and complex patterns).

| Script | Description | Input / Dependencies | Output |
|---|---|---|---|
| `import_tecplot_and_interpolate.py` | Interpolates 2-D CFD Tecplot flowfields into 3-D grids and computes plasma ($\omega_p, \gamma$) values. | Tecplot `.plt`, `scipy`, `pandas` | `optional_geom_bulk.bin`, `material_properties.npy`, diagnostic plots |
| `import_plasma_properties.py` | Formats `material_properties.npy` data into multi-pole plasma definitions for the solver. | `material_properties.npy` | Solver plasma input streams |
| `convert_to_touchstone.py` | Converts S-parameter CSVs into `.sNp` files with angle-dependent $Z_0$ normalization to $50\,\Omega$. | S-parameter `.csv`, `data.dat`, `skrf` | `.sNp` Touchstone files |
| `parallel_kmax.py` | Manages parallel execution of multiple `kmax` solver runs across compute threads. | `master.py`, binaries | Subdirectories with sweep runs |
| `spiral.py` | Generates 2-D Archimedean spiral coordinates for planar antenna creation. | NumPy | Coordinate arrays |
| `plot_kmax_f_vs_angle.py` | Plots frequency vs. incident/observation angles for `kmax` datasets. | `kmax` sweep datasets | Plot files (`.png`, `.pdf`) |
| `plot_kmax_angle_vs_angle.py` | Plots angle-1 vs. angle-2 variation mappings for `kmax` datasets. | `kmax` sweep datasets | Plot files (`.png`, `.pdf`) |
