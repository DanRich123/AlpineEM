# Utility Scripts Folder

This folder contains the key pre- and post-processing scripts for the main FDTD solver. Below is a list of the scripts with a brief description of what each one does. See the [`examples/`](../examples) folder for examples of how some of the primary scripts are used.

## Main Utility Scripts

| Script | Description | Output |
|---|---|---|
| `post_processor.py` | Takes a few user inputs (top of file only) and reads `.dat` files produced by `master.py`. | S-parameters, realized gain, etc. as `.csv` files |
| `post_processor_kmax.py` | Same as above, for the `kmax` variant. | S-parameters, realized gain, etc. as `.csv` files |
| `python_fields_viewer.py` | Takes a few user inputs (top of file only) and reads 2-D slices of binary (`.bin`) field data produced by `master.py` (if that option was selected). | 2-D time- or frequency-domain videos of the simulation |
| `paraview_fields_viewer.py` | Takes a few user inputs (SETUP section) and reads 3-D slices of binary (`.bin`) field data produced by `master.py` (if that option was selected). | 3-D time- or frequency-domain ParaView files — view still frames or videos in ParaView |
| `make_optional_geom_bulk.py` | Creates an optional geometry `.bin` file. Drawing objects (blocks, cylinders, spheres) conventionally in `master.py` can be tedious; this script makes it easier to build a geometry as a NumPy array and save it to the binary format required by the solver. The entire FDTD grid must be accounted for in the binary — any object drawn in `master.py` will overwrite sections of it. | A `.bin` file used by `master.py` for geometries that are otherwise hard to draw |
| `make_optional_geom_and_conform_bulk.py` | Same as above, with the added ability to draw more detailed geometries (finer than the FDTD grid itself) and map them to the grid using conformal averaging with anisotropic materials. Useful for dielectrics and lossy materials, including Drude materials. Requires `conformal_builder.py`. | A `.bin` file for hard-to-draw geometries, plus a NumPy array of conformally-averaged materials data for use in `master.py` |
| `conformal_builder.py` | Class-based script required by `make_optional_geom_and_conform_bulk.py`. | N/A |


## Other Utility Scripts

| Script | Description | Output |
|---|---|---|
| `spiral.py` | Highly customizable script for creating a 2-D Archimedean-spiral-like pattern. | NumPy arrays that can be imported directly into `master.py` to create spiral antennas (or similar) from thin sheets |
| `plot_kmax_f_vs_angle.py` | Highly customizable script for plotting frequency vs. angle results from `kmax`. | Image files |
| `plot_kmax_angle_vs_angle.py` | Highly customizable script for plotting angle-1 vs. angle-2 results from `kmax`. | Image files |
| `parallel_kmax.py` | Highly customizable script for setting up and running many `kmax` simulations in parallel. | Multiple results folders and files |
| `import_tecplot_and_interpolate.py` | Highly customizable script for importing, interpolating, and preparing FDTD geometries from custom CFD simulations on structured grids. | FDTD-ready optional geometry files (with appropriate grid sizes printed out) and a materials-info NumPy array |
| `import_plasma_properties.py` | Highly customizable script showing how to use the materials results from `import_tecplot_and_interpolate.py` in `master.py`. | — |
| `convert_to_touchstone.py` | Highly customizable script for creating Touchstone files and renormalized results from `kmax`. | Touchstone files |
