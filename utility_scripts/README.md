# Utility Scripts Folder (In progress..)

This folder contains key pre- and post-processing scripts for the main FDTD solver. Below is the list of the scripts with a brief description of what each one is intended to do. Also see the examples folder [`examples/`](../examples) for examples where some of the primary scripts are used.

| Script                      | Description |  Output  |
|-------------------------------|---------------------|----------------------|
| `post_processor.py`                | This script takes a few inputs from the user (top of file only) and reads in `.dat` files produced from `master.py`        | S-parameters, realized gain, etc. as `.csv` files |
| `post_processor_kmax.py`                | This script takes a few inputs from the user (top of file only) and reads in `.dat` files produced from `master.py`  for the `kmax` variant      | S-parameters, realized gain, etc. as `.csv` files |
| `python_fields_viewer.py` | This script takes a few inputs from the user (top of file only) and reads 2-D slices of binary (`.bin`) field data produced by `master.py` if that option was chosen by the user | 2-D time or frequency domain videos of the simulation |
| `paraview_fields_viewer.py` | This script takes a few inputs from the user (SETUP section) and reads 3-D slices of binary (`.bin`) field data produced by `master.py` if that option was chosen by the user | 3-D time or frequency domain Paraview files of the simulation - can view still frames and videos using Paraview software |
| `make_optional_geom_bulk.py` | This script takes inputs from the user to create an optional geometry `.bin` file, if desired. The conventional drawing of objects (blocks, cylinders, spheres) in `master.py` can be tedious and sometimes it is easier to create a geometry in a Numpy array, then save it to the binary file format required by the main FDTD solver. This script aids in doing this. The entire FDTD grid must be accounted for in the binary, and any object drawn in `master.py` will overwrite sections of this optional geometry. | A binary (`.bin`) file used by `master.py` to create a harder to draw geometry in the main FDTD solver |
| `make_optional_geom_and_conform_bulk.py` | Identical to the above script with the added benefit of drawing very detailed geometries (more refined than the FDTD grid itself) and using conformal averaging techniques to map it to the actual FDTD grid using anisotropic materials. Requires `conformal_builder.py`. This technique is useful for dielectrics and lossy materials, including Drude materials. | A binary (`.bin`) file used by `master.py` to create a harder to draw geometry in the main FDTD solver, along with a Numpy array of materials information created under conformal averaging for using in `master.py` |
| `conformal_builder.py` | Class based script required by the above `make_optional_geom_and_conform_bulk.py` | N/A |
| xx|xx | xx|
| `plot_kmax_f_vs_angle.py`|xx | xx|
| `plot_kmax_angle_vs_angle.py`|xx | xx|
| xx|xx | xx|
| xx|xx | xx|
| xx|xx | xx|
| `spiral.py`    | A highly customizable script designed to create a 2-D Archimedean spiral-like pattern         | Numpy arrays that can be directly imported to `master.py` for creating spiral antenna patterns, or similar, from thin sheets |
