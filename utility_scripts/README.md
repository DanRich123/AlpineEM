# Utility Scripts Folder (In progress..)

This folder contains key pre- and post-processing scripts for the main FDTD solver. Below is the list of the scripts with a brief description of what each one is intended to do. Also see the examples folder [`examples/`](../examples) for examples where some of them are used.

| Script                      | Description |  Output  |
|-------------------------------|---------------------|----------------------|
| `post_processor.py`                | This script takes a few inputs from the user (top of file) and reads in `.dat` files produced from `master.py`        | S-parameters, realized gain, etc. as `.csv` files |
| `post_processor_kmax.py`                | This script takes a few inputs from the user (top of file) and reads in `.dat` files produced from `master.py`  for the `kmax` variant      | S-parameters, realized gain, etc. as `.csv` files |
| `spiral.py`    | A highly customizable script designed to create a 2-D Archimedean spiral-like pattern         | Numpy arrays that can be directly imported to `master.py` for creating spiral antenna patterns, or similar, from thin sheets |
| xx|xx | xx|
