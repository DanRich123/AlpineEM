# Transmission and Reflection Through a Multilayered Quasi-1D Infinite Slab of Dielectrics and Thin Sheets (AlpineEM FDTD)

This example uses AlpineEM to perform FDTD simulations that calculate the transmission and reflection through a quasi-1D infinite multilayered slab of materials. The slab of material is infinite in the y and z directions, and finite in the x direction. The infinite condition is created using periodic boudnary conditions. The multilayered stack up is as follows (also see Paraview image below):

- Layer 1: 50 ohms/sq thin sheet
- Layer 2: 6mm of dielectric 2.1
- Layer 3: 25 ohms/sq thin sheet
- Layer 4: 6mm of dielectric 3.5

It outputs time-domain data for post-processing, along with binary geometry files for the **object** (multilayer present), **clear** (multilayer absent), and **metal** (multilayer replaced with metal) cases.

This example is designed to teach and validate the core mechanisms of the solver including normal incidence plane waves with periodic boundary conditions, dielectric materials, the thin sheets method, and the three builder options (single-threaded CPU, multi-threaded CPU (OpenMP), and the GPU (OpenACC)).

## Workflow

### 1. Compile the FDTD solver

Choose one of three build options depending on the resources available to you:

- **Single-threaded CPU**
- **Multi-threaded CPU** (OpenMP)
- **GPU** (OpenACC)

There are batch scripts included for building all versions (`compile.sh`, `compile_mp.sh`, and `compile_acc.sh`).

### 2. Run the object case — `master.py`

Configures and runs the simulation with the object present.

- Generates a text file of inputs, then executes the binary compiled in Step 1.
- You must set the solver name in `master.py` (the single-threaded version is selected by default).
- Produces several output files used for post-processing and geometry viewing.

  > **Note** Run for zero time steps if viewing the geometry (see step 5) is desired before running a full simulation.

### 3. Run the clear case — `master_clear.py`

Performs the same steps as `master.py`, but without the object present, to establish the reference (background) fields.

- You must set the solver name in `master_clear.py` (the single-threaded version is selected by default).
- Produces several output files used for post-processing and geometry viewing.

### 4. Run the metal case — `master_metal.py`

Performs the same steps as `master.py`, but the object has been replaced with metal, to finalize the last piece needed for TRL calibration.

- You must set the solver name in `master_metal.py` (the single-threaded version is selected by default).
- Produces several output files used for post-processing and geometry viewing.
- Note: set `use_metal=True` to `use_metal=False` within `post_processor.py` if the user wishes to bypass using the metal case.
- Using the metal case provides more accurate results, but it shift the phase center away from the wave port to the surface of the metal.

### 5. Post-process — `post_process.py`

Combines the object, clear, and metal case outputs to generate transmission and reflection data as a `.csv` file.

- Example Slurm batch scripts are included and can be adapted to your cluster environment.

### 6. (Optional) View the geometry

To visualize the simulation geometry:

1. Run `fdtd_geometry_maker.py` to generate ParaView files.
2. Open the resulting **single** ParaView file directly in ParaView — it references an accompanying folder of associated files, so leave that folder in place and do not open its contents individually.
3. For an easier setup, load the included macro, `fdtd_macro.py`, into ParaView (**Macros** tab) to automatically configure common viewing filters.

<p align="center">
  <img src="./geometry.png" alt="Model geom" width="50%" />
</p>

### 7. Validation

The `plot.py` file can be used to compare the FDTD result with the analytic result:

![Model Accuracy](./Comparison_plot.png)

## Performance reference

Approximate per-simulation runtimes measured on the author's hardware:

| Solver                       | Time per simulation |
|-------------------------------|---------------------|
| OpenACC (GPU)                 | ~6.5 seconds         |
| OpenMP (multi-threaded CPU)   | ~1.8 seconds          |
| Single-threaded (default)     | ~3.8 seconds          |

> **Note:** These timings depend heavily on hardware, problem size, and system load. Use them only as a rough point of reference, not a direct benchmark against other software. OpenACC (GPU) performed worse due to the size and scaling of the problem. Unless the simulation is large enough, the overhead required for GPU calculations will dominate. I.e., this problem was not large enough in size for the GPU to dominate. See the monstatic scattering of a sphere example where the GPU handedly dominates the simulation time.
