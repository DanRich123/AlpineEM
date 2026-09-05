# A Non-Foster Recieving Monopole Antenna (AlpineEM FDTD) (Work in progress....)

This example uses AlpineEM to perform FDTD simulations that calculate the effective aperture area of a non-Foster loaded monopole antenna over an infinite ground plane. A uniform plane wave incident on the antenna is received at the non-Foster loaded port (square coaxial port), and the effective aperture area is calculated from the received voltage and the incident plane wave information. 

It outputs time-domain data for post-processing, along with a binary geometry file. Because SPICE is used, post processor produced values that sometimes need correcting. The post processing uses 50 ohms (constant over frequency) for the impedance by default when performing calculations algebraic calculations. This value is purely for post processing - the correct port impedance was used in the simulation, but the post processor doesn't natively know this value. Thus, the user must manually correct for the actual impedance.

Though a square coaxial-like port is used, there is another script in [`statics_solver/`](./statics_solver) that supports circular coaxial ports as well - `coax_example.py`. The user can define any port shape they want via the static solver, but square and circle coax examples are available since these are most common.

## Workflow

### 1. Compile the FDTD solver

Choose one of three build options depending on the resources available to you:

- **Single-threaded CPU**
- **Multi-threaded CPU** (OpenMP)
- **GPU** (OpenACC)

### 2. Run the static solver to produce the TEM E and H field weightings — `square_coax_example.py`

This python script calls `EM2Dsolver.py` to solve the statics problem. The weights also apply to TEM modes if used appropriately.

- Generates a binary file for E and H field weightings that will be directly read into and executed via the binary compiled in Step 1.
- User must set the permittivity and permeability of the 2D geometry, it must match the shape intended in the `master.py` script.

### 4. Run the object case — `master.py`

Configures and runs the simulation.

- Generates a text file of inputs, then executes the binary compiled in Step 1.
- You must set the solver name in `master.py` (the multi-threaded CPU version (OpenMP) is currently selected by default).
- When entering gridded feed information, manually enter the name of the binary created in Step 3.
- Produces several output files used for post-processing and geometry viewing.

### 5. Post-process — `post_process.py`

Uses the outputs to generate the effective aperture area data as a `.csv` file.

- Example Slurm batch scripts are included and can be adapted to your cluster environment.

### 6. (Optional) View the geometry

To visualize the simulation geometry:

1. Run `fdtd_geometry_maker.py` to generate ParaView files.
2. Open the resulting **single** ParaView file directly in ParaView — it references an accompanying folder of associated files, so leave that folder in place and do not open its contents individually.
3. For an easier setup, load the included macro, `fdtd_macro.py`, into ParaView (**Macros** tab) to automatically configure common viewing filters.

## Performance reference

Approximate per-simulation runtimes measured on the author's hardware:

| Solver                       | Time per simulation |
|-------------------------------|---------------------|
| OpenACC (GPU)                 | ~2.3 minutes          |
| OpenMP (multi-threaded CPU)   | ~2.1 minutes          |
| Single-threaded (default)     | ~2.9 minutes          |

> **Note:** SPICE is the larger time consumption in this specific example, and it does not run natively on the GPU - leading OpenMP to produce a faster runtime in this specific example. See the metal sphere monostatic scattering example for more drastic differences in runtime where OpenACC dominates. These timings depend heavily on hardware, problem size, and system load. Use them only as a rough point of reference, not a direct benchmark against other software. 
