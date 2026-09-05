# A Non-Foster Receiving Monopole Antenna (AlpineEM FDTD) (Work in progress...)

This example uses AlpineEM to run FDTD simulations that calculate the effective aperture area of a non-Foster loaded monopole antenna over an infinite ground plane. A uniform plane wave incident on the antenna is received at the non-Foster loaded port (a square coaxial port), and the effective aperture area is calculated from the received voltage and the incident plane wave information.

The simulation outputs time-domain data for post-processing, along with a binary geometry file. Because SPICE is used, some post-processed values need correcting: by default, the post processor performs its algebraic calculations using a constant 50 ohm impedance. This value is used for post-processing only — the correct port impedance is used in the simulation itself — but since the post processor has no way of knowing that value, the user must manually correct for the actual impedance.

This example uses a square coaxial-like port, but [`statics_solver/`](./statics_solver) also includes `coax_example.py`, which supports circular coaxial ports. The static solver can accommodate any port shape, but square and circular coax examples are provided since they're the most common.

## Workflow

### 1. Compile the FDTD solver
Choose one of three build options depending on the resources available to you:
- **Single-threaded CPU**
- **Multi-threaded CPU** (OpenMP)
- **GPU** (OpenACC)

### 2. Run the static solver to produce the TEM E- and H-field weightings — `square_coax_example.py`
This Python script calls `EM2Dsolver.py` to solve the statics problem. The resulting weights also apply to TEM modes when used appropriately.
- Generates a binary file of E- and H-field weightings, which is read directly into and executed by the binary compiled in Step 1.
- The user must set the permittivity and permeability of the 2D geometry to match the shape used in `master.py`.

### 3. Run the object case — `master.py`
Configures and runs the simulation.
- Generates a text file of inputs, then executes the binary compiled in Step 1.
- You must set the solver name in `master.py` (the multi-threaded CPU/OpenMP version is selected by default).
- When entering gridded feed information, manually enter the name of the binary created in Step 2.
- Produces several output files used for post-processing and geometry viewing.

### 4. Post-process — `post_process.py`
Uses the simulation outputs to generate the effective aperture area data as a `.csv` file.
- Example Slurm batch scripts are included and can be adapted to your cluster environment.

### 5. (Optional) View the geometry
To visualize the simulation geometry:
1. Run `fdtd_geometry_maker.py` to generate ParaView files.
2. Open the resulting **single** ParaView file directly in ParaView — it references an accompanying folder of associated files, so leave that folder in place and don't open its contents individually.
3. For easier setup, load the included macro, `fdtd_macro.py`, into ParaView (**Macros** tab) to automatically configure common viewing filters.

## Performance reference
Approximate per-simulation runtimes measured on the author's hardware:

| Solver                     | Time per simulation |
|-----------------------------|---------------------|
| OpenACC (GPU)               | ~2.3 minutes         |
| OpenMP (multi-threaded CPU) | ~2.1 minutes         |
| Single-threaded (default)   | ~2.9 minutes         |

> **Note:** SPICE accounts for most of the runtime in this particular example, due to a small SPICE time step, and it doesn't run natively on the GPU — which is why OpenMP outperforms OpenACC here. See the metal sphere monostatic scattering example for a case with more drastic runtime differences, where OpenACC dominates. These timings depend heavily on hardware, problem size, and system load — use them only as a rough point of reference, not a direct benchmark against other software.
