# A Non-Foster Recieving Monopole Antenna (AlpineEM FDTD) (Work in progress....)

This example uses AlpineEM to perform FDTD simulations that calculate the effective aperture area of a non-Foster loaded monopole antenna over an infinite ground plane. A uniform plane wave incident on the antenna is received at the non-Foster loaded port (square coaxial port), and the effective aperture area is calculated from the received voltage and the incident plane wave information. 

It outputs time-domain data for post-processing, along with a binary geometry file. Because SPICE is used, the port impedance needs to be corrected after post processing - the post processing uses 50 ohms (constant over frequency) for the impedance instead of the correct impedance value. This value is purely for post processing - the correct port impedance was used in the simulation, but the post processor doesn't natively know this value. The file YY... corrects for this...

Though a square coaxial-like port is used, there is another script in ZZ that supports circular ports as well. The user can define any port shape they want via the static solver, but square and circle coax examples are available...

## Workflow

### 1. Compile the FDTD solver

Choose one of three build options depending on the resources available to you:

- **Single-threaded CPU**
- **Multi-threaded CPU** (OpenMP)
- **GPU** (OpenACC)

### 2. Run the static solver to produce the feed pattern... — `....py`

Configures and runs...

- Generates... used for...
- You must set...
- Produces several output files used for...

### 4. Run the object case — `master.py`

Configures and runs the simulation.

- Generates a text file of inputs, then executes the binary compiled in Step 1.
- You must set the solver name in `master.py` (the GPU/OpenACC version is selected by default).
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
| OpenACC (GPU)                 | ~xx seconds         |
| OpenMP (multi-threaded CPU)   | ~xx minutes          |
| Single-threaded (default)     | ~xx minutes          |

> **Note:** These timings depend heavily on hardware, problem size, and system load. Use them only as a rough point of reference, not a direct benchmark against other software. Note why OpenMP outperforms here...
