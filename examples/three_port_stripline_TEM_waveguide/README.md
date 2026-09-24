# A Three Port Stripline-like TEM Waveguide (AlpineEM FDTD) (In progres...)

This example uses AlpineEM to run FDTD simulations that calculate the S-parameters of a 3-port stripline-like waveguide that is TEM below around 250 MHz or so. Each of the 3 ports is a square coaxial-like gridded port, matched to 50-ohms, though in one case there is a non-Foster circuit before the 50-ohm load of port 3. 

This example builds on a prior example given in `/examples/non-Foster_recieving_monopole` and utilizes the method developed at GTRI, published in [1], for determining the realized gain of an electrically small antenna using a transmission line measurement within a waveguide to compare the waveguide extracted realized gain with the free space (w/ infinite ground plane) determined realized gain.

The simulation outputs time-domain data for post-processing, along with a binary geometry file. Because SPICE is used, some post-processed values need correcting: by default, the post processor performs its algebraic calculations using a constant 1 ohm impedance. This value is used for post-processing only — the correct port impedance is used in the simulation itself — but since the post processor has no way of directly knowing that value, the user must sometimes manually correct for the actual impedance in the `.csv` output files (there is a note within the header if normalization is needed).

This example uses two different square coaxial-like ports, but [`statics_solver/`](./statics_solver) also includes `coax_example.py`, which supports circular coaxial ports. The static solver can accommodate any port shape, but square and circular coax examples are provided since they're the most common.

This example is primarily designed to demonstrate how to use SPICE, gridded feeds, and extract S-parameters. Additionally, the self-consistency and accuracy between free space (w/ infinite ground plane) determined realized gain and transmission line waveguide produced realized gain emphasizes the accuracy of the solver given two very different simulations and processes to arrive at the same result. Note the final mismatch in data can very likely be attributed to cell size and resolution given the 0.5dB shift in free space (w/ infinite ground plane) model results with a 1/8 cell size reduction (1/2 size reduction in each direction).

## Workflow

### 1. Compile the FDTD solver
Choose one of three build options depending on the resources available to you:
- **Single-threaded CPU**
- **Multi-threaded CPU** (OpenMP)
- **GPU** (OpenACC)

There are batch scripts included for building all versions (`compile_spice.sh`, `compile_spice_mp.sh`, and `compile_spice_acc.sh`). **However, due to the size of the simulation, using resources other than a GPU (OpenACC) will results in several hours of computation time.**

> **Note** SPICE must be installed. There are two example batch scripts included in this folder, `build_my_ngspice.sh` and `build_my_ngspice_acc.sh`, for installing and exporting paths for the different build options. They will likely need to be modified to include the user's specific information.

### 2. Run the static solver to produce the TEM E- and H-field weightings — `square_coax_example.py`
This Python script calls `EM2Dsolver.py` to solve the statics problem. The resulting weights also apply to TEM modes when used appropriately.
- Generates a binary file of E- and H-field weightings, which is read directly into and executed by the binary compiled in Step 1.
- The user must set the permittivity and permeability of the 2D geometry to match the shape used in `master.py`.
- For reference, the expected binary output named `gridded_feed_12.bin` is included in this tutorial.
- Additionally, the E and H field plots are included with voltage and Az (normal component) as well:
<p align="center">
  <img src="./electric fields_12.png" alt="Model E" width="49%" />
  <img src="./magnetic fields_12.png" alt="Model H" width="49%" />
</p>

Anywhere the weightings are zero within the mask (outside the coaxial cable), the FDTD solver will ignore those cells as part of the port.

**The gridded feed file created in `/examples/non-Foster_recieving_monopole` was also copied over and used here as well.**

> **Note** Z-direction is always normal in this statics solver. This mask can still be used, as is done in this case, when a different direction is actually normal in the FDTD solver. The FDTD solver accounts for this information correctly as long as the intended FDTD direction in `master.py` is selected.

### 3. Create the waveguide through the optional geometry method — `make_optional_geom_bulk.py`
Creates a geometry binary file to be read in and used by the `master.py` script. Any geometry created in `master.py` after the import will overwrite relevant sections of the optional geometry immport. For example, the coaxial feeds will be manually drawn to overwrite portions of the waveguide design. `make_optional_geom_bulk.py` will import a `.npy` file that can be created using `build_waveguide.py` that was created by an AI for extracting the approximate waveguide design from [1].

### 4. Run the object case — `master.py`
Configures and runs the simulation.
- Generates a text file of inputs, then executes the binary compiled in Step 1.
- You must set the solver name in `master.py` (the GPU version using OpenACC is selected by default here for drastic speed improvements).
- When entering gridded feed information, manually enter the name of the binary created in Step 2, along with the other binary created from `/examples/non-Foster_recieving_monopole`.
- Produces several output files used for post-processing and geometry viewing.

  > **Note** Run for zero time steps if viewing the geometry (see step 6) is desired before running a full simulation.

### 5. Post-process — `post_processor.py`
Uses the simulation outputs to generate the S-parameter data as a `.csv` file.
- Example Slurm batch scripts are included and can be adapted to your cluster environment.
- A normalization is needed to account for the correct impedance of each port location.

### 6. Convert aperture area to realized gain and plot both — `plot.py`
Imports the `.csv` file for the S-parameters, correct them, determine the realized gain using [1], and then plot against data from `/examples/non-Foster_recieving_monopole` using the small cell size variant:

<p align="center">
  <img src="./Realized gain.png" alt="Model G" width="50%" />
</p>

### 7. (Optional) View the geometry
To visualize the simulation geometry:
1. Run `fdtd_geometry_maker.py` to generate ParaView files.
2. Open the resulting **single** ParaView file directly in ParaView — it references an accompanying folder of associated files, so leave that folder in place and don't open its contents individually.
3. For easier setup, load the included macro, `fdtd_macro.py`, into ParaView (**Macros** tab) to automatically configure common viewing filters.

The Paraview rendered image is:

<p align="center">
  <img src="./geometry.jpg" alt="Model geom" width="80%" />
</p>

## Performance reference
Approximate per-simulation runtimes measured on the author's hardware:

| Solver                     | Time per simulation |
|-----------------------------|---------------------|
| OpenACC (GPU)               | ~22 minutes         |
| OpenMP (multi-threaded CPU) | --         |
| Single-threaded (default)   | --       |

> **Note:** Due to the size and number of time steps, OpenACC was over 15x faster than the default version. These timings depend heavily on hardware, problem size, and system load — use them only as a rough point of reference, not a direct benchmark against other software.

## References
[1] Richardson et al....
