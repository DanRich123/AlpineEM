# Time and Frequency Electromagnetic Field Visualization (AlpineEM FDTD)

This example uses AlpineEM to perform FDTD simulations of a plane wave incident on a dielectric cube in free space. It outputs electromagnetic (EM) field data for plotting, along with a binary geometry file. Unlike other examples, it does not include post-processing to determine derived quantities such as S-parameters or realized gain.

The sole purpose of this example is to demonstrate the capabilities of the EM field viewer options. A few of the capabilities are shown explicitly, while many others are only highlighted. When simulating EM problems, it is often useful to view the fields in time (or in the frequency domain) to better understand what is physically happening within the simulation, even if the derived post-processed quantities are the values of actual interest.

There are two options for viewing fields:

- **Python viewer**: creates videos or still frames of 2D (cell-centered) slices.
- **ParaView viewer**: creates videos or still frames of the full 3D fields. This option still uses Python to create the ParaView files, similar to the geometry viewer. It is the most general option, and 2D slices can also be created in ParaView, but the memory required is much larger.

## Workflow

### 1. Compile the FDTD solver

Choose one of three build options depending on the resources available to you:

- **Single-threaded CPU**
- **Multi-threaded CPU** (OpenMP)
- **GPU** (OpenACC)

Batch scripts are included for building each version (`compile.sh`, `compile_mp.sh`, and `compile_acc.sh`).

For field viewing, it is often desirable to run shorter simulations, even if longer simulation times are required for accurate post-processing. Storing all EM field values is costly in both memory and time, and viewing fields over shorter time frames can still give a good understanding of the physics. Because of this, there is no substantial advantage in choosing the higher-performance options. Note that OpenMP and OpenACC are not necessarily faster for every case, depending on how the EM field memory is distributed and copied back to the host CPU for writing to binary files.

### 2. Run the object case — `master.py`

Configures and runs the simulation with the cube present.

- Generates a text file of inputs, then executes the binary compiled in Step 1.
- You must set the solver name in `master.py` (the default version is selected).
- Produces several output files used for field and geometry viewing.

> **Note:** Run for zero time steps if you want to view the geometry (see Step 4) before running a full simulation.

### 3. (Optional) Run the clear case — `master_clear.py`

Performs the same steps as `master.py`, but without the cube present, to establish the reference (background) fields.

- You must set the solver name in `master_clear.py` (the default version is selected).
- Produces several output files used for field and geometry viewing.

Example Slurm batch scripts are included and can be adapted to your cluster environment.

### 4. (Optional) View the geometry

To visualize the simulation geometry:

1. Run `fdtd_geometry_maker.py` to generate the ParaView files. It reads the FDTD inputs file and the corresponding geometry binary.
2. Open the resulting **single** ParaView file directly in ParaView. It references an accompanying folder of associated files, so leave that folder in place and do not open its contents individually.
3. For an easier setup, load the included macro, `fdtd_macro.py`, into ParaView (**Macros** tab) to automatically configure common viewing filters.

<p align="center">
  <img src="./geometry.png" alt="Model geometry in ParaView" width="50%" />
</p>

### 5. Create 2D field slices with the Python viewer — `python_fields_viewer.py`

The top of the script has many user-configurable options. Both time-domain and frequency-domain (still frame and phasor) views are available, with or without the geometry displayed on the 2D Yee cell slice. Additional settings allow viewing the incident or scattered fields only, if the clear case was also run. There are also many options for changing colors and scaling.

All expected E and H field components are saved as `.mp4` files for reference. An example frame from the time-domain `Ex.mp4` file is shown below.

<p align="center">
  <img src="./Ex.jpg" alt="Time-domain Ex field slice from the Python viewer" width="50%" />
</p>

### 6. Create 3D field files with the ParaView viewer — `paraview_fields_viewer.py`

The top of the script has a few user-configurable options. Both time-domain and frequency-domain outputs are available. There are fewer settings here because ParaView itself offers extensive options for controlling what information is viewed and how.

The most convenient workflow is to open either a new ParaView instance or the geometry file (if it is already open and viewable), then load all time- or frequency-domain files as one grouped file (the default option). From there, you can easily add slices or volume renderings, change colors and scales, compare magnitude versus individual components, and more. See the ParaView documentation for best practices. Examples are shown below.

<p align="center">
  <img src="./example2.png" alt="Time-domain 3D fields in ParaView (view 1)" width="49%" />
  <img src="./example2-2.png" alt="Time-domain 3D fields in ParaView (view 2)" width="49%" />
</p>

> **Caution:** ParaView often applies unexpected interpolation schemes or creates grids that look uniform but are actually nonuniform. For example, symmetric fields can appear asymmetric if you are not careful. ParaView is a very powerful tool, but it does a lot under the hood without explicitly stating what it is doing.

## Performance reference

Simulation times are on the order of 10–30 seconds per simulation or script. The full 3D fields total approximately 3 GB for all six components. The 2D field slices are substantially smaller.

> **Note:** These timings depend heavily on hardware, problem size, and system load. Use them only as a rough point of reference, not as a direct benchmark against other software.
