# In progress...
# Oblique Angle TM Transmission and Reflection Through a Quasi-1D Infinite Slab of Dielectric Material (AlpineEM FDTD)

This example uses AlpineEM to perform FDTD simulations that calculate the oblique angle transmission and reflection through a quasi-1D infinite slab of material. The slab of material is infinite in the y and z directions, and finite in the x direction. The infinite condition is created using periodic boundary conditions. Because the angles are oblique with respect to the periodic boundaries, the `kmax` variant is required.

It outputs time-domain data for post-processing, along with binary geometry files for the **object** (material present), **clear** (material absent), and **metal** (material replaced with metal) cases.

This example is designed to teach and validate the core mechanisms of the solver including oblique incidence plane waves with periodic boundary conditions, dielectric materials, and the utiltity scripts associated with processing `kmax` data.

This FDTD variant utilizes what is commonly referred to as the constant k vector method. Each frequency component of the transmission and reflection extracted from the time domain simulations corresponds to a different angle of incidence, requiring additional steps after post processing to understand and visualize the data. Similarly other S-parameter and far field quantities (even w/ Bloch phase removed) can be hard to meaningfully interpret so great care must be taken when using the `kmax` variant.

## Workflow

### 1. Compile the FDTD solver and Prep the Batch Script

Choose one of three build options depending on the resources available to you:

- **Single-threaded CPU**
- **Multi-threaded CPU** (OpenMP)
- **GPU** (OpenACC)

There are batch scripts included for building all versions (`compile_kmax.sh`, `compile_kmax_mp.sh`, and `compile_kmax_acc.sh`).

The example batch script (`run_parallel.sh`) included here was designed for running many simulations in parallel on an HPC system. Due to HPC resource availability, only the single-threaded option is demonstrated - each simulation is small, and 1 core is sufficient for each. The other versions can be utilized here, if resources allow.

There are two examples demonstrated here. The only difference between the examples is the number of cases, which k vectors are chosen for each case, and then subsequently which plotting script should be used to understand the results. To perform either example, only the `run_parallel.sh` needs to be modified. Example 2 is currently set as default. In example 1, 400 cases (3 sims per case - object, metal, clear) were run within the `run_parallel.sh` script where the line for `parallel_kmax_f_vs_angle.py` can be uncommented and `parallel_kmax_angle_vs_angle.py` should then be commented out. Likewise, to run example 2 where 225 cases were simulated, the number of simulations can be set to 225 and the `parallel_kmax_angle_vs_angle.py` can be uncommented and `parallel_kmax_f_vs_angle.py` can be commented out.

The batch script will run all simulations but the breakdown, like in other examples, is recorded here.

### 2. Run the object case — `master.py`

Configures and runs the simulation with the object present.

- Generates a text file of inputs, then executes the binary compiled in Step 1.
- You must set the solver name in `master.py` (the single-threaded version is selected by default).
- Produces several output files used for post-processing and geometry viewing.

  > **Note** Run for zero time steps if viewing the geometry (see step 6) is desired before running a full simulation.

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

### 5. Post-process — `post_processor_kmax.py`

Combines the object, clear, and metal case outputs to generate transmission and reflection data as a `.csv` file.

### 6. (Optional) View the geometry

To visualize the simulation geometry (same for each k vector combination chosen):

1. Run `fdtd_geometry_maker.py` to generate ParaView files.
2. Open the resulting **single** ParaView file directly in ParaView — it references an accompanying folder of associated files, so leave that folder in place and do not open its contents individually.
3. For an easier setup, load the included macro, `fdtd_macro.py`, into ParaView (**Macros** tab) to automatically configure common viewing filters.

The Paraview rendered image is shown below:

<p align="center">
  <img src="./geometry.png" alt="Model geom" width="50%" />
</p>

### 7a. Validation - example 1

The `plot_kmax_f_vs_angle.py` file can be used to compare the FDTD result with the analytic result.

<p align="center">
  <img src="./ex1.png" alt="ex1" width="100%" />
</p>

If you want to go to higher angles for the same frequency you need to lower the pulse parameter frequency and then go to higher k values, or move beyond the 20 dB rule in the post processor. Moving beyond the 20 dB rule means data produced and saved beyond that limit is more qualitative and less trustworthy. However, you can get more data this way and it's sometime not that bad of a result.

### 7b. Validation - example 2

The `plot_kmax_angle_vs_angle.py` file can be used to plot the results.

<p align="center">
  <img src="./ex2.png" alt="ex2" width="100%" />
</p>

Though the definition defined angles for ky,kz outside the data shown, the frequency parameter using the 20dB rule doesn't allow for those angle data to be trusted so it is not saved in the csv file - hence the white space in the plots. As discussed above in 7a., it that information is desired, use beyond than the 20 dB rule in the post processor or lower the frequency parameter with higher k values.
