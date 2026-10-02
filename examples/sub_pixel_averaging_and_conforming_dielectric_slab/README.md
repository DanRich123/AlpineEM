# In progress...

# Sub-pixel Averaging and Conforming Media: Transmission and Reflection Through a Quasi-1D Infinite Slab of Dielectric Material (AlpineEM FDTD)

This example uses AlpineEM to perform FDTD simulations that calculate the transmission and reflection through a quasi-1D infinite slab of dielectric material. The slab is infinite in the y and z directions and finite in the x direction. The infinite condition is created using periodic boundary conditions.

Using the `conformal_builder.py` and `make_optional_geom_and_conform_bulk.py` scripts, objects can be drawn into a much larger (more cells) and more refined (smaller cells) Yee cell grid. Through anisotropic harmonic and arithmetic sub-pixel averaging schemes, the materials can then be conformed to fit the actual Yee cell grid of the intended simulation. Dispersive media properties are not directly supported yet, and this method is most accurate for low-loss dielectric materials. Users will notice minimal improvement when conforming high-conductivity materials.

This example considers and compares two cases:

1. The dielectric slab fits neatly within the Yee cell grid.
2. The dielectric slab is shifted by 1/2 cell in space, so it does not fit neatly within the Yee cell grid and must be conformed through sub-pixel averaging.

It outputs time-domain data for post-processing, along with binary geometry files for the **object** (dielectric present, both cases above), **clear** (dielectric absent), and **metal** (dielectric replaced with metal) cases.

This example is designed to teach and validate the core mechanisms of the solver, including normal-incidence plane waves with periodic boundary conditions, dielectric materials, sub-pixel and conformal building, and the three builder options (single-threaded CPU, multi-threaded CPU (OpenMP), and GPU (OpenACC)).

## Workflow

### 1. Compile the FDTD solver

Choose one of three build options depending on the resources available to you:

- **Single-threaded CPU**
- **Multi-threaded CPU** (OpenMP)
- **GPU** (OpenACC)

Batch scripts are included for building each version (`compile.sh`, `compile_mp.sh`, and `compile_acc.sh`). Example Slurm batch scripts are also included for running steps 3-6 and can be adapted to your cluster environment.

### 2. Perform sub-pixel averaging and conform the media — `make_optional_geom_and_conform_bulk.py`

This custom script uses the class found in `conformal_builder.py`. It performs harmonic and arithmetic averaging of the dielectric properties, as discussed above.

It outputs two files, `materials_id_opfile.npy` and `optional_geom_bulk.bin`. `materials_id_opfile.npy` is a NumPy array of material property information that can be loaded and used in `master.py`. `optional_geom_bulk.bin` is a geometry binary that `master.py` reads in. Anything drawn in `master.py` after the import overwrites the corresponding regions of the imported geometry.

### 3. Run the object case — `master.py` or `master_perfect_alignment.py`

Configures and runs the simulation with the object present.

- Generates a text file of inputs, then executes the binary compiled in Step 1.
- You must set the solver name in `master.py` (the single-threaded version is selected by default).
- Produces several output files used for post-processing and geometry viewing.

  > **Note:** To view the geometry (see Step 7) before running a full simulation, run for zero time steps.

### 4. Run the clear case — `master_clear.py`

Performs the same steps as `master.py`, but without the object present, to establish the reference (background) fields.

- You must set the solver name in `master_clear.py` (the single-threaded version is selected by default).
- Produces several output files used for post-processing and geometry viewing.

### 5. Run the metal case — `master_metal.py`

Performs the same steps as `master.py`, but with the object replaced by metal, to provide the last piece needed for TRL calibration.

- You must set the solver name in `master_metal.py` (the single-threaded version is selected by default).
- Produces several output files used for post-processing and geometry viewing.
- To bypass the metal case, set `use_metal=False` in `post_processor.py`.
- Using the metal case provides more accurate results, but it shifts the phase center away from the wave port to the surface of the metal.

### 6. Post-process — `post_processor.py`

Combines the object, clear, and metal case outputs to generate transmission and reflection data as a `.csv` file.

### 7. (Optional) View the geometry

To visualize the simulation geometry:

1. Run `fdtd_geometry_maker.py` to generate ParaView files. It reads the FDTD inputs file and the corresponding geometry binary to create them.
2. Open the resulting **single** ParaView file directly in ParaView. It references an accompanying folder of associated files, so leave that folder in place and do not open its contents individually.
3. For an easier setup, load the included macro, `fdtd_macro.py`, into ParaView (**Macros** tab) to automatically configure common viewing filters.

The ParaView-rendered image is shown below:

<p align="center">
  <img src="./geometry.png" alt="Model geom" width="50%" />
</p>

### 8. Validation

The `plot.py` file can be used to compare the two FDTD results. There is excellent agreement between the amplitude and phase of both transmission and reflection.

Because the conformally averaged case is shifted 1/2 cell relative to the perfectly aligned case, a 1-cell phase correction must be applied manually to the reflection phase to compare apples to apples. The transmission phase is unaffected by this placement shift, as expected.

At high frequencies there is a growing error, most notable in the transmission plot. This is expected: sub-pixel averaging is not perfect, and failures to conform show up first at the higher frequencies.

![Model Accuracy](./Comparison_plot.png)

## Performance reference

Approximate per-simulation runtimes measured on the author's hardware:

| Solver                      | Time per simulation |
|-----------------------------|---------------------|
| OpenACC (GPU)               | ~3.8 seconds        |
| OpenMP (multi-threaded CPU) | ~0.9 seconds        |
| Single-threaded (default)   | ~1.8 seconds        |

> **Note:** OpenACC (GPU) performed worse here because of the size and scaling of the problem. Unless the simulation is large enough, the overhead of GPU calculations dominates. This problem uses a severely non-cubic grid (75x30x30) and is small in both total cells and time steps, so the GPU did not come out ahead. See the monostatic scattering of a sphere example, where the GPU handily dominates. These timings depend heavily on hardware, problem size, and system load; use them only as a rough point of reference, not as a direct benchmark against other software.
