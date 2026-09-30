# utility_scripts

Pre- and post-processing tools for the main FDTD solver. See [`examples/`](../examples) for demonstration workflows and [`main_fdtd/`](../main_fdtd) for the solver itself.

## How these fit together

```
geometry helpers ──► master.py ──► solver ──► data.dat (+ clear.dat, metal.dat)
(optional)                            │              │
                                      │              └─► post_processor*.py ──► .csv results
                                      │                                          └─► plot / touchstone scripts
                                      └─► *_video_*.bin / *_full_*.bin ──► field viewers
```

Conventions shared by most scripts:

- **Copy the script into your simulation folder** and run it there. Most scripts read and write files in the current directory.
- **Edit the setup block at the top**, not the body. Setup blocks are marked with comment banners.
- The post-processors and viewers read header lines of the solver's `.dat` file by position, so use them with output from the matching solver version.
- Several scripts are deliberately project-specific (marked below). Treat them as templates to adapt.

## Standard workflow scripts

| Script | Description | Input | Output |
|---|---|---|---|
| `post_processor.py` | Computes S-parameters, scattering far field, and realized antenna gain for the standard (non-`kmax`) solver. Supports plane-wave and antenna excitation, optional SPICE-port excitation, and optional metal / infinite-ground-plane reference cases | `data.dat`; plus `clear.dat`, `metal.dat`, `IGP_full_clear.dat` as needed | `S_parameters.csv`, `Scattering_Far_Field.csv`, `Realized_antenna_gain.csv` (which are written depends on the simulation; names are configurable) |
| `post_processor_kmax.py` | Same outputs for the `kmax` solver. Has a `dB_rule` setting that trades trustworthy bandwidth for qualitative data | `data-k.dat`, `clear-k.dat`, `metal-k.dat` | `S_parameters_k.csv`, `Scattering_Far_Field_k.csv`, `Realized_antenna_gain_k.csv` |
| `python_fields_viewer.py` | Animations (time domain) or still/phase-sweep images (frequency domain) of the **2D slice** field export. Can show total, incident, or scattered fields and overlay the geometry | `<name>_video_*.bin`, `<name>_geometry.bin`, `data.dat` (and `clear.dat` for incident/scattered) | One `.mp4` per field (or `.png` for frequency-domain stills) |
| `paraview_fields_viewer.py` | Converts the **full-volume** field export to VTK image files, in time or frequency (magnitude) | `<name>_full_*.bin` | `Output/E_Field_###.vti`, `Output/H_Field_###.vti` |
| `make_optional_geom_bulk.py` | Minimal example of drawing a bulk material array directly in NumPy | none | `optional_geom_bulk.bin` |
| `make_optional_geom_and_conform_bulk.py` | Draws geometry on a fine grid and conformally averages it down to the FDTD grid | `conformal_builder.py` | `optional_geom_bulk.bin`, `materials_id_opfile.npy` |
| `conformal_builder.py` | The `ConformalGeometry` class used above: harmonic/arithmetic averaging of ε and σ, with new blended materials created as needed. No Drude or permeability support yet | NumPy (Matplotlib for `plot_slice`) | Imported, not run directly |

### Notes on the workflow scripts

**Post-processors.** `padding` zero-pads the time series to interpolate in frequency. It can cause ripples if the fields have not converged, so check the time series first. Port impedances for SPICE ports default to 50 Ω, and the CSV headers say so; edit the script if yours differ.

**Field viewers.** Two different exports feed two different tools: the *slice* export (`video`) goes to `python_fields_viewer.py`, the *full-volume* export (`full`) goes to `paraview_fields_viewer.py`. `python_fields_viewer.py` needs `ffmpeg` on your path for video output. In `paraview_fields_viewer.py`, enter the grid size, time step, and step count by hand in the setup block (it will fail if they don't match the files), and set `kmax=True` for complex-valued `kmax` data. In Paraview, open the `Output/` series as a file group, and be aware that Paraview's interpolation can make symmetric fields look asymmetric on slices.

**Optional bulk geometry files.** These are flat `float32` arrays in Fortran (column-major) order. The solver expects exactly this format. Zero means "ignored", and material IDs must be registered in `master.py`. Python indexing is 0-based while `master.py` locations are 1-based; `make_optional_geom_bulk.py` shows the conversion. For the conformal version, the fine grid must be `conform_num` times the FDTD grid in each direction, and the `.npy` table (columns: id, εx, εy, εz, σx, σy, σz) makes registering the materials, including any newly blended ones, easier.

## Specialized and non-routine scripts

| Script | Description | Input / dependencies | Output |
|---|---|---|---|
| `import_tecplot_and_interpolate.py` | Imports a 2D CFD flowfield in a specific Tecplot layout, interpolates it onto a 3D cylindrical-symmetry grid, and computes plasma parameters (ωp, γ) per species. Options for a quarter model, radome layers, and material binning | Tecplot `.plt` file, `scipy`, `pandas`, `matplotlib` | `optional_geom_bulk.bin`, `material_properties.npy`, diagnostic PNGs |
| `import_plasma_properties.py` | **Code snippet, not a standalone script.** Paste it into the materials section of `master.py` (it writes to `master.py`'s open file `f`). Turns `material_properties.npy` into 6-pole Drude materials, then appends radome, metal, sheet, and air materials (IDs numbered after the plasmas) | `material_properties.npy` | Material lines in the solver input file |
| `convert_to_touchstone.py` | **Custom script.** Builds a 4-port `.s4p` (MA format) from two lumped-port and two wave-port (TE/TM) S-parameter CSVs, using angle-dependent wave-port impedances and renormalizing to 50 Ω. Hard-coded for that port layout and x-normal incidence | `S_parameters_TE.csv`, `S_parameters_TM.csv`, `S_parameters_antenna1.csv`, `S_parameters_antenna2.csv`, `data.dat`, `scikit-rf` | `my_sparams_MA.s4p` (name configurable) |
| `parallel_kmax.py` | Runs **one** point of a `kmax` sweep per invocation: `python parallel_kmax.py <index>`. Makes `working/sim_<ky>_<kz>/`, copies the current folder into it, then runs `master.py`, `master_clear.py`, and `post_processor_kmax.py` with that ky, kz. Launch many indices in parallel with a Slurm job array or a shell loop | `master.py` and `master_clear.py` that accept `ky kz` as command-line arguments, solver binary, `post_processor_kmax.py` | One sub-folder per sweep point |
| `plot_kmax_f_vs_angle.py` | Frequency vs. angle (and vs. kz) heat maps of TM reflection/transmission from a sweep | `working/sim_*/S_parameters_k.csv` | `refl_f_vs_angle.png`, `refl_f_vs_kz.png`, `trans_f_vs_angle.png`, `trans_f_vs_kz.png` |
| `plot_kmax_angle_vs_angle.py` | θ vs. φ heat maps of TM reflection/transmission at one frequency, interpolated from a sweep | `working/sim_*/S_parameters_k.csv` | `reflection_transmission_<freq>GHz_angle_vs_angle.png` |
| `spiral.py` | Generates a discretized two-arm Archimedean spiral as a 0/1 mask, e.g., for placing planar-antenna geometry | NumPy, Matplotlib | `mask.npy`, `spiral.png` |

### Notes on the specialized scripts

- The two `plot_kmax_*` scripts must be edited so their ky/kz arrays match the ones in `parallel_kmax.py`, since they locate results by folder name. The plotted CSV columns (TM reflection and transmission) and the angle formula are hard-coded; adjust them for your geometry and port ordering.
- `convert_to_touchstone.py` and the Tecplot scripts were written for specific projects. They are included as working references for similar tasks.

## Dependencies

| Package / tool | Needed by |
|---|---|
| `numpy`, `matplotlib` | Nearly everything |
| `scipy` | `paraview_fields_viewer.py`, `import_tecplot_and_interpolate.py` |
| `pandas` | `import_tecplot_and_interpolate.py` |
| `pyevtk` | `paraview_fields_viewer.py` |
| `scikit-rf` | `convert_to_touchstone.py` |
| `ffmpeg` | `python_fields_viewer.py` (video output) |
