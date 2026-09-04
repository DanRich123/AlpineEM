# AlpineEM

**An FDTD–SPICE electromagnetic simulation suite**, coupling a full-wave 3D Finite-Difference Time-Domain (FDTD) solver with SPICE circuit co-simulation via ngspice. Built for a wide variety of electromagnetics problems.

Developed by Daniel Richardson at the Center for National Security Initiatives (NSI), University of Colorado Boulder and as an independent.

> **Status:** Active development. APIs, file formats, and folder structure may still change — see [`ITEMS TO ADD.txt`](./ITEMS%20TO%20ADD.txt) for the current roadmap.

---

## Features

- **3D FDTD solver** on a standard cubic Yee grid, written in Fortran (`.f90`), with OpenMP / OpenACC build variants for multithreaded CPU or GPU execution
- **SPICE co-simulation** — couple FDTD lumped ports directly to an ngspice circuit netlist
- **Convolutional PML** boundary option (ADE formulation)
- **Periodic boundary support** via a standard periodic boundary option and a dedicated `kmax` solver variant for oblique-incidence plane waves (constant k-vector method)
- **Infinite ground plane** boundary option
- **Two excitation types** — Gaussian and normalized differentiated Gaussian — usable as lumped-port or TF/SF plane-wave sources
- **Adaptive on-the-fly** time-domain far-field (at select angles) and S-parameter extraction, both producing broadband information
- **Sub-cell thin-sheet modeling** (Smith–Mahoney method), including zero-impedance approximations for PEC sheets
- **Dispersive media support** via an auxiliary differential equation (ADE) approach — currently only Drude (plasma) media
- **Anisotropic media** — diagonal elements for both permittivity and electrical conductivity
- **Vacuum permeability** — only vacuum permeability is currently supported
- **Statics solver** for generating non-dispersive E/H field patterns used in gridded lumped ports (TEM-mode focused)
- **Optimization examples** — adjoint (gradient-based) and genetic-algorithm workflows
- **Paraview integration** for geometry visualization, including a ready-to-import macro
- **Slurm submission scripts** for HPC/cluster runs (especially useful for OpenMP/OpenACC builds)
- An experimental, incomplete **pure-Python (PyTorch/TensorFlow) trial version** of the solver

## Repository structure

| Path | Description |
|---|---|
| [`main_fdtd/`](./main_fdtd) | Core FDTD solver (`.f90`) and Python execution script; includes SPICE-linked build files |
| [`statics_solver/`](./statics_solver) | Statics solver for generating TEM-mode E/H patterns for gridded lumped ports |
| [`compile_scripts/`](./compile_scripts) | Scripts for compiling the FDTD solver and (optionally) building SPICE |
| [`slurm_scripts/`](./slurm_scripts) | Slurm job submission scripts, primarily for OpenMP/OpenACC runs |
| [`optimization_scripts/`](./optimization_scripts) | Example adjoint (gradient) and genetic-algorithm optimization workflows |
| [`utility_scripts/`](./utility_scripts) | Post-processing utilities, including field plotting |
| [`paraview/`](./paraview) | Geometry-viewing macro for Paraview |
| [`fdtd_python_trial_versions/`](./fdtd_python_trial_versions) | Incomplete, lightly-tested pure-Python (PyTorch/TensorFlow) version of the solver for educational purposes |
| [`examples/`](./examples) | Example simulation files, updated periodically (may occasionally lag behind the current version) |
| [`ITEMS TO ADD.txt`](./ITEMS%20TO%20ADD.txt) | Planned fixes and features |
| [`LICENSE.txt`](./LICENSE.txt) | License and third-party acknowledgments |

## Solver variants

The FDTD core is built in four flavors, combined from two axes:

| | Standard incidence | Oblique incidence (periodic, `kmax`) |
|---|---|---|
| **No SPICE** | ✅ | ✅ (1 source type only) |
| **With SPICE** | ✅ | ✅ (1 source type only) |

`kmax` refers to a solver variant specifically for oblique-angle plane waves incident on periodic boundary conditions.

## Requirements

- Fortran compiler with OpenMP / OpenACC support (e.g. `ifx/ifort`, `gfortran`, `nvfortran`)
- [ngspice](https://ngspice.sourceforge.io/) (only required for SPICE-coupled builds)
- C compiler with OpenMP / OpenACC support (only required for SPICE-coupled builds, e.g. `icx/icc`, `gcc`, `nvc`)
- Python 3 (recommended) for the execution wrapper script and post-processing utilities
- [Paraview](https://www.paraview.org/) (optional, for geometry visualization)
- Slurm (optional, only needed for cluster job submission)

## Building

Compile scripts for each solver variant live in [`compile_scripts/`](./compile_scripts). At a high level:

```bash
# Example — replace with your target variant's script
./compile.sh        # standard build
./compile_spice.sh  # SPICE-coupled build (requires ngspice)
```

The compilation process is generally straightforward since there are few dependencies (other than ngspice), so it's easy to modify if needed.

There are also two bash scripts for building ngspice, if the user has not already done so.

## Quick start

1. Create a new folder and choose a solver variant (standard vs. `kmax`, with or without SPICE) to compile using the appropriate script in [`compile_scripts/`](./compile_scripts).
2. Copy any needed files for your selected variant from [`main_fdtd/`](./main_fdtd) to this new folder (e.g. `fdtd_solver.f90`, `circuit.f90`, etc.) and compile the solver.
3. Copy `master.py` from [`main_fdtd/`](./main_fdtd) and any utility files needed from [`utility_scripts/`](./utility_scripts), such as `post_processor.py`.
4. Set up your geometry and source parameters in `master.py` (see [`examples/`](./examples) for reference input files).
5. If using gridded lumped ports, run the [`statics_solver/`](./statics_solver) first to generate the TEM-mode E/H field pattern.
6. Run the solver via the Python execution script `master.py` that was copied from [`main_fdtd/`](./main_fdtd) (locally, or submit via [`slurm_scripts/`](./slurm_scripts) on a cluster).
7. You may need to run a "clear" case depending on what you are simulating (e.g. a clear case is the identical simulation, typically with no geometry present).
8. Post-process results with the tools in [`utility_scripts/`](./utility_scripts), and view geometry in Paraview using the macro in [`paraview/`](./paraview).

## Additional information and known limitations

- **CFL-reducible time step:** Yee-grid FDTD on a standard rectangular coordinate system with the standard CFL condition, but it's user-reducible.
- **Thin-sheet filtering:** Vacuum sheets are intentionally skipped (used as the existence-filter mechanism); PEC sheets are approximated as high-conductivity with finite thickness rather than zero-impedance.
- **Lumped port information:** Internal (non-SPICE) ports intentionally omit the FDTD cell capacitance; SPICE-linked ports include it. SPICE FDTD locations must use non-dispersive permittivity. Gridded lumped ports accept arbitrary geometry but require known non-dispersive E/H field coefficients — currently generated via the statics solver for TEM modes.
- **Adaptive far-field time-domain calculations** do not parallelize efficiently via OpenMP and OpenACC due to atomics; a fully parallel version would be RAM-intensive, though. A future release will include a broadband-angle-at-select-frequencies option; this is expected to parallelize well with no issues.
- **Non-TEM mode generation** for lumped or wave ports requires a 2D Helmholtz-equation solver, which does not yet exist, along with modifications to the main code to support dispersive port behavior.
- **Non-LTI (narrow-band/CW) source excitation** is not currently available but planned.
- **Pure-Python trial version** is incomplete and not well tested — this is primarily for teaching purposes.

See [`ITEMS TO ADD.txt`](./ITEMS%20TO%20ADD.txt) for the full list of planned additions, including further source types, a Helmholtz solver, and dispersive lumped-port coefficients.

## Citing this work

If you use AlpineEM in academic work, please cite it, e.g.:

```bibtex
@software{alpineem,
  author  = {Richardson, Daniel},
  title   = {AlpineEM: An FDTD-SPICE Software Suite},
  year    = {2026},
  url     = {https://github.com/DanRich123/AlpineEM}
}
```

## Acknowledgments

Portions of this software are derived from open-source work by others:

- FDTD/CPML implementation derived from code by Jamesina J. Simpson (University of New Mexico) — [CPML implementation](https://github.com/cvarin/FDTD/blob/master/Taflove/fdtd3D_CPML.f90)
- ngspice Fortran interface derived from code by Alberto Gascon ([OpenSEMBA/ngspice_fortran_interface](https://github.com/OpenSEMBA/ngspice_fortran_interface))

Full attribution is in [`LICENSE.txt`](./LICENSE.txt).

## License

MIT License — © 2025–2026 The Regents of the University of Colorado, a body corporate; © 2025–2026 Daniel Richardson. See [`LICENSE.txt`](./LICENSE.txt) for full terms.

