# AlpineEM

**An FDTD–SPICE electromagnetic simulation suite**, coupling a full-wave 3D Finite-Difference Time-Domain (FDTD) solver with SPICE circuit co-simulation via ngspice. Built for antenna, periodic-structure, and lumped-port circuit problems that need both field-level accuracy and circuit-level fidelity.

![License](https://img.shields.io/github/license/DanRich123/AlpineEM)
![Language](https://img.shields.io/github/languages/top/DanRich123/AlpineEM)
![Last commit](https://img.shields.io/github/last-commit/DanRich123/AlpineEM)
![Status](https://img.shields.io/badge/status-active%20development-yellow)

Developed by Daniel Richardson at the Center for National Security Initiatives (NSI), University of Colorado Boulder.

> **Status:** Active development. APIs, file formats, and folder structure may still change — see [`ITEMS TO ADD.txt`](./ITEMS%20TO%20ADD.txt) for the current roadmap.

---

## Features

- **3D FDTD solver** on a standard Yee grid, written in Fortran (`.f90`), with OpenMP / OpenACC build variants for multithreaded CPU or GPU execution
- **SPICE co-simulation** — couple FDTD lumped ports directly to an ngspice circuit netlist
- **Convolutional PML** boundaries (ADE formulation)
- **Periodic boundary support** via a dedicated `kmax` solver variant for oblique-incidence plane waves
- **Two excitation types** — Gaussian and normalized differentiated Gaussian — usable as lumped-port or TF/SF plane-wave sources
- **Adaptive on-the-fly** time-domain far-field and S-parameter extraction (broadband, at selected angles)
- **Sub-cell thin-sheet modeling** (Smith–Mahoney method), including high-conductivity approximations for PEC sheets
- **Dispersive media support** via an auxiliary differential equation (ADE) approach — currently Drude (plasma) media
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
| [`fdtd_python_trial_versions/`](./fdtd_python_trial_versions) | Incomplete, lightly-tested pure-Python (PyTorch/TensorFlow) version of the solver |
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

- Fortran compiler with OpenMP / OpenACC support (e.g. `gfortran`, `nvfortran`)
- [ngspice](https://ngspice.sourceforge.io/) (only required for SPICE-coupled builds)
- Python 3 (for the execution wrapper script and post-processing utilities)
- [Paraview](https://www.paraview.org/) (optional, for geometry visualization)
- Slurm (optional, only needed for cluster job submission)

## Building

Compile scripts for each solver variant live in [`compile_scripts/`](./compile_scripts). At a high level:

```bash
# Example — replace with your target variant's script
cd compile_scripts
./compile_fdtd.sh        # standard build
./compile_fdtd_spice.sh  # SPICE-coupled build (requires ngspice)
```

> Update this section with the exact script names/flags from `compile_scripts/` — see note above.

## Quick start

1. Choose a solver variant (standard vs. `kmax`, with or without SPICE) and compile it using the appropriate script in `compile_scripts/`.
2. Set up your geometry and source parameters (see [`examples/`](./examples) for reference input files).
3. If using gridded lumped ports, run the [`statics_solver/`](./statics_solver) first to generate the TEM-mode E/H field pattern.
4. Run the solver via the Python execution script in [`main_fdtd/`](./main_fdtd) (locally, or submit via [`slurm_scripts/`](./slurm_scripts) on a cluster).
5. Post-process results with the tools in [`utility_scripts/`](./utility_scripts), and view geometry in Paraview using the macro in [`paraview/`](./paraview).

## Methods

- **Field update:** Yee-grid FDTD on a standard rectangular coordinate system, standard CFL condition (user-reducible)
- **Boundaries:** Convolutional PML (ADE approach)
- **Dispersive media:** ADE approach, currently limited to Drude (plasma) media
- **Thin sheets:** Smith–Mahoney sub-cell method. Vacuum sheets are intentionally skipped (used as the existence-filter mechanism); PEC sheets are approximated as high-conductivity with finite thickness rather than zero-impedance
- **Lumped ports:** Internal (non-SPICE) ports intentionally omit the FDTD cell capacitance; SPICE-linked ports include it. SPICE FDTD locations must use non-dispersive permittivity. Gridded lumped ports accept arbitrary geometry but require known non-dispersive E/H field coefficients — currently generated via the statics solver for TEM modes
- **Far field / S-parameters:** Adaptive on-the-fly time-domain method for broadband data at selected angles (note: computationally harder to parallelize well due to atomics on both CPU/OpenMP and GPU/OpenACC)

### Known limitations

- Non-TEM mode generation requires a Helmholtz-equation solver, which does not yet exist
- Non-LTI (narrow-band/CW) source excitation is not currently available
- The pure-Python trial version is incomplete and not well tested
- Far-field and adaptive time-domain calculations do not parallelize efficiently; a fully parallel version would be RAM-intensive

See [`ITEMS TO ADD.txt`](./ITEMS%20TO%20ADD.txt) for the full list of planned additions, including further source types, a Helmholtz solver, and dispersive lumped-port coefficients.

## Citing this work

If you use AlpineEM in academic work, please cite it, e.g.:

```bibtex
@software{alpineem,
  author  = {Richardson, Daniel},
  title   = {AlpineEM: An FDTD-SPICE Software Suite},
  year    = {2026},
  url     = {https://github.com/DanRich123/AlpineEM},
  note    = {Center for National Security Initiatives, University of Colorado Boulder}
}
```

> If you archive a release on [Zenodo](https://zenodo.org/) (free, and gives you a DOI), swap the `url` field for a versioned DOI so citations point to the exact release someone used — Zenodo can auto-generate the BibTeX entry for you.

## Acknowledgments

Portions of this software are derived from open-source work by others:

- FDTD/CPML implementation derived from code by Jamesina J. Simpson (University of New Mexico)
- ngspice Fortran interface derived from code by Alberto Gascon ([OpenSEMBA/ngspice_fortran_interface](https://github.com/OpenSEMBA/ngspice_fortran_interface))

Full attribution is in [`LICENSE.txt`](./LICENSE.txt).

## License

MIT License — © 2025–2026 The Regents of the University of Colorado, a body corporate; © 2025–2026 Daniel Richardson. See [`LICENSE.txt`](./LICENSE.txt) for full terms.
