# Statics Solver & Gridded Feed Generator

2-D finite-difference electrostatic and magnetostatic Poisson solver designed to compute static field profiles for initializing wave sources in the main FDTD solver. Currently, this workflow supports **TEM modes only**.

The solver extracts $E$-field and $H$-field weightings directly mapped to Yee-cell grids, exports binary feed distribution files (`gridded_feed.bin`), and generates visualization field plots. See [`examples/`](../examples) for demonstration workflows.

---

## Files Overview

| Script / File | Description | Input / Dependencies | Output |
|---|---|---|---|
| `EM2Dsolver.py` | Core class (`EM2DSolver`) solving 2-D Laplace/Poisson equations via sparse linear systems (`scipy.sparse.linalg.spsolve`). Computes electrostatic potential $V$, vector potential $A_z$, edge $E/H$ fields, and integrated surface currents. | `numpy`, `scipy` | Object class dependency |
| `coax_example.py` | Example script for circular coaxial geometries. Computes $E$ and $H$ static fields, symmetry checks, characteristic impedance ($Z_0$), per-unit-length $L$ and $C$, and exports binary grid feeds. | `EM2Dsolver.py`, `matplotlib` | `gridded_feed.bin`, `electric fields.png`, `magnetic fields.png` |
| `square_coax_example.py` | Example script setup specifically for square-profile coaxial cross-sections. Computes mode profiles, per-unit-length params ($L, C, Z_0$), and exports binary feeds. | `EM2Dsolver.py`, `matplotlib` | `gridded_feed.bin`, `electric fields.png`, `magnetic fields.png` |

---

## Technical Notes

* **Formulation:** Electrostatics ($\nabla \cdot (\epsilon_r \nabla V) = 0$) and Magnetostatics ($\nabla \cdot (\mu_r^{-1} \nabla A_z) = 0$) share identical differential operator mechanics. $V$ and $A_z$ are defined at grid nodes, while material parameters ($\epsilon_r, \mu_r^{-1}$) are defined across cells.
* **Current Normalization:** $H$-fields are normalized by total simulated loop current $I_{\text{total}}$ to produce standardized unit per-ampere field distribution arrays.
* **FDTD Binary Format:** `gridded_feed.bin` is written out as a single flattened float32 binary array storing 4 grid layers ($E_x$, $E_y$, $H_x$, $H_y$) ordered in Fortran memory layout (`order='F'`) matching the primary FDTD solver Yee grid setup.
