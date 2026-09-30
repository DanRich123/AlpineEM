# FDTD Solver Build & Compilation Scripts

Build utilities and Bash compilation scripts for compiling the FDTD solver across various hardware target architectures (CPU OpenMP, NVIDIA GPU OpenACC) and solver variants (Standard, KMAX periodicity, NGSPICE circuit co-simulation). The specifics might vary for a given user. Example: Loading Autotools was required on CU's Alpine HPC system user nodes, but was not required on my personal Linux machine.

---

## Files Overview

### NGSPICE Library Installation
| Script | Description | Environment / Compiler Dependencies | Output |
|---|---|---|---|
| `build_my_ngspice.sh` | Downloads, configures, and builds the NGSPICE shared library (`libngspice`) from source with CIDER support enabled for CPU-bound circuit co-simulations. | Intel compiler suite (`intel/2024.2.1`), `autotools` | `libngspice` under `~/fdtd-fortran-ngspice/` |
| `build_my_ngspice_acc.sh` | Builds the shared `libngspice` library configured with `-fPIC` and OpenMP disabled to prevent thread contention during GPU OpenACC co-simulations. | `gcc/11.2.0`, `autotools` | `libngspice` under `~/fdtd-fortran-ngspice-acc/` |

---

### Solver Build Scripts
| Script | Description | Feature Preprocessor Flags | Target Executable |
|---|---|---|---|
| `compile.sh` | Basic serial FDTD solver build. | None (`-fpp`) | `FDTD` |
| `compile_acc.sh` | GPU-accelerated FDTD solver using OpenACC. | `-acc=gpu -Mpreprocess` | `FDTD-ACC` |
| `compile_mp.sh` | Multi-core CPU FDTD solver using OpenMP. | `-qopenmp -fpp` | `FDTD-MP` |
| `compile_kmax.sh` | Serial FDTD solver with KMAX periodicity enabled. | `-Duse_kmax_version` | `FDTD-KMAX` |
| `compile_kmax_acc.sh` | GPU-accelerated FDTD solver with KMAX periodicity. | `-acc=gpu -Duse_kmax_version` | `FDTD-KMAX-ACC` |
| `compile_kmax_mp.sh` | Multi-core OpenMP FDTD solver with KMAX periodicity. | `-qopenmp -Duse_kmax_version` | `FDTD-KMAX-MP` |
| `compile_spice.sh` | Serial FDTD solver linked with NGSPICE C/Fortran interfaces. | `-Duse_spice_version` | `FDTD-SPICE` |
| `compile_spice_acc.sh` | GPU-accelerated OpenACC FDTD solver with NGSPICE co-simulation. | `-acc=gpu -Duse_spice_version` | `FDTD-SPICE-ACC` |
| `compile_spice_mp.sh` | Multi-core OpenMP FDTD solver with NGSPICE co-simulation. | `-qopenmp -Duse_spice_version` | `FDTD-SPICE-MP` |
| `compile_kmax_spice.sh` | Serial FDTD solver combining KMAX periodicity and NGSPICE co-simulation. | `-Duse_kmax_version -Duse_spice_version` | `FDTD-KMAX-SPICE` |
| `compile_kmax_spice_acc.sh` | GPU OpenACC FDTD solver combining KMAX periodicity and NGSPICE. | `-acc=gpu -Duse_kmax_version -Duse_spice_version` | `FDTD-KMAX-SPICE-ACC` |
| `compile_kmax_spice_mp.sh` | Multi-core OpenMP FDTD solver combining KMAX periodicity and NGSPICE. | `-qopenmp -Duse_kmax_version -Duse_spice_version` | `FDTD-KMAX-SPICE-MP` |

---

## Technical Notes

* **Toolchain Requirements:**
  * **Intel Workloads:** Uses `ifx` (Fortran compiler) and `icx` (C compiler for NGSPICE wrapper interface).
  * **NVIDIA GPU Workloads:** Uses `nvfortran` and `nvc` from the NVIDIA HPC SDK (`nvhpc_sdk`) alongside CUDA modules.
* **SPICE Pre-compilation Dependency:** Any `*_spice*.sh` build requires building the `ngspice` shared library target first using either `build_my_ngspice.sh` (for CPU/MP execution) or `build_my_ngspice_acc.sh` (for GPU ACC execution).
