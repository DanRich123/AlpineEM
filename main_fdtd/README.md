# main_fdtd


The core AlpineEM solver: a 3D FDTD code in Fortran, optionally coupled to ngspice for circuit co-simulation, plus the Python script that drives it. For the big picture and quick-start steps, see the [top-level README](../README.md). Worked examples are in [`examples/`](../examples).

## Files

| File | Used by | Purpose |
|---|---|---|
| `fdtd_solver.f90` | All variants | The FDTD solver (single program; variants are selected at compile time) |
| `master.py` | All variants | Writes the solver input file, optionally writes a SPICE netlist, then runs the solver. Heavily commented; it doubles as the input-format reference |
| `circuit.F90` | SPICE only | Fortran circuit module: steps ngspice, exchanges port voltages and currents with the FDTD grid |
| `ngspice_interface.F90` | SPICE only | Fortran (`iso_c_binding`) interface to the C wrapper |
| `ngspice_interfaces.c`, `ngspice_interfaces.h` | SPICE only | C wrapper around the ngspice shared-library API |

Non-SPICE builds need only `fdtd_solver.f90`.

## Choosing a variant

The solver is one source file with two preprocessor flags:

| Flag | Effect |
|---|---|
| `use_spice_version` | Enables ngspice co-simulation (needs the four SPICE files above and ngspice) |
| `use_kmax_version` | Enables the oblique-incidence periodic (`kmax`) formulation |

Use neither, either, or both. OpenMP / OpenACC are separate build options. The scripts in [`compile_scripts/`](../compile_scripts) set the flags for you.

## How a run works

1. `master.py` writes a plain-text input file (default `inputs.txt`).
2. It then calls the compiled solver with that file as its only argument: `./<solver> inputs.txt`.
3. The solver writes its results next to the input file.

The input file is **read in order, line by line**, so each line's meaning depends on its position and on the choices made earlier (e.g., excitation type changes which lines follow). Edit the values in `master.py` rather than reordering lines, and use its commented-out blocks as templates for each option.

### Output files

The output name you choose in `master.py` (e.g., `data.dat`) is the base for the rest. Use a three-character extension: the solver strips the last four characters when it builds the other names.

| File | When |
|---|---|
| `<name>.dat` | Always: S-parameter and far-field data (plus port/SPICE data if used) for the post-processor |
| `<name>_geometry.bin` | Always: geometry for the Paraview macro |
| `<name>_video_{Ex,Ey,Ez,Hx,Hy,Hz}.bin` | If you request a 2D field slice |
| `<name>_full_{Ex,Ey,Ez,Hx,Hy,Hz}.bin` | If you request full-volume fields (slower; H is half a time step earlier) |

## Practical guidelines

These are the rules of thumb documented in `master.py`; see the comments there for details.

- Keep geometry about 20 cells away from the PML (which is 10 cells thick).
- Use at least 10 grid points per wavelength (20 recommended), accounting for dielectrics and oblique angles.
- 1000-4000 time steps is typical; check the output time series to confirm the fields have decayed.
- The trustworthy frequency band is set by the pulse type and the frequency parameter (`kmax`: also by the k-vector). The formulas are in `master.py`.
- For `kmax`, avoid k values near integer multiples of 2π/d (a known numerical resonance at the Brillouin-zone center).

## SPICE co-simulation notes

- SPICE-linked ports must sit in non-dispersive media.
- Ports are matched to netlist elements by name, so the names in `master.py` and the netlist must agree. The capacitor named for each basic port is overwritten by the solver with the FDTD cell capacitance.
- `kmax` + SPICE: SPICE cannot run transient simulations on complex data, so each port must be entered **twice**, one after the other, with different names and identical duplicate circuits (one copy for the real part, one for the imaginary). The post-processor recombines them into a single port.
- Extra SPICE nodes (e.g., a load voltage) can be recorded as additional "ports" with no connection to the grid.
- Installation guides for ngspice are provided separately; see [`compile_scripts/`](../compile_scripts).
