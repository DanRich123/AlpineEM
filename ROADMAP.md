# AlpineEM Roadmap
 
Planned fixes and features for [AlpineEM](./README.md). Items are grouped by area and are not listed in priority order. This is a working list, so plans may change.
 
## Solver core and materials
 
- [ ] Permeability beyond the free-space value for bulk materials
- [ ] More plasma terms (thermal pressure, cyclotron motion, etc.) in the multipole Drude model (currently up to 6 poles)
- [ ] Optional double-precision arithmetic in Fortran
- [ ] Thin sub-cell wires
- [ ] Additional dispersive and dissipative materials (Debye, Lorentz, etc.)
- [ ] Optional x, y, z sheet binary read-ins, similar to the optional bulk material array
- [ ] Output of E and H fields at cell surfaces, not just at cell centers
## Sources and far-field
 
- [ ] Additional internal (non-SPICE generated) source types (e.g., CW, raised-cosine envelope)
- [ ] Additional far-field method: choose a frequency and compute all angles at that frequency - better for many angles and parallelization in general for openMP and openACC
## Ports and feeds
 
- [ ] 2D Helmholtz solver for frequency-dependent (dispersive or non-dispersive) mode generation for gridded feeds
- [ ] Support for dispersive field coefficients in gridded feeds
- [ ] Allow lumped ports at the first cube (currently an error is thrown to warn the user if this is selected accidentally)
## Oblique incidence (`kmax`)
 
- [ ] Improve CPML and incident plane-wave configurations for `kmax` (oblique-angle periodic boundary conditions)
## Performance and scaling
 
- [ ] MPI support to spread memory across multiple nodes (multiple computers)
## Visualization
 
- [ ] ParaView and Python viewer renderings of wave-port locations for unit-cell designs
## Education
 
- [ ] Improve and validate the pure-Python version of the solver (intended for student education, not performance)
