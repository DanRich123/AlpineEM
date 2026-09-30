# Examples

These examples are chosen to validate the FDTD+SPICE solver and its features, and/or demonstrate how to model various scenarios of interest. Each folder within the examples folder is considered self-contained.

Current examples include:

- The monostatic scattering of an incident plane wave off of a metal sphere
- The aperture area and realized gain of a receiving non-Foster loaded monopole-like antenna (gridded feed/port type) mounted over an infinite ground plane
- Oblique angle plane wave incidence with periodic boundary conditions - a `kmax` variant required example 
- A radiating dipole antenna with the basic port type (parallel plate-like gaps) 
- A radiating dipole antenna with the basic port type (parallel plate-like gaps) that uses a SPICE source instead of a local FDTD generated source
- A 3-port stripline wave guide design that operates in a quasi-TEM mode below around 250 MHz that builds on the non-Foster antenna example
- Time and frequency domain electromagnetic field visualization using 2D cell centered field slices (python) or full 3D cell centered fields (ParaView)
- Transmission and reflection through a quasi-1D infinite (periodic boundary conditions) single slab of a 3 pole Drude (plasma) material
- Transmission and reflection through a quasi-1D infinite (periodic boundary conditions) multilayered slab made of dielectrics and thin sheets
  
Forthcoming examples include:
- The monostatic scattering of an incident plane wave off of a metal half sphere over an infinite ground plane - validation of an infinite ground plane usage (incident and NFFF)
- A dielectric conformal averaging (permttivity averaging) design
- An Archimedean spiral antenna design
- Simple bandstop FSS design
- Optimization examples (adjoint-based gradient descent, genetic algorithm w/ and w/out and HPC system, and Bayesian statistical optimization)
