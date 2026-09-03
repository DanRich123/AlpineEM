LICENSE file is included in this folder and applies to all files within this directory as stated in the license.

INFORMATION ON INCLUDED FOLDERS AND FILES:
There is a main fdtd folder with the main .f90 file and a python script for execution, along with spice files if compiling to use with spice. 
There is a utility scripts folder with scripts for things such as post processing of data and creating plots for fields.
There is a folder with scripts for compiling FDTD and building spice if needed.
There is a folder with slurm submission scripts, this is especially useful for openMP or ACC compilations.
There is also a folder with a statics solver that is primarily used for TEM modes for gridded lumped port generation in the fdtd solver but it has other uses.
There is a folder for geometry viewing. Paraview is used for geometry viewing, a macro that can be uploaded to Paraview for easy viewing is also available.
There are several optimization script examples given as well in the optimization folder - adjoint (gradient) optimization and genetic algorithms.
There is a full python version (pytorch and tensorflow) that is trial. It is incomplete and not well tested.
There is a file that lists (most) updates since the last version was published.
There is also a folder of examples that will be updated periodically. Some example files might become slightly outdated over time so use caution.
Lastly, there is a list of items I am planning to fix and/or add.

METHODS USED IN THE SOLVER:
Yee method for FDTD using a standard rectangular coordinate system.
Standard CFL with option for user to reduce further.
Convolutional PML is used via the ADE approach.
There are 4 versions - normal and kmax, w/ and w/out spice.
'normal' refers to the bulk of the solver, while kmax refers to a specific alteration designed for oblique angle incident plane waves into periodic boundary conditions.
2 source types are available for most simulations - Gaussian and a normalized differentiated Gaussian.
kmax only has 1 incident wave type.
The source types can be used for lumped port excitation or plane wave excitation via the TF/SF approach.
Current far field method utilizes the adaptive on the fly time domain version for broadband frequency information at selective angles.
Similarly, adaptive on the fly approach for S parameter calculations as well.
Sheet method is the sub cell method by Smith and Mahoney.
Auxiliary differential equation approach for dispersive media - only plasmas (Drude) currently allowed.
My code doesn't do vacuum sheets, they are skipped over intentionally as the mechanism for filtering if sheets exist or not.
My code doesn't do pec for sheets at zero impedance, it does a high sigma w/ thickness in mind.
I specifically, on purpose, don't use the FDTD cell capacitance for internal (non-spice) lumped ports, I only do this for spice lumped ports. 
The spice FDTD locations must have non-dispersive FDTD capacitance (aka non-dispersive permittivity in the FDTD grid).
No specific 'method' for gridded lumped ports. These ports can be arbitrary geometries but expect E,H field non-dispersive coefficients to know what the pattern is.
I currently use the statics python solver to generate E,H patterns for TEM modes. 

FINAL NOTES:
There are several key notes related to 1) dispersive modes and 2) mode types in gridded lumped ports:
If operating in the LTI regime, one can still use a non-dispersive mode excitation in a dispersive context, but only keep applicable frequencies from the solver.
If in the non-LTI regime, a narrow banded source (or CW) could be used, but are not currently available. 
A Helmholtz equation solver is needed for non-TEM mode generation.
There are future plans for additional source types, a Helmholtz solver, and to allow for dispersive coefficients.
Far field calculations are hard to parallelize well on multithreaded CPU (openMP) or GPU (openACC) due to atomics so far field angles can slow it down.
The adaptive on the fly time domain version is not well suited to parallelization but a fully parallelized version would be RAM extensive.