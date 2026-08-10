#!/bin/bash

export INSTALL_DIR=/home/$USER/fdtd-fortran-ngspice-acc
# if using acc for gpu usage, I use an optional flag in case we want to omit something but that would be for testing purposes
OPT_FLAGS="-O3 -acc=gpu -Mpreprocess -Duse_kmax_version -Duse_spice_version"

# load packages
ml purge
ml cuda
ml nvhpc_sdk
ml gcc/11.2.0
ml autotools/2.71

# Now add NGSPICE executables and libraries to PATH
export PATH=$INSTALL_DIR/ngspice-sf-mirror-ngspice-44.2/local/bin:$PATH
export LD_LIBRARY_PATH=$INSTALL_DIR/ngspice-sf-mirror-ngspice-44.2/local/lib:$LD_LIBRARY_PATH

# Now compile my own code - move to the correct folder manually before doing so
nvc -I${INSTALL_DIR}/ngspice-sf-mirror-ngspice-44.2/src/include/ngspice -c ngspice_interfaces.c
nvfortran -c ngspice_interface.F90
nvfortran -c circuit.F90
nvfortran ${OPT_FLAGS} -c fdtd_solver.f90
nvfortran ${OPT_FLAGS} fdtd_solver.o circuit.o ngspice_interface.o ngspice_interfaces.o -L${INSTALL_DIR}/ngspice-sf-mirror-ngspice-44.2/local/lib -lngspice -o FDTD-KMAX-SPICE-ACC