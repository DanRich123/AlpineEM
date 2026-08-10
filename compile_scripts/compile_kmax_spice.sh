#!/bin/bash

export INSTALL_DIR=/home/$USER/fdtd-fortran-ngspice
# normal optional flag
OPT_FLAGS="-O3 -fpp -Duse_kmax_version -Duse_spice_version"

# load packages
ml purge
ml intel/2024.2.1
ml autotools/2.71

# Now add NGSPICE executables and libraries to PATH
export PATH=$INSTALL_DIR/ngspice-sf-mirror-ngspice-44.2/local/bin:$PATH
export LD_LIBRARY_PATH=$INSTALL_DIR/ngspice-sf-mirror-ngspice-44.2/local/lib:$LD_LIBRARY_PATH

# Now compile my own code - move to the correct folder manually before doing so
icc -I${INSTALL_DIR}/ngspice-sf-mirror-ngspice-44.2/src/include/ngspice -c ngspice_interfaces.c
ifx -c ngspice_interface.F90
ifx -c circuit.F90
ifx ${OPT_FLAGS} -c fdtd_solver.f90
ifx ${OPT_FLAGS} fdtd_solver.o circuit.o ngspice_interface.o ngspice_interfaces.o -L${INSTALL_DIR}/ngspice-sf-mirror-ngspice-44.2/local/lib -lngspice -o FDTD-KMAX-SPICE