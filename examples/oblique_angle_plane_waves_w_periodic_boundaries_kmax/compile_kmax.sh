#!/bin/bash

# optional flags
OPT_FLAGS="-O3 -fpp -Duse_kmax_version"

# load modules
ml purge
ml intel

# Now compile my code
ifx ${OPT_FLAGS} fdtd_solver.f90 -o FDTD-KMAX