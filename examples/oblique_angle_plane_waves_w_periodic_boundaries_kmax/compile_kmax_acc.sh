#!/bin/bash

# optional flags
OPT_FLAGS="-O3 -acc=gpu -Mpreprocess -Duse_kmax_version"

# load modules
module purge
module load cuda
module load nvhpc_sdk

# Now compile my code
nvfortran ${OPT_FLAGS} fdtd_solver.f90 -o FDTD-KMAX-ACC