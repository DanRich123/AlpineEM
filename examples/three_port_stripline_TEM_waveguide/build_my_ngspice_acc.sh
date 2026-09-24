#!/bin/bash

export INSTALL_DIR=/home/$USER/fdtd-fortran-ngspice-acc

# Load packages
ml purge
ml gcc/11.2.0
ml autotools/2.71

# Build NGSPICE from source
mkdir -p $INSTALL_DIR
cd $INSTALL_DIR
wget https://github.com/danchitnis/ngspice-sf-mirror/archive/refs/tags/ngspice-44.2.tar.gz
tar -xf ngspice-44.2.tar.gz
cd ngspice-sf-mirror-ngspice-44.2/

./autogen.sh

./configure \
    --with-x \
    --enable-ngshared \
    --with-ngshared \
    --enable-cider \
    --disable-openmp \
    --prefix="$PWD/local" \
    --libdir="$PWD/local/lib" \
    CFLAGS="-O2 -fPIC" \
    CXXFLAGS="-O2 -fPIC" \
    FFLAGS="-O2 -fPIC"

make -j 8
make install