#!/bin/bash

export INSTALL_DIR=/home/$USER/fdtd-fortran-ngspice

# load packages
ml purge
ml intel/2024.2.1
ml autotools/2.71

# build NGSPICE from source:
mkdir -p $INSTALL_DIR
cd $INSTALL_DIR
wget https://github.com/danchitnis/ngspice-sf-mirror/archive/refs/tags/ngspice-44.2.tar.gz
tar -xf ngspice-44.2.tar.gz
cd ngspice-sf-mirror-ngspice-44.2/
./autogen.sh
./configure --with-x --with-ngshared --enable-cider --prefix="$PWD/local" --libdir="$PWD/local/lib" CFLAGS="-m64 -O2" LDFLAGS="-m64 -s"
make -j 8
make install