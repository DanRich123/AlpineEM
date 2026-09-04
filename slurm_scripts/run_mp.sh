#!/bin/bash
#SBATCH --job-name=mp_fdtd
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
# #SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=mp_%j.out
#SBATCH --error=mp_%j.err

# Print job info
echo "Job started on $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"

# Load modules
module purge
module load anaconda
module load intel
# if using spice, then use these (use this intel version and add autotools)
# module load intel/2024.2.1
# module load autotools/2.71

# Activate environment
conda activate sandbox

# Additionally if spice, make sure spice libraries are accessible to slurm
# export INSTALL_DIR=/home/$USER/fdtd-fortran-ngspice
# export PATH=$INSTALL_DIR/ngspice-sf-mirror-ngspice-44.2/local/bin:$PATH
# export LD_LIBRARY_PATH=$INSTALL_DIR/ngspice-sf-mirror-ngspice-44.2/local/lib:$LD_LIBRARY_PATH

# OpenMP settings
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close      
export OMP_PLACES=cores

# Run simulation
python master.py
# python master_clear.py
# python master_metal.py
python post_processor.py
