#!/bin/bash
#SBATCH --job-name=fdtd
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# #SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --array=0-399
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err

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

# establish some directories for placement
if [ $SLURM_ARRAY_TASK_ID -eq 0 ]; then
    mkdir -p "working"
    mkdir -p "touchstones"
    mkdir -p "logs"
fi

# Run simulation
python parallel_kmax.py $SLURM_ARRAY_TASK_ID
