#!/bin/bash
#SBATCH --job-name=acc_fdtd
#SBATCH --partition=aa100
#SBATCH --qos=gpu-normal #or gpu-testing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1       
#SBATCH --gres=gpu:1
# #SBATCH --mem=16G  # I think this is cpu memory only
#SBATCH --time=01:00:00
#SBATCH --output=acc_%j.out
#SBATCH --error=acc_%j.err

# Print job info
echo "Job started on $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

# Load modules
module purge
module load cuda
module load nvhpc_sdk
module load anaconda
# if using spice
# module load gcc/11.2.0
# module load autotools/2.71

# Activate environment
conda activate sandbox

# Additionally if spice, make sure spice libraries are accessible to slurm
# export INSTALL_DIR=/home/$USER/fdtd-fortran-ngspice-acc
# export PATH=$INSTALL_DIR/ngspice-sf-mirror-ngspice-44.2/local/bin:$PATH
# export LD_LIBRARY_PATH=$INSTALL_DIR/ngspice-sf-mirror-ngspice-44.2/local/lib:$LD_LIBRARY_PATH

# Run simulation
python master.py
# python master_clear.py
# python master_metal.py
python post_processor.py
