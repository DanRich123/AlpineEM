#!/bin/bash
#SBATCH --job-name=fdtd
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --mem=8G ## indicates slurm will skip it
#SBATCH --time=06:00:00
#SBATCH --output=cpu_%j.out
#SBATCH --error=cpu_%j.err

# Set type to either tensorflow or pytorch
export type=tensorflow

# Set filename
export filename=inputs.txt

# Print job info
echo "Job started on $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"

# Load modules
module purge
module load anaconda

# Activate environment
conda activate sandbox

# Run simulation
python run.py $filename
python fdtd_$type.py $filename
