#!/bin/bash
#SBATCH --job-name=fdtd
#SBATCH --partition=aa100
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1       
#SBATCH --gres=gpu:1
##SBATCH --mem=16G               ## ignored by slurm, I think this is cpu memory only
#SBATCH --time=01:00:00
#SBATCH --output=gpu_%j.out
#SBATCH --error=gpu_%j.err

# Set type to either tensorflow or pytorch
export type=tensorflow

# Set filename
export filename=inputs.txt

# Print job info
echo "Job started on $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

# Load modules
module purge
module load cuda
module load anaconda

# Activate environment
conda activate sandbox

# Run simulation
python run.py $filename
python check_gpu_$type.py
python fdtd_$type.py $filename
