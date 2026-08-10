#!/bin/bash
#SBATCH --job-name=ga_eval
#SBATCH --qos=normal
#SBATCH --time=06:00:00
#SBATCH --partition=amilan
#SBATCH --no-requeue
#SBATCH --output=%x/%A_%a.out
#SBATCH --error=%x/%A_%a.err

# ============================================================
# SLURM Array Job for GA Individual Evaluation
# 
# This script is called by main_ga.py with:
#   sbatch --array=0-N run_generation.sh generation_X
#
# Each array task evaluates one individual
# ============================================================

echo "======================================"
echo "SLURM Job Information"
echo "======================================"
echo "Job ID:           $SLURM_JOB_ID"
echo "Array Job ID:     $SLURM_ARRAY_JOB_ID"
echo "Array Task ID:    $SLURM_ARRAY_TASK_ID"
echo "Node:             $SLURM_NODELIST"
echo "CPUs:             $SLURM_CPUS_PER_TASK"
echo "Memory:           $SLURM_MEM_PER_NODE MB"
echo "Start time:       $(date)"
echo "======================================"

# Parse arguments
if [ $# -eq 0 ]; then
    echo "ERROR: No generation folder provided"
    echo "Usage: sbatch --array=0-N run_generation.sh <gen_folder>"
    exit 1
fi

GEN_FOLDER=$1
IDX=$SLURM_ARRAY_TASK_ID

echo ""
echo "Generation folder: $GEN_FOLDER"
echo "Individual index:  $IDX"
echo ""

# Ensure generation folder exists
if [ ! -d "$GEN_FOLDER" ]; then
    echo "ERROR: Generation folder does not exist: $GEN_FOLDER"
    exit 1
fi

# Create logs subdirectory in generation folder
LOGS_DIR="${GEN_FOLDER}/logs"
mkdir -p "$LOGS_DIR"

# Redirect output to generation-specific log files
OUT_FILE="${LOGS_DIR}/ind_${IDX}.out"
ERR_FILE="${LOGS_DIR}/ind_${IDX}.err"

# Function to log both to console and file
log_output() {
    # This captures all further output
    exec > >(tee -a "$OUT_FILE") 2> >(tee -a "$ERR_FILE" >&2)
}

# Start logging
log_output

echo "======================================"
echo "Environment Setup"
echo "======================================"

# Load required modules
module purge  # Clean slate
module load anaconda 2>/dev/null || echo "Warning: anaconda module not found"
module load intel 2>/dev/null || echo "Warning: intel module not found"

# Activate conda environment
# IMPORTANT: Use 'source activate' not 'conda activate' in SLURM scripts
echo "Activating conda environment: sandbox"
source activate sandbox

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to activate conda environment 'sandbox'"
    echo "Available environments:"
    conda env list
    exit 1
fi

echo "Conda environment activated successfully"
echo "Python location: $(which python)"
echo "Python version: $(python --version)"

# Print current working directory
echo ""
echo "Working directory: $(pwd)"
echo ""

echo ""
echo "======================================"
echo "Running Evaluation"
echo "======================================"

# Run the Python evaluation script
python evaluate_offspring.py "$GEN_FOLDER" "$IDX"

EXIT_CODE=$?

echo ""
echo "======================================"
echo "Job Complete"
echo "======================================"
echo "Exit code:  $EXIT_CODE"
echo "End time:   $(date)"
echo "======================================"

# Check if fitness file was created
FITNESS_FILE="${GEN_FOLDER}/fit_${IDX}.npy"
if [ -f "$FITNESS_FILE" ]; then
    echo "✓ Fitness file created: $FITNESS_FILE"
else
    echo "✗ WARNING: Fitness file not found: $FITNESS_FILE"
    EXIT_CODE=1
fi

echo ""

exit $EXIT_CODE