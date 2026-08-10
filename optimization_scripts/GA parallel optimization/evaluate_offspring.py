import numpy as np
import sys
import os
import shutil
import traceback
from pathlib import Path

# Import the fitness function from main GA script
from main_ga import run_program, NEEDED_FILES

def main():
    if len(sys.argv) != 3:
        print("Usage: python evaluate_offspring.py <gen_folder> <individual_index>", flush=True)
        sys.exit(1)
    
    gen_folder = sys.argv[1]
    idx = int(sys.argv[2])
    
    print("="*60, flush=True)
    print(f"Evaluating Individual {idx} in {gen_folder}", flush=True)
    print("="*60, flush=True)
    
    # Load the offspring population
    offspring_file = os.path.join(gen_folder, "offspring.npy")
    if not os.path.exists(offspring_file):
        print(f"ERROR: Offspring file not found: {offspring_file}", flush=True)
        sys.exit(1)
    
    offspring = np.load(offspring_file, allow_pickle=True)
    
    if idx >= len(offspring):
        print(f"ERROR: Individual index {idx} out of range (max: {len(offspring)-1})", flush=True)
        sys.exit(1)
    
    individual = offspring[idx]
    print(f"Individual parameters: {list(individual)}", flush=True)
    
    # Create subfolder for this individual
    ind_folder = os.path.join(gen_folder, f"ind_{idx}")
    os.makedirs(ind_folder, exist_ok=True)
    print(f"Working directory: {ind_folder}", flush=True)
    
    # Copy necessary simulation files
    print("\nCopying simulation files...", flush=True)
    files_copied = 0
    files_missing = []
    
    for f in NEEDED_FILES:
        if os.path.exists(f):
            try:
                # Handle both files and directories
                if os.path.isdir(f):
                    shutil.copytree(f, os.path.join(ind_folder, os.path.basename(f)), 
                                  dirs_exist_ok=True)
                else:
                    shutil.copy2(f, ind_folder)
                files_copied += 1
                print(f"  ✓ {f}", flush=True)
            except Exception as e:
                print(f"  ✗ Failed to copy {f}: {e}", flush=True)
                files_missing.append(f)
        else:
            print(f"  ⚠ Not found: {f}")
            files_missing.append(f)
    
    print(f"\nCopied {files_copied}/{len(NEEDED_FILES)} files", flush=True)
    
    if files_missing:
        print(f"WARNING: Missing files: {files_missing}", flush=True)
        print("Proceeding anyway - simulation may fail if these are required.", flush=True)
    
    # Run the fitness evaluation
    print("\n" + "="*60, flush=True)
    print("Fitness Evaluation", flush=True)
    print("="*60, flush=True)
    
    try:
        fitness = run_program(individual, work_dir=ind_folder)
        print("\n" + "="*60)
        print(f"Simulation Complete")
        print(f"Fitness: {fitness[0]:.6f}")
        print("="*60)
        
        # Save fitness result
        fitness_file = os.path.join(gen_folder, f"fit_{idx}.npy")
        np.save(fitness_file, fitness)
        print(f"\nSaved fitness to: {fitness_file}")
        
        # Create a summary file
        summary_file = os.path.join(ind_folder, "summary.txt")
        with open(summary_file, "w") as f:
            f.write(f"Individual {idx}\n")
            f.write(f"Parameters: {list(individual)}\n")
            f.write(f"Fitness: {fitness[0]:.6f}\n")
        
        print(f"SUCCESS: Individual {idx} evaluated successfully")
        return 0
        
    except Exception as e:
        print("\n" + "="*60)
        print("ERROR: Simulation Failed")
        print("="*60)
        print(f"Exception: {e}")
        print("\nTraceback:")
        traceback.print_exc()
        
        # Save worst fitness on error
        worst_fitness = (-1e10,)
        fitness_file = os.path.join(gen_folder, f"fit_{idx}.npy")
        np.save(fitness_file, worst_fitness)
        print(f"\nSaved worst fitness ({worst_fitness[0]}) to: {fitness_file}")
        
        # Save error log
        error_file = os.path.join(ind_folder, "error.txt")
        with open(error_file, "w") as f:
            f.write(f"Individual {idx} FAILED\n")
            f.write(f"Parameters: {list(individual)}\n")
            f.write(f"Error: {str(e)}\n")
            f.write(f"\nTraceback:\n")
            traceback.print_exc(file=f)
        
        print(f"FAILED: Individual {idx} evaluation failed")
        return 1

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)