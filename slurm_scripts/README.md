# SLURM Submission Scripts

SLURM batch execution scripts for submitting FDTD solver runs, post-processing pipelines, and parameter sweeps across HPC cluster environments.

---

## Submission Scripts Overview

| Script | Description | Hardware / Target Partition | Workload / Execution Flow |
|---|---|---|---|
| `run.sh` | Single-core serial CPU batch execution script. Loads basic Intel compiler modules and Conda environment. | `acpu` partition (1 node, 1 task, 1 CPU core) | Runs `master.py` followed automatically by `post_processor.py`. |
| `run_mp.sh` | OpenMP multi-threaded CPU execution script with core-binding optimizations (`OMP_PROC_BIND=close`, `OMP_PLACES=cores`). | `acpu` partition (1 node, 1 task, 8 CPU cores) | Configures OpenMP environment threads and runs `master.py` followed by `post_processor.py`. |
| `run_acc.sh` | OpenACC GPU-accelerated batch execution script. | `aa100` partition (1 node, 1 task, 1 NVIDIA A100 GPU) | Loads CUDA and NVIDIA HPC SDK (`nvhpc_sdk`) modules, running `master.py` followed by `post_processor.py`. |
| `run_parallel.sh` | High-throughput SLURM job array script for executing large parameter sweeps in parallel. | `acpu` partition (Job Array: tasks 0 to 399) | Initializes directory structures (`working/`, `touchstones/`, `logs/`) and dispatches `parallel_kmax.py $SLURM_ARRAY_TASK_ID`. |

---

## Technical & Cluster Setup Notes

* **Conda Environment:** All job scripts automatically activate the customized `sandbox` Conda environment prior to execution.
* **ngspice / SPICE Integration:** Commented configuration blocks are provided across all scripts to load `autotools` and export dynamic library paths (`LD_LIBRARY_PATH`) when coupling circuit co-simulations via ngspice.
