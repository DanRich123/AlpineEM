# GPU Diagnostics and PyTorch FDTD Electromagnetic Simulator

This repository contains Python scripts for validating GPU hardware acceleration (via PyTorch and TensorFlow) and performing 3D Finite-Difference Time-Domain (FDTD) electromagnetic field simulations using PyTorch.

> **⚠️ Notice / Disclaimer**  
> These scripts are part of an **experimental / trial codebase** designed primarily for teaching, learning, and educational purposes. The implementation contains **known logic and physical errors**—including incomplete periodic boundary conditions and unhandled material averaging at edges—and should **not** be used for production-level electromagnetic simulations without significant revision.

---

## Overview

* **`check_gpu_pytorch.py`**  
  Diagnoses PyTorch GPU support by inspecting CUDA availability, listing detected NVIDIA GPU hardware properties, validating driver accessibility through `nvidia-smi`, and checking for local CUDA toolkit installations (`nvcc`).

* **`check_gpu_tensorflow.py`**  
  Performs GPU environment diagnostics for TensorFlow. It checks detected physical GPU devices, runs a small tensor matrix multiplication on `/GPU:0` to test operational usability, and reports driver and build information.

* **`fdtd_pytorch.py`**  
  An FDTD electromagnetic wave simulation tool built on PyTorch GPU tensors. Features include:
  * **GPU Profiling:** Custom `GPUProfiler` context manager for tracking GPU memory allocation and runtime per simulation phase.
  * **Input Parsing:** Reads grid dimensions, simulation timing, excitation sources (antennas or plane waves), thin sheet resistive parameters, geometric primitives (blocks, spheres, cylinders), and far-field angle requests.
  * **PML Setup:** Initializes Perfectly Matched Layer (PML) boundary conditions and pre-calculates field absorption/update coefficients on the GPU.
  * **Geometry Voxelization:** Maps volumetric materials and 2D thin conductive sheets onto a 3D Yee grid.
