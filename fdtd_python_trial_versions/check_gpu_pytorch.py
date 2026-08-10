import torch as th
import subprocess
import sys

print("="*60)
print("CUDA/GPU DIAGNOSTIC")
print("="*60)

# Check PyTorch installation
print(f"\nPyTorch version: {th.__version__}")
print(f"CUDA available: {th.cuda.is_available()}")
print(f"CUDA built with PyTorch: {th.version.cuda}")

if th.cuda.is_available():
    print(f"\nNumber of GPUs: {th.cuda.device_count()}")
    for i in range(th.cuda.device_count()):
        print(f"\nGPU {i}:")
        print(f"  Name: {th.cuda.get_device_name(i)}")
        props = th.cuda.get_device_properties(i)
        print(f"  Compute capability: {props.major}.{props.minor}")
        print(f"  Total memory: {props.total_memory / 1e9:.2f} GB")
        print(f"  Multi-processors: {props.multi_processor_count}")
else:
    print("\n❌ CUDA is NOT available in PyTorch")
    print("\nPossible reasons:")
    print("1. PyTorch was installed without CUDA support (CPU-only version)")
    print("2. No NVIDIA GPU in system")
    print("3. NVIDIA drivers not installed or outdated")
    print("4. CUDA toolkit version mismatch")

# Check for nvidia-smi
print("\n" + "="*60)
print("NVIDIA DRIVER CHECK")
print("="*60)
try:
    result = subprocess.run(['nvidia-smi'], capture_output=True, text=True)
    if result.returncode == 0:
        print("\n✓ nvidia-smi found - GPU drivers are installed")
        print("\nGPU Status:")
        print(result.stdout)
    else:
        print("\n❌ nvidia-smi failed")
except FileNotFoundError:
    print("\n❌ nvidia-smi not found - NVIDIA drivers may not be installed")

# Check CUDA installation
print("\n" + "="*60)
print("CUDA TOOLKIT CHECK")
print("="*60)
try:
    result = subprocess.run(['nvcc', '--version'], capture_output=True, text=True)
    if result.returncode == 0:
        print("\n✓ CUDA toolkit found:")
        print(result.stdout)
    else:
        print("\n❌ nvcc not found")
except FileNotFoundError:
    print("\n⚠ nvcc not found - CUDA toolkit may not be in PATH")
    print("  (Note: PyTorch includes its own CUDA runtime, so this is optional)")

# Recommendations
print("\n" + "="*60)
print("RECOMMENDATIONS")
print("="*60)

if not th.cuda.is_available():
    print("\nTo enable GPU support:")
    print("\n1. Check if you have an NVIDIA GPU:")
    print("   lspci | grep -i nvidia")
    
    print("\n2. Install NVIDIA drivers if missing:")
    print("   sudo ubuntu-drivers autoinstall  # Ubuntu/Debian")
    print("   # or download from: https://www.nvidia.com/Download/index.aspx")
    
    print("\n3. Reinstall PyTorch with CUDA support:")
    print("   # For CUDA 11.8:")
    print("   pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118")
    print("\n   # For CUDA 12.1:")
    print("   pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121")
    print("\n   # Check available versions at: https://pytorch.org/get-started/locally/")
    
    print("\n4. Verify installation:")
    print("   python3 -c 'import torch; print(torch.cuda.is_available())'")
else:
    print("\n✓ GPU is properly configured and ready to use!")

print("\n" + "="*60)