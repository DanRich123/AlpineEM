import tensorflow as tf
import subprocess
import sys

print("="*60)
print("CUDA/GPU DIAGNOSTIC FOR TENSORFLOW")
print("="*60)

# Check TensorFlow installation
print(f"\nTensorFlow version: {tf.__version__}")

# Get list of physical devices
gpus = tf.config.list_physical_devices('GPU')
print(f"GPU devices detected: {len(gpus)}")

if gpus:
    print(f"\n✓ TensorFlow GPU support is available")
    
    for i, gpu in enumerate(gpus):
        print(f"\nGPU {i}:")
        print(f"  Device: {gpu.name}")
        print(f"  Type: {gpu.device_type}")
        
        # Try to get device details (may not work on all systems)
        try:
            details = tf.config.experimental.get_device_details(gpu)
            if details:
                print(f"  Details: {details}")
        except:
            print(f"  Details: Not available through TensorFlow API")
    
    # Check if GPU is actually usable
    print("\nTesting GPU computation:")
    try:
        with tf.device('/GPU:0'):
            a = tf.constant([[1.0, 2.0], [3.0, 4.0]])
            b = tf.constant([[1.0, 1.0], [0.0, 1.0]])
            c = tf.matmul(a, b)
        print("✓ GPU computation successful!")
        print(f"  Test result: \n{c.numpy()}")
    except Exception as e:
        print(f"❌ GPU computation failed: {e}")
    
    # Show memory growth setting
    print("\nMemory growth settings:")
    for i, gpu in enumerate(gpus):
        try:
            growth = tf.config.experimental.get_memory_growth(gpu)
            print(f"  GPU {i}: Memory growth = {growth}")
        except:
            print(f"  GPU {i}: Cannot query memory growth")
            
else:
    print("\n❌ GPU is NOT available in TensorFlow")
    print("\nPossible reasons:")
    print("1. TensorFlow was installed without CUDA support (CPU-only version)")
    print("2. No NVIDIA GPU in system")
    print("3. NVIDIA drivers not installed or outdated")
    print("4. CUDA/cuDNN version mismatch")

# Check CUDA build info
print("\n" + "="*60)
print("TENSORFLOW BUILD INFO")
print("="*60)
print(f"\nBuilt with CUDA: {tf.test.is_built_with_cuda()}")
if tf.test.is_built_with_cuda():
    print("✓ TensorFlow was built with CUDA support")
else:
    print("❌ TensorFlow was built WITHOUT CUDA support (CPU-only)")

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
    print("  (Note: TensorFlow includes its own CUDA runtime, so this is optional)")

# Check cuDNN
print("\n" + "="*60)
print("cuDNN CHECK")
print("="*60)
print(f"Built with cuDNN: {tf.test.is_built_with_cuda()}")  # cuDNN is included with CUDA build

# Recommendations
print("\n" + "="*60)
print("RECOMMENDATIONS")
print("="*60)

if not gpus or not tf.test.is_built_with_cuda():
    print("\nTo enable GPU support:")
    print("\n1. Check if you have an NVIDIA GPU:")
    print("   lspci | grep -i nvidia")
    
    print("\n2. Install NVIDIA drivers if missing:")
    print("   sudo ubuntu-drivers autoinstall  # Ubuntu/Debian")
    print("   # or download from: https://www.nvidia.com/Download/index.aspx")
    
    print("\n3. Reinstall TensorFlow with GPU support:")
    print("   # For TensorFlow 2.15+ (includes CUDA/cuDNN):")
    print("   pip install tensorflow[and-cuda]")
    print("\n   # For older versions, you may need:")
    print("   pip install tensorflow-gpu")
    print("\n   # Or use conda:")
    print("   conda install tensorflow-gpu")
    
    print("\n4. Verify CUDA compatibility:")
    print("   Check: https://www.tensorflow.org/install/source#gpu")
    print("   TensorFlow 2.15+: CUDA 12.2, cuDNN 8.9")
    print("   TensorFlow 2.10-2.14: CUDA 11.8, cuDNN 8.6")
    
    print("\n5. Verify installation:")
    print("   python3 -c 'import tensorflow as tf; print(len(tf.config.list_physical_devices(\"GPU\")))'")
else:
    print("\n✓ GPU is properly configured and ready to use!")
    print("\nOptional: Enable memory growth to prevent TensorFlow from allocating all GPU memory:")
    print("   gpus = tf.config.list_physical_devices('GPU')")
    print("   for gpu in gpus:")
    print("       tf.config.experimental.set_memory_growth(gpu, True)")

print("\n" + "="*60)

# Additional TensorFlow-specific info
print("\nADDITIONAL TENSORFLOW INFO")
print("="*60)
print(f"Eager execution: {tf.executing_eagerly()}")
print(f"XLA available: {tf.test.is_built_with_xla()}")
print(f"MKL support: {'MKL' in tf.__config__.show()}")

# List all devices
print("\nAll available devices:")
devices = tf.config.list_physical_devices()
for device in devices:
    print(f"  {device.device_type}: {device.name}")

print("\n" + "="*60)