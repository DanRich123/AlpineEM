import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

# this script needs to be customized throughout based on:
# k vectors (kx,ky,kz)
# angle(s) desired for plotting
# file names
# what traces to plot up
# etc.

def generate_ky_kz_arrays(n_theta, n_phi, k_magnitude):
    # Create uniform spacing in angles
    theta = np.linspace(np.pi/6, np.pi/6+np.pi*2/3, n_theta)  # 30 to 150 degrees
    phi = np.linspace(-np.pi/3, np.pi/3, n_phi)  # -60 to 60 degrees
    # Create meshgrid
    THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
    # Convert to ky, kz (with theta from Z axis)
    ky = k_magnitude * np.sin(THETA) * np.sin(PHI)
    kz = k_magnitude * np.cos(THETA)
    # Flatten for indexing
    ky_flat = ky.flatten()
    kz_flat = kz.flatten()
    
    # Also return theta and phi arrays for plotting
    theta_flat = THETA.flatten()
    phi_flat = PHI.flatten()
    
    return ky_flat, kz_flat, theta_flat, phi_flat

# Target frequency
target_freq = 9.0  # GHz

# set num angles for each
n_theta = 15
n_phi = 15

# Generate arrays
ky, kz, theta_array, phi_array = generate_ky_kz_arrays(n_theta=n_theta, n_phi=n_phi, k_magnitude=2*np.pi*target_freq*1E9/3E8)

# set number used by loops below
num_k = len(ky)

# First pass — collect all frequency ranges
all_freqs = []

for i in range(num_k):
    k_val1 = round(ky[i], 2)
    k_val2 = round(kz[i], 2)
    file_path = f'working/sim_{k_val1}_{k_val2}/S_parameters_k.csv'
    
    if not os.path.exists(file_path):
        print(f"[WARNING] Missing file: {file_path}")
        continue
    
    data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
    all_freqs.append(data[:, 0])

# Build common frequency grid
all_freqs_flat = np.concatenate(all_freqs)
freq_common = np.linspace(np.min(all_freqs_flat), np.max(all_freqs_flat), 500)
num_freq = len(freq_common)

# Initialize arrays
refl = np.full(num_k, np.nan)
trans = np.full(num_k, np.nan)

# Second pass — load and interpolate at 9 GHz
for i in range(num_k):
    k_val1 = round(ky[i], 2)
    k_val2 = round(kz[i], 2)
    file_path = f'working/sim_{k_val1}_{k_val2}/S_parameters_k.csv'
    
    if not os.path.exists(file_path):
        continue
    
    data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
    f = data[:, 0]  # Frequency
    r = data[:, 3]  # Reflection TM (dB)
    t = data[:, 7]  # Transmission TM (dB)

    # Interpolate at target frequency
    try:
        r_interp = interp1d(f, r, bounds_error=False, fill_value=np.nan)
        refl[i] = r_interp(target_freq)
    except Exception as e:
        print(f"[ERROR] Interpolation failed at k={k_val1},{k_val2}: {e}")
        continue
    
    try:
        t_interp = interp1d(f, t, bounds_error=False, fill_value=np.nan)
        trans[i] = t_interp(target_freq)
    except Exception as e:
        print(f"[ERROR] Interpolation failed at k={k_val1},{k_val2}: {e}")
        continue

# Convert theta and phi to degrees for plotting
theta_deg = theta_array * (180 / np.pi)
phi_deg = phi_array * (180 / np.pi)

# reshape for plotting
theta_2d = theta_deg.reshape(n_theta, n_phi)
phi_2d = phi_deg.reshape(n_theta, n_phi)
refl_2d = refl.reshape(n_theta, n_phi)
trans_2d = trans.reshape(n_theta, n_phi)

# Plot Reflection at target frequency vs theta and phi
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
mesh1 = plt.pcolormesh(phi_2d, theta_2d, refl_2d, shading='gouraud', cmap='rainbow', vmin=-25, vmax=0)
plt.colorbar(mesh1, label='Reflection (dB)')
plt.xlabel('Phi (degrees)')
plt.ylabel('Theta (degrees)')
plt.title('TM Metasurface Reflection at {} GHz vs Theta and Phi'.format(target_freq))
#plt.xlim(-45,45)
#plt.ylim(30,150)
plt.grid()
plt.tight_layout()

# Plot Transmission at target frequency vs theta and phi
plt.subplot(1, 2, 2)
mesh2 = plt.pcolormesh(phi_2d, theta_2d, trans_2d, shading='gouraud', cmap='rainbow', vmin=-25, vmax=0)
plt.colorbar(mesh2, label='Transmission (dB)')
plt.xlabel('Phi (degrees)')
plt.ylabel('Theta (degrees)')
plt.title('TM Metasurface Transmission at {} GHz vs Theta and Phi'.format(target_freq))
#plt.xlim(-50,50)
#plt.ylim(30,150)
plt.grid()
plt.tight_layout()

plt.savefig('reflection_transmission_{}GHz_angle_vs_angle.png'.format(target_freq), dpi=300, bbox_inches='tight')