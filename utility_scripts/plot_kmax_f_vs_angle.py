import numpy as np
import matplotlib.pyplot as plt
import os

# this script needs to be customized throughout based on:
# k vectors (kx,ky,kz)
# angle desired for plotting
# file names
# what traces to plot up
# etc.

# Generate arrays
index_length = 400
ky = np.linspace(0, 0, index_length)
kz = np.linspace(0, 399, index_length)

# set number of simulations we did (ky and kz are the same length, even if ky or kz is all zeros for example)
num_k = len(ky)

# Setup the figures before the loop so we can "individually" add points to them
fig_refl, ax_refl = plt.subplots(figsize=(10, 6))
fig_trans, ax_trans = plt.subplots(figsize=(10, 6))

# Lists to manually keep track of the max/min data plotted for the colorbar limits
all_r, all_t, all_f, all_theta = [], [], [], []

for i in range(num_k):
    k_val1 = round(ky[i], 2)
    k_val2 = round(kz[i], 2)
    file_path = f'working/sim_{k_val1}_{k_val2}/S_parameters_k.csv'
    
    if not os.path.exists(file_path):
        continue
    
    # Load raw data directly - no artificial grids
    data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
    f = data[:, 0]  # Raw Frequency
    r = data[:, 3]  # Raw Reflection TM (dB)
    t = data[:, 7]  # Raw Transmission TM (dB)

    # Compute angle - what the angle represents depends on k vectors and formula used
    arg = (kz[i] * 3e8) / (2 * np.pi * (f + 1E-9) * 1e9)
    arg = np.clip(arg, 0, 1)
    theta_raw = np.arcsin(arg) * (180 / np.pi)
    
    # Save arrays to a list to use for scatter plotting later
    all_f.append(f)
    all_theta.append(theta_raw)
    all_r.append(r)
    all_t.append(t)

# Flatten lists for unified plotting
all_f = np.concatenate(all_f)
all_theta = np.concatenate(all_theta)
all_r = np.concatenate(all_r)
all_t = np.concatenate(all_t)

# --- Plot Reflection ---
sc_refl = ax_refl.scatter(all_f, all_theta, c=all_r, cmap='rainbow', vmin=-40, vmax=0, s=2, marker='s')
fig_refl.colorbar(sc_refl, ax=ax_refl, label='Reflection (dB)')
ax_refl.set_xlabel('Frequency (GHz)')
ax_refl.set_ylabel('Theta (degrees)')
ax_refl.set_title('Reflection vs Frequency and Theta')
ax_refl.set_ylim(0, 90)
fig_refl.tight_layout()
fig_refl.savefig('refl_f_vs_angle.png')

# --- Plot Transmission ---
sc_trans = ax_trans.scatter(all_f, all_theta, c=all_t, cmap='rainbow', vmin=-2, vmax=0, s=2, marker='s')
fig_trans.colorbar(sc_trans, ax=ax_trans, label='Transmission (dB)')
ax_trans.set_xlabel('Frequency (GHz)')
ax_trans.set_ylabel('Theta (degrees)')
ax_trans.set_title('Transmission vs Frequency and Theta')
ax_trans.set_ylim(0, 90)
fig_trans.tight_layout()
fig_trans.savefig('trans_f_vs_angle.png')

plt.show()