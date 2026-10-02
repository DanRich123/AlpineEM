import numpy as np
import os
from matplotlib import pyplot as plt

# move to the current directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# import the conformed media S-parameter data
S = np.loadtxt('S_parameters.csv', delimiter=',', skiprows=1)
S = np.transpose(S)

# import the perfectly aligned media S-parameter data
S_pf = np.loadtxt('S_parameters_perfect_alignment.csv', delimiter=',', skiprows=1)
S_pf = np.transpose(S_pf)

# The reflection phase will need a 1 cell total phase shift due to 1/2 cell placement difference - transmission doesn't need one
S[2] = S[2] + 2*np.pi*S[0]*1E9/3E8*0.5E-3

# Create 2x2 grid
fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True)

# 1. Reflection - Amplitude (Top Left)
axs[0, 0].plot(S[0], S[1], label='Conformal')
axs[0, 0].plot(S_pf[0], S_pf[1], '--', color='k', label='Perfectly Aligned')
axs[0, 0].set_ylabel('Amplitude (dB)')
axs[0, 0].set_title('Reflection Amplitude')
axs[0, 0].grid(True)
axs[0, 0].legend()

# 2. Transmission - Amplitude (Top Right)
axs[0, 1].plot(S[0], S[3], label='Conformal')
axs[0, 1].plot(S_pf[0], S_pf[3], '--', color='k', label='Perfectly Aligned')
axs[0, 1].set_ylabel('Amplitude (dB)')
axs[0, 1].set_title('Transmission Amplitude')
axs[0, 1].grid(True)
axs[0, 1].legend()

# 3. Reflection - Phase (Bottom Left)
axs[1, 0].plot(S[0], np.unwrap(S[2]), label='Conformal')
axs[1, 0].plot(S_pf[0], np.unwrap(S_pf[2]), '--', color='k', label='Perfectly Aligned')
axs[1, 0].set_xlabel('Frequency (GHz)')
axs[1, 0].set_ylabel('Phase (Radians)')
axs[1, 0].set_title('Reflection Phase')
axs[1, 0].grid(True)
axs[1, 0].legend()

# 4. Transmission - Phase (Bottom Right)
axs[1, 1].plot(S[0], np.unwrap(S[4]), label='Conformal')
axs[1, 1].plot(S_pf[0], np.unwrap(S_pf[4]), '--', color='k', label='Perfectly Aligned')
axs[1, 1].set_xlabel('Frequency (GHz)')
axs[1, 1].set_ylabel('Phase (Radians)')
axs[1, 1].set_title('Transmission Phase')
axs[1, 1].grid(True)
axs[1, 1].legend()

plt.tight_layout()
plt.savefig('Combined_S_Parameters_Plot.png')
plt.show()

