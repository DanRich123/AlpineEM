import numpy as np
from matplotlib import pyplot as plt
import os

# move to the current directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# import S parameter file for old data - no spice source
S_data_old = np.loadtxt('S_parameters_old.csv', delimiter = ',', skiprows = 1)
S_data_old = np.transpose(S_data_old)
freq_old = S_data_old[0]
S_11_amp_old = S_data_old[1]
S_11_phase_old = S_data_old[2]
# import S parameter file - for new data - spice source
S_data = np.loadtxt('S_parameters.csv', delimiter = ',', skiprows = 1)
S_data = np.transpose(S_data)
freq = S_data[0]
S_11_amp = S_data[1]
S_11_phase = S_data[2]

# import realized gain file for old data - no spice source
G_data_old = np.loadtxt('Realized_antenna_gain_old.csv', delimiter = ',', skiprows = 1)
G_data_old = np.transpose(G_data_old)
# only import theta polarization (main pol of interest polarization)
G_amp_ang1_old = G_data_old[1]
G_phase_ang1_old = G_data_old[2]
# import realized gain file for new data - spice source
G_data = np.loadtxt('Realized_antenna_gain.csv', delimiter = ',', skiprows = 1)
G_data = np.transpose(G_data)
# only import theta polarization (main pol of interest polarization) and correct known impedance as needed
G_amp_ang1 = G_data[1] + 10*np.log10(50)
G_phase_ang1 = G_data[2]

plt.figure()
plt.plot(freq_old, S_11_amp_old, color = 'k', label = 'S11 amplitude - no spice source')
plt.plot(freq_old, S_11_phase_old, color='k', label = 'S11 phase - no spice source')
plt.plot(freq, S_11_amp, '--', color = 'r', label = 'S11 amplitude - spice source')
plt.plot(freq, S_11_phase, '--', color='r', label = 'S11 phase - spice source')
plt.grid()
plt.legend()
plt.xlim(2.5,22)
plt.title('S11')
plt.ylabel('Amplitude (dB)')
plt.xlabel('Frequency (GHz)')
plt.savefig('S_parameter_compare')

plt.figure()
plt.plot(freq_old, G_amp_ang1_old, color = 'k', label = 'Theta Polarization: Theta=90, Phi=0, Amplitude - no spice source')
plt.plot(freq_old, G_phase_ang1_old, color = 'k', label = 'Theta Polarization: Theta=90, Phi=0, Phase - no spice source')
plt.plot(freq, G_amp_ang1, '--', color = 'r', label = 'Theta Polarization: Theta=90, Phi=0, Amplitude - spice source')
plt.plot(freq, G_phase_ang1, '--', color = 'r', label = 'Theta Polarization: Theta=90, Phi=0, Phase - spice source')
plt.ylabel('Amplitude (dB)')
plt.xlabel('Frequency (GHz)')
plt.xlim(2.5,22)
plt.ylim(-20,5)
plt.title('Realized Gain')
plt.grid()
plt.legend()
plt.savefig('Realized_gain_plot_compare')

