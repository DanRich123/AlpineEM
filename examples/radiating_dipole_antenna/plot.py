import numpy as np
from matplotlib import pyplot as plt
import os

# move to the current directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# import S parameter file
S_data = np.loadtxt('S_parameters.csv', delimiter = ',', skiprows = 1)
S_data = np.transpose(S_data)
freq = S_data[0]
S_11_amp = S_data[1]
S_11_phase = S_data[2]

# import realized gain file
G_data = np.loadtxt('Realized_antenna_gain.csv', delimiter = ',', skiprows = 1)
G_data = np.transpose(G_data)
# only import theta polarization (main pol of interest polarization)
G_amp_ang1 = G_data[1]
G_phase_ang1 = G_data[2]
G_amp_ang2 = G_data[5]
G_phase_ang2 = G_data[6]
G_amp_ang3 = G_data[9]
G_phase_ang3 = G_data[10]

plt.figure()
plt.plot(freq, S_11_amp, color = 'k', label = 'S11 amplitude')
#plt.plot(freq, S_11_phase, label = 'S11 phase')
plt.grid()
plt.legend()
plt.xlim(2.5,22)
plt.title('S11')
plt.ylabel('Amplitude (dB)')
plt.xlabel('Frequency (GHz)')
plt.savefig('S_parameter')

plt.figure()
plt.plot(freq, G_amp_ang1, color = 'r', label = 'Theta Polarization: Theta=90, Phi=0, Amplitude')
#plt.plot(freq, G_phase_ang2, color = 'r', label = 'Theta Polarization: Theta=90, Phi=0, Phase')
plt.plot(freq, G_amp_ang2, 'x', markersize='5', markevery=50, color = 'b', label = 'Theta Polarization: Theta=90, Phi=180, Amplitude')
#plt.plot(freq, G_phase_ang2, 'x', markersize='5', markevery=50, color = 'b', label = 'Theta Polarization: Theta=90, Phi=180, Phase')
plt.plot(freq, G_amp_ang3, color = 'g', label = 'Theta Polarization: Theta=45, Phi=0, Amplitude')
#plt.plot(freq, G_phase_ang3, color = 'g', label = 'Theta Polarization: Theta=45, Phi=0, Phase')
plt.ylabel('Amplitude (dB)')
plt.xlabel('Frequency (GHz)')
plt.xlim(2.5,22)
plt.ylim(-20,5)
plt.title('Realized Gain')
plt.grid()
plt.legend()
plt.savefig('Realized_gain_plot')

