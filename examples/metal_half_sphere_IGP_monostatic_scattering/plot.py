import numpy as np
import os
from matplotlib import pyplot as plt

# change working directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# import the analytic/theoretical data saved in a NumPy array
analytic_RCS = np.load('analytic_RCS.npy')
analytic_freq = np.load('analytic_freq.npy')

# import fdtd data
fdtd_data = np.loadtxt('Scattering_Far_Field.csv', delimiter=',', skiprows=1)
fdtd_data = np.transpose(fdtd_data)

plt.figure()
plt.plot(fdtd_data[0], fdtd_data[1], 'x', label = 'FDTD')
plt.plot(analytic_freq, analytic_RCS, color='k', label = 'Analytic')
plt.ylabel('Amplitude (dB)')
plt.xlabel('Frequency (GHz)')
plt.title('RCS Comparison - FDTD vs. Analytic')
plt.tight_layout()
plt.legend()
plt.xscale('log')
plt.grid()
#plt.show()
plt.savefig('Comparison_plot')
