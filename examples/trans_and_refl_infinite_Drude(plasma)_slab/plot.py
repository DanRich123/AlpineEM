import numpy as np
from matplotlib import pyplot as plt

# import analytic solution
analytic_s11 = np.load('s11.npy')
analytic_s21 = np.load('s21.npy')
analytic_freq = np.load('freq.npy')

# import fdtd data
data = np.loadtxt('S_parameters.csv', delimiter=',', skiprows=1)
data = np.transpose(data)
#data_mp = np.loadtxt('S_parameters_MP.csv', delimiter=',', skiprows=1)
#data_mp = np.transpose(data_mp)
#data_acc = np.loadtxt('S_parameters_ACC.csv', delimiter=',', skiprows=1)
#data_acc = np.transpose(data_acc)

# make figures
plt.figure(figsize=(6, 4))

# Left plot: S11 comparison
plt.plot(data[0], data[1], 'x', markersize='2', label='FDTD S11')
#plt.plot(data_mp[0], data_mp[1], label='FDTD S11 - MP')
#plt.plot(data_acc[0], data_acc[1], label='FDTD S11 - ACC')
plt.plot(analytic_freq, analytic_s11, color='k', label='Analytic S11')
plt.grid()
plt.legend()
plt.xlabel('Frequency (GHz)')
plt.ylabel('Amplitude (dB)')

# Right plot: S21 comparison
plt.plot(data[0], data[3], 'x', markersize='2', label='FDTD S21')
#plt.plot(data_mp[0], data_mp[3], label='FDTD S21 - MP')
#plt.plot(data_acc[0], data_acc[3], label='FDTD S21 - ACC')
plt.plot(analytic_freq, analytic_s21, color='k', label='Analytic S21')
plt.grid()
plt.legend()
plt.xlabel('Frequency (GHz)')
plt.ylabel('Amplitude (dB)')
plt.title('S-parameter Comparison')


plt.tight_layout()
plt.savefig('Comparison_plot')