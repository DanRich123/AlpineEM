import numpy as np
from matplotlib import pyplot as plt

# load the non foster aperture area data
data = np.loadtxt('Aperture_area.csv', delimiter=',', skiprows=1)
data = np.transpose(data)

# load the aperture area where the non foster element is not used
data2 = np.loadtxt('Aperture_area_50.csv', delimiter=',', skiprows=1)
data2 = np.transpose(data2)

# output is A_e * 50 in a log scale since the spice load of interest here is a simple 50 ohm load - true for both simulations
# need to divide by 50 to get A_e but since in log scale, we subtract a log

A_e_corrected_mag = data[3] - 10*np.log10(50)
A_e_corrected_mag2 = data2[3] - 10*np.log10(50)

plt.figure()
plt.plot(data[0]*1000, A_e_corrected_mag, label = 'Aperture Area - NF Circuit w/ 50 ohm load')
plt.plot(data2[0]*1000, A_e_corrected_mag2, label = 'Aperture Area - 50 ohm load only')
plt.xlabel('Frequency (MHz)')
plt.ylabel('Amplitude (dB)')
plt.ylim(-90,-55)
plt.grid()
plt.title('Aperture Area Comparison')
plt.savefig('Aperture area comparison')

lam = 3E8 / (data[0] * 1E9) #same for data2 as well

G_r = A_e_corrected_mag + 10 * np.log10(4*np.pi/(lam**2))
G_r2 = A_e_corrected_mag2 + 10 * np.log10(4*np.pi/(lam**2))

plt.figure()
plt.plot(data[0]*1000, G_r, label = 'Realized Gain - NF Circuit w/ 50 ohm load')
plt.plot(data2[0]*1000, G_r2, label = 'Realized Gain - 50 ohm load only')
plt.xlabel('Frequency (MHz)')
plt.ylabel('Amplitude (dB)')
plt.ylim(-100,-35)
plt.grid()
plt.title('Realized Gain Comparison')
plt.savefig('Realized gain comparison')