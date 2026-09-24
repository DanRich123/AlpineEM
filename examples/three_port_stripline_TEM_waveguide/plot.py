import numpy as np
from matplotlib import pyplot as plt
import os 

# move to working directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# load the non foster aperture area data
data = np.loadtxt('Aperture_area_NF.csv', delimiter=',', skiprows=1)
data = np.transpose(data)

# load the aperture area where the non foster element is not used
data2 = np.loadtxt('Aperture_area.csv', delimiter=',', skiprows=1)
data2 = np.transpose(data2)

data3 = np.loadtxt('Aperture_area_3mm_NF.csv', delimiter=',', skiprows=1)
data3 = np.transpose(data3)

data4 = np.loadtxt('Aperture_area_3mm_50.csv', delimiter=',', skiprows=1)
data4 = np.transpose(data4)

# output is A_e * 50 in a log scale since the spice load of interest here is a simple 50 ohm load - true for both simulations
# need to divide by 50 to get A_e but since in log scale, we subtract a log
A_e_corrected_mag = data[3] - 10*np.log10(50)
A_e_corrected_mag2 = data2[3] - 10*np.log10(50)
A_e_corrected_mag3 = data3[3] - 10*np.log10(50)
A_e_corrected_mag4 = data4[3] - 10*np.log10(50)
lam = 3E8 / (data[0] * 1E9) #same for data2 as well
lam2 = 3E8 / (data2[0] * 1E9) #same for data2 as well
lam3 = 3E8 / (data3[0] * 1E9) #same for data2 as well
lam4 = 3E8 / (data4[0] * 1E9) #same for data2 as well

G_r = A_e_corrected_mag + 10 * np.log10(4*np.pi/(lam**2))
G_r2 = A_e_corrected_mag2 + 10 * np.log10(4*np.pi/(lam2**2))
G_r3 = A_e_corrected_mag3 + 10 * np.log10(4*np.pi/(lam3**2))
G_r4 = A_e_corrected_mag4 + 10 * np.log10(4*np.pi/(lam4**2))

test_50=np.loadtxt('S_parameters_50.csv', delimiter=',', skiprows=1)
test_50=np.transpose(test_50)

test_3=np.loadtxt('S_parameters_3mm_NF.csv', delimiter=',', skiprows=1)
test_3=np.transpose(test_3)

test_4=np.loadtxt('S_parameters_3mm_50.csv', delimiter=',', skiprows=1)
test_4=np.transpose(test_4)

test_nf=np.loadtxt('S_parameters_NF.csv', delimiter=',', skiprows=1)
test_nf=np.transpose(test_nf)

freq_test_50=test_50[0]
lam_test_50 = 3E8 / (freq_test_50 * 1E9)
data_test_50=test_50[7] - 10 * np.log10(50)

freq_test_nf=test_nf[0]
lam_test_nf = 3E8 / (freq_test_nf * 1E9)
data_test_nf=test_nf[7] - 10 * np.log10(50)

freq_test_3=test_3[0]
lam_test_3 = 3E8 / (freq_test_3 * 1E9)
data_test_3=test_3[7] - 10 * np.log10(50)

freq_test_4=test_4[0]
lam_test_4 = 3E8 / (freq_test_4 * 1E9)
data_test_4=test_4[7] - 10 * np.log10(50)

A1=0.0155
A2=0.0155

#note value for mode area chosen based on 50 ohm case that matched well, I didn't calibrate it the more accurate way.
#just a proof of concept so it should be fine.
# need one more correction for C(k) - the same in both scenarios (50 and NF)
G_r_test_50 = data_test_50 + 10*np.log10(4/np.pi*A1*(2*np.pi/lam_test_50)**2) + test_50[3]/2 
G_r_test_3 = data_test_3 + 10*np.log10(4/np.pi*A2*(2*np.pi/lam_test_3)**2) + test_3[3]/2 
G_r_test_nf = data_test_nf + 10*np.log10(4/np.pi*A1*(2*np.pi/lam_test_nf)**2) + test_nf[3]/2
G_r_test_4 = data_test_4 + 10*np.log10(4/np.pi*A2*(2*np.pi/lam_test_4)**2) + test_4[3]/2 

plt.plot(data[0]*1000, G_r, 'x',color ='r', label = 'NF Circuit w/ 50 ohm load - free space')
plt.plot(freq_test_nf*1000, G_r_test_nf, '--', color='r', label = 'NF Circuit w/ 50 ohm load - waveguide extraction')
plt.plot(data2[0]*1000, G_r2, 'x',color='b' , label = '50 ohm load only - free space')
plt.plot(freq_test_50*1000, G_r_test_50, '--', color='b', label = '50 ohm load only - Waveguide extraction')
#plt.plot(data3[0]*1000, G_r3, 'x',color='g' , label = '3mm NF - free space')
#plt.plot(freq_test_3*1000, G_r_test_3, '--', color='g', label = '3mm NF - Waveguide extraction')
#plt.plot(data4[0]*1000, G_r4, 'x',color='k' , label = '3mm 50 - free space')
#plt.plot(freq_test_4*1000, G_r_test_4, '--', color='k', label = '3mm 50 - Waveguide extraction')
plt.xlabel('Frequency (MHz)')
plt.ylabel('Amplitude (dB)')
plt.ylim(-120,-35)
plt.xlim(0,320)
plt.grid()
plt.legend()
plt.title('Realized Gain Comparison')
plt.show()
