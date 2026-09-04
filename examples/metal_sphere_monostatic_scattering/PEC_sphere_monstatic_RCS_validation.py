import numpy as np
from scipy.special import spherical_jn,spherical_yn,hankel2
#from mpmath import hankel2
import matplotlib.pyplot as plt
import os

# change working directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# setup analytic calcultions
radius=12.5E-3 #meters
freq=np.linspace(0.7E9,50E9,1000)
c=299792458
lam=c/freq
k=2*np.pi/lam
RCS_PEC_Sphere=np.zeros(len(freq), dtype=np.complex128)
interation_number=50 #too many causes demoninator to diverge numerically and become unstable - 1/infinity. 50 should be enough.
z=k*radius
for i in range(interation_number):
    h2 = np.sqrt(np.pi*z/2) * hankel2(i+1 + 0.5, z)
    h2_nminus1 = np.sqrt(np.pi*z/2) * hankel2(i + 0.5, z)
    h2_nplus1 = np.sqrt(np.pi*z/2) * hankel2(i+2 + 0.5, z)
    h2_prime = 1/2*(h2_nminus1-(h2_nplus1)+h2/z)
    test=np.array(1/(h2*h2_prime),dtype=np.complex128)
    RCS_PEC_Sphere=RCS_PEC_Sphere+((-1)**(i+1))*(2*(i+1)+1)*test
RCS_PEC_Sphere=(lam**2/(4*np.pi))*np.abs(RCS_PEC_Sphere)**2

#import and noramlize/process the fdtd data
FDTD=np.loadtxt('Scattering_Far_Field.csv', delimiter=',', skiprows=1)
FDTD=np.transpose(FDTD)
lam_FDTD=(c/(FDTD[0]*1E9))
fdtd_rcs=10**(FDTD[1]/10)/np.pi/radius**2

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))

# Normalized RCS vs Circumference Per Wavelength
ax1.plot(2*np.pi*radius/lam_FDTD, fdtd_rcs, 'x', label='FDTD Monostatic RCS')
ax1.plot(2*np.pi*radius/lam, RCS_PEC_Sphere/np.pi/radius**2, color='k', label='Analytic Monostatic RCS')
ax1.set_xlabel('Circumference Per Wavelength')
ax1.set_ylabel('Normalized RCS (Projected Area Asymptote)')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_title('PEC Sphere of Radius {} Meters'.format(radius))
ax1.grid(True)
ax1.legend()

# RCS (dB) vs Frequency (GHz)
ax2.plot(FDTD[0], FDTD[1], 'x', label='FDTD Monostatic RCS')
ax2.plot(c/lam/1E9, 10*np.log10(RCS_PEC_Sphere), color='k', label='Analytic Monostatic RCS')
ax2.set_xlabel('Frequency (GHz)')
ax2.set_ylabel('RCS (dB)')
ax2.set_xscale('log')
ax2.grid(True)
ax2.set_title('PEC Sphere of Radius {} Meters'.format(radius))
ax2.legend()

plt.tight_layout()
plt.show()