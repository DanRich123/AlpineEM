import numpy as np
import skrf as rf # Import scikit-rf for Network handling
from typing import List
import os

#os.chdir('/Users/dari6475/Desktop/patchgrid/me_test')

###SETUP#####################################################################
#this script takes in the S parameter data from fdtd simulations and formats them for Ryan's algorithm into a touchstone file
#currently assumes 2 waveports and 2 lumped ports and formats them as he desires
#it's also currenlty formatted for x normal but that can be adjusted easily as well
#so this is really a custom script
S_TE_name='S_parameters_TE.csv'
S_TM_name='S_parameters_TM.csv'
S_lumped1_name='S_parameters_antenna1.csv'
S_lumped2_name='S_parameters_antenna2.csv'
#os.chdir('/Users/dari6475/Desktop/patchgrid/me_test')
#name of file from fdtd with output data for pulling simulation info
data_filename='data.dat'
#name of the file you want to name the touchstone file
output_file_name='my_sparams_MA.s4p'
#############################################################################

cells_size_y=int(np.loadtxt(data_filename,usecols=3,skiprows=1,max_rows=1)) #cells
cells_size_z=int(np.loadtxt(data_filename,usecols=3,skiprows=2,max_rows=1)) #cells
step_size_y=np.loadtxt(data_filename,usecols=4,skiprows=3,max_rows=1) #meters
step_size_z=np.loadtxt(data_filename,usecols=5, skiprows=3,max_rows=1) #meters
size_y=cells_size_y*step_size_y #meters
size_z=cells_size_z*step_size_z #meters
ky_var=np.loadtxt(data_filename,usecols=0,skiprows=17,max_rows=1) #rad/meter
kz_var=np.loadtxt(data_filename,usecols=1,skiprows=17,max_rows=1) #rad/meter

data_TE=np.loadtxt(S_TE_name, skiprows=1, delimiter=',')
data_TE=np.transpose(data_TE)
data_TM=np.loadtxt(S_TM_name, skiprows=1, delimiter=',')
data_TM=np.transpose(data_TM)
data_ant1=np.loadtxt(S_lumped1_name, skiprows=1, delimiter=',')
data_ant1=np.transpose(data_ant1)
data_ant2=np.loadtxt(S_lumped2_name, skiprows=1, delimiter=',')
data_ant2=np.transpose(data_ant2)

freq = data_TE[0]*1E9
s = np.zeros((len(freq),4,4), dtype=complex)

#impedance of ports w/ incident angle
k=2*np.pi*freq/3E8
if (ky_var==0) and (kz_var==0):
    cos_angle=1.0+0*freq
if (ky_var!=0) or (kz_var!=0):
    cos_angle=np.sqrt(k**2-ky_var**2-kz_var**2)/(k)
z0_array = np.column_stack([
    np.full_like(freq, 50.0),
    np.full_like(freq, 50.0),
    376.3 * cos_angle,
    376.3 / cos_angle
])

#TE incident data is S3
s[:,2,2]=10**(data_TE[1]/20)*np.exp(1j*data_TE[2])
s[:,3,2]=10**(data_TE[3]/20)*np.exp(1j*data_TE[4])
s[:,0,2]=10**(data_TE[9]/20)*np.exp(1j*data_TE[10])
s[:,1,2]=10**(data_TE[11]/20)*np.exp(1j*data_TE[12])
#TM incident data is S4
s[:,2,3]=10**(data_TM[1]/20)*np.exp(1j*data_TM[2])
s[:,3,3]=10**(data_TM[3]/20)*np.exp(1j*data_TM[4])
s[:,0,3]=10**(data_TM[9]/20)*np.exp(1j*data_TM[10])
s[:,1,3]=10**(data_TM[11]/20)*np.exp(1j*data_TM[12])
#port 1 incident data is S1
s[:,2,0]=10**(data_ant1[1]/20)*np.exp(1j*data_ant1[2])
s[:,3,0]=10**(data_ant1[3]/20)*np.exp(1j*data_ant1[4])
s[:,0,0]=10**(data_ant1[9]/20)*np.exp(1j*data_ant1[10])
s[:,1,0]=10**(data_ant1[11]/20)*np.exp(1j*data_ant1[12])
#port 2 incident data is S2
s[:,2,1]=10**(data_ant2[1]/20)*np.exp(1j*data_ant2[2])
s[:,3,1]=10**(data_ant2[3]/20)*np.exp(1j*data_ant2[4])
s[:,0,1]=10**(data_ant2[9]/20)*np.exp(1j*data_ant2[10])
s[:,1,1]=10**(data_ant2[11]/20)*np.exp(1j*data_ant2[12])


def write_touchstone_ma_matrix(filename, freq, s, comments=None):

    n_freqs, n_ports, _ = s.shape
    with open(filename, 'w') as f:
        # Optional comment lines
        if comments:
            for line in comments:
                f.write(f"! {line}\n")

        # Touchstone header
        f.write(f"# GHZ S MA 50\n")

        # Write data in matrix-style blocks
        for i in range(n_freqs):
            f.write(f"{freq[i]/1e9:.6f}\n")
            for m in range(n_ports):
                for n in range(n_ports):
                    mag = (np.abs(s[i, m, n]))
                    ang = np.angle(s[i, m, n], deg=True)
                    f.write(f"  {mag:.6e} {ang:.6f}")
                f.write("\n")  # new row per port

# Create scikit-rf Network
ntwk = rf.Network()
ntwk.frequency = rf.Frequency.from_f(freq, unit='hz')
ntwk.s = s

# Impedance per port, per frequency (important!)
#imped=[50,50,377,377]
#z0_array = np.tile(imped, (len(freq), 1))  # shape (n_freqs, n_ports)

ntwk.z0 = z0_array

# In-place renormalize
ntwk.renormalize(50)

# Set name and extract renormalized S-matrix
ntwk.name = output_file_name
s = ntwk.s  # IMPORTANT: now using renormalized S-params

comments = [
    "Format: Magnitude/Angle in degrees",
    "ky={}".format(ky_var),
    "dy={}".format(size_y),
    "kydy={}".format(ky_var*size_y),
    "kz={}".format(kz_var),
    "dz={}".format(size_z),
    "kzdz={}".format(kz_var*size_z),
    "ky/kz in rad/meter and dy/dz in meters",
    "1 - lumped port 1",
    "2 - lumped port 2",
    "3 - TE port",
    "4 - TM port",
    "Lastly, note that each ky,kz combination has a different frequency range"
]

write_touchstone_ma_matrix(output_file_name, freq, s, comments)