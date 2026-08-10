import os
import sys
import shutil
import numpy as np

def generate_ky_kz_arrays(n_theta, n_phi, k_magnitude):
    # Create uniform spacing in angles
    theta = np.linspace(np.pi/6, np.pi/6+np.pi*2/3, n_theta)  # 30 to 150 degrees
    phi = np.linspace(-np.pi/3, np.pi/3, n_phi)  # -60 to 60 degrees
    # Create meshgrid
    THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
    # Convert to ky, kz (with theta from Z axis)
    ky = k_magnitude * np.sin(THETA) * np.sin(PHI)
    kz = k_magnitude * np.cos(THETA)
    # Flatten for indexing
    ky_flat = ky.flatten()
    kz_flat = kz.flatten()
    return ky_flat, kz_flat

# Generate arrays via definition (use this for angle vs angle heat maps at single frequency):
# Target frequency
#target_freq = 9.0  # GHz
# set num angles for each
#n_theta = 15
#n_phi = 15
# total length of ky and kz will be n_theta*n_phi =  num of simulations =  len ky and kz
#ky, kz = generate_ky_kz_arrays(n_theta=n_theta, n_phi=n_phi, k_magnitude=2*np.pi*target_freq*1E9/3E8)

# other method for direct control of k values - more useful for angle vs frequency heat maps at a single other angle:
index_length = 400
ky = np.linspace(0, 0, index_length)
kz = np.linspace(0, 399, index_length)

current_index_str = sys.argv[1]
index = int(current_index_str)
this_ky = round(ky[index], 2)
this_kz = round(kz[index], 2)

sim_dir = f'working/sim_{this_ky}_{this_kz}'

# Create directory if it doesn't exist
os.makedirs(sim_dir, exist_ok=True)

# Copy all needed files to the directory
for file in os.listdir('.'):
    if os.path.isfile(file):
        shutil.copy(file, sim_dir)

# Change to this directory
os.chdir(sim_dir)

#run all scripts needed here
os.system('python master.py {} {}'.format(this_ky,this_kz))
os.system('python master_clear.py {} {}'.format(this_ky,this_kz))
os.system('python post_processor_kmax.py {} {}'.format(this_ky,this_kz))

#os.system('python convert_to_touchstone.py {} {}'.format(this_ky,this_kz))

#then copy the touchstone only to the touchstone folder for convenience
#os.chdir('../..')
#os.system('cp {}/my_sparams_MA_{}_{}.s4p ./touchstones'.format(sim_dir,this_ky,this_kz))
