import numpy as np
from pyevtk.hl import imageToVTK
import os
from scipy.io import FortranFile

# LICENSE FILE IS INCLUDED IN THE PARENT FOLDER OF THIS FILE.
# This will load the E,H cell centered, leave them staggered in time by 1/2 time step, then put fields into numpy arrays.
# It will then make a paraview object for easy viewing of the fields.
# Open paraview first w/ option to load all files for E,H simultanteously in time/freq.
# Note that the interpolation in Paraview for data can be very odd, be careful when looking at fields as they will be interpolated and might appear assymetric/irregular at locations.
# Ex. when plotting a slice, the squares across the axis might not be evenly spaced, so symmetric fields might appear assymetric.

# move to the current directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

#SETUP
#################################################################################
# Simulation parameters - input your simulation parameters, will error out if the sizes are wrong.
Nx = 65
Ny = 30
Nz = 30
del_t = 9.4365835E-13
time_steps = 500
# Use H fields, and if so, what filenames
Use_H=True
Hx_field_name = "data_full_Hx.bin"
Hy_field_name = "data_full_Hy.bin"
Hz_field_name = "data_full_Hz.bin"
# Use E fields, and if so, what filenames
Use_E=True
Ex_field_name = "data_full_Ex.bin"
Ey_field_name = "data_full_Ey.bin"
Ez_field_name = "data_full_Ez.bin"
# Output file info
output_folder = "./Output"
# not setup yet but kmax or not kmax for data types that come in and what we do with them
kmax=False
# Select info on still frame, freq, time, etc.
time_or_freq="freq" #'time' or 'freq'
# Select time or freq frame, [start,stop] frames, both won't be used (only time or freq used, not both)
time_frame=[0,99] # Ex. if you wanted the first 100 frames you would use [0,99]
freq_frame=[0,19] # At end of program is will print what these frequencies are in GHz for reference
##################################################################################

# Output folder for VTK files
os.makedirs(output_folder, exist_ok=True)

# some helper definitions
def read_fortran_records(file_path, shape):
    frames = []
    with FortranFile(file_path, 'r') as f:
        while True:
            try:
                if kmax==False:
                    rec = f.read_reals(dtype=np.float32)
                if kmax==True:
                    rec = f.read_reals(dtype=np.complex64)
            except Exception:
                break
            frames.append(rec.reshape(shape, order='F'))
    return np.array(frames)

def load_field_data(field_name):
    file_path = '{}'.format(field_name)
    return read_fortran_records(file_path, (Nx,Ny,Nz))

# now import and process the fields - E,H as desired

if Use_H==True:
    # import h fields and process
    hx = load_field_data(Hx_field_name)
    hy = load_field_data(Hy_field_name)
    hz = load_field_data(Hz_field_name)

    if time_or_freq=='time':
        # Export h field w/ time steps for ParaView
        for t in range(time_frame[0],time_frame[1]+1):

            hx_t = hx[t, :, :, :]
            hy_t = hy[t, :, :, :]
            hz_t = hz[t, :, :, :]

            # Save each timestep with a padded index format
            file_path = os.path.join(output_folder, f"H_Field_{t:03d}")

            imageToVTK(
                file_path,
                spacing=(1.0, 1.0, 1.0),
                origin=(0.5, 0.5, 0.5),
                cellData={"H_Field": (hx_t, hy_t, hz_t)},
            )

    if time_or_freq=='freq':
        #fft the data
        hx_f = np.fft.fft(hx,axis=0)
        hy_f = np.fft.fft(hy,axis=0)
        hz_f = np.fft.fft(hz,axis=0)
        freq = np.fft.fftfreq(time_steps,del_t)
        # Export h field at freqs desired for ParaView
        for t in range(freq_frame[0],freq_frame[1]+1):

            hx_t = np.abs(hx_f[t, :, :, :])
            hy_t = np.abs(hy_f[t, :, :, :])
            hz_t = np.abs(hz_f[t, :, :, :])

            # Save each timestep with a padded index format
            file_path = os.path.join(output_folder, f"H_Field_{t:03d}")

            imageToVTK(
                file_path,
                spacing=(1.0, 1.0, 1.0),
                origin=(0.5, 0.5, 0.5),
                cellData={"H_Field": (hx_t, hy_t, hz_t)},
            )


if Use_E==True:
    # import e fields and process
    ex = load_field_data(Ex_field_name)
    ey = load_field_data(Ey_field_name)
    ez = load_field_data(Ez_field_name)

    if time_or_freq=='time':
        # Export e field w/ time steps for ParaView
        for t in range(time_frame[0],time_frame[1]+1):

            ex_t = ex[t, :, :, :]
            ey_t = ey[t, :, :, :]
            ez_t = ez[t, :, :, :]

            # Save each timestep with a padded index format
            file_path = os.path.join(output_folder, f"E_Field_{t:03d}")

            imageToVTK(
                file_path,
                spacing=(1.0, 1.0, 1.0),
                origin=(0.5, 0.5, 0.5),
                cellData={"E_Field": (ex_t, ey_t, ez_t)},
            )

    if time_or_freq=='freq':
        #fft the data
        ex_f = np.fft.fft(ex,axis=0)
        ey_f = np.fft.fft(ey,axis=0)
        ez_f = np.fft.fft(ez,axis=0)
        freq = np.fft.fftfreq(time_steps,del_t)
        # Export h field at freqs desired for ParaView
        for t in range(freq_frame[0],freq_frame[1]+1):

            ex_t = np.abs(ex_f[t, :, :, :])
            ey_t = np.abs(ey_f[t, :, :, :])
            ez_t = np.abs(ez_f[t, :, :, :])

            # Save each timestep with a padded index format
            file_path = os.path.join(output_folder, f"E_Field_{t:03d}")

            imageToVTK(
                file_path,
                spacing=(1.0, 1.0, 1.0),
                origin=(0.5, 0.5, 0.5),
                cellData={"E_Field": (ex_t, ey_t, ez_t)},
            )

if time_or_freq=='freq':
    print('User request frequencies (GHz) ', freq[freq_frame[0]]/1E9, ' through ', freq[freq_frame[1]]/1E9)