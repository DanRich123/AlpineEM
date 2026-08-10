import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree

#this script imports data provided by Nick from CFD
#it exepcts a specific format and imports, then interpolates, and can create an optional geometry input file
#it also outputs a npy file with material properties - there is an example script showing how to use these materials

##################################
#SETUP############################
##################################
file_path='outputCC_Park_71km.plt'
#del_r should be set to del_x,del_y value for fdtd
del_z, del_r = 0.01, 0.01 
num_materials = 15
#if use all bins the num_materials will be overwritten below with the actual number of max IDs
use_all_bins=True
#can fill interior with radome if we want to, otherwise it just lines the interior surface
fill_interior=False
num_layers=1
#option to turn on making the fdtd output txt file we need - takes a while so can turn off if testing
fdtd_prep=True
#old version swapped x and z so use True to swap them.
swap_axis=True 
#do we want to use all of the species in addition to the electrons
use_all_materials=True
# Define your z range for interior radome (adjust these values as needed) for placement
z_min_boundary = 0.02  # minimum z value
z_max_boundary = 1.28  # maximum z value
#set quarter = true for quarter or false for full
#quarter is nice for visualizion what's going on in 3D paraview, hard to really see the full version patterns well
quarter=False
#if we want to round down always we set this to true. Otherwise it will round normally for 2D to 3D mapping using radius.
round_down=False
##################################
#END SETUP########################
##################################

##################################
#READ IN DATA AND INTEPROLATE#####
##################################
all_records = []
# Using lists to store everything zone-by-zone to ensure matching lengths
z_coords_all = []
r_coords_all = []
indices_all = []
rho_e_all_faces = [] # This will store rho_e values duplicated for each triangle
#this will do the same for ions and temperatures if we select all species
if use_all_materials==True:
    rho_N2_all_faces  = []
    rho_O2_all_faces  = []
    rho_NO_all_faces  = []
    rho_N_all_faces   = []
    rho_O_all_faces   = []
    rho_NOp_all_faces = []
    rho_N2p_all_faces = []
    rho_O2p_all_faces = []
    rho_Np_all_faces  = []
    rho_Op_all_faces  = []
    temperature_electrons_all_faces = []
    temperature_ions_all_faces = []
node_offset = 0

# --- 1. Parse File ---
with open(file_path, 'r') as f:
    header = f.readline()
    variables = re.findall(r'"([^"]*)"', header)
    
    def get_tokens(file_handle):
        for line in file_handle:
            line_strip = line.strip()
            if line_strip.startswith("ZONE"):
                yield "NEW_ZONE"; yield line; continue
            for word in line.split(): yield word
    
    tokens = get_tokens(f)

    while True:
        try:
            token = next(tokens)
            if token == "NEW_ZONE": zone_line = next(tokens)
            elif "ZONE" in token: zone_line = token
            else: continue

            N = int(re.search(r'N=(\d+)', zone_line).group(1))
            E = int(re.search(r'E=(\d+)', zone_line).group(1))
            
            zone_data = {}
            for i, var_name in enumerate(variables):
                count = N if (i+1 <= 2 or i+1 == 231) else E
                zone_data[var_name] = [float(next(tokens)) for _ in range(count)]

            z_coords = zone_data['Z[m]']
            r_coords = zone_data['r[m]']
            rho_e_vals = zone_data['rho_e']
            if use_all_materials==True:
                #first ions
                rho_NOp_vals = zone_data['rho_NOp']
                rho_N2p_vals = zone_data['rho_N2p']
                rho_O2p_vals = zone_data['rho_O2p']
                rho_Np_vals = zone_data['rho_Np']
                rho_Op_vals = zone_data['rho_Op']
                #then neutrals
                rho_NO_vals = zone_data['rho_NO']
                rho_N2_vals = zone_data['rho_N2']
                rho_O2_vals = zone_data['rho_O2']
                rho_N_vals = zone_data['rho_N']
                rho_O_vals = zone_data['rho_O']
                #the temperatures
                temperature_electrons_vals = zone_data['Tv']
                temperature_ions_vals = zone_data['T']

            z_coords_all.extend(z_coords)
            r_coords_all.extend(r_coords)

            for i in range(E):
                nodes = [int(next(tokens)) for _ in range(4)]
                idx0 = [n - 1 for n in nodes]
                
                # Centroid for interpolation
                z_avg = sum(z_coords[idx] for idx in idx0) / 4
                r_avg = sum(r_coords[idx] for idx in idx0) / 4
                cell_size = max([np.sqrt((z_coords[idx]-z_avg)**2 + (r_coords[idx]-r_avg)**2) for idx in idx0])
                
                if use_all_materials==False:
                    all_records.append({'Z': z_avg, 'R': r_avg, 'rho_e': rho_e_vals[i], 'size': cell_size})
                if use_all_materials==True:
                    all_records.append({
                        #first all old data - to keep same format as without all materials
                        'Z': z_avg, 
                        'R': r_avg, 
                        'rho_e': rho_e_vals[i], 
                        'size': cell_size,
                        #now do the ions
                        'rho_NOp': rho_NOp_vals[i],
                        'rho_N2p': rho_N2p_vals[i],
                        'rho_O2p': rho_O2p_vals[i],
                        'rho_Np': rho_Np_vals[i],
                        'rho_Op': rho_Op_vals[i],
                        #now do the neutral ones
                        'rho_NO': rho_NO_vals[i],
                        'rho_N2': rho_N2_vals[i],
                        'rho_O2': rho_O2_vals[i],
                        'rho_N': rho_N_vals[i],
                        'rho_O': rho_O_vals[i],
                        #and finally the temperatures
                        'Tv': temperature_electrons_vals[i],
                        'T': temperature_ions_vals[i]
                    })

                # Connectivity: Split quad into two triangles
                abs_idx = [n - 1 + node_offset for n in nodes]
                indices_all.append([abs_idx[0], abs_idx[1], abs_idx[2]])
                indices_all.append([abs_idx[0], abs_idx[2], abs_idx[3]])
                
                # For every 1 quad, we have 2 triangles. Both need the same rho_e color.
                rho_e_all_faces.extend([rho_e_vals[i], rho_e_vals[i]])
                #same for others if we use all species/materials
                if use_all_materials==True:
                    #first the ions
                    rho_NOp_all_faces.extend([rho_NOp_vals[i], rho_NOp_vals[i]])
                    rho_N2p_all_faces.extend([rho_N2p_vals[i], rho_N2p_vals[i]])
                    rho_O2p_all_faces.extend([rho_O2p_vals[i], rho_O2p_vals[i]])
                    rho_Np_all_faces.extend([rho_Np_vals[i], rho_Np_vals[i]])
                    rho_Op_all_faces.extend([rho_Op_vals[i], rho_Op_vals[i]])
                    #then the neutral ones
                    rho_NO_all_faces.extend([rho_NO_vals[i], rho_NO_vals[i]])
                    rho_N2_all_faces.extend([rho_N2_vals[i], rho_N2_vals[i]])
                    rho_O2_all_faces.extend([rho_O2_vals[i], rho_O2_vals[i]])
                    rho_N_all_faces.extend([rho_N_vals[i], rho_N_vals[i]])
                    rho_O_all_faces.extend([rho_O_vals[i], rho_O_vals[i]])
                    #and finally the temperatures
                    temperature_electrons_all_faces.extend([temperature_electrons_vals[i], temperature_electrons_vals[i]])
                    temperature_ions_all_faces.extend([temperature_ions_vals[i], temperature_ions_vals[i]])
            node_offset += N
        except StopIteration: break

# --- 2. Grid and Material Settings ---
z_arr, r_arr = np.array(z_coords_all), np.array(r_coords_all)
indices_arr = np.array(indices_all)
face_colors_arr = np.array(rho_e_all_faces)
#do the same as electrons to others if needed:
if use_all_materials==True:
    #first the ions
    face_colors_NOp = np.array(rho_NOp_all_faces)
    face_colors_N2p = np.array(rho_N2p_all_faces)
    face_colors_O2p = np.array(rho_O2p_all_faces)
    face_colors_Np = np.array(rho_Np_all_faces)
    face_colors_Op = np.array(rho_Op_all_faces)
    #then the neutral ones
    face_colors_NO = np.array(rho_NO_all_faces)
    face_colors_N2 = np.array(rho_N2_all_faces)
    face_colors_O2 = np.array(rho_O2_all_faces)
    face_colors_N = np.array(rho_N_all_faces)
    face_colors_O = np.array(rho_O_all_faces)
    #and finally the temperatures
    face_colors_temperature_electrons = np.array(temperature_electrons_all_faces)
    face_colors_temperature_ions = np.array(temperature_ions_all_faces)

z_min, z_max = z_arr.min(), z_arr.max()
r_min, r_max = r_arr.min(), r_arr.max()

z_vec = np.arange(z_min, z_max + del_z, del_z)
r_vec = np.arange(r_min, r_max + del_r, del_r)
grid_z, grid_r = np.meshgrid(z_vec, r_vec, indexing='ij')

# --- 3. Interpolation and Masking ---
points_orig = np.array([[rec['Z'], rec['R']] for rec in all_records])
sizes_orig = np.array([rec['size'] for rec in all_records])

values_orig = np.array([rec['rho_e'] for rec in all_records])
interp_func = LinearNDInterpolator(points_orig, values_orig)
grid_rho = interp_func(grid_z, grid_r)

#now we will do the rest of them if needed
if use_all_materials==True:
    species_keys = [
        'rho_NOp', 'rho_N2p', 'rho_O2p', 'rho_Np', 'rho_Op', 
        'rho_NO', 'rho_N2', 'rho_O2', 'rho_N', 'rho_O', 'Tv', 'T'
    ]
    grid_results = {}
    for key in species_keys:
        values_orig_others = np.array([rec[key] for rec in all_records])
        interp_func_others = LinearNDInterpolator(points_orig, values_orig_others)
        grid_results[key] = interp_func_others(grid_z, grid_r)

    #first the ions
    grid_rho_NOp = grid_results['rho_NOp']
    grid_rho_N2p = grid_results['rho_N2p']
    grid_rho_O2p = grid_results['rho_O2p']
    grid_rho_Np  = grid_results['rho_Np']
    grid_rho_Op  = grid_results['rho_Op']
    #then the neutral ones
    grid_rho_NO  = grid_results['rho_NO']
    grid_rho_N2  = grid_results['rho_N2']
    grid_rho_O2  = grid_results['rho_O2']
    grid_rho_N   = grid_results['rho_N']
    grid_rho_O   = grid_results['rho_O']
    #and finally the temperatures
    grid_Tv = grid_results['Tv']
    grid_T = grid_results['T']

tree = cKDTree(points_orig)
dist, idx = tree.query(np.c_[grid_z.ravel(), grid_r.ravel()])
dist = dist.reshape(grid_z.shape)
mask = dist > (sizes_orig[idx].reshape(grid_z.shape) * 1.1)
grid_rho[mask] = np.nan
#now the rest of them if neeed
if use_all_materials==True:
    #first the ions
    grid_rho_NOp[mask] = np.nan
    grid_rho_N2p[mask] = np.nan
    grid_rho_O2p[mask] = np.nan
    grid_rho_Np[mask]  = np.nan
    grid_rho_Op[mask]  = np.nan
    #then the neutral ones
    grid_rho_NO[mask]  = np.nan
    grid_rho_N2[mask]  = np.nan
    grid_rho_O2[mask]  = np.nan
    grid_rho_N[mask]   = np.nan
    grid_rho_O[mask]   = np.nan
    #and finally the temperatures
    grid_Tv[mask] = np.nan
    grid_T[mask]  = np.nan

valid_vals = values_orig[values_orig > 0]
#old content saved for historical reasons - might want again someday
#bins = np.logspace(np.log10(valid_vals.min()), np.log10(valid_vals.max()), num_materials + 1)
#grid_material = np.digitize(grid_rho, bins)
#grid_material = grid_material.astype(float)
#grid_material[np.isnan(grid_rho)] = np.nan
#if use_all_bins==True:
    #num_materials = grid_rho.size
    #print('Max number of IDs based on grid size is set to {}'.format(num_materials))
    ##testing
    #count_test = np.sum((grid_rho > 0) & (~np.isnan(grid_rho)))
    #print(f"Count: {count_test}")
##################################
#END READ IN DATA AND INTEPROLATE#
##################################

##################################
# VISUALIZE TO CONFIRM INTERP ####
##################################

#used in all plots throughout for making them look nice goes here
plt.rcParams.update({
    'axes.titlesize': 25,    # Title size
    'axes.labelsize': 20,    # X and Y label size
    'xtick.labelsize': 18,   # X tick numbers
    'ytick.labelsize': 18,   # Y tick numbers
    'font.size': 18          # Default (affects colorbar labels/numbers)
})

# Added layout='constrained' to prevent colorbar shrinking the plots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 8), sharey=True, layout='constrained')

cmap_cont = plt.colormaps['nipy_spectral'].copy()
cmap_cont.set_bad(color='white', alpha=0)

norm_log = LogNorm(vmin=valid_vals.min(), vmax=valid_vals.max())

# Plot 1: Unstructured
tp = ax1.tripcolor(z_arr, r_arr, indices_arr, 
                   facecolors=face_colors_arr, cmap=cmap_cont, norm=norm_log)
ax1.set_title("1. Original Structured Mesh")
ax1.set_ylabel("R [m]")
ax1.set_aspect('auto') # Forces it to fill the box

# Plot 2: Continuous
im2 = ax2.imshow(grid_rho.T, extent=(z_min, z_max, r_min, r_max), 
                 origin='lower', aspect='auto', cmap=cmap_cont, norm=norm_log)
ax2.set_title(fr"2. Continuous Grid ($\Delta_Z={del_z}$,$\Delta_R={del_r}$)")

# Colorbars (layout='constrained' handles the sizing automatically now)
fig.colorbar(tp, ax=ax1, label='Electron Mass Density')
fig.colorbar(im2, ax=ax2, label='Electron Mass Density')

for ax in [ax1, ax2]: ax.set_xlabel("Z [m]")

# plt.tight_layout() # Removed in favor of constrained_layout above
plt.savefig('rho_e_material_analysis.png', dpi=300)
print("Success! Comparison plot saved with large fonts.")
plt.show()

#if we use all materials then we care about temperature too
#I won't plot the other rho values (mass densities) but I will plot the two temperatures
if use_all_materials==True:
    #now do the same for electron temperature
    # Added layout='constrained' to prevent colorbar shrinking the plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 8), sharey=True, layout='constrained')

    cmap_cont = plt.colormaps['nipy_spectral'].copy()
    cmap_cont.set_bad(color='white', alpha=0)
    norm_log = LogNorm(vmin=np.nanmin(face_colors_temperature_electrons),vmax=np.nanmax(face_colors_temperature_electrons))

    # Plot 1: Unstructured
    tp = ax1.tripcolor(z_arr, r_arr, indices_arr, 
                    facecolors=face_colors_temperature_electrons, cmap=cmap_cont, norm=norm_log)
    ax1.set_title("1. Original Structured Mesh")
    ax1.set_ylabel("R [m]")
    ax1.set_aspect('auto') # Forces it to fill the box

    # Plot 2: Continuous
    im2 = ax2.imshow(grid_Tv.T, extent=(z_min, z_max, r_min, r_max), 
                    origin='lower', aspect='auto', cmap=cmap_cont, norm=norm_log)
    ax2.set_title(fr"2. Continuous Grid ($\Delta_z={del_z}$)")

    # Colorbars (layout='constrained' handles the sizing automatically now)
    fig.colorbar(tp, ax=ax1, label='Electron Temperature (K)')
    fig.colorbar(im2, ax=ax2, label='Electron Temperature (K)')

    for ax in [ax1, ax2]: ax.set_xlabel("Z [m]")

    # plt.tight_layout() # Removed in favor of constrained_layout above
    plt.savefig('Tv_electrons_analysis.png', dpi=300)
    print("Success! Comparison plot saved with large fonts.")
    plt.show()

    #now do the same for ion temperature
    # Added layout='constrained' to prevent colorbar shrinking the plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 8), sharey=True, layout='constrained')

    cmap_cont = plt.colormaps['nipy_spectral'].copy()
    cmap_cont.set_bad(color='white', alpha=0)

    #set to electron max for comparison
    norm_log = LogNorm(vmin=np.nanmin(face_colors_temperature_electrons),vmax=np.nanmax(face_colors_temperature_electrons))

    # Plot 1: Unstructured
    tp = ax1.tripcolor(z_arr, r_arr, indices_arr, 
                    facecolors=face_colors_temperature_ions, cmap=cmap_cont, norm=norm_log)
    ax1.set_title("1. Original Structured Mesh")
    ax1.set_ylabel("R [m]")
    ax1.set_aspect('auto') # Forces it to fill the box

    # Plot 2: Continuous
    im2 = ax2.imshow(grid_T.T, extent=(z_min, z_max, r_min, r_max), 
                    origin='lower', aspect='auto', cmap=cmap_cont, norm=norm_log)
    ax2.set_title(fr"2. Continuous Grid ($\Delta_z={del_z}$)")

    # Colorbars (layout='constrained' handles the sizing automatically now)
    fig.colorbar(tp, ax=ax1, label='Ion &  Neutral Temperature (K)')
    fig.colorbar(im2, ax=ax2, label='Ion & Neutral Temperature (K)')

    for ax in [ax1, ax2]: ax.set_xlabel("Z [m]")

    # plt.tight_layout() # Removed in favor of constrained_layout above
    plt.savefig('T_ions_analysis.png', dpi=300)
    print("Success! Comparison plot saved with large fonts.")
    plt.show()

##################################
#END VISUALIZE TO CONFIRM INTERP##
##################################

##################################
#CREATE 3D########################
##################################
interpolated_3d = {}
print(f"\nRevolving 2D axisymmetric data into 3D Cartesian grid...")

if quarter==True:
    x_size=len(grid_rho.T)
    y_size=len(grid_rho.T)
    z_size=len(grid_rho.T[0])
    x_offset, y_offset, z_offset = 0, 0, 0

    del_x=del_r
    del_y=del_r

    # Create 3D Cartesian coordinate arrays
    x_coords = (np.arange(x_size) - x_offset) * del_x
    y_coords = (np.arange(y_size) - y_offset) * del_y
    z_coords = (np.arange(z_size) - z_offset) * del_z

    # Initialize 3D array
    grid_3d = np.zeros((x_size, y_size, z_size))
    #initialize the others if using all of them
    if use_all_materials==True:
        #first the ions
        grid_3d_NOp = np.zeros((x_size, y_size, z_size))
        grid_3d_N2p = np.zeros((x_size, y_size, z_size))
        grid_3d_O2p = np.zeros((x_size, y_size, z_size))
        grid_3d_Np  = np.zeros((x_size, y_size, z_size))
        grid_3d_Op  = np.zeros((x_size, y_size, z_size))
        #then the neutral ones
        grid_3d_NO  = np.zeros((x_size, y_size, z_size))
        grid_3d_N2  = np.zeros((x_size, y_size, z_size))
        grid_3d_O2  = np.zeros((x_size, y_size, z_size))
        grid_3d_N   = np.zeros((x_size, y_size, z_size))
        grid_3d_O   = np.zeros((x_size, y_size, z_size))
        #then the temperatures
        grid_3d_Tv = np.zeros((x_size, y_size, z_size))
        grid_3d_T  = np.zeros((x_size, y_size, z_size))

    # Fill in the 3D grid
    for ix, x in enumerate(x_coords):
        for iy, y in enumerate(y_coords):
            # Calculate r for this (x,y) position
            r = np.sqrt(x**2 + y**2)

            # Only fill if r is within the data range (keep it circular)
            if r <= r_max:
                # Find closest r index in our grid
                ir = np.argmin(np.abs(r_vec - r))

                if round_down==True:
                    ir = int(np.floor(r/del_r))+1
                
                # Copy entire z column
                grid_3d[ix, iy, :] = grid_rho[:, ir]
                # if using all of the others
                if use_all_materials==True:
                    #first the ions
                    grid_3d_NOp[ix, iy, :] = grid_rho_NOp[:, ir]
                    grid_3d_N2p[ix, iy, :] = grid_rho_N2p[:, ir]
                    grid_3d_O2p[ix, iy, :] = grid_rho_O2p[:, ir]
                    grid_3d_Np[ix, iy, :]  = grid_rho_Np[:, ir]
                    grid_3d_Op[ix, iy, :]  = grid_rho_Op[:, ir]
                    #then the neutral ones
                    grid_3d_NO[ix, iy, :]  = grid_rho_NO[:, ir]
                    grid_3d_N2[ix, iy, :]  = grid_rho_N2[:, ir]
                    grid_3d_O2[ix, iy, :]  = grid_rho_O2[:, ir]
                    grid_3d_N[ix, iy, :]   = grid_rho_N[:, ir]
                    grid_3d_O[ix, iy, :]   = grid_rho_O[:, ir]
                    #then the temperatures
                    grid_3d_Tv[ix, iy, :] = grid_Tv[:, ir]
                    grid_3d_T[ix, iy, :]  = grid_T[:, ir]

if quarter==False:
    # First, determine size for full revolution (centered at origin)
    # We need to span from -r_max to +r_max in both x and y
    x_size = int(2 * r_max / del_r) + 1  # +1 to include both -r_max and +r_max
    y_size = int(2 * r_max / del_r) + 1
    z_size = len(grid_rho.T[0])

    # Create grid centered at origin (spans negative and positive)
    x_offset = (x_size - 1) / 2  # Center offset
    y_offset = (y_size - 1) / 2
    z_offset = 0

    del_x = del_r
    del_y = del_r

    # Create 3D Cartesian coordinate arrays (centered at origin)
    x_coords = (np.arange(x_size) - x_offset) * del_x
    y_coords = (np.arange(y_size) - y_offset) * del_y
    z_coords = (np.arange(z_size) - z_offset) * del_z

    # Initialize 3D array
    grid_3d = np.zeros((x_size, y_size, z_size))
    #initialize the others if using all of them
    if use_all_materials==True:
        #first the ions
        grid_3d_NOp = np.zeros((x_size, y_size, z_size))
        grid_3d_N2p = np.zeros((x_size, y_size, z_size))
        grid_3d_O2p = np.zeros((x_size, y_size, z_size))
        grid_3d_Np  = np.zeros((x_size, y_size, z_size))
        grid_3d_Op  = np.zeros((x_size, y_size, z_size))
        #then the neutral ones
        grid_3d_NO  = np.zeros((x_size, y_size, z_size))
        grid_3d_N2  = np.zeros((x_size, y_size, z_size))
        grid_3d_O2  = np.zeros((x_size, y_size, z_size))
        grid_3d_N   = np.zeros((x_size, y_size, z_size))
        grid_3d_O   = np.zeros((x_size, y_size, z_size))
        #then the temperatures
        grid_3d_Tv = np.zeros((x_size, y_size, z_size))
        grid_3d_T  = np.zeros((x_size, y_size, z_size))

    # Fill in the 3D grid (full 360 degrees)
    for ix, x in enumerate(x_coords):
        for iy, y in enumerate(y_coords):
            # Calculate r for this (x,y) position
            r = np.sqrt(x**2 + y**2)

            # Only fill if r is within the data range (keep it circular)
            if r <= r_max:
                # Find closest r index in our grid
                ir = np.argmin(np.abs(r_vec - r))

                if round_down==True:
                    ir = int(np.floor(r/del_r))+1
                
                # Copy entire z column
                grid_3d[ix, iy, :] = grid_rho[:, ir]
                # if using all of the others
                if use_all_materials==True:
                    #first the ions
                    grid_3d_NOp[ix, iy, :] = grid_rho_NOp[:, ir]
                    grid_3d_N2p[ix, iy, :] = grid_rho_N2p[:, ir]
                    grid_3d_O2p[ix, iy, :] = grid_rho_O2p[:, ir]
                    grid_3d_Np[ix, iy, :]  = grid_rho_Np[:, ir]
                    grid_3d_Op[ix, iy, :]  = grid_rho_Op[:, ir]
                    #then the neutral ones
                    grid_3d_NO[ix, iy, :]  = grid_rho_NO[:, ir]
                    grid_3d_N2[ix, iy, :]  = grid_rho_N2[:, ir]
                    grid_3d_O2[ix, iy, :]  = grid_rho_O2[:, ir]
                    grid_3d_N[ix, iy, :]   = grid_rho_N[:, ir]
                    grid_3d_O[ix, iy, :]   = grid_rho_O[:, ir]
                    #then the temperatures
                    grid_3d_Tv[ix, iy, :] = grid_Tv[:, ir]
                    grid_3d_T[ix, iy, :]  = grid_T[:, ir]

    print(f"3D grid shape (centered): {grid_3d.shape}")
    print(f"X range: [{x_coords[0]:.3f}, {x_coords[-1]:.3f}]")
    print(f"Y range: [{y_coords[0]:.3f}, {y_coords[-1]:.3f}]")
    print(f"Z range: [{z_coords[0]:.3f}, {z_coords[-1]:.3f}]")

##################################
#END CREATE 3D####################
##################################

##################################
#Calculate parameters and plot####
##################################

def plot_omega_p(array_to_plot,name_to_save):

    fig, axes = plt.subplots(1, 1, figsize=(14, 6))
    ax2 = axes

    slice_data = array_to_plot[int(x_offset), :, :]

    im = ax2.pcolormesh(
        z_coords,
        x_coords,
        slice_data,
        cmap='nipy_spectral',
        shading='auto', norm=LogNorm()
    )

    ax2.set_xlabel('Z [m]')
    ax2.set_ylabel('R [m]')
    ax2.set_title(f'Interpolated \u03C9$_p$')
    ax2.set_aspect('equal')
    #ax2.set_ylim(0, 0.6)
    cbar = plt.colorbar(im, ax=ax2)
    cbar.set_label('rad/s')

    plt.tight_layout()
    plt.savefig('comparison_omega_p_{}.png'.format(name_to_save), dpi=150, bbox_inches='tight')
    plt.show()

def plot_gamma(array_to_plot,name_to_save):
    fig, axes = plt.subplots(1, 1, figsize=(14, 6))
    ax2 = axes

    slice_data = array_to_plot[int(x_offset), :, :]

    im = ax2.pcolormesh(
        z_coords,
        x_coords,
        slice_data,
        cmap='nipy_spectral',
        shading='auto', norm=LogNorm()
    )

    ax2.set_xlabel('Z [m]')
    ax2.set_ylabel('R [m]')
    ax2.set_title(f'Interpolated gamma total')
    ax2.set_aspect('equal')
    #ax2.set_ylim(0, 0.6)
    cbar = plt.colorbar(im, ax=ax2)
    cbar.set_label('rad/s')

    plt.tight_layout()
    plt.savefig('comparison_gamma_{}.png'.format(name_to_save), dpi=150, bbox_inches='tight')
    plt.show()

#these are for plotting omega plasma, several values provided by Nick and Gemini
#could replace Z_charges with an array if needed but it's 1 for all for now
Z_charges=1
e_charge=1.60217663e-19
mass_e=9.10938356E-31
mass_ion=1.66053907e-27
ep_0=8.85418782e-12
kappa_function=0.5
#e-,NOp,N2p,O2p,Np,Op is the ordering
#term refers to mass density to number density conversion
term=[1/0.00054860*1e3*(6.0221409e23),1/29.9994514*1e3*(6.0221409e23),1/27.9994514*1e3*(6.0221409e23),1/31.9994514*1e3*(6.0221409e23),1/13.9994514*1e3*(6.0221409e23),1/15.9994514*1e3*(6.0221409e23)]
mass=[mass_e,mass_ion*30.006,mass_ion*28.014,mass_ion*31.998,mass_ion*14.007,mass_ion*15.999]
omega_convert=np.sqrt(np.array(term))*Z_charges*e_charge/np.sqrt(ep_0*np.array(mass))
omega_p=[]
omega_p.append(np.sqrt(grid_3d)*omega_convert[0])
if use_all_materials==True:
    omega_p.append(np.sqrt(grid_3d_NOp)*omega_convert[1])
    omega_p.append(np.sqrt(grid_3d_N2p)*omega_convert[2])
    omega_p.append(np.sqrt(grid_3d_O2p)*omega_convert[3])
    omega_p.append(np.sqrt(grid_3d_Np)*omega_convert[4])
    omega_p.append(np.sqrt(grid_3d_Op)*omega_convert[5])

plot_omega_p(omega_p[0],'e')

#if using others do this to get all gammas and plot both gammas and omega_p values
if use_all_materials==True:

    #now we are ready for the linear collisional frequency term
    #refer to ppt for where these values come from
    #first term for electrons is electron-ion, second term is electron-neutral. They all get summed together
    #order for ions is similar ion-electron is first, the ion ion for non same paris, and then ion-neutral.
    #They all get summed together too.

    #these are for all of them
    kb=1.380649E-23
    gamma_total=[]
    max_distance=(np.sqrt((ep_0*kb/(e_charge**2))/((grid_3d*np.array(term[0])/grid_3d_Tv)+
        Z_charges*(grid_3d_NOp*np.array(term[1])/grid_3d_T)+
        Z_charges*(grid_3d_N2p*np.array(term[2])/grid_3d_T)+
        Z_charges*(grid_3d_O2p*np.array(term[3])/grid_3d_T)+
        Z_charges*(grid_3d_Np*np.array(term[4])/grid_3d_T)+
        Z_charges*(grid_3d_Op*np.array(term[5])/grid_3d_T)
        )))

    #now electrons up first, specifically electron-ion scattering
    #I'm storing in arrays here and below instead outright summing because I might want to look at them individually at some point
    gamma_electron_ion=[]
    for i in range(5):
        #these two change each electron/ion
        mass_reduced=((mass_e*mass[1+i])/(mass_e+mass[1+i]))
        v_term=(kappa_function*np.sqrt((8*kb*grid_3d_Tv)/(np.pi*mass_e)+(8*kb*grid_3d_T)/(np.pi*mass[1+i])))
        #these equations are the same each electron/ion
        min_distance=(Z_charges**2*e_charge**2/(4*np.pi*ep_0*mass_reduced*v_term**2))
        ln_A=(np.log(np.array(max_distance)/min_distance))
        gamma_electron_ion.append(Z_charges**2*e_charge**4*ln_A/(4*np.pi*ep_0**2*mass_reduced**2*v_term**3))

    #these will be summed together shortly
    #the general form is the same each time with some slight differences based on which ones we are using
    #we will add these 5 plus contributions from the electron neutral scattering
    gamma_electron_ion[0]=gamma_electron_ion[0]*np.array(term[1])*grid_3d_NOp
    gamma_electron_ion[1]=gamma_electron_ion[1]*np.array(term[2])*grid_3d_N2p
    gamma_electron_ion[2]=gamma_electron_ion[2]*np.array(term[3])*grid_3d_O2p
    gamma_electron_ion[3]=gamma_electron_ion[3]*np.array(term[4])*grid_3d_Np
    gamma_electron_ion[4]=gamma_electron_ion[4]*np.array(term[5])*grid_3d_Op

    #now for electron-neutral species scattering terms
    #these formulas are taking from look up tables
    #the 'term' is used for mass to number density conversion.
    #I am using ion term conversion in place of a neutral equivalent. So it will be off by ~0.1% of actual.
    gamma_electron_neutral=[]
    gamma_electron_neutral.append(np.array(term[1])*grid_3d_NO*(4.3E-17*np.sqrt(grid_3d_Tv)))
    gamma_electron_neutral.append(np.array(term[2])*grid_3d_N2*(2.82e-16*np.sqrt(grid_3d_Tv)*(1+0.036*np.sqrt(grid_3d_Tv))))
    gamma_electron_neutral.append(np.array(term[3])*grid_3d_O2*(1.8E-16*np.sqrt(grid_3d_Tv)*(1+0.036*np.sqrt(grid_3d_Tv))))
    gamma_electron_neutral.append(np.array(term[4])*grid_3d_N*(4.5E-17*np.sqrt(grid_3d_Tv)))
    gamma_electron_neutral.append(np.array(term[5])*grid_3d_O*(8.9E-17*np.sqrt(grid_3d_Tv)*(1+1.4E-4*(grid_3d_Tv))))

    #final for electron gamma term which will be gamma_total[0]
    gamma_total.append(gamma_electron_ion[0]+gamma_electron_ion[1]+gamma_electron_ion[2]+gamma_electron_ion[3]+gamma_electron_ion[4]+gamma_electron_neutral[0]+gamma_electron_neutral[1]+gamma_electron_neutral[2]+gamma_electron_neutral[3]+gamma_electron_neutral[4])

    #next up is ions. I will do something a little differnt here.
    #I will do all ions-charged scattering, then do all the ions-neutral scattering, then add them all together
    #next up is NOp, specifically NOp scattering with electrons and other ions
    gamma_NOp_charged=[]
    for i in range(6):
        if i!=1:
            #these two change each electron/ion
            mass_reduced=((mass[1]*mass[i])/(mass[1]+mass[i]))
            if i==0:
                v_term=(kappa_function*np.sqrt((8*kb*grid_3d_T)/(np.pi*mass[1])+(8*kb*grid_3d_Tv)/(np.pi*mass[0])))
            if i!=0:
                v_term=(kappa_function*np.sqrt((8*kb*grid_3d_T)/(np.pi*mass[1])+(8*kb*grid_3d_T)/(np.pi*mass[i])))
            #these equations are the same each electron/ion
            min_distance=(Z_charges**2*e_charge**2/(4*np.pi*ep_0*mass_reduced*v_term**2))
            ln_A=(np.log(np.array(max_distance)/min_distance))
            gamma_NOp_charged.append(Z_charges**2*e_charge**4*ln_A/(4*np.pi*ep_0**2*mass_reduced**2*v_term**3))

    #these will be summed together shortly
    #the general form is the same each time with some slight differences based on which ones we are using
    #we will add these 5 plus contributions from the electron neutral scattering
    gamma_NOp_charged[0]=gamma_NOp_charged[0]*np.array(term[0])*grid_3d
    gamma_NOp_charged[1]=gamma_NOp_charged[1]*np.array(term[2])*grid_3d_N2p
    gamma_NOp_charged[2]=gamma_NOp_charged[2]*np.array(term[3])*grid_3d_O2p
    gamma_NOp_charged[3]=gamma_NOp_charged[3]*np.array(term[4])*grid_3d_Np
    gamma_NOp_charged[4]=gamma_NOp_charged[4]*np.array(term[5])*grid_3d_Op

    #next up is N2p, specifically N2p scattering with electrons and other ions
    gamma_N2p_charged=[]
    for i in range(6):
        if i!=2:
            #these two change each electron/ion
            mass_reduced=((mass[2]*mass[i])/(mass[2]+mass[i]))
            if i==0:
                v_term=(kappa_function*np.sqrt((8*kb*grid_3d_T)/(np.pi*mass[2])+(8*kb*grid_3d_Tv)/(np.pi*mass[0])))
            if i!=0:
                v_term=(kappa_function*np.sqrt((8*kb*grid_3d_T)/(np.pi*mass[2])+(8*kb*grid_3d_T)/(np.pi*mass[i])))
            #these equations are the same each electron/ion
            min_distance=(Z_charges**2*e_charge**2/(4*np.pi*ep_0*mass_reduced*v_term**2))
            ln_A=(np.log(np.array(max_distance)/min_distance))
            gamma_N2p_charged.append(Z_charges**2*e_charge**4*ln_A/(4*np.pi*ep_0**2*mass_reduced**2*v_term**3))

    #these will be summed together shortly
    #the general form is the same each time with some slight differences based on which ones we are using
    #we will add these 5 plus contributions from the electron neutral scattering
    gamma_N2p_charged[0]=gamma_N2p_charged[0]*np.array(term[0])*grid_3d
    gamma_N2p_charged[1]=gamma_N2p_charged[1]*np.array(term[1])*grid_3d_NOp
    gamma_N2p_charged[2]=gamma_N2p_charged[2]*np.array(term[3])*grid_3d_O2p
    gamma_N2p_charged[3]=gamma_N2p_charged[3]*np.array(term[4])*grid_3d_Np
    gamma_N2p_charged[4]=gamma_N2p_charged[4]*np.array(term[5])*grid_3d_Op

    #next up is O2p, specifically O2p scattering with electrons and other ions
    gamma_O2p_charged=[]
    for i in range(6):
        if i!=3:
            #these two change each electron/ion
            mass_reduced=((mass[3]*mass[i])/(mass[3]+mass[i]))
            if i==0:
                v_term=(kappa_function*np.sqrt((8*kb*grid_3d_T)/(np.pi*mass[3])+(8*kb*grid_3d_Tv)/(np.pi*mass[0])))
            if i!=0:
                v_term=(kappa_function*np.sqrt((8*kb*grid_3d_T)/(np.pi*mass[3])+(8*kb*grid_3d_T)/(np.pi*mass[i])))
            #these equations are the same each electron/ion
            min_distance=(Z_charges**2*e_charge**2/(4*np.pi*ep_0*mass_reduced*v_term**2))
            ln_A=(np.log(np.array(max_distance)/min_distance))
            gamma_O2p_charged.append(Z_charges**2*e_charge**4*ln_A/(4*np.pi*ep_0**2*mass_reduced**2*v_term**3))

    #these will be summed together shortly
    #the general form is the same each time with some slight differences based on which ones we are using
    #we will add these 5 plus contributions from the electron neutral scattering
    gamma_O2p_charged[0]=gamma_O2p_charged[0]*np.array(term[0])*grid_3d
    gamma_O2p_charged[1]=gamma_O2p_charged[1]*np.array(term[1])*grid_3d_NOp
    gamma_O2p_charged[2]=gamma_O2p_charged[2]*np.array(term[2])*grid_3d_N2p
    gamma_O2p_charged[3]=gamma_O2p_charged[3]*np.array(term[4])*grid_3d_Np
    gamma_O2p_charged[4]=gamma_O2p_charged[4]*np.array(term[5])*grid_3d_Op

    #next up is Np, specifically Np scattering with electrons and other ions
    gamma_Np_charged=[]
    for i in range(6):
        if i!=4:
            #these two change each electron/ion
            mass_reduced=((mass[4]*mass[i])/(mass[4]+mass[i]))
            if i==0:
                v_term=(kappa_function*np.sqrt((8*kb*grid_3d_T)/(np.pi*mass[4])+(8*kb*grid_3d_Tv)/(np.pi*mass[0])))
            if i!=0:
                v_term=(kappa_function*np.sqrt((8*kb*grid_3d_T)/(np.pi*mass[4])+(8*kb*grid_3d_T)/(np.pi*mass[i])))
            #these equations are the same each electron/ion
            min_distance=(Z_charges**2*e_charge**2/(4*np.pi*ep_0*mass_reduced*v_term**2))
            ln_A=(np.log(np.array(max_distance)/min_distance))
            gamma_Np_charged.append(Z_charges**2*e_charge**4*ln_A/(4*np.pi*ep_0**2*mass_reduced**2*v_term**3))

    #these will be summed together shortly
    #the general form is the same each time with some slight differences based on which ones we are using
    #we will add these 5 plus contributions from the electron neutral scattering
    gamma_Np_charged[0]=gamma_Np_charged[0]*np.array(term[0])*grid_3d
    gamma_Np_charged[1]=gamma_Np_charged[1]*np.array(term[1])*grid_3d_NOp
    gamma_Np_charged[2]=gamma_Np_charged[2]*np.array(term[2])*grid_3d_N2p
    gamma_Np_charged[3]=gamma_Np_charged[3]*np.array(term[3])*grid_3d_O2p
    gamma_Np_charged[4]=gamma_Np_charged[4]*np.array(term[5])*grid_3d_Op

    #Last up is Op, specifically Op scattering with electrons and other ions
    gamma_Op_charged=[]
    for i in range(6):
        if i!=5:
            #these two change each electron/ion
            mass_reduced=((mass[5]*mass[i])/(mass[5]+mass[i]))
            if i==0:
                v_term=(kappa_function*np.sqrt((8*kb*grid_3d_T)/(np.pi*mass[5])+(8*kb*grid_3d_Tv)/(np.pi*mass[0])))
            if i!=0:
                v_term=(kappa_function*np.sqrt((8*kb*grid_3d_T)/(np.pi*mass[5])+(8*kb*grid_3d_T)/(np.pi*mass[i])))
            #these equations are the same each electron/ion
            min_distance=(Z_charges**2*e_charge**2/(4*np.pi*ep_0*mass_reduced*v_term**2))
            ln_A=(np.log(np.array(max_distance)/min_distance))
            gamma_Op_charged.append(Z_charges**2*e_charge**4*ln_A/(4*np.pi*ep_0**2*mass_reduced**2*v_term**3))

    #these will be summed together shortly
    #the general form is the same each time with some slight differences based on which ones we are using
    #we will add these 5 plus contributions from the electron neutral scattering
    gamma_Op_charged[0]=gamma_Op_charged[0]*np.array(term[0])*grid_3d
    gamma_Op_charged[1]=gamma_Op_charged[1]*np.array(term[1])*grid_3d_NOp
    gamma_Op_charged[2]=gamma_Op_charged[2]*np.array(term[2])*grid_3d_N2p
    gamma_Op_charged[3]=gamma_Op_charged[3]*np.array(term[3])*grid_3d_O2p
    gamma_Op_charged[4]=gamma_Op_charged[4]*np.array(term[4])*grid_3d_Np

    #Now we have all of the ions - charged scattering terms
    #Now we will do all the ions - neutral scattering terms before summing together appropriately
    #the ions - neutral are much more straightforward as only the same to same has a complex formula
    #first term will be a formula followed by all the other terms with a very similar formula

    #need divid by mass of ion
    #also need ot use densities...
    gamma_ion_neutral=[]
    alpha=[1.7,1.74,1.58,1.1,0.8]
    densities=[term[1]*grid_3d_NO,term[2]*grid_3d_N2,term[3]*grid_3d_O2,term[4]*grid_3d_N,term[5]*grid_3d_O]
    term_con=0
    for i in range(5):
        if i!=0:
            term_con+=2.6E-15*np.sqrt((alpha[i])/((mass[i+1]*mass[1]/(mass[i+1]+mass[1]))/mass_ion))*densities[i]
    gamma_ion_neutral.append(8.4E-16*densities[0]+term_con)
    term_con=0
    for i in range(5):
        if i!=1:
            term_con+=2.6E-15*np.sqrt((alpha[i])/((mass[i+1]*mass[2]/(mass[i+1]+mass[2]))/mass_ion))*densities[i]
    gamma_ion_neutral.append(densities[1]*5.14E-17*np.sqrt(grid_3d_T)*(1-0.069*np.log10(grid_3d_T))**2+term_con)
    term_con=0
    for i in range(5):
        if i!=2:
            term_con+=2.6E-15*np.sqrt((alpha[i])/((mass[i+1]*mass[3]/(mass[i+1]+mass[3]))/mass_ion))*densities[i]
    gamma_ion_neutral.append(densities[2]*4.59E-17*np.sqrt(grid_3d_T)*(1-0.073*np.log10(grid_3d_T))**2+term_con)
    term_con=0
    for i in range(5):
        if i!=3:
            term_con+=2.6E-15*np.sqrt((alpha[i])/((mass[i+1]*mass[4]/(mass[i+1]+mass[4]))/mass_ion))*densities[i]
    gamma_ion_neutral.append(densities[3]*3.84E-17*np.sqrt(grid_3d_T)*(1-0.076*np.log10(grid_3d_T))**2+term_con)
    term_con=0
    for i in range(5):
        if i!=4:
            term_con+=2.6E-15*np.sqrt((alpha[i])/((mass[i+1]*mass[5]/(mass[i+1]+mass[5]))/mass_ion))*densities[i]
    gamma_ion_neutral.append(densities[4]*3.63E-17*np.sqrt(grid_3d_T)*(1-0.076*np.log10(grid_3d_T))**2+term_con)

    #now sum all totals
    gamma_total.append(gamma_NOp_charged[0]+gamma_NOp_charged[1]+gamma_NOp_charged[2]+gamma_NOp_charged[3]+gamma_NOp_charged[4]+gamma_ion_neutral[0])
    gamma_total.append(gamma_N2p_charged[0]+gamma_N2p_charged[1]+gamma_N2p_charged[2]+gamma_N2p_charged[3]+gamma_N2p_charged[4]+gamma_ion_neutral[1])
    gamma_total.append(gamma_O2p_charged[0]+gamma_O2p_charged[1]+gamma_O2p_charged[2]+gamma_O2p_charged[3]+gamma_O2p_charged[4]+gamma_ion_neutral[2])
    gamma_total.append(gamma_Np_charged[0]+gamma_Np_charged[1]+gamma_Np_charged[2]+gamma_Np_charged[3]+gamma_Np_charged[4]+gamma_ion_neutral[3])
    gamma_total.append(gamma_Op_charged[0]+gamma_Op_charged[1]+gamma_Op_charged[2]+gamma_Op_charged[3]+gamma_Op_charged[4]+gamma_ion_neutral[4])

    #now make all plots for visualizing

    #omega_p values
    plot_omega_p(omega_p[1],'NOp')
    plot_omega_p(omega_p[2],'N2p')
    plot_omega_p(omega_p[3],'O2p')
    plot_omega_p(omega_p[4],'Np')
    plot_omega_p(omega_p[5],'Op')
    #gammas
    plot_gamma(gamma_total[0],'e')
    plot_gamma(gamma_total[1],'NOp')
    plot_gamma(gamma_total[2],'N2p')
    plot_gamma(gamma_total[3],'O2p')
    plot_gamma(gamma_total[4],'Np')
    plot_gamma(gamma_total[5],'Op')

##################################
#End Calculate parameters and plot
##################################

##################################
#CREATE MATERIAL ID###############
##################################

material_ids_3d = {}
print(f"\nCreating material IDs...")
mat_ids = np.zeros_like(grid_3d, dtype=int)
nonzero_mask = grid_3d > 0
if nonzero_mask.any():
    # Get the data for processing
    valid_data = grid_3d[nonzero_mask]
    if use_all_bins==True:
        # OPTION A: Every unique value gets its own ID
        unique_vals, inverse_indices = np.unique(valid_data, return_inverse=True)
        # Shift indices by 1 so ID 0 remains reserved for "empty/zero"
        mat_ids[nonzero_mask] = inverse_indices + 1
        actual_num_used = len(unique_vals)
        print(f"    Mode: Full Resolution ({actual_num_used} unique values detected)")
    else:
        # OPTION B: Binning into num_materials
        vmin, vmax = valid_data.min(), valid_data.max()
        # Create logarithmic bins
        bins = np.logspace(np.log10(vmin), np.log10(vmax), num_materials + 1)
        # Digitize and clip to prevent the "vmax overflow" bin
        raw_indices = np.digitize(valid_data, bins)
        mat_ids[nonzero_mask] = np.clip(raw_indices, 1, num_materials)
        actual_num_used = num_materials
        print(f"Mode: Log-spaced binning (Max {num_materials} IDs)")
        print(f"Bins from {vmin:.6e} to {vmax:.6e}")
else:
    print("Warning: No non-zero values found in grid! There is an issue here.")
    actual_num_used = 0

# Report
unique_ids = np.unique(mat_ids)
print(f"    Unique material IDs: {unique_ids}")
print(f"    Number of Unique material IDs: {len(unique_ids)}")
print(f"    ID 0 (zero values): {(mat_ids == 0).sum()} points")

# Dynamically set boundary ID based on what was actually used
boundary_material_id = int(mat_ids.max() + 1)
print(f"\nAdding inner boundary material (ID {boundary_material_id})...")

boundary_material_id_arr = np.linspace(boundary_material_id, boundary_material_id + num_layers - 1, num_layers).astype(int)
if num_layers==0:
    boundary_material_id_arr=[boundary_material_id]
# Precompute r for all (ix, iy) positions
X, Y = np.meshgrid(x_coords, y_coords, indexing='ij')  # shape: (nx, ny)
R = np.sqrt(X**2 + Y**2)  # shape: (nx, ny)
# Precompute z mask
Z = z_coords  # 1D
z_mask = (Z >= z_min_boundary) & (Z <= z_max_boundary)  # shape: (nz,)

if fill_interior == False:
    for i in range(num_layers):
        bid = boundary_material_id_arr[i]
        # Occupied voxels: mat_ids > 0 and not already this boundary layer
        occupied = (mat_ids > 0) & (mat_ids != bid)  # (nx, ny, nz)
        # Apply z mask
        occupied &= z_mask[np.newaxis, np.newaxis, :]
        # For each of the 4 XY neighbors and the +Z neighbor,
        # mark if: neighbor is empty AND neighbor r < current r
        # Neighbor ix-1: exists at [1:], source at [1:]
        mask = occupied[1:, :, :] & (mat_ids[:-1, :, :] == 0)
        mask &= (R[:-1, :] < R[1:, :])[:, :, np.newaxis]
        mat_ids[:-1, :, :][mask] = bid
        # Neighbor ix+1: exists at [:-1], source at [:-1]
        mask = occupied[:-1, :, :] & (mat_ids[1:, :, :] == 0)
        mask &= (R[1:, :] < R[:-1, :])[:, :, np.newaxis]
        mat_ids[1:, :, :][mask] = bid
        # Neighbor iy-1
        mask = occupied[:, 1:, :] & (mat_ids[:, :-1, :] == 0)
        mask &= (R[:, :-1] < R[:, 1:])[:, :, np.newaxis]
        mat_ids[:, :-1, :][mask] = bid
        # Neighbor iy+1
        mask = occupied[:, :-1, :] & (mat_ids[:, 1:, :] == 0)
        mask &= (R[:, 1:] < R[:, :-1])[:, :, np.newaxis]
        mat_ids[:, 1:, :][mask] = bid
        # Neighbor iz+1
        mask = occupied[:, :, :-1] & (mat_ids[:, :, 1:] == 0)
        mat_ids[:, :, 1:][mask] = bid
    # Collapse all boundary layers back to base ID
    for i in range(num_layers - 1):
        mat_ids[mat_ids == boundary_material_id_arr[i + 1]] = boundary_material_id_arr[0]

if fill_interior == True:
    bid = boundary_material_id_arr[0]
    while True:
        occupied = (mat_ids > 0)
        occupied &= z_mask[np.newaxis, np.newaxis, :]
        prev_count = (mat_ids == bid).sum()
        mask = occupied[1:, :, :] & (mat_ids[:-1, :, :] == 0) & (R[:-1, :] < R[1:, :])[:, :, np.newaxis]
        mat_ids[:-1, :, :][mask] = bid
        mask = occupied[:-1, :, :] & (mat_ids[1:, :, :] == 0) & (R[1:, :] < R[:-1, :])[:, :, np.newaxis]
        mat_ids[1:, :, :][mask] = bid
        mask = occupied[:, 1:, :] & (mat_ids[:, :-1, :] == 0) & (R[:, :-1] < R[:, 1:])[:, :, np.newaxis]
        mat_ids[:, :-1, :][mask] = bid
        mask = occupied[:, :-1, :] & (mat_ids[:, 1:, :] == 0) & (R[:, 1:] < R[:, :-1])[:, :, np.newaxis]
        mat_ids[:, 1:, :][mask] = bid
        mask = occupied[:, :, :-1] & (mat_ids[:, :, 1:] == 0)
        mat_ids[:, :, 1:][mask] = bid
        if (mat_ids == bid).sum() == prev_count:
            break  # No new voxels added, done

boundary_material_id = boundary_material_id_arr[0]
boundary_count = int((mat_ids == boundary_material_id).sum())

"""
#this is the old way kept just in case we decide we want it back.

# Dynamically set boundary ID based on what was actually used
boundary_material_id = int(mat_ids.max() + 1)
print(f"\nAdding inner boundary material (ID {boundary_material_id})...")

# Mark boundaries at r-1 and z+1 where data exists
boundary_count = 0

# Loop through all points in 3D Cartesian grid
for ix in range(mat_ids.shape[0]):  # Loop over X
    for iy in range(mat_ids.shape[1]):  # Loop over Y
        for iz in range(mat_ids.shape[2]):  # Loop over Z
            
            # Calculate r for this (x,y) position
            x = x_coords[ix]
            y = y_coords[iy]
            z = z_coords[iz]  # Get the z coordinate
            r = np.sqrt(x**2 + y**2)
            
            # Check if z is within the specified range
            if z < z_min_boundary or z > z_max_boundary:
                continue  # Skip this point if z is out of range

            if fill_interior == True:
                if mat_ids[ix, iy, iz] > 0:
                    
                    # Mark neighbors moving inward (toward r=0)
                    # Check all 4 neighbors in X-Y plane and mark if they reduce r
                    
                    # Neighbor at ix-1
                    if ix > 0:
                        r_neighbor = np.sqrt(x_coords[ix-1]**2 + y**2)
                        if r_neighbor < r and mat_ids[ix-1, iy, iz] == 0:
                            mat_ids[ix-1, iy, iz] = boundary_material_id
                            boundary_count += 1
                    
                    # Neighbor at ix+1
                    if ix < mat_ids.shape[0]-1:
                        r_neighbor = np.sqrt(x_coords[ix+1]**2 + y**2)
                        if r_neighbor < r and mat_ids[ix+1, iy, iz] == 0:
                            mat_ids[ix+1, iy, iz] = boundary_material_id
                            boundary_count += 1
                    
                    # Neighbor at iy-1
                    if iy > 0:
                        r_neighbor = np.sqrt(x**2 + y_coords[iy-1]**2)
                        if r_neighbor < r and mat_ids[ix, iy-1, iz] == 0:
                            mat_ids[ix, iy-1, iz] = boundary_material_id
                            boundary_count += 1
                    
                    # Neighbor at iy+1
                    if iy < mat_ids.shape[1]-1:
                        r_neighbor = np.sqrt(x**2 + y_coords[iy+1]**2)
                        if r_neighbor < r and mat_ids[ix, iy+1, iz] == 0:
                            mat_ids[ix, iy+1, iz] = boundary_material_id
                            boundary_count += 1
                    
                    # Mark z+1 if it exists and is currently zero
                    if iz < mat_ids.shape[2]-1 and mat_ids[ix, iy, iz+1] == 0:
                        mat_ids[ix, iy, iz+1] = boundary_material_id
                        boundary_count += 1
            
            if fill_interior == False:
                if mat_ids[ix, iy, iz] > 0 and mat_ids[ix, iy, iz] != boundary_material_id:
                    
                    # Check all 4 neighbors in X-Y plane, mark those moving toward r=0
                    
                    # Neighbor at ix-1
                    if ix > 0:
                        r_neighbor = np.sqrt(x_coords[ix-1]**2 + y**2)
                        if r_neighbor < r and mat_ids[ix-1, iy, iz] == 0:
                            mat_ids[ix-1, iy, iz] = boundary_material_id
                            boundary_count += 1
                    
                    # Neighbor at ix+1
                    if ix < mat_ids.shape[0]-1:
                        r_neighbor = np.sqrt(x_coords[ix+1]**2 + y**2)
                        if r_neighbor < r and mat_ids[ix+1, iy, iz] == 0:
                            mat_ids[ix+1, iy, iz] = boundary_material_id
                            boundary_count += 1
                    
                    # Neighbor at iy-1
                    if iy > 0:
                        r_neighbor = np.sqrt(x**2 + y_coords[iy-1]**2)
                        if r_neighbor < r and mat_ids[ix, iy-1, iz] == 0:
                            mat_ids[ix, iy-1, iz] = boundary_material_id
                            boundary_count += 1
                    
                    # Neighbor at iy+1
                    if iy < mat_ids.shape[1]-1:
                        r_neighbor = np.sqrt(x**2 + y_coords[iy+1]**2)
                        if r_neighbor < r and mat_ids[ix, iy+1, iz] == 0:
                            mat_ids[ix, iy+1, iz] = boundary_material_id
                            boundary_count += 1
                    
                    # Mark z+1 if it exists and is currently zero
                    if iz < mat_ids.shape[2]-1 and mat_ids[ix, iy, iz+1] == 0:
                        mat_ids[ix, iy, iz+1] = boundary_material_id
                        boundary_count += 1
"""
##################################
#END CREATE MATERIAL ID###########
##################################

##################################
#PLOT MATERIAL ID#################
##################################
min_col = 0
max_col = boundary_material_id

fig, axes = plt.subplots(1, 1, figsize=(14, 6))

# --- Interpolated Slice ---
ax2 = axes

slice_data = mat_ids[int(x_offset), :, :]

im = ax2.pcolormesh(
    z_coords,
    x_coords,
    slice_data,
    cmap='nipy_spectral',
    shading='auto', vmin=min_col, vmax=max_col
)

ax2.set_xlabel('Z [m]')
ax2.set_ylabel('R [m]')
ax2.set_title(f'Interpolated')
ax2.set_aspect('equal')
#ax2.set_ylim(0, 0.6)
cbar = plt.colorbar(im, ax=ax2)
cbar.set_label('ID Number')

plt.tight_layout()
plt.savefig(f'comparison_ID.png', dpi=150, bbox_inches='tight')
plt.show()

fig, axes = plt.subplots(1, 1, figsize=(14, 6))
print(f"\nSaved comparison plot for rho_e")
##################################
#END PLOT MATERIAL ID#############
##################################


#### PREP FOR FDTD #######
###########################
############################
#############################
if fdtd_prep==True:
    print(f"\nGenerating binaries from numpy arrays for FDTD read-in...")

    # work with mat_ids, make sure data type is correct
    # Get the material IDs
    mat_ids_export=mat_ids
    if swap_axis==True:
        mat_ids_export = np.transpose(mat_ids_export, (2, 1, 0))
    mat_ids_export = np.array(mat_ids_export, dtype=np.float32)
    
    #create the array

    #use these for buffer on all sides +- of the object
    x_shift_fdtd=60
    y_shift_fdtd=60
    z_shift_fdtd=60

    #Calculate the ranges based on the total size then we want
    x_size_fdtd = 2*x_shift_fdtd + mat_ids_export.shape[0]
    y_size_fdtd = 2*y_shift_fdtd + mat_ids_export.shape[1]
    z_size_fdtd = 2*z_shift_fdtd + mat_ids_export.shape[2]

    #zeros will be ignored by fdtd - specific to this optional read in section
    data=np.zeros((x_size_fdtd,y_size_fdtd,z_size_fdtd), dtype=np.float32)

    #print(mat_ids_export.shape[0],mat_ids_export.shape[1],mat_ids_export.shape[2])

    #Perform the assignment
    data[x_shift_fdtd:x_size_fdtd-x_shift_fdtd, y_shift_fdtd:y_size_fdtd-y_shift_fdtd, z_shift_fdtd:z_size_fdtd-z_shift_fdtd] = mat_ids_export

    print(f"\nPrinting the size of the array for use in FDTD...")
    print(x_size_fdtd,y_size_fdtd,z_size_fdtd)

    # save
    # mat ids of zero will be ignored on fdtd import by if statement selection
    data.flatten(order='F').tofile('optional_geom_bulk.bin')

    print(f"\nFinished binary generation. Now prepping material ID parameter arrays for export...")

    """
    # Also save material properties lookup
    material_properties = []

    # Get unique material IDs (excluding 0)
    unique_mat_ids = np.unique(mat_ids[mat_ids > 0])

    for mat_id in unique_mat_ids:
        mask = mat_ids == mat_id
        if use_all_materials==False:
            #just omega_p for electrons
            avg_value = omega_p[0][mask].mean() if mask.any() else 0.0
            material_properties.append([avg_value])
        if use_all_materials==True:
            #first plasma frequencies
            avg_value  = omega_p[0][mask].mean() if mask.any() else 0.0
            avg_value_NOp  = omega_p[1][mask].mean() if mask.any() else 0.0
            avg_value_N2p  = omega_p[2][mask].mean() if mask.any() else 0.0
            avg_value_O2p  = omega_p[3][mask].mean() if mask.any() else 0.0
            avg_value_Np   = omega_p[4][mask].mean() if mask.any() else 0.0
            avg_value_Op   = omega_p[5][mask].mean() if mask.any() else 0.0
            #now collision frequencies
            avg_value_g  = gamma_total[0][mask].mean() if mask.any() else 0.0
            avg_value_NOp_g  = gamma_total[1][mask].mean() if mask.any() else 0.0
            avg_value_N2p_g  = gamma_total[2][mask].mean() if mask.any() else 0.0
            avg_value_O2p_g  = gamma_total[3][mask].mean() if mask.any() else 0.0
            avg_value_Np_g   = gamma_total[4][mask].mean() if mask.any() else 0.0
            avg_value_Op_g   = gamma_total[5][mask].mean() if mask.any() else 0.0

            material_properties.append([
                avg_value,avg_value_g,
                avg_value_NOp,avg_value_NOp_g,
                avg_value_N2p,avg_value_N2p_g,
                avg_value_O2p,avg_value_O2p_g,
                avg_value_Np,avg_value_Np_g,
                avg_value_Op,avg_value_Op_g
            ])
    """
    # this is the new section from above - it is slightly different numbers due to machine precsion
    # but it was tested and the difference was very, very tiny so seem sokay
    # it is way faster
    # Also save material properties lookup
    material_properties = []

    # Get unique material IDs (excluding 0)
    unique_mat_ids = np.unique(mat_ids[mat_ids > 0])

    flat_ids = mat_ids.ravel()
    counts = np.bincount(flat_ids)

    if use_all_materials == False:
        sums = np.bincount(flat_ids, weights=omega_p[0].ravel())
        means = np.divide(sums, counts, out=np.zeros_like(sums), where=counts>0)

        material_properties = means[unique_mat_ids][:, None].tolist()

    if use_all_materials == True:

        # plasma frequencies
        wp0 = np.divide(np.bincount(flat_ids, weights=omega_p[0].ravel()), counts, out=np.zeros_like(counts, dtype=float), where=counts>0)
        wp1 = np.divide(np.bincount(flat_ids, weights=omega_p[1].ravel()), counts, out=np.zeros_like(counts, dtype=float), where=counts>0)
        wp2 = np.divide(np.bincount(flat_ids, weights=omega_p[2].ravel()), counts, out=np.zeros_like(counts, dtype=float), where=counts>0)
        wp3 = np.divide(np.bincount(flat_ids, weights=omega_p[3].ravel()), counts, out=np.zeros_like(counts, dtype=float), where=counts>0)
        wp4 = np.divide(np.bincount(flat_ids, weights=omega_p[4].ravel()), counts, out=np.zeros_like(counts, dtype=float), where=counts>0)
        wp5 = np.divide(np.bincount(flat_ids, weights=omega_p[5].ravel()), counts, out=np.zeros_like(counts, dtype=float), where=counts>0)

        # collision frequencies
        g0 = np.divide(np.bincount(flat_ids, weights=gamma_total[0].ravel()), counts, out=np.zeros_like(counts, dtype=float), where=counts>0)
        g1 = np.divide(np.bincount(flat_ids, weights=gamma_total[1].ravel()), counts, out=np.zeros_like(counts, dtype=float), where=counts>0)
        g2 = np.divide(np.bincount(flat_ids, weights=gamma_total[2].ravel()), counts, out=np.zeros_like(counts, dtype=float), where=counts>0)
        g3 = np.divide(np.bincount(flat_ids, weights=gamma_total[3].ravel()), counts, out=np.zeros_like(counts, dtype=float), where=counts>0)
        g4 = np.divide(np.bincount(flat_ids, weights=gamma_total[4].ravel()), counts, out=np.zeros_like(counts, dtype=float), where=counts>0)
        g5 = np.divide(np.bincount(flat_ids, weights=gamma_total[5].ravel()), counts, out=np.zeros_like(counts, dtype=float), where=counts>0)

        material_properties = np.stack([
            wp0, g0,
            wp1, g1,
            wp2, g2,
            wp3, g3,
            wp4, g4,
            wp5, g5
        ], axis=1)[unique_mat_ids].tolist()

    material_properties = np.array(material_properties)
    np.save('material_properties.npy', material_properties)
    print(f"\nSaved material_properties.npy with {len(material_properties)} materials")