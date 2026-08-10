import time as tp
from matplotlib import pyplot as plt
import torch as th
import math
import sys

import os
os.chdir('C:/Users/dari6475/Desktop')

##########################################################
# GPU SETUP AND VERIFICATION
##########################################################
def setup_gpu():
    """Setup and verify GPU availability"""
    if th.cuda.is_available():
        device = th.device('cuda')
        print(f"✓ GPU detected: {th.cuda.get_device_name(0)}")
        print(f"✓ CUDA version: {th.version.cuda}")
        print(f"✓ GPU memory: {th.cuda.get_device_properties(0).total_memory / 1e9:.2f} GB")
        
        # Set some optimizations
        th.backends.cudnn.benchmark = True  # Auto-tune kernels
        th.cuda.empty_cache()  # Clear cache
        
        return device
    else:
        print("⚠ GPU not available, using CPU (will be much slower)")
        return th.device('cpu')

device = setup_gpu()
dtype = th.float32
compile_functions=False

##########################################################
# PROFILING UTILITIES
##########################################################
class GPUProfiler:
    """Context manager for profiling GPU operations"""
    def __init__(self, name, device):
        self.name = name
        self.device = device
        
    def __enter__(self):
        if self.device.type == 'cuda':
            th.cuda.synchronize()  # Wait for GPU to finish
        self.start = tp.time()
        return self
        
    def __exit__(self, *args):
        if self.device.type == 'cuda':
            th.cuda.synchronize()  # Wait for GPU to finish
        self.elapsed = tp.time() - self.start
        print(f"{self.name}: {self.elapsed:.4f} seconds")
        if self.device.type == 'cuda':
            print(f"  GPU memory allocated: {th.cuda.memory_allocated() / 1e9:.2f} GB")
            print(f"  GPU memory reserved: {th.cuda.memory_reserved() / 1e9:.2f} GB")

##########################################################
# SETUP SECTION
##########################################################
with GPUProfiler("Simulation Setup", device):

    # initalize several constants - some of which will be overwritten by read in file
    # FORTRAN required intializing all variables, but python doesn't inherently require this. 
    # However some still need initializing in some instances for the specific way it's used here.
    # When in doubt, just initalize it. It can't hurt.

    # used for sheet setup section under ga,gb
    sheet_ep_average_x = 0.0 
    sheet_sig_average_x = 0.0
    sheet_ep_average_y = 0.0
    sheet_sig_average_y = 0.0
    sheet_ep_average_z = 0.0
    sheet_sig_average_z = 0.0

    # sheet thickness - will be set below - unclear if it needs initializing
    sheet_thickness = 0.0

    # used for determining if sheets are present for ga,gb setup section
    counter_x_sheet = 0      
    counter_y_sheet = 0
    counter_z_sheet = 0

    # used for sheets in main fdtd looping for speed improvements
    count_unique_sheets_x = 0
    count_unique_sheets_y = 0
    count_unique_sheets_z = 0

    #related to mirrored plane waves for IGP conditions
    xlow_wall = 1.0
    ylow_wall = 1.0
    zlow_wall = 1.0
    xhigh_wall = 1.0
    yhigh_wall = 1.0
    zhigh_wall = 1.0
    Mx_mirror = 1.0
    My_mirror = 1.0
    Mz_mirror = 1.0
    Jx_mirror = 1.0
    Jy_mirror = 1.0
    Jz_mirror = 1.0
    x_mirror_offset = 0.0
    y_mirror_offset = 0.0
    z_mirror_offset = 0.0
    i_mirror = 0.0
    j_mirror = 0.0
    k_mirror = 0.0           
    use_x_mirror = 0.0
    use_y_mirror = 0.0
    use_z_mirror = 0.0
    is_mirror=0
    mirror_type=0
    theta_mirror = phi_mirror = mirror_height = 0

    # plane wave incidence intialization
    phi=theta=0

    # used for videos - saved fields in time for viewing later
    vid_size1 = 0            
    vid_size2 = 0
    slice_type=0
    slice_location=0
    video_on=0

    # far field array starting place - starts at 1 and gets shifted by algorithm as needed
    data_out_phi = 1
    data_out_theta = 1
    
    # wave amplitudes - unclear if they need initializing or not
    plane_wave_amp=0
    antenna_amp=0

    # number of items - defaulting to zero until read in section
    num_ports=0
    num_sheets=0
    num_far_field_angles=0
    num_blocks=0
    num_spheres=0
    num_cylinders=0
    num_sheets_x=0
    num_sheets_y=0
    num_sheets_z=0
    num_materials=0
    num_sheet_materials=0

    # if pbc version is used, there are numerious places where we shift things by 1 to accomodate 1 size now being a duplicate.
    pbc_shift=0

    # Constants
    c = 299792458.0
    mu_0 = 1.25663706E-6
    ep_0 = 8.85418782E-12
    pi = 3.14159265358979323846

    # Simulation boundary variable - seperates TF from SF
    buffer = 3

    # PML variables
    nxPML_1 = nxPML_2 = nyPML_1 = nyPML_2 = nzPML_1 = nzPML_2 = 10
    m = 3
    ma = 1 

    fdtd_input_file_name = sys.argv[1] 
    with open(fdtd_input_file_name, 'r') as f:
        x_size, y_size, z_size = map(int, f.readline().split(','))
        del_x, del_y, del_z = map(float, f.readline().split(','))
        time_steps = int(f.readline())
        f_center = float(f.readline())
        pulse_type = int(f.readline())
        type_sim=f.readline().strip()
        if type_sim=='uc':
            type_sim=0
        if type_sim=='pml':
            type_sim=1
        excitation_type=f.readline().strip()
        if excitation_type=='antenna':
            antenna_amp, antenna_time_delay = map(float, f.readline().split(','))
            port_excited = int(f.readline())
        if excitation_type=='plane wave':
            plane_wave_amp, plane_wave_time_delay = map(float, f.readline().split(','))
            theta, phi, pol = map(float, f.readline().split(','))

        num_far_field_angles = int(f.readline())
        far_field_angles = []
        for _ in range(num_far_field_angles):
            angle_theta, angle_phi = map(float, f.readline().split(','))
            angle_theta_rad = pi / 180.0 * angle_theta
            angle_phi_rad = pi / 180.0 * angle_phi
            far_field_angles.append([angle_theta_rad, angle_phi_rad])
        
        char_is_mirror = f.readline()
        if char_is_mirror.strip() == 'IGP':
            is_mirror = 1
            char_mirror_type = f.readline()
            if char_mirror_type == '-x':
                mirror_type = 1
                theta_mirror = theta
                phi_mirror = 180.0 - phi
            elif char_mirror_type == '+x':
                mirror_type = 2
                theta_mirror = theta
                phi_mirror = 180.0 - phi
            elif char_mirror_type == '-y':
                mirror_type = 3
                theta_mirror = theta
                phi_mirror = 360.0 - phi
            elif char_mirror_type == '+y':
                mirror_type = 4
                theta_mirror = theta
                phi_mirror = 360.0 - phi
            elif char_mirror_type == '-z':
                mirror_type = 5
                theta_mirror = 180.0 - theta
                phi_mirror = phi
            elif char_mirror_type == '+z':
                mirror_type = 6
                theta_mirror = 180.0 - theta
                phi_mirror = phi
            mirror_height = f.readline()

        video_on_char = f.readline()
        if video_on_char.strip() == 'yes':
            video_on = 1
            char_slice = f.readline()
            if char_slice.strip() == 'x':
                v_slice = 0
            elif char_slice.strip() == 'y':
                v_slice = 1
            elif char_slice.strip() == 'z':
                v_slice = 2
            slice_location = f.readline()

        filename = f.readline()

        num_total_materials = int(f.readline())
        i = 0
        j = 0
        if num_total_materials > 0:
            # max size it could be, will resize below
            materials_properties = th.zeros((num_total_materials, 3, 3), device=device)
            sheet_properties = th.zeros((num_total_materials, 4, 3), device=device)
        for ii in range(num_total_materials):
            material_type = f.readline()
            if material_type.strip() == 'volume':
                num_materials = num_materials + 1
                materials_properties[i,0,0] = int(f.readline())
                materials_properties[i,1,0], materials_properties[i,1,1], materials_properties[i,1,2] = list(map(float, f.readline().split(',')))
                materials_properties[i,2,0], materials_properties[i,2,1], materials_properties[i,2,2] = list(map(float, f.readline().split(',')))
                i = i + 1
            elif material_type.strip() == 'sheet':
                num_sheet_materials = num_sheet_materials + 1
                # old method imported thickness and ep,sig anisotropic for each and was a 'true' way to do it.
                # but this meant having to fix a thickness that would numerically work while then choosing a sigma,ep that matched manufacterer reporting impedance.
                # there was also a constraint for managing the thickness of each layer correctly. I had defaulted to requiring a fixed thickness per elevation.
                # new method here simplifies but restricts the user inputs - should work better in long run I think.
                # I left the bones the same below so that it can be modified easily for future use if we desire a permittivity modification by the user.
                sheet_properties[j,0,0] = int(f.readline())
                imped_value = float(f.readline())
                sheet_thickness = min(del_x, del_y, del_z) / 5000.0
                sheet_properties[j,2,0] = 1.0
                sheet_properties[j,2,1] = 1.0
                sheet_properties[j,2,2] = 1.0
                if imped_value == 0:
                    # should account for min 1E-12 thickness or ~1E-9 cell size with copper plating
                    sheet_properties[j,3,0] = 1E20
                    sheet_properties[j,3,1] = 1E20
                    sheet_properties[j,3,2] = 1E20
                if imped_value != 0:
                    sheet_properties[j,3,0] = 1.0 / (imped_value * sheet_thickness)
                    sheet_properties[j,3,1] = 1.0 / (imped_value * sheet_thickness)
                    sheet_properties[j,3,2] = 1.0 / (imped_value * sheet_thickness)
                j = j + 1
        if num_materials > 0 and num_materials < num_total_materials:
            materials_properties = materials_properties[0:num_materials, :, :]
        elif num_materials == 0 and materials_properties is not None:
            del materials_properties
        if num_sheet_materials < num_total_materials and num_sheet_materials > 0:
            sheet_properties = sheet_properties[0:num_sheet_materials, :, :]
        elif num_sheet_materials == 0 and sheet_properties is not None:
            del sheet_properties

        num_objects = int(f.readline())
        if num_objects > 0:
            # max size it could be, will resize below
            object_type = [''] * num_objects
            spheres = th.zeros((num_objects, 3, 3), device=device)
            blocks = th.zeros((num_objects, 3, 3), device=device)
            cylinders = th.zeros((num_objects, 4, 3), device=device)
        else:
            blocks=[]
            spheres=[]
            cylinders=[]
        i = 0
        j = 0
        k = 0
        for ii in range(num_objects):
            object_type[ii] = f.readline()
            if object_type[ii].strip() == 'block':
                num_blocks = num_blocks + 1
                blocks[i,0,0] = int(f.readline())
                blocks[i,1,0], blocks[i,1,1], blocks[i,1,2] = list(map(int, f.readline().split(',')))
                blocks[i,2,0], blocks[i,2,1], blocks[i,2,2] = list(map(int, f.readline().split(',')))
                i = i + 1
            elif object_type[ii].strip() == 'sphere':
                num_spheres = num_spheres + 1
                spheres[j,0,0] = int(f.readline())
                spheres[j,1,0], spheres[j,1,1], spheres[j,1,2] = list(map(int, f.readline().split(',')))
                spheres[j,2,0] = int(f.readline())
                j = j + 1
            elif object_type[ii].strip() == 'cylinder':
                num_cylinders = num_cylinders + 1
                cylinders[k,0,0] = int(f.readline())
                char_cylinders = f.readline()
                if char_cylinders.strip() == 'x':
                    cylinders[k,1,0] = 0
                elif char_cylinders.strip() == 'y':
                    cylinders[k,1,0] = 1
                elif char_cylinders.strip() == 'z':
                    cylinders[k,1,0] = 2
                cylinders[k,2,0], cylinders[k,2,1], cylinders[k,2,2] = list(map(int, f.readline().split(',')))
                cylinders[k,3,0], cylinders[k,3,1] = list(map(int, f.readline().split(',')))
                k = k + 1
        if num_blocks > 0 and num_blocks < num_objects:
            blocks = blocks[0:num_blocks, :, :]
        elif num_blocks == 0 and blocks is not None:
            del blocks
        if num_spheres > 0 and num_spheres < num_objects:
            spheres = spheres[0:num_spheres, :, :]
        elif num_spheres == 0 and spheres is not None:
            del spheres
        if num_cylinders > 0 and num_cylinders < num_objects:
            cylinders = cylinders[0:num_cylinders, :, :]
        elif num_cylinders == 0 and cylinders is not None:
            del cylinders

        num_sheets = int(f.readline())
        if num_sheets > 0:
            sheets_x = th.zeros((num_sheets, 4, 2), device=device)
            sheets_y = th.zeros((num_sheets, 4, 2), device=device)
            sheets_z = th.zeros((num_sheets, 4, 2), device=device)
        else:
            sheets_x=[]
            sheets_y=[]
            sheets_z=[]
        i = 0
        j = 0
        k = 0
        for ii in range(num_sheets):
            sheet_position_index = f.readline()
            if sheet_position_index.strip() == 'x':
                num_sheets_x = num_sheets_x + 1
                sheets_x[i,0,0] = int(f.readline())
                sheets_x[i,1,0] = int(f.readline())
                sheets_x[i,2,0], sheets_x[i,2,1] = list(map(int, f.readline().split(',')))
                sheets_x[i,3,0], sheets_x[i,3,1] = list(map(int, f.readline().split(',')))
                i = i + 1
            elif sheet_position_index.strip() == 'y':
                num_sheets_y = num_sheets_y + 1
                sheets_y[j,0,0] = int(f.readline())
                sheets_y[j,1,0] = int(f.readline())
                sheets_y[j,2,0], sheets_y[j,2,1] = list(map(int, f.readline().split(',')))
                sheets_y[j,3,0], sheets_y[j,3,1] = list(map(int, f.readline().split(',')))
                j = j + 1
            elif sheet_position_index.strip() == 'z':
                num_sheets_z = num_sheets_z + 1
                sheets_z[k,0,0] = int(f.readline())
                sheets_z[k,1,0] = int(f.readline())
                sheets_z[k,2,0], sheets_z[k,2,1] = list(map(int, f.readline().split(',')))
                sheets_z[k,3,0], sheets_z[k,3,1] = list(map(int, f.readline().split(',')))
                k = k + 1
        
        if num_sheets_x > 0 and num_sheets_x < num_sheets:
            sheets_x = sheets_x[0:num_sheets_x, :, :]
        elif num_sheets_x == 0 and sheets_x is not None:
            del sheets_x
        if num_sheets_y > 0 and num_sheets_y < num_sheets:
            sheets_y = sheets_y[0:num_sheets_y, :, :]
        elif num_sheets_y == 0 and sheets_y is not None:
            del sheets_y
        if num_sheets_z > 0 and num_sheets_z < num_sheets:
            sheets_z = sheets_z[0:num_sheets_z, :, :]
        elif num_sheets_z == 0 and sheets_z is not None:
            del sheets_z

        num_ports = int(f.readline())
        if num_ports > 0:
            ports = th.zeros((num_ports, 4, 4), device=device)
        for i in range(num_ports):
            char_ports = f.readline()
            if char_ports.strip() == 'x':
                ports[i,0,0] = 0
            elif char_ports.strip() == 'y':
                ports[i,0,0] = 1
            elif char_ports.strip() == 'z':
                ports[i,0,0] = 2
            imped_port = float(f.readline())
            speed_port = c
            # L is 2 and C is 3 - 1 and 4 are unused from old notation in case we ever add back R and G
            ports[i,1,1] = imped_port / speed_port
            ports[i,1,2] = 1.0 / (imped_port * speed_port)
            ports[i,2,0], ports[i,2,1], ports[i,2,2] = list(map(int, f.readline().split(',')))
            ports[i,3,0], ports[i,3,1], ports[i,3,2] = list(map(int, f.readline().split(',')))
            if char_ports.strip() == 'x':
                ports[i,3,1] = ports[i,3,1] + 1
                ports[i,3,2] = ports[i,3,2] + 1
            elif char_ports.strip() == 'y':
                ports[i,3,0] = ports[i,3,0] + 1
                ports[i,3,2] = ports[i,3,2] + 1
            elif char_ports.strip() == 'z':
                ports[i,3,0] = ports[i,3,0] + 1
                ports[i,3,1] = ports[i,3,1] + 1

    #simulation time constant based on CFL
    del_t = 0.98 / (c * ((1.0/del_x)**2 + (1.0/del_y)**2 + (1.0/del_z)**2)**0.5)

    # PML equation parameters
    sig_x_max = 0.75 * (0.8*(m+1)/(del_x*(mu_0/ep_0*1.0)**0.5))
    sig_y_max = 0.75 * (0.8*(m+1)/(del_y*(mu_0/ep_0*1.0)**0.5))
    sig_z_max = 0.75 * (0.8*(m+1)/(del_z*(mu_0/ep_0*1.0)**0.5))
    alpha_x_max = alpha_y_max = alpha_z_max = 0.1
    kappa_x_max = kappa_y_max = kappa_z_max = 1

    # Far field distance - only used in time delay of array element placement, so actally distance doesn't have to be the actual far field
    # it only needs to exceed a certain value in the far field time shifts associated with simulation distances in the interior
    r_for_time_relay=1.0

    # shift sizes up by 1 because the actual simulation only goes to size-1 in all directions
    x_size=x_size+1
    if (type_sim==0):
        pbc_shift=1
        y_size=y_size+2
        z_size=z_size+2
    if (type_sim==1):
        y_size=y_size+1
        z_size=z_size+1

    # this is slight overkill of condition - it can technically be a little smaller on lHS of >
    # additionally, the rhs of the = sign could be less too...
    # this updates if 1.0 above is too small
    if ((del_x*x_size+del_y*y_size+del_z*z_size)>r_for_time_relay):
        r_for_time_relay=(del_x*x_size+del_y*y_size+del_z*z_size)*2

    # This is slightly larger than it needs to be, I think, but is fine - not a huge time or memory constraint
    # this is because not all array elements are updated - only those that are used. 
    # will revisit this later to improve code maybe
    len_far_field_arrays=int(time_steps+(r_for_time_relay+del_x*x_size+del_y*y_size+del_z*z_size)/(c*del_t))
        
    # size needed for no signal coming back into the simultion + a 100 buffer zone
    # based on speed in a vacuum to get end and back, so should be good indefinitely.
    port_array_size=int(time_steps*del_t*c/min(del_x,del_y,del_z)/2.0+100)

    #this will need a replacement - this worked for FORTRAN but here I'll need to just omit those PML updates in the fdtd looping with an if statement or similar.
    #if (type_sim==0):
    #    nyPML_1=0
    #    nyPML_2=0
    #    nzPML_1=0
    #    nzPML_2=0

    # for videos if they are asked for
    if video_on==1:
        if (v_slice==0):
                vid_size1=y_size
                vid_size2=z_size
        if (v_slice==1):
            vid_size1=x_size
            vid_size2=z_size
            slice_location=slice_location+pbc_shift
        if (v_slice==2):
            vid_size1=x_size
            vid_size2=y_size
            slice_location=slice_location+pbc_shift

    # Setup pulse information
    if pulse_type == 1:
        spread = 1 / (2 * 3.141 * f_center)
        t_spread = 5.0 * spread
    elif pulse_type == 2:
        spread = 1 / (2 * 3.141 * f_center)
        t_spread = 6.0 * spread

    # now begin any calculations of variables
    # convert to radians
    theta=pi/180.0*(theta)
    phi=pi/180.0*(phi)
    pol=pi/180.0*(pol)

    # mirror conversion
    theta_mirror=pi/180.0*(theta_mirror)
    phi_mirror=pi/180.0*(phi_mirror)
    
    # convert them to tensors now (both mirror and regular)
    theta = th.tensor(theta, dtype=dtype, device=device)
    phi = th.tensor(phi, dtype=dtype, device=device)
    pol = th.tensor(pol, dtype=dtype, device=device)

    theta_mirror = th.tensor(theta_mirror, dtype=dtype, device=device)
    phi_mirror = th.tensor(phi_mirror, dtype=dtype, device=device)

    # far field data uses coordinate system centered in the middle of the simulation space
    # currently only used in the PML version, not pbc version - pbc not present here anyway
    ic=int((x_size-1)/2.0) 
    jc=int((y_size-1)/2.0)
    kc=int((z_size-1)/2.0)

    # These are used in a number of places - mostly far fields
    if (type_sim==0):
        xlow=nxPML_1+buffer
        xhigh=x_size-nxPML_2-buffer
    if (type_sim==1):
        xlow=nxPML_1+buffer
        xhigh=x_size-nxPML_2-buffer
        ylow=nyPML_1+buffer
        yhigh=y_size-nyPML_2-buffer   
        zlow=nzPML_1+buffer
        zhigh=z_size-nzPML_2-buffer

    # if mirror is used: there are several parameters I use that are now ready to be set before moving on.
    # based on direction of which plane +-x,y,z
    # recall mirror types are 1,2,3,4,5,6 for -x,+x,-y,+y,-z,+z
    if (is_mirror==1):
        if mirror_type==1:
            xlow=mirror_height
            xlow_wall=0
            x_mirror_offset = 2.0 * (xlow)
            ic=xlow
            nxPML_1=0
            use_x_mirror=1.0
        if mirror_type==2:
            xhigh=mirror_height
            xhigh_wall=0
            x_mirror_offset = 2.0 * (xhigh)
            ic=xhigh
            nxPML_2=0
            use_x_mirror=1.0
        if mirror_type==3:
            ylow=mirror_height
            ylow_wall=0
            y_mirror_offset = 2.0 * (ylow)
            jc=ylow
            nyPML_1=0
            use_y_mirror=1.0
        if mirror_type==4:
            yhigh=mirror_height
            yhigh_wall=0
            y_mirror_offset = 2.0 * (yhigh)
            jc=yhigh
            nyPML_2=0
            use_y_mirror=1.0
        if mirror_type==5:
            zlow=mirror_height
            zlow_wall=0
            z_mirror_offset = 2.0 * (zlow)
            kc=zlow
            nzPML_1=0
            use_z_mirror=1.0
        if mirror_type==6:
            zhigh=mirror_height
            zhigh_wall=0
            z_mirror_offset = 2.0 * (zhigh)
            kc=zhigh
            nzPML_2=0
            use_z_mirror=1.0

    # delays for plane waves
    # varying with incident angles, so this is a lookup table
    if ((theta>=0.0) and (theta<=(pi/2.0))):
        if ((phi>=0.0) and (phi<=(pi/2.0))):
            x_delay=xlow
            y_delay=ylow
            z_delay=zlow
        if ((phi>(pi/2.0)) and (phi<=pi)):
            x_delay=xhigh
            y_delay=ylow
            z_delay=zlow
        if ((phi>pi) and (phi<=(3.0*pi/2.0))):
            x_delay=xhigh
            y_delay=yhigh
            z_delay=zlow
        if ((phi>(3.0*pi/2.0)) and (phi<(2.0*pi))):
            x_delay=xlow
            y_delay=yhigh
            z_delay=zlow
    if ((theta>(pi/2.0)) and (theta<=pi)):
        if ((phi>=0) and (phi<=(pi/2.0))):
            x_delay=xlow
            y_delay=ylow
            z_delay=zhigh
        if ((phi>(pi/2.0)) and (phi<=pi)):
            x_delay=xhigh
            y_delay=ylow
            z_delay=zhigh
        if ((phi>pi) and (phi<=(3.0*pi/2.0))):
            x_delay=xhigh
            y_delay=yhigh
            z_delay=zhigh
        if ((phi>(3.0*pi/2.0)) and (phi<(2.0*pi))):
            x_delay=xlow
            y_delay=yhigh
            z_delay=zhigh

    # weights for plane waves
    # default is 1 V/m but user can modify
    WHx=plane_wave_amp*(th.sin(pol)*th.sin(phi)+th.cos(pol)*th.cos(theta)*th.cos(phi))*(ep_0/mu_0)**(0.5)
    WHy=plane_wave_amp*(-1*th.sin(pol)*th.cos(phi)+th.cos(pol)*th.cos(theta)*th.sin(phi))*(ep_0/mu_0)**(0.5)
    WHz=plane_wave_amp*(-1*th.cos(pol)*th.sin(theta))*(ep_0/mu_0)**(0.5)
    WEx=plane_wave_amp*(th.cos(pol)*th.sin(phi)-th.sin(pol)*th.cos(theta)*th.cos(phi))
    WEy=plane_wave_amp*(-1*th.cos(pol)*th.cos(phi)-th.sin(pol)*th.cos(theta)*th.sin(phi))
    WEz=plane_wave_amp*(th.sin(pol)*th.sin(theta))

    # delay_mirrors for mirrored plane waves
    # varying with incident angles, so this is a lookup table
    if ((theta_mirror>=0.0) and (theta_mirror<=(pi/2.0))):
        if ((phi_mirror>=0.0) and (phi_mirror<=(pi/2.0))):
            x_delay_mirror=xlow
            y_delay_mirror=ylow
            z_delay_mirror=zlow
        if ((phi_mirror>(pi/2.0)) and (phi_mirror<=pi)):
            x_delay_mirror=xhigh
            y_delay_mirror=ylow
            z_delay_mirror=zlow
        if ((phi_mirror>pi) and (phi_mirror<=(3.0*pi/2.0))):
            x_delay_mirror=xhigh
            y_delay_mirror=yhigh
            z_delay_mirror=zlow
        if ((phi_mirror>(3.0*pi/2.0)) and (phi_mirror<(2.0*pi))):
            x_delay_mirror=xlow
            y_delay_mirror=yhigh
            z_delay_mirror=zlow
    if ((theta_mirror>(pi/2.0)) and (theta_mirror<=pi)):
        if ((phi_mirror>=0) and (phi_mirror<=(pi/2.0))):
            x_delay_mirror=xlow
            y_delay_mirror=ylow
            z_delay_mirror=zhigh
        if ((phi_mirror>(pi/2.0)) and (phi_mirror<=pi)):
            x_delay_mirror=xhigh
            y_delay_mirror=ylow
            z_delay_mirror=zhigh
        if ((phi_mirror>pi) and (phi_mirror<=(3.0*pi/2.0))):
            x_delay_mirror=xhigh
            y_delay_mirror=yhigh
            z_delay_mirror=zhigh
        if ((phi_mirror>(3.0*pi/2.0)) and (phi_mirror<(2.0*pi))):
            x_delay_mirror=xlow
            y_delay_mirror=yhigh
            z_delay_mirror=zhigh

    # specific to mirror delays only.
    # we are now ready to define more mirror variables if needed.
    # need to shift the delay so origin makes PEC plane halfway.
    # recall mirror types are 1,2,3,4,5,6 for -x,+x,-y,+y,-z,+z
    if (is_mirror==1):
        if mirror_type==1:
            x_delay_mirror=mirror_height+(mirror_height-xhigh)
        if mirror_type==2:
            x_delay_mirror=mirror_height+(mirror_height-xlow)
        if mirror_type==3:
            y_delay_mirror=mirror_height+(mirror_height-yhigh)
        if mirror_type==4:
            y_delay_mirror=mirror_height+(mirror_height-ylow)
        if mirror_type==5:
            z_delay_mirror=mirror_height+(mirror_height-zhigh)
        if mirror_type==6:
            z_delay_mirror=mirror_height+(mirror_height-zlow)

    # weights for mirrored plane waves if needed
    # recall mirror types are 1,2,3,4,5,6 for -x,+x,-y,+y,-z,+z
    if (is_mirror==1):
        if ((mirror_type==1) or (mirror_type==2)): #yz-plane is ground, x is normal

            WEx_mirror = WEx         # Normal to ground - not flipped for PEC
            WEy_mirror = -1.0 * WEy  # Tangential - flip for PEC
            WEz_mirror = -1.0 * WEz  # Tangential - flip for PEC
            WHx_mirror = -1.0 * WHx  # Normal - flip naturally to cancel
            WHy_mirror = WHy         # Tangential - don't flip
            WHz_mirror = WHz         # Tangential - don't flip

            #Opposite from E,H because of cross product
            Jy_mirror=-1.0
            Jz_mirror=-1.0
            Mx_mirror=-1.0

        if ((mirror_type==3) or (mirror_type==4)): #xz-plane is ground, y is normal

            WEx_mirror = -1.0 * WEx  # Tangential - flip
            WEy_mirror = WEy         # Normal - not flipped
            WEz_mirror = -1.0 * WEz  # Tangential - flip
            WHx_mirror = WHx         # Tangential - don't flip
            WHy_mirror = -1.0 * WHy  # Normal - flip naturally
            WHz_mirror = WHz         # Tangential - don't flip

            #Opposite from E,H because of cross product
            Jx_mirror=-1.0
            Jz_mirror=-1.0
            My_mirror=-1.0

        if ((mirror_type==5) or (mirror_type==6)): #xy-plane is ground, z is normal

            WEx_mirror = -1.0 * WEx  # Tangential - flip
            WEy_mirror = -1.0 * WEy  # Tangential - flip
            WEz_mirror = WEz         # Normal - not flipped
            WHx_mirror = WHx         # Tangential - don't flip
            WHy_mirror = WHy         # Tangential - don't flip
            WHz_mirror = -1.0 * WHz  # Normal - flip naturally

            #Opposite from E,H because of cross product
            Jx_mirror=-1.0
            Jy_mirror=-1.0
            Mz_mirror=-1.0

    #Set far field looping bounds
    ff_xlow=xlow-1
    ff_xhigh=xhigh+2
    ff_ylow=ylow-1
    ff_yhigh=yhigh+2
    ff_zlow=zlow-1
    ff_zhigh=zhigh+2

    # if mirror: we are changing how we handle one of the walls
    # this ensures the 4 perpendicular walls will only loop until the pec plane, not on or through it.
    # recall mirror types are 1,2,3,4,5,6 for -x,+x,-y,+y,-z,+z
    if (mirror_type==1):
        ff_xlow=xlow+1
    if (mirror_type==2):
        ff_xhigh=xhigh-1
    if (mirror_type==3):
        ff_ylow=ylow+1
    if (mirror_type==4):
        ff_yhigh=yhigh-1
    if (mirror_type==5):
        ff_zlow=zlow+1
    if (mirror_type==6):
        ff_zhigh=zhigh-1

    # Establish minimum steps required - but likely needs way more
    min_steps=2*t_spread/del_t+th.sqrt(((th.sin(theta)*th.cos(phi)*(x_size)*del_x)**2 \
    +(th.sin(theta)*th.sin(phi)*(y_size)*del_y)**2+(th.cos(theta)*(z_size)*del_z)**2))/(c*del_t)

##########################################################
# FDTD INITIALIZATION
##########################################################
with GPUProfiler("FDTD Initialization", device):

    # Initialize all arrays ON GPU with proper dtype
    
    # Used for sizing several options below
    Field_shape = (x_size, y_size, z_size)

    # Main E,H Fields on GPU
    Ex = th.zeros(Field_shape, dtype=dtype, device=device)
    Ey = th.zeros(Field_shape, dtype=dtype, device=device)
    Ez = th.zeros(Field_shape, dtype=dtype, device=device)
    Hx = th.zeros(Field_shape, dtype=dtype, device=device)
    Hy = th.zeros(Field_shape, dtype=dtype, device=device)
    Hz = th.zeros(Field_shape, dtype=dtype, device=device)

    # Material cells (relative_ep_X_cell, sigma_X_cell)
    relative_ep_x_cell = th.ones(Field_shape, dtype=dtype, device=device)
    relative_ep_y_cell = th.ones(Field_shape, dtype=dtype, device=device)
    relative_ep_z_cell = th.ones(Field_shape, dtype=dtype, device=device)
    sigma_x_cell = th.zeros(Field_shape, dtype=dtype, device=device)
    sigma_y_cell = th.zeros(Field_shape, dtype=dtype, device=device)
    sigma_z_cell = th.zeros(Field_shape, dtype=dtype, device=device)

    # Grid points (relative_ep_X, sigma_X)
    relative_ep_x = th.ones(Field_shape, dtype=dtype, device=device)
    relative_ep_y = th.ones(Field_shape, dtype=dtype, device=device)
    relative_ep_z = th.ones(Field_shape, dtype=dtype, device=device)
    sigma_x = th.zeros(Field_shape, dtype=dtype, device=device)
    sigma_y = th.zeros(Field_shape, dtype=dtype, device=device)
    sigma_z = th.zeros(Field_shape, dtype=dtype, device=device)

    # For permeability (scalars for now, but will modify at some point)
    da = 1.0
    db = del_t / mu_0
    
    # Field aux for permittivity and electrical conductivity on GPU
    gax = th.zeros(Field_shape, dtype=dtype, device=device)
    gbx = th.zeros(Field_shape, dtype=dtype, device=device)
    gay = th.zeros(Field_shape, dtype=dtype, device=device)
    gby = th.zeros(Field_shape, dtype=dtype, device=device)
    gaz = th.zeros(Field_shape, dtype=dtype, device=device)
    gbz = th.zeros(Field_shape, dtype=dtype, device=device)
    
    # Initialize psi arrays on GPU
    psi_Ezx_1 = th.zeros((nxPML_1, y_size, z_size), dtype=dtype, device=device)
    psi_Ezx_2 = th.zeros((nxPML_2, y_size, z_size), dtype=dtype, device=device)
    psi_Hyx_1 = th.zeros((nxPML_1-1, y_size, z_size), dtype=dtype, device=device)
    psi_Hyx_2 = th.zeros((nxPML_2-1, y_size, z_size), dtype=dtype, device=device)
    psi_Ezy_1 = th.zeros((x_size, nyPML_1, z_size), dtype=dtype, device=device)
    psi_Ezy_2 = th.zeros((x_size, nyPML_2, z_size), dtype=dtype, device=device)
    psi_Hxy_1 = th.zeros((x_size, nyPML_1-1, z_size), dtype=dtype, device=device)
    psi_Hxy_2 = th.zeros((x_size, nyPML_2-1, z_size), dtype=dtype, device=device)
    psi_Hxz_1 = th.zeros((x_size, y_size-1, nzPML_1-1), dtype=dtype, device=device)
    psi_Hxz_2 = th.zeros((x_size, y_size-1, nzPML_2-1), dtype=dtype, device=device)
    psi_Hyz_1 = th.zeros((x_size-1, y_size, nzPML_1-1), dtype=dtype, device=device)
    psi_Hyz_2 = th.zeros((x_size-1, y_size, nzPML_2-1), dtype=dtype, device=device)
    psi_Exz_1 = th.zeros((x_size-1, y_size, nzPML_1), dtype=dtype, device=device)
    psi_Exz_2 = th.zeros((x_size-1, y_size, nzPML_2), dtype=dtype, device=device)
    psi_Eyz_1 = th.zeros((x_size, y_size-1, nzPML_1), dtype=dtype, device=device)
    psi_Eyz_2 = th.zeros((x_size, y_size-1, nzPML_2), dtype=dtype, device=device)
    psi_Hzx_1 = th.zeros((nxPML_1-1, y_size-1, z_size), dtype=dtype, device=device)
    psi_Eyx_1 = th.zeros((nxPML_1, y_size-1, z_size), dtype=dtype, device=device)
    psi_Hzx_2 = th.zeros((nxPML_2-1, y_size-1, z_size), dtype=dtype, device=device)
    psi_Eyx_2 = th.zeros((nxPML_2, y_size-1, z_size), dtype=dtype, device=device)
    psi_Hzy_1 = th.zeros((x_size-1, nyPML_1-1, z_size), dtype=dtype, device=device)
    psi_Exy_1 = th.zeros((x_size-1, nyPML_1, z_size), dtype=dtype, device=device)
    psi_Hzy_2 = th.zeros((x_size-1, nyPML_2-1, z_size), dtype=dtype, device=device)
    psi_Exy_2 = th.zeros((x_size-1, nyPML_2, z_size), dtype=dtype, device=device)
    
    # Initialize PML parameter arrays on GPU
    alphae_x_PML_1 = th.zeros(nxPML_1, dtype=dtype, device=device)
    sige_x_PML_1 = th.zeros(nxPML_1, dtype=dtype, device=device)
    kappae_x_PML_1 = th.zeros(nxPML_1, dtype=dtype, device=device)
    alphah_x_PML_1 = th.zeros(nxPML_1-1, dtype=dtype, device=device)
    sigh_x_PML_1 = th.zeros(nxPML_1-1, dtype=dtype, device=device)
    kappah_x_PML_1 = th.zeros(nxPML_1-1, dtype=dtype, device=device)
    
    alphae_x_PML_2 = th.zeros(nxPML_2, dtype=dtype, device=device)
    sige_x_PML_2 = th.zeros(nxPML_2, dtype=dtype, device=device)
    kappae_x_PML_2 = th.zeros(nxPML_2, dtype=dtype, device=device)
    alphah_x_PML_2 = th.zeros(nxPML_2-1, dtype=dtype, device=device)
    sigh_x_PML_2 = th.zeros(nxPML_2-1, dtype=dtype, device=device)
    kappah_x_PML_2 = th.zeros(nxPML_2-1, dtype=dtype, device=device)
    
    alphae_y_PML_1 = th.zeros(nyPML_1, dtype=dtype, device=device)
    sige_y_PML_1 = th.zeros(nyPML_1, dtype=dtype, device=device)
    kappae_y_PML_1 = th.zeros(nyPML_1, dtype=dtype, device=device)
    alphah_y_PML_1 = th.zeros(nyPML_1-1, dtype=dtype, device=device)
    sigh_y_PML_1 = th.zeros(nyPML_1-1, dtype=dtype, device=device)
    kappah_y_PML_1 = th.zeros(nyPML_1-1, dtype=dtype, device=device)
    
    alphae_y_PML_2 = th.zeros(nyPML_2, dtype=dtype, device=device)
    sige_y_PML_2 = th.zeros(nyPML_2, dtype=dtype, device=device)
    kappae_y_PML_2 = th.zeros(nyPML_2, dtype=dtype, device=device)
    alphah_y_PML_2 = th.zeros(nyPML_2-1, dtype=dtype, device=device)
    sigh_y_PML_2 = th.zeros(nyPML_2-1, dtype=dtype, device=device)
    kappah_y_PML_2 = th.zeros(nyPML_2-1, dtype=dtype, device=device)
    
    alphae_z_PML_1 = th.zeros(nzPML_1, dtype=dtype, device=device)
    sige_z_PML_1 = th.zeros(nzPML_1, dtype=dtype, device=device)
    kappae_z_PML_1 = th.zeros(nzPML_1, dtype=dtype, device=device)
    alphah_z_PML_1 = th.zeros(nzPML_1-1, dtype=dtype, device=device)
    sigh_z_PML_1 = th.zeros(nzPML_1-1, dtype=dtype, device=device)
    kappah_z_PML_1 = th.zeros(nzPML_1-1, dtype=dtype, device=device)
    
    alphae_z_PML_2 = th.zeros(nzPML_2, dtype=dtype, device=device)
    sige_z_PML_2 = th.zeros(nzPML_2, dtype=dtype, device=device)
    kappae_z_PML_2 = th.zeros(nzPML_2, dtype=dtype, device=device)
    alphah_z_PML_2 = th.zeros(nzPML_2-1, dtype=dtype, device=device)
    sigh_z_PML_2 = th.zeros(nzPML_2-1, dtype=dtype, device=device)
    kappah_z_PML_2 = th.zeros(nzPML_2-1, dtype=dtype, device=device)
    
    # Initialize b and c coefficient arrays on GPU
    be_x_1 = th.zeros((nxPML_1,y_size,z_size), dtype=dtype, device=device)
    ce_x_1 = th.zeros((nxPML_1,y_size,z_size), dtype=dtype, device=device)
    bh_x_1 = th.zeros((nxPML_1-1,y_size,z_size), dtype=dtype, device=device)
    ch_x_1 = th.zeros((nxPML_1-1,y_size,z_size), dtype=dtype, device=device)

    be_x_2 = th.zeros((nxPML_2,y_size,z_size), dtype=dtype, device=device)
    ce_x_2 = th.zeros((nxPML_2,y_size,z_size), dtype=dtype, device=device)
    bh_x_2 = th.zeros((nxPML_2-1,y_size,z_size), dtype=dtype, device=device)
    ch_x_2 = th.zeros((nxPML_2-1,y_size,z_size), dtype=dtype, device=device)

    be_y_1 = th.zeros((x_size,nyPML_1,z_size), dtype=dtype, device=device)
    ce_y_1 = th.zeros((x_size,nyPML_1,z_size), dtype=dtype, device=device)
    bh_y_1 = th.zeros((x_size,nyPML_1-1,z_size), dtype=dtype, device=device)
    ch_y_1 = th.zeros((x_size,nyPML_1-1,z_size), dtype=dtype, device=device)

    be_y_2 = th.zeros((x_size,nyPML_2,z_size), dtype=dtype, device=device)
    ce_y_2 = th.zeros((x_size,nyPML_2,z_size), dtype=dtype, device=device)
    bh_y_2 = th.zeros((x_size,nyPML_2-1,z_size), dtype=dtype, device=device)
    ch_y_2 = th.zeros((x_size,nyPML_2-1,z_size), dtype=dtype, device=device)

    be_z_1 = th.zeros((x_size,y_size,nzPML_1), dtype=dtype, device=device)
    ce_z_1 = th.zeros((x_size,y_size,nzPML_1), dtype=dtype, device=device)
    bh_z_1 = th.zeros((x_size,y_size,nzPML_1-1), dtype=dtype, device=device)
    ch_z_1 = th.zeros((x_size,y_size,nzPML_1-1), dtype=dtype, device=device)

    be_z_2 = th.zeros((x_size,y_size,nzPML_2), dtype=dtype, device=device)
    ce_z_2 = th.zeros((x_size,y_size,nzPML_2), dtype=dtype, device=device)
    bh_z_2 = th.zeros((x_size,y_size,nzPML_2-1), dtype=dtype, device=device)
    ch_z_2 = th.zeros((x_size,y_size,nzPML_2-1), dtype=dtype, device=device)

    # Denominators for update equations (pre-compute on GPU)
    den_ex = th.full((x_size-1,), 1/del_x, dtype=dtype, device=device)
    den_ey = th.full((y_size-1,), 1/del_y, dtype=dtype, device=device)
    den_ez = th.full((z_size-1,), 1/del_z, dtype=dtype, device=device)
    den_hx = th.full((x_size-1,), 1/del_x, dtype=dtype, device=device)
    den_hy = th.full((y_size-1,), 1/del_y, dtype=dtype, device=device)
    den_hz = th.full((z_size-1,), 1/del_z, dtype=dtype, device=device)

    # Voltage and Current Aux Arrays if needed
    if num_ports>0:
        Voltage = th.zeros((num_ports, port_array_size), dtype=dtype, device=device)
        Current = th.zeros((num_ports, port_array_size), dtype=dtype, device=device)
        Voltage_out = th.zeros((num_ports, time_steps), dtype=dtype, device=device)
        Current_out = th.zeros((num_ports, time_steps), dtype=dtype, device=device)

    # Video Arrays (3D Fields) of needed
    if video_on=='yes':
        Ex_video = th.zeros((vid_size1, vid_size2, time_steps), dtype=dtype, device=device)
        Hx_video = th.zeros((vid_size1, vid_size2, time_steps), dtype=dtype, device=device)
        Ey_video = th.zeros((vid_size1, vid_size2, time_steps), dtype=dtype, device=device)
        Hy_video = th.zeros((vid_size1, vid_size2, time_steps), dtype=dtype, device=device)
        Ez_video = th.zeros((vid_size1, vid_size2, time_steps), dtype=dtype, device=device)
        Hz_video = th.zeros((vid_size1, vid_size2, time_steps), dtype=dtype, device=device)

    # Incident and Output Field Arrays (1D) ---
    incident = th.zeros((time_steps,), dtype=dtype, device=device)
    E_reflected = th.zeros((time_steps,), dtype=dtype, device=device)
    E_transmitted = th.zeros((time_steps,), dtype=dtype, device=device)

    # Far Field Post Processing Arrays (2D Boundary Surfaces) if needed
    if num_far_field_angles > 0:
        # x-low/x-high faces (depend on y_size, z_size)
        My_xlow = th.zeros((y_size, z_size), dtype=dtype, device=device)
        Mz_xlow = th.zeros((y_size, z_size), dtype=dtype, device=device)
        Jy_xlow = th.zeros((y_size, z_size), dtype=dtype, device=device)
        Jz_xlow = th.zeros((y_size, z_size), dtype=dtype, device=device)
        My_xlow_oldt = th.zeros((y_size, z_size), dtype=dtype, device=device)
        Mz_xlow_oldt = th.zeros((y_size, z_size), dtype=dtype, device=device)
        Jy_xlow_oldt = th.zeros((y_size, z_size), dtype=dtype, device=device)
        Jz_xlow_oldt = th.zeros((y_size, z_size), dtype=dtype, device=device)

        My_xhigh = th.zeros((y_size, z_size), dtype=dtype, device=device)
        Mz_xhigh = th.zeros((y_size, z_size), dtype=dtype, device=device)
        Jy_xhigh = th.zeros((y_size, z_size), dtype=dtype, device=device)
        Jz_xhigh = th.zeros((y_size, z_size), dtype=dtype, device=device)
        My_xhigh_oldt = th.zeros((y_size, z_size), dtype=dtype, device=device)
        Mz_xhigh_oldt = th.zeros((y_size, z_size), dtype=dtype, device=device)
        Jy_xhigh_oldt = th.zeros((y_size, z_size), dtype=dtype, device=device)
        Jz_xhigh_oldt = th.zeros((y_size, z_size), dtype=dtype, device=device)

        # y-low/y-high faces (depend on x_size, z_size)
        Mx_ylow = th.zeros((x_size, z_size), dtype=dtype, device=device)
        Mz_ylow = th.zeros((x_size, z_size), dtype=dtype, device=device)
        Jx_ylow = th.zeros((x_size, z_size), dtype=dtype, device=device)
        Jz_ylow = th.zeros((x_size, z_size), dtype=dtype, device=device)
        Mx_ylow_oldt = th.zeros((x_size, z_size), dtype=dtype, device=device)
        Mz_ylow_oldt = th.zeros((x_size, z_size), dtype=dtype, device=device)
        Jx_ylow_oldt = th.zeros((x_size, z_size), dtype=dtype, device=device)
        Jz_ylow_oldt = th.zeros((x_size, z_size), dtype=dtype, device=device)

        Mx_yhigh = th.zeros((x_size, z_size), dtype=dtype, device=device)
        Mz_yhigh = th.zeros((x_size, z_size), dtype=dtype, device=device)
        Jx_yhigh = th.zeros((x_size, z_size), dtype=dtype, device=device)
        Jz_yhigh = th.zeros((x_size, z_size), dtype=dtype, device=device)
        Mx_yhigh_oldt = th.zeros((x_size, z_size), dtype=dtype, device=device)
        Mz_yhigh_oldt = th.zeros((x_size, z_size), dtype=dtype, device=device)
        Jx_yhigh_oldt = th.zeros((x_size, z_size), dtype=dtype, device=device)
        Jz_yhigh_oldt = th.zeros((x_size, z_size), dtype=dtype, device=device)

        # z-low/z-high faces (depend on x_size, y_size)
        Mx_zlow = th.zeros((x_size, y_size), dtype=dtype, device=device)
        My_zlow = th.zeros((x_size, y_size), dtype=dtype, device=device)
        Jx_zlow = th.zeros((x_size, y_size), dtype=dtype, device=device)
        Jy_zlow = th.zeros((x_size, y_size), dtype=dtype, device=device)
        Mx_zlow_oldt = th.zeros((x_size, y_size), dtype=dtype, device=device)
        My_zlow_oldt = th.zeros((x_size, y_size), dtype=dtype, device=device)
        Jx_zlow_oldt = th.zeros((x_size, y_size), dtype=dtype, device=device)
        Jy_zlow_oldt = th.zeros((x_size, y_size), dtype=dtype, device=device)

        Mx_zhigh = th.zeros((x_size, y_size), dtype=dtype, device=device)
        My_zhigh = th.zeros((x_size, y_size), dtype=dtype, device=device)
        Jx_zhigh = th.zeros((x_size, y_size), dtype=dtype, device=device)
        Jy_zhigh = th.zeros((x_size, y_size), dtype=dtype, device=device)
        Mx_zhigh_oldt = th.zeros((x_size, y_size), dtype=dtype, device=device)
        My_zhigh_oldt = th.zeros((x_size, y_size), dtype=dtype, device=device)
        Jx_zhigh_oldt = th.zeros((x_size, y_size), dtype=dtype, device=device)
        Jy_zhigh_oldt = th.zeros((x_size, y_size), dtype=dtype, device=device)

        # Far Field Summation and Output Tensors (W, U, E)
        W_shape = (num_far_field_angles, len_far_field_arrays)
        E_out_shape = (num_far_field_angles, time_steps)

        Wx = th.zeros(W_shape, dtype=dtype, device=device)
        Wy = th.zeros(W_shape, dtype=dtype, device=device)
        Wz = th.zeros(W_shape, dtype=dtype, device=device)
        Ux = th.zeros(W_shape, dtype=dtype, device=device)
        Uy = th.zeros(W_shape, dtype=dtype, device=device)
        Uz = th.zeros(W_shape, dtype=dtype, device=device)
        W_theta = th.zeros(W_shape, dtype=dtype, device=device)
        W_phi = th.zeros(W_shape, dtype=dtype, device=device)
        U_theta = th.zeros(W_shape, dtype=dtype, device=device)
        U_phi = th.zeros(W_shape, dtype=dtype, device=device)

        # Combined outputs
        E_theta = th.zeros(W_shape, dtype=dtype, device=device)
        E_phi = th.zeros(W_shape, dtype=dtype, device=device)
        E_theta_out = th.zeros(E_out_shape, dtype=dtype, device=device)
        E_phi_out = th.zeros(E_out_shape, dtype=dtype, device=device)
        data_out_time = th.zeros((num_far_field_angles,), dtype=dtype, device=device)

    # Special Field and Sheet List Arrays
    Ex_special = th.zeros(Field_shape, dtype=dtype, device=device)
    Ey_special = th.zeros(Field_shape, dtype=dtype, device=device)
    Ez_special = th.zeros(Field_shape, dtype=dtype, device=device)

    x_sheet_list = th.zeros((x_size - 1,), dtype=dtype, device=device)
    y_sheet_list = th.zeros((y_size - 1,), dtype=dtype, device=device)
    z_sheet_list = th.zeros((z_size - 1,), dtype=dtype, device=device)

    # Sheet Material Arrays (x-normal)
    sheet_ep_x_cell_x = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_ep_y_cell_x = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_ep_z_cell_x = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_sig_x_cell_x = th.zeros(Field_shape, dtype=dtype, device=device)
    sheet_sig_y_cell_x = th.zeros(Field_shape, dtype=dtype, device=device)
    sheet_sig_z_cell_x = th.zeros(Field_shape, dtype=dtype, device=device)

    sheet_ep_x_x = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_ep_y_x = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_ep_z_x = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_sig_x_x = th.zeros(Field_shape, dtype=dtype, device=device)
    sheet_sig_y_x = th.zeros(Field_shape, dtype=dtype, device=device)
    sheet_sig_z_x = th.zeros(Field_shape, dtype=dtype, device=device)

    # Sheet Material Arrays (y-normal)
    sheet_ep_x_cell_y = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_ep_y_cell_y = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_ep_z_cell_y = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_sig_x_cell_y = th.zeros(Field_shape, dtype=dtype, device=device)
    sheet_sig_y_cell_y = th.zeros(Field_shape, dtype=dtype, device=device)
    sheet_sig_z_cell_y = th.zeros(Field_shape, dtype=dtype, device=device)

    sheet_ep_x_y = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_ep_y_y = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_ep_z_y = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_sig_x_y = th.zeros(Field_shape, dtype=dtype, device=device)
    sheet_sig_y_y = th.zeros(Field_shape, dtype=dtype, device=device)
    sheet_sig_z_y = th.zeros(Field_shape, dtype=dtype, device=device)

    # Sheet Material Arrays (z-normal)
    sheet_ep_x_cell_z = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_ep_y_cell_z = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_ep_z_cell_z = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_sig_x_cell_z = th.zeros(Field_shape, dtype=dtype, device=device)
    sheet_sig_y_cell_z = th.zeros(Field_shape, dtype=dtype, device=device)
    sheet_sig_z_cell_z = th.zeros(Field_shape, dtype=dtype, device=device)

    sheet_ep_x_z = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_ep_y_z = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_ep_z_z = th.ones(Field_shape, dtype=dtype, device=device)
    sheet_sig_x_z = th.zeros(Field_shape, dtype=dtype, device=device)
    sheet_sig_y_z = th.zeros(Field_shape, dtype=dtype, device=device)
    sheet_sig_z_z = th.zeros(Field_shape, dtype=dtype, device=device)


    # Initialize arrays for saving data
    input_data = th.zeros(time_steps, dtype=dtype, device=device)
    
    # Pre-compute pulse time values
    time_values = th.arange(time_steps, dtype=dtype, device=device) * del_t

##########################################################
# PML SETUP
##########################################################
with GPUProfiler("PML Calculations", device):
    #setup the parameters - probably a faster/better way to do this long term but copied and slightly modified from fortran
    #the PML are only 10 in size so it's pretty fast regardless of methods and cpu/gpu settings
    for i in range(nxPML_1):
        sige_x_PML_1[i] = sig_x_max * ( (nxPML_1 - (i+1)) / (nxPML_1 - 1.0) )**m
        alphae_x_PML_1[i] = alpha_x_max*(((i+1)-1)/(nxPML_1-1.0))**ma
        kappae_x_PML_1[i] = 1.0+(kappa_x_max-1.0)*((nxPML_1 - (i+1)) / (nxPML_1 - 1.0))**m
        be_x_1[i,:,:] = math.exp(-(sige_x_PML_1[i] / kappae_x_PML_1[i] + alphae_x_PML_1[i])*del_t/ep_0)
        if ((sige_x_PML_1[i] == 0.0) and (alphae_x_PML_1[i] == 0.0) and (i == nxPML_1-1)):
            ce_x_1[i,:,:] = 0.0
        else:
            ce_x_1[i,:,:] = sige_x_PML_1[i]*(be_x_1[i,:,:]-1.0)/(sige_x_PML_1[i]+kappae_x_PML_1[i]*alphae_x_PML_1[i]) / kappae_x_PML_1[i]

    for i in range(nxPML_1-1):
        sigh_x_PML_1[i] = sig_x_max * ( (nxPML_1 - (i+1) - 0.5)/(nxPML_1-1.0))**m
        alphah_x_PML_1[i] = alpha_x_max*(((i+1)-0.5)/(nxPML_1-1.0))**ma
        kappah_x_PML_1[i] = 1.0+(kappa_x_max-1.0)*((nxPML_1 - (i+1) - 0.5) / (nxPML_1 - 1.0))**m
        bh_x_1[i,:,:] = math.exp(-(sigh_x_PML_1[i] / kappah_x_PML_1[i] + alphah_x_PML_1[i])*del_t/ep_0)
        ch_x_1[i,:,:] = sigh_x_PML_1[i]*(bh_x_1[i,:,:]-1.0)/(sigh_x_PML_1[i]+kappah_x_PML_1[i]*alphah_x_PML_1[i])/kappah_x_PML_1[i]

    for i in range(nxPML_2):
        sige_x_PML_2[i] = sig_x_max * ( (nxPML_2 - (i+1)) / (nxPML_2 - 1.0) )**m
        alphae_x_PML_2[i] = alpha_x_max*(((i+1)-1)/(nxPML_2-1.0))**ma
        kappae_x_PML_2[i] = 1.0+(kappa_x_max-1.0)*((nxPML_2 - (i+1)) / (nxPML_2 - 1.0))**m
        be_x_2[i,:,:] = math.exp(-(sige_x_PML_2[i] / kappae_x_PML_2[i] +alphae_x_PML_2[i])*del_t/ep_0)
        if ((sige_x_PML_2[i] == 0.0) and (alphae_x_PML_2[i] == 0.0) and (i == nxPML_2-1)):
            ce_x_2[i,:,:] = 0.0
        else:
            ce_x_2[i,:,:] = sige_x_PML_2[i]*(be_x_2[i,:,:]-1.0)/(sige_x_PML_2[i]+kappae_x_PML_2[i]*alphae_x_PML_2[i]) / kappae_x_PML_2[i]

    for i in range(nxPML_2-1):
        sigh_x_PML_2[i] = sig_x_max * ( (nxPML_2 - (i+1) - 0.5)/(nxPML_2-1.0))**m
        alphah_x_PML_2[i] = alpha_x_max*(((i+1)-0.5)/(nxPML_2-1.0))**ma
        kappah_x_PML_2[i] = 1.0+(kappa_x_max-1.0)*((nxPML_2 - (i+1) - 0.5) / (nxPML_2 - 1.0))**m
        bh_x_2[i,:,:] = math.exp(-(sigh_x_PML_2[i] / kappah_x_PML_2[i] +alphah_x_PML_2[i])*del_t/ep_0)
        ch_x_2[i,:,:] = sigh_x_PML_2[i]*(bh_x_2[i,:,:]-1.0)/(sigh_x_PML_2[i]+kappah_x_PML_2[i]*alphah_x_PML_2[i])/kappah_x_PML_2[i]

    for j in range(nyPML_1):
        sige_y_PML_1[j] = sig_y_max * ( (nyPML_1 - (j+1)) / (nyPML_1 - 1.0) )**m
        alphae_y_PML_1[j] = alpha_y_max*(((j+1)-1)/(nyPML_1-1.0))**ma
        kappae_y_PML_1[j] = 1.0+(kappa_y_max-1.0)*((nyPML_1 - (j+1)) / (nyPML_1 - 1.0))**m
        be_y_1[:,j,:] = math.exp(-(sige_y_PML_1[j] / kappae_y_PML_1[j] +alphae_y_PML_1[j])*del_t/ep_0)
        if ((sige_y_PML_1[j] == 0.0) and (alphae_y_PML_1[j] == 0.0) and (j == nyPML_1-1)):
            ce_y_1[:,j,:] = 0.0
        else:
            ce_y_1[:,j,:] = sige_y_PML_1[j]*(be_y_1[:,j,:]-1.0)/(sige_y_PML_1[j]+kappae_y_PML_1[j]*alphae_y_PML_1[j]) / kappae_y_PML_1[j]

    for j in range(nyPML_1-1):
        sigh_y_PML_1[j] = sig_y_max * ( (nyPML_1 - (j+1) - 0.5)/(nyPML_1-1.0))**m
        alphah_y_PML_1[j] = alpha_y_max*(((j+1)-0.5)/(nyPML_1-1.0))**ma
        kappah_y_PML_1[j] = 1.0+(kappa_y_max-1.0)*((nyPML_1 - (j+1) - 0.5) / (nyPML_1 - 1.0))**m
        bh_y_1[:,j,:] = math.exp(-(sigh_y_PML_1[j] / kappah_y_PML_1[j] +alphah_y_PML_1[j])*del_t/ep_0)
        ch_y_1[:,j,:] = sigh_y_PML_1[j]*(bh_y_1[:,j,:]-1.0)/(sigh_y_PML_1[j]+kappah_y_PML_1[j]*alphah_y_PML_1[j])/ kappah_y_PML_1[j]

    for j in range(nyPML_2):
        sige_y_PML_2[j] = sig_y_max * ( (nyPML_2 - (j+1)) / (nyPML_2 - 1.0) )**m
        alphae_y_PML_2[j] = alpha_y_max*(((j+1)-1)/(nyPML_2-1.0))**ma
        kappae_y_PML_2[j] = 1.0+(kappa_y_max-1.0)*((nyPML_2 - (j+1)) / (nyPML_2 - 1.0))**m
        be_y_2[:,j,:] = math.exp(-(sige_y_PML_2[j] / kappae_y_PML_2[j] +alphae_y_PML_2[j])*del_t/ep_0)
        if ((sige_y_PML_2[j] == 0.0) and(alphae_y_PML_2[j] == 0.0) and(j == nyPML_2-1)):
            ce_y_2[:,j,:] = 0.0
        else:
            ce_y_2[:,j,:] = sige_y_PML_2[j]*(be_y_2[:,j,:]-1.0)/(sige_y_PML_2[j]+kappae_y_PML_2[j]*alphae_y_PML_2[j])/kappae_y_PML_2[j]

    for j in range(nyPML_2-1):
        sigh_y_PML_2[j] = sig_y_max * ( (nyPML_2 - (j+1) - 0.5)/(nyPML_2-1.0))**m
        alphah_y_PML_2[j] = alpha_y_max*(((j+1)-0.5)/(nyPML_2-1.0))**ma
        kappah_y_PML_2[j] = 1.0+(kappa_y_max-1.0)*((nyPML_2 - (j+1) - 0.5) / (nyPML_2 - 1.0))**m
        bh_y_2[:,j,:] = math.exp(-(sigh_y_PML_2[j] / kappah_y_PML_2[j] +alphah_y_PML_2[j])*del_t/ep_0)
        ch_y_2[:,j,:] = sigh_y_PML_2[j]*(bh_y_2[:,j,:]-1.0)/(sigh_y_PML_2[j]+kappah_y_PML_2[j]*alphah_y_PML_2[j])/kappah_y_PML_2[j]

    for k in range(nzPML_1):
        sige_z_PML_1[k] = sig_z_max * ( (nzPML_1 - (k+1)) / (nzPML_1 - 1.0) )**m
        alphae_z_PML_1[k] = alpha_z_max*(((k+1)-1)/(nzPML_1-1.0))**ma
        kappae_z_PML_1[k] = 1.0+(kappa_z_max-1.0)*((nzPML_1 - (k+1)) / (nzPML_1 - 1.0))**m
        be_z_1[:,:,k] = math.exp(-(sige_z_PML_1[k] / kappae_z_PML_1[k] +alphae_z_PML_1[k])*del_t/ep_0)
        if ((sige_z_PML_1[k] == 0.0) and (alphae_z_PML_1[k] == 0.0) and(k == nzPML_1-1)):
            ce_z_1[:,:,k] = 0.0
        else:
            ce_z_1[:,:,k] = sige_z_PML_1[k]*(be_z_1[:,:,k]-1.0)/(sige_z_PML_1[k]+kappae_z_PML_1[k]*alphae_z_PML_1[k]) / kappae_z_PML_1[k]

    for k in range(nzPML_1-1):
        sigh_z_PML_1[k] = sig_z_max * ( (nzPML_1 - (k+1) - 0.5)/(nzPML_1-1.0))**m
        alphah_z_PML_1[k] = alpha_z_max*(((k+1)-0.5)/(nzPML_1-1.0))**ma
        kappah_z_PML_1[k] = 1.0+(kappa_z_max-1.0)*((nzPML_1 - (k+1) - 0.5) / (nzPML_1 - 1.0))**m
        bh_z_1[:,:,k] = math.exp(-(sigh_z_PML_1[k] / kappah_z_PML_1[k] +alphah_z_PML_1[k])*del_t/ep_0)
        ch_z_1[:,:,k] = sigh_z_PML_1[k]*(bh_z_1[:,:,k]-1.0)/(sigh_z_PML_1[k]+kappah_z_PML_1[k]*alphah_z_PML_1[k])/ kappah_z_PML_1[k]

    for k in range(nzPML_2):
        sige_z_PML_2[k] = sig_z_max * ( (nzPML_2 - (k+1)) / (nzPML_2 - 1.0) )**m
        alphae_z_PML_2[k] = alpha_z_max*(((k+1)-1)/(nzPML_2-1.0))**ma
        kappae_z_PML_2[k] = 1.0+(kappa_z_max-1.0)*((nzPML_2 - (k+1)) / (nzPML_2 - 1.0))**m
        be_z_2[:,:,k] = math.exp(-(sige_z_PML_2[k] / kappae_z_PML_2[k] +alphae_z_PML_2[k])*del_t/ep_0)
        if ((sige_z_PML_2[k] == 0.0) and (alphae_z_PML_2[k] == 0.0) and(k == nzPML_2-1)):
            ce_z_2[:,:,k] = 0.0
        else:
            ce_z_2[:,:,k] = sige_z_PML_2[k]*(be_z_2[:,:,k]-1.0)/(sige_z_PML_2[k]+kappae_z_PML_2[k]*alphae_z_PML_2[k])/ kappae_z_PML_2[k]

    for k in range(nzPML_2-1):
        sigh_z_PML_2[k] = sig_z_max * ( (nzPML_2 - (k+1) - 0.5)/(nzPML_2-1.0))**m
        alphah_z_PML_2[k] = alpha_z_max*(((k+1)-0.5)/(nzPML_2-1.0))**ma
        kappah_z_PML_2[k] = 1.0+(kappa_z_max-1.0)*((nzPML_2 - (k+1) - 0.5) / (nzPML_2 - 1.0))**m
        bh_z_2[:,:,k] = math.exp(-(sigh_z_PML_2[k] / kappah_z_PML_2[k] + alphah_z_PML_2[k])*del_t/ep_0)
        ch_z_2[:,:,k] = sigh_z_PML_2[k]*(bh_z_2[:,:,k]-1.0)/(sigh_z_PML_2[k]+kappah_z_PML_2[k]*alphah_z_PML_2[k]) / kappah_z_PML_2[k]

        # Indexing for main FDTD loop used for pml updates
        i_indices_h = th.arange(x_size - nxPML_2, x_size - 1, device=device)
        ii_indices_h = th.arange(nxPML_2 - 2, -1, -1, device=device)
        i_plus_1_indices_h = i_indices_h + 1

        j_indices_h = th.arange(y_size - nyPML_2, y_size - 1, device=device) 
        jj_indices_h = th.arange(nyPML_2 - 2, -1, -1, device=device) 
        j_plus_1_indices_h = j_indices_h + 1

        k_indices_h = th.arange(z_size - nzPML_2, z_size - 1, device=device)
        kk_indices_h = th.arange(nzPML_2 - 2, -1, -1, device=device)
        k_plus_1_indices_h = k_indices_h + 1

        i_indices_e = th.arange(x_size - nxPML_2, x_size - 1, device=device)
        ii_indices_e = th.arange(nxPML_2 - 1, 0, -1, device=device)
        i_minus_1_indices_e = i_indices_e - 1
        
        j_indices_e = th.arange(y_size - nyPML_2, y_size - 1, device=device)
        jj_indices_e = th.arange(nyPML_2 - 1, 0, -1, device=device)
        j_minus_1_indices_e = j_indices_e - 1

        k_indices_e = th.arange(z_size - nzPML_2, z_size - 1, device=device)
        kk_indices_e = th.arange(nzPML_2 - 1, 0, -1, device=device)
        k_minus_1_indices_e = k_indices_e - 1

##########################################################
# Denominator SETUP
##########################################################
with GPUProfiler("Denominator Calculations", device):

    den_hx[0 : nxPML_1 - 1] = 1.0 / (kappah_x_PML_1[:nxPML_1 - 1] * del_x)
    den_hx[x_size - nxPML_2 : x_size - 1] = 1.0 / (kappah_x_PML_2[:nxPML_2 - 1].flip(dims=[0]) * del_x)

    den_hy[0 : nyPML_1 - 1] = 1.0 / (kappah_y_PML_1[:nyPML_1 - 1] * del_y)
    den_hy[y_size - nyPML_2 : y_size - 1] = 1.0 / (kappah_y_PML_2[:nyPML_2 - 1].flip(dims=[0]) * del_y)

    den_hz[0 : nzPML_1 - 1] = 1.0 / (kappah_z_PML_1[:nzPML_1 - 1] * del_z)
    den_hz[z_size - nzPML_2 : z_size - 1] = 1.0 / (kappah_z_PML_2[:nzPML_2 - 1].flip(dims=[0]) * del_z)

    den_ex[0 : nxPML_1] = 1.0 / (kappae_x_PML_1[:nxPML_1] * del_x)
    den_ex[x_size - nxPML_2 : x_size - 1] = 1.0 / (kappae_x_PML_2[1:nxPML_2].flip(dims=[0]) * del_x)

    den_ey[0 : nyPML_1] = 1.0 / (kappae_y_PML_1[:nyPML_1] * del_y)
    den_ey[y_size - nyPML_2 : y_size - 1] = 1.0 / (kappae_y_PML_2[1:nyPML_2].flip(dims=[0]) * del_y)

    den_ez[0 : nzPML_1] = 1.0 / (kappae_z_PML_1[:nzPML_1] * del_z)
    den_ez[z_size - nzPML_2 : z_size - 1] = 1.0 / (kappae_z_PML_2[1:nzPML_2].flip(dims=[0]) * del_z)

##########################################################
# Geometry SETUP
##########################################################
with GPUProfiler("Geometry Calculations", device):
    """
    This section is long and messy
    Key items we keep when done with this section are:
    ga,gb variables (x,y,z)
    masks for main fdtd loop for sheets
    ep,sig matrices for sheets (xx,yy,zz) for special subcell update in E
    """

    # First we need to put all blocks,spheres,cylinders, and sheets into material cells
    ep_sig_tensors = [
        relative_ep_x_cell, relative_ep_y_cell, relative_ep_z_cell,
        sigma_x_cell, sigma_y_cell, sigma_z_cell
    ]
    gi = th.arange(relative_ep_x_cell.shape[0], device=device).view(-1, 1, 1)
    gj = th.arange(relative_ep_x_cell.shape[1], device=device).view(1, -1, 1)
    gk = th.arange(relative_ep_x_cell.shape[2], device=device).view(1, 1, -1)
    ii, jj, kk = 0, 0, 0
    for counter in range(num_objects):
        obj_type = object_type[counter].strip()
        
        # --- BLOCK ---
        if obj_type == 'block':
            m_id = (materials_properties[:, 0, 0] == blocks[ii, 0, 0]).nonzero(as_tuple=True)[0][0]
            
            # Define ranges
            i_s, i_e = int(blocks[ii, 1, 0]), int(blocks[ii, 1, 0] + blocks[ii, 2, 0])
            j_s, j_e = int(blocks[ii, 1, 1] + pbc_shift), int(blocks[ii, 1, 1] + blocks[ii, 2, 1] + pbc_shift)
            k_s, k_e = int(blocks[ii, 1, 2] + pbc_shift), int(blocks[ii, 1, 2] + blocks[ii, 2, 2] + pbc_shift)
            
            # Property values from materials_properties: ep(x,y,z) is row 1, sig(x,y,z) is row 2
            vals = th.cat([materials_properties[m_id, 1, :], materials_properties[m_id, 2, :]])
            
            # Apply slices to all 6 tensors
            for idx, tensor in enumerate(ep_sig_tensors):
                tensor[i_s:i_e, j_s:j_e, k_s:k_e] = vals[idx]
            
            ii += 1

        # --- SPHERE ---
        elif obj_type == 'sphere':
            m_id = (materials_properties[:, 0, 0] == spheres[jj, 0, 0]).nonzero(as_tuple=True)[0][0]
            rad = int(spheres[jj, 2, 0])
            ci, cj, ck = int(spheres[jj, 1, 0]), int(spheres[jj, 1, 1]), int(spheres[jj, 1, 2])
            
            # Bounding box slices
            s_i = slice(ci - rad, ci + rad)
            s_j = slice(cj - rad + pbc_shift, cj + rad + pbc_shift)
            s_k = slice(ck - rad + pbc_shift, ck + rad + pbc_shift)
            
            # Geometric mask within bounding box
            dist_sq = (gi[s_i] - ci)**2 + (gj[:, s_j, :] - cj)**2 + (gk[:, :, s_k] - ck)**2
            mask = dist_sq <= rad**2
            
            vals = th.cat([materials_properties[m_id, 1, :], materials_properties[m_id, 2, :]])
            
            # Apply mask to all 6 tensors
            for idx, tensor in enumerate(ep_sig_tensors):
                tensor[s_i, s_j, s_k][mask] = vals[idx]
                
            jj += 1

        elif obj_type == 'cylinder':
            m_id = (materials_properties[:, 0, 0] == cylinders[kk, 0, 0]).nonzero(as_tuple=True)[0][0]
            axis, length, rad = int(cylinders[kk, 1, 0]), int(cylinders[kk, 3, 0]), int(cylinders[kk, 3, 1])
            bi, bj, bk = cylinders[kk, 2, 0], cylinders[kk, 2, 1], cylinders[kk, 2, 2]

            D, H, W = relative_ep_x_cell.shape

            # 1. Calculate integer bounds with clamping
            if axis == 0: # X-Axis
                s_i = slice(max(0, int(bi)), min(D, int(bi + length)))
                s_j = slice(max(0, int(bj - rad + pbc_shift)), min(H, int(bj + rad + pbc_shift)))
                s_k = slice(max(0, int(bk - rad + pbc_shift)), min(W, int(bk + rad + pbc_shift)))
                # Cross-section is in J and K
                mj = gj[:, s_j, :].expand(-1, -1, -1) # Ensure 3D
                mk = gk[:, :, s_k].expand(-1, -1, -1)
                mask = (mj - bj)**2 + (mk - bk)**2 <= rad**2
            elif axis == 1: # Y-Axis
                s_i = slice(max(0, int(bi - rad)), min(D, int(bi + rad)))
                s_j = slice(max(0, int(bj + pbc_shift)), min(H, int(bj + length + pbc_shift)))
                s_k = slice(max(0, int(bk - rad + pbc_shift)), min(W, int(bk + rad + pbc_shift)))
                # Cross-section is in I and K
                mi = gi[s_i].expand(-1, -1, -1)
                mk = gk[:, :, s_k].expand(-1, -1, -1)
                mask = (mi - bi)**2 + (mk - bk)**2 <= rad**2
            else: # Z-Axis
                s_i = slice(max(0, int(bi - rad)), min(D, int(bi + rad)))
                s_j = slice(max(0, int(bj - rad + pbc_shift)), min(H, int(bj + rad + pbc_shift)))
                s_k = slice(max(0, int(bk + pbc_shift)), min(W, int(bk + length + pbc_shift)))
                # Cross-section is in I and J
                mi = gi[s_i].expand(-1, -1, -1)
                mj = gj[:, s_j, :].expand(-1, -1, -1)
                mask = (mi - bi)**2 + (mj - bj)**2 <= rad**2

            # 2. Prepare properties
            vals = th.cat([materials_properties[m_id, 1, :], materials_properties[m_id, 2, :]])
            
            # 3. Apply to tensors
            for idx, tensor in enumerate(ep_sig_tensors):
                sub_volume = tensor[s_i, s_j, s_k]
                
                # CRITICAL FIX: Explicitly expand the mask to match the sub_volume shape
                # This handles cases where one dimension of the mask is 1 (like your 1x10x10 case)
                current_mask = mask.expand(sub_volume.shape)
                
                sub_volume[current_mask] = vals[idx]
                
            kk += 1
  
    # Create a helper to find material IDs for any sheet tensor
    def get_material_ids(sheets, properties):
        # Broadcast comparison: (num_sheets, 1) == (1, num_materials)
        matches = sheets[:, 0:1, 0] == properties[:, 0, 0]
        return matches.float().argmax(dim=1)

    # Pre-calculate for all three orientations
    if num_sheets_x>0:
        m_ids_x = get_material_ids(sheets_x, sheet_properties)
    if num_sheets_y>0:
        m_ids_y = get_material_ids(sheets_y, sheet_properties)
    if num_sheets_z>0:
        m_ids_z = get_material_ids(sheets_z, sheet_properties)

    # --- X-Oriented Sheets ---
    for ii in range(num_sheets_x):
        i = int(sheets_x[ii, 1, 0])
        j_s, j_e = int(sheets_x[ii, 2, 0] + pbc_shift), int(sheets_x[ii, 2, 0] + sheets_x[ii, 3, 0] + pbc_shift)
        k_s, k_e = int(sheets_x[ii, 2, 1] + pbc_shift), int(sheets_x[ii, 2, 1] + sheets_x[ii, 3, 1] + pbc_shift)
        
        mid = m_ids_x[ii]
        # Slicing the whole material property block at once
        sheet_ep_x_cell_x[i, j_s:j_e, k_s:k_e] = sheet_properties[mid, 2, 0]
        sheet_ep_y_cell_x[i, j_s:j_e, k_s:k_e] = sheet_properties[mid, 2, 1]
        sheet_ep_z_cell_x[i, j_s:j_e, k_s:k_e] = sheet_properties[mid, 2, 2]
        sheet_sig_x_cell_x[i, j_s:j_e, k_s:k_e] = sheet_properties[mid, 3, 0]
        sheet_sig_y_cell_x[i, j_s:j_e, k_s:k_e] = sheet_properties[mid, 3, 1]
        sheet_sig_z_cell_x[i, j_s:j_e, k_s:k_e] = sheet_properties[mid, 3, 2]

    # --- Y-Oriented Sheets ---
    for ii in range(num_sheets_y):
        j = int(sheets_y[ii, 1, 0] + pbc_shift)
        i_s, i_e = int(sheets_y[ii, 2, 0]), int(sheets_y[ii, 2, 0] + sheets_y[ii, 3, 0])
        k_s, k_e = int(sheets_y[ii, 2, 1] + pbc_shift), int(sheets_y[ii, 2, 1] + sheets_y[ii, 3, 1] + pbc_shift)
        
        mid = m_ids_y[ii]
        sheet_ep_x_cell_y[i_s:i_e, j, k_s:k_e] = sheet_properties[mid, 2, 0]
        sheet_ep_y_cell_y[i_s:i_e, j, k_s:k_e] = sheet_properties[mid, 2, 1]
        sheet_ep_z_cell_y[i_s:i_e, j, k_s:k_e] = sheet_properties[mid, 2, 2]
        sheet_sig_x_cell_y[i_s:i_e, j, k_s:k_e] = sheet_properties[mid, 3, 0]
        sheet_sig_y_cell_y[i_s:i_e, j, k_s:k_e] = sheet_properties[mid, 3, 1]
        sheet_sig_z_cell_y[i_s:i_e, j, k_s:k_e] = sheet_properties[mid, 3, 2]

    # --- Z-Oriented Sheets ---
    for ii in range(num_sheets_z):
        k = int(sheets_z[ii, 1, 0] + pbc_shift)
        i_s, i_e = int(sheets_z[ii, 2, 0]), int(sheets_z[ii, 2, 0] + sheets_z[ii, 3, 0])
        j_s, j_e = int(sheets_z[ii, 2, 1] + pbc_shift), int(sheets_z[ii, 2, 1] + sheets_z[ii, 3, 1] + pbc_shift)
        
        mid = m_ids_z[ii]
        sheet_ep_x_cell_z[i_s:i_e, j_s:j_e, k] = sheet_properties[mid, 2, 0]
        sheet_ep_y_cell_z[i_s:i_e, j_s:j_e, k] = sheet_properties[mid, 2, 1]
        sheet_ep_z_cell_z[i_s:i_e, j_s:j_e, k] = sheet_properties[mid, 2, 2]
        sheet_sig_x_cell_z[i_s:i_e, j_s:j_e, k] = sheet_properties[mid, 3, 0]
        sheet_sig_y_cell_z[i_s:i_e, j_s:j_e, k] = sheet_properties[mid, 3, 1]
        sheet_sig_z_cell_z[i_s:i_e, j_s:j_e, k] = sheet_properties[mid, 3, 2]

    # now if pbc we need to account for this in material cells before moving to yee cell creation
    if type_sim == 0:
        # z-direction PBC for i,j
        relative_ep_x_cell[:,:,0] = relative_ep_x_cell[:,:,z_size-2]
        relative_ep_y_cell[:,:,0] = relative_ep_y_cell[:,:,z_size-2]
        relative_ep_z_cell[:,:,0] = relative_ep_z_cell[:,:,z_size-2]
        sigma_x_cell[:,:,0] = sigma_x_cell[:,:,z_size-2]
        sigma_y_cell[:,:,0] = sigma_y_cell[:,:,z_size-2]
        sigma_z_cell[:,:,0] = sigma_z_cell[:,:,z_size-2]

        sheet_ep_x_cell_x[:,:,0] = sheet_ep_x_cell_x[:,:,z_size-2]
        sheet_ep_y_cell_x[:,:,0] = sheet_ep_y_cell_x[:,:,z_size-2]
        sheet_ep_z_cell_x[:,:,0] = sheet_ep_z_cell_x[:,:,z_size-2]
        sheet_sig_x_cell_x[:,:,0] = sheet_sig_x_cell_x[:,:,z_size-2]
        sheet_sig_y_cell_x[:,:,0] = sheet_sig_y_cell_x[:,:,z_size-2]
        sheet_sig_z_cell_x[:,:,0] = sheet_sig_z_cell_x[:,:,z_size-2]

        sheet_ep_x_cell_y[:,:,0] = sheet_ep_x_cell_y[:,:,z_size-2]
        sheet_ep_y_cell_y[:,:,0] = sheet_ep_y_cell_y[:,:,z_size-2]
        sheet_ep_z_cell_y[:,:,0] = sheet_ep_z_cell_y[:,:,z_size-2]
        sheet_sig_x_cell_y[:,:,0] = sheet_sig_x_cell_y[:,:,z_size-2]
        sheet_sig_y_cell_y[:,:,0] = sheet_sig_y_cell_y[:,:,z_size-2]
        sheet_sig_z_cell_y[:,:,0] = sheet_sig_z_cell_y[:,:,z_size-2]

        sheet_ep_x_cell_z[:,:,0] = sheet_ep_x_cell_z[:,:,z_size-2]
        sheet_ep_y_cell_z[:,:,0] = sheet_ep_y_cell_z[:,:,z_size-2]
        sheet_ep_z_cell_z[:,:,0] = sheet_ep_z_cell_z[:,:,z_size-2]
        sheet_sig_x_cell_z[:,:,0] = sheet_sig_x_cell_z[:,:,z_size-2]
        sheet_sig_y_cell_z[:,:,0] = sheet_sig_y_cell_z[:,:,z_size-2]
        sheet_sig_z_cell_z[:,:,0] = sheet_sig_z_cell_z[:,:,z_size-2]

        # y-direction PBC for i,k
        relative_ep_x_cell[:,0,:] = relative_ep_x_cell[:,y_size-2,:]
        relative_ep_y_cell[:,0,:] = relative_ep_y_cell[:,y_size-2,:]
        relative_ep_z_cell[:,0,:] = relative_ep_z_cell[:,y_size-2,:]
        sigma_x_cell[:,0,:] = sigma_x_cell[:,y_size-2,:]
        sigma_y_cell[:,0,:] = sigma_y_cell[:,y_size-2,:]
        sigma_z_cell[:,0,:] = sigma_z_cell[:,y_size-2,:]

        sheet_ep_x_cell_x[:,0,:] = sheet_ep_x_cell_x[:,y_size-2,:]
        sheet_ep_y_cell_x[:,0,:] = sheet_ep_y_cell_x[:,y_size-2,:]
        sheet_ep_z_cell_x[:,0,:] = sheet_ep_z_cell_x[:,y_size-2,:]
        sheet_sig_x_cell_x[:,0,:] = sheet_sig_x_cell_x[:,y_size-2,:]
        sheet_sig_y_cell_x[:,0,:] = sheet_sig_y_cell_x[:,y_size-2,:]
        sheet_sig_z_cell_x[:,0,:] = sheet_sig_z_cell_x[:,y_size-2,:]

        sheet_ep_x_cell_y[:,0,:] = sheet_ep_x_cell_y[:,y_size-2,:]
        sheet_ep_y_cell_y[:,0,:] = sheet_ep_y_cell_y[:,y_size-2,:]
        sheet_ep_z_cell_y[:,0,:] = sheet_ep_z_cell_y[:,y_size-2,:]
        sheet_sig_x_cell_y[:,0,:] = sheet_sig_x_cell_y[:,y_size-2,:]
        sheet_sig_y_cell_y[:,0,:] = sheet_sig_y_cell_y[:,y_size-2,:]
        sheet_sig_z_cell_y[:,0,:] = sheet_sig_z_cell_y[:,y_size-2,:]

        sheet_ep_x_cell_z[:,0,:] = sheet_ep_x_cell_z[:,y_size-2,:]
        sheet_ep_y_cell_z[:,0,:] = sheet_ep_y_cell_z[:,y_size-2,:]
        sheet_ep_z_cell_z[:,0,:] = sheet_ep_z_cell_z[:,y_size-2,:]
        sheet_sig_x_cell_z[:,0,:] = sheet_sig_x_cell_z[:,y_size-2,:]
        sheet_sig_y_cell_z[:,0,:] = sheet_sig_y_cell_z[:,y_size-2,:]
        sheet_sig_z_cell_z[:,0,:] = sheet_sig_z_cell_z[:,y_size-2,:]

    # now material cells to yee cell creation
    # Bulk materials
    relative_ep_x[1:-1,1:-1,1:-1] = (
        relative_ep_x_cell[1:-1,1:-1,1:-1] + 
        relative_ep_x_cell[1:-1,0:-2,1:-1] + 
        relative_ep_x_cell[1:-1,1:-1,0:-2] + 
        relative_ep_x_cell[1:-1,0:-2,0:-2]
    )/4.0

    relative_ep_y[1:-1,1:-1,1:-1] = (
        relative_ep_y_cell[1:-1,1:-1,1:-1] +
        relative_ep_y_cell[0:-2,1:-1,1:-1] +
        relative_ep_y_cell[1:-1,1:-1,0:-2] +
        relative_ep_y_cell[0:-2,1:-1,0:-2]
    )/4.0

    relative_ep_z[1:-1,1:-1,1:-1] = (
        relative_ep_z_cell[1:-1,1:-1,1:-1] +
        relative_ep_z_cell[0:-2,1:-1,1:-1] +
        relative_ep_z_cell[1:-1,0:-2,1:-1] +
        relative_ep_z_cell[0:-2,0:-2,1:-1]
    )/4.0

    sigma_x[1:-1,1:-1,1:-1] = (
        sigma_x_cell[1:-1,1:-1,1:-1] +
        sigma_x_cell[1:-1,0:-2,1:-1] +
        sigma_x_cell[1:-1,1:-1,0:-2] +
        sigma_x_cell[1:-1,0:-2,0:-2]
    )/4.0

    sigma_y[1:-1,1:-1,1:-1] = (
        sigma_y_cell[1:-1,1:-1,1:-1] +
        sigma_y_cell[0:-2,1:-1,1:-1] +
        sigma_y_cell[1:-1,1:-1,0:-2] +
        sigma_y_cell[0:-2,1:-1,0:-2]
    )/4.0

    sigma_z[1:-1,1:-1,1:-1] = (
        sigma_z_cell[1:-1,1:-1,1:-1] +
        sigma_z_cell[0:-2,1:-1,1:-1] +
        sigma_z_cell[1:-1,0:-2,1:-1] +
        sigma_z_cell[0:-2,0:-2,1:-1]
    )/4.0

    # These next 3 if statements just save time for material to yee if nothing it's 1 and 0s anyway in all places
    # Now put sheets into yee cells from material cells
    if num_sheets_x>0:
        # x-normal sheets
        sheet_ep_x_x[1:-1, 1:-1, 1:-1] = (
            sheet_ep_x_cell_x[1:-1, 1:-1, 1:-1] +
            sheet_ep_x_cell_x[1:-1, 0:-2, 1:-1] +
            sheet_ep_x_cell_x[1:-1, 1:-1, 0:-2] +
            sheet_ep_x_cell_x[1:-1, 0:-2, 0:-2]
        ) / 4.0

        sheet_ep_y_x[1:-1, 1:-1, 1:-1] = (
            sheet_ep_y_cell_x[1:-1, 1:-1, 1:-1] +
            sheet_ep_y_cell_x[1:-1, 1:-1, 0:-2] +
            sheet_ep_y_cell_x[1:-1, 1:-1, 1:-1] +
            sheet_ep_y_cell_x[1:-1, 1:-1, 0:-2]
        ) / 4.0

        sheet_ep_z_x[1:-1, 1:-1, 1:-1] = (
            sheet_ep_z_cell_x[1:-1, 1:-1, 1:-1] +
            sheet_ep_z_cell_x[1:-1, 0:-2, 1:-1] +
            sheet_ep_z_cell_x[1:-1, 1:-1, 1:-1] +
            sheet_ep_z_cell_x[1:-1, 0:-2, 1:-1]
        ) / 4.0

        sheet_sig_x_x[1:-1, 1:-1, 1:-1] = (
            sheet_sig_x_cell_x[1:-1, 1:-1, 1:-1] +
            sheet_sig_x_cell_x[1:-1, 0:-2, 1:-1] +
            sheet_sig_x_cell_x[1:-1, 1:-1, 0:-2] +
            sheet_sig_x_cell_x[1:-1, 0:-2, 0:-2]
        ) / 4.0

        sheet_sig_y_x[1:-1, 1:-1, 1:-1] = (
            sheet_sig_y_cell_x[1:-1, 1:-1, 1:-1] +
            sheet_sig_y_cell_x[1:-1, 1:-1, 0:-2] +
            sheet_sig_y_cell_x[1:-1, 1:-1, 1:-1] +
            sheet_sig_y_cell_x[1:-1, 1:-1, 0:-2]
        ) / 4.0

        sheet_sig_z_x[1:-1, 1:-1, 1:-1] = (
            sheet_sig_z_cell_x[1:-1, 1:-1, 1:-1] +
            sheet_sig_z_cell_x[1:-1, 0:-2, 1:-1] +
            sheet_sig_z_cell_x[1:-1, 1:-1, 1:-1] +
            sheet_sig_z_cell_x[1:-1, 0:-2, 1:-1]
        ) / 4.0

    if num_sheets_y>0:
        # y-normal sheets
        sheet_ep_x_y[1:-1, 1:-1, 1:-1] = (
            sheet_ep_x_cell_y[1:-1, 1:-1, 1:-1] +
            sheet_ep_x_cell_y[1:-1, 1:-1, 0:-2] +
            sheet_ep_x_cell_y[1:-1, 1:-1, 1:-1] +
            sheet_ep_x_cell_y[1:-1, 1:-1, 0:-2]
        ) / 4.0

        sheet_ep_y_y[1:-1, 1:-1, 1:-1] = (
            sheet_ep_y_cell_y[1:-1, 1:-1, 1:-1] +
            sheet_ep_y_cell_y[1:-1, 1:-1, 0:-2] +
            sheet_ep_y_cell_y[0:-2, 1:-1, 1:-1] +
            sheet_ep_y_cell_y[0:-2, 1:-1, 0:-2]
        ) / 4.0

        sheet_ep_z_y[1:-1, 1:-1, 1:-1] = (
            sheet_ep_z_cell_y[1:-1, 1:-1, 1:-1] +
            sheet_ep_z_cell_y[0:-2, 1:-1, 1:-1] +
            sheet_ep_z_cell_y[1:-1, 1:-1, 1:-1] +
            sheet_ep_z_cell_y[0:-2, 1:-1, 1:-1]
        ) / 4.0

        sheet_sig_x_y[1:-1, 1:-1, 1:-1] = (
            sheet_sig_x_cell_y[1:-1, 1:-1, 1:-1] +
            sheet_sig_x_cell_y[1:-1, 1:-1, 0:-2] +
            sheet_sig_x_cell_y[1:-1, 1:-1, 1:-1] +
            sheet_sig_x_cell_y[1:-1, 1:-1, 0:-2]
        ) / 4.0

        sheet_sig_y_y[1:-1, 1:-1, 1:-1] = (
            sheet_sig_y_cell_y[1:-1, 1:-1, 1:-1] +
            sheet_sig_y_cell_y[1:-1, 1:-1, 0:-2] +
            sheet_sig_y_cell_y[0:-2, 1:-1, 1:-1] +
            sheet_sig_y_cell_y[0:-2, 1:-1, 0:-2]
        ) / 4.0

        sheet_sig_z_y[1:-1, 1:-1, 1:-1] = (
            sheet_sig_z_cell_y[1:-1, 1:-1, 1:-1] +
            sheet_sig_z_cell_y[0:-2, 1:-1, 1:-1] +
            sheet_sig_z_cell_y[1:-1, 1:-1, 1:-1] +
            sheet_sig_z_cell_y[0:-2, 1:-1, 1:-1]
        ) / 4.0

    if num_sheets_z>0:
        # z-normal sheets
        sheet_ep_x_z[1:-1, 1:-1, 1:-1] = (
            sheet_ep_x_cell_z[1:-1, 1:-1, 1:-1] +
            sheet_ep_x_cell_z[1:-1, 0:-2, 1:-1] +
            sheet_ep_x_cell_z[1:-1, 1:-1, 1:-1] +
            sheet_ep_x_cell_z[1:-1, 0:-2, 1:-1]
        ) / 4.0

        sheet_ep_y_z[1:-1, 1:-1, 1:-1] = (
            sheet_ep_y_cell_z[1:-1, 1:-1, 1:-1] +
            sheet_ep_y_cell_z[0:-2, 1:-1, 1:-1] +
            sheet_ep_y_cell_z[1:-1, 1:-1, 1:-1] +
            sheet_ep_y_cell_z[0:-2, 1:-1, 1:-1]
        ) / 4.0

        sheet_ep_z_z[1:-1, 1:-1, 1:-1] = (
            sheet_ep_z_cell_z[1:-1, 1:-1, 1:-1] +
            sheet_ep_z_cell_z[0:-2, 1:-1, 1:-1] +
            sheet_ep_z_cell_z[1:-1, 0:-2, 1:-1] +
            sheet_ep_z_cell_z[0:-2, 0:-2, 1:-1]
        ) / 4.0

        sheet_sig_x_z[1:-1, 1:-1, 1:-1] = (
            sheet_sig_x_cell_z[1:-1, 1:-1, 1:-1] +
            sheet_sig_x_cell_z[1:-1, 0:-2, 1:-1] +
            sheet_sig_x_cell_z[1:-1, 1:-1, 1:-1] +
            sheet_sig_x_cell_z[1:-1, 0:-2, 1:-1]
        ) / 4.0

        sheet_sig_y_z[1:-1, 1:-1, 1:-1] = (
            sheet_sig_y_cell_z[1:-1, 1:-1, 1:-1] +
            sheet_sig_y_cell_z[0:-2, 1:-1, 1:-1] +
            sheet_sig_y_cell_z[1:-1, 1:-1, 1:-1] +
            sheet_sig_y_cell_z[0:-2, 1:-1, 1:-1]
        ) / 4.0

        sheet_sig_z_z[1:-1, 1:-1, 1:-1] = (
            sheet_sig_z_cell_z[1:-1, 1:-1, 1:-1] +
            sheet_sig_z_cell_z[0:-2, 1:-1, 1:-1] +
            sheet_sig_z_cell_z[1:-1, 0:-2, 1:-1] +
            sheet_sig_z_cell_z[0:-2, 0:-2, 1:-1]
        ) / 4.0

    # now if pbc exists we need to apply the same condition like we did for material cells
    if type_sim == 0:
        # z-boundary
        relative_ep_x[:, :, 0] = relative_ep_x[:, :, z_size-2]
        relative_ep_y[:, :, 0] = relative_ep_y[:, :, z_size-2]
        relative_ep_z[:, :, 0] = relative_ep_z[:, :, z_size-2]
        sigma_x[:, :, 0] = sigma_x[:, :, z_size-2]
        sigma_y[:, :, 0] = sigma_y[:, :, z_size-2]
        sigma_z[:, :, 0] = sigma_z[:, :, z_size-2]

        # y-boundary
        relative_ep_x[:, 0, :] = relative_ep_x[:, y_size-2, :]
        relative_ep_y[:, 0, :] = relative_ep_y[:, y_size-2, :]
        relative_ep_z[:, 0, :] = relative_ep_z[:, y_size-2, :]
        sigma_x[:, 0, :] = sigma_x[:, y_size-2, :]
        sigma_y[:, 0, :] = sigma_y[:, y_size-2, :]
        sigma_z[:, 0, :] = sigma_z[:, y_size-2, :]

        # similar to above, these 3 if statements are just to save time - 1 and 0s if nothing added so the vectors wouldn't do anything
        if num_sheets_x>0:
            # z-boundary
            sheet_ep_x_x[:, :, 0] = sheet_ep_x_x[:, :, z_size-2]
            sheet_ep_y_x[:, :, 0] = sheet_ep_y_x[:, :, z_size-2]
            sheet_ep_z_x[:, :, 0] = sheet_ep_z_x[:, :, z_size-2]
            sheet_sig_x_x[:, :, 0] = sheet_sig_x_x[:, :, z_size-2]
            sheet_sig_y_x[:, :, 0] = sheet_sig_y_x[:, :, z_size-2]
            sheet_sig_z_x[:, :, 0] = sheet_sig_z_x[:, :, z_size-2]

            # y-boundary
            sheet_ep_x_x[:, 0, :] = sheet_ep_x_x[:, y_size-2, :]
            sheet_ep_y_x[:, 0, :] = sheet_ep_y_x[:, y_size-2, :]
            sheet_ep_z_x[:, 0, :] = sheet_ep_z_x[:, y_size-2, :]
            sheet_sig_x_x[:, 0, :] = sheet_sig_x_x[:, y_size-2, :]
            sheet_sig_y_x[:, 0, :] = sheet_sig_y_x[:, y_size-2, :]
            sheet_sig_z_x[:, 0, :] = sheet_sig_z_x[:, y_size-2, :]

        if num_sheets_y>0:
            # z-boundary
            sheet_ep_x_y[:, :, 0] = sheet_ep_x_y[:, :, z_size-2]
            sheet_ep_y_y[:, :, 0] = sheet_ep_y_y[:, :, z_size-2]
            sheet_ep_z_y[:, :, 0] = sheet_ep_z_y[:, :, z_size-2]
            sheet_sig_x_y[:, :, 0] = sheet_sig_x_y[:, :, z_size-2]
            sheet_sig_y_y[:, :, 0] = sheet_sig_y_y[:, :, z_size-2]
            sheet_sig_z_y[:, :, 0] = sheet_sig_z_y[:, :, z_size-2]

            # y-boundary
            sheet_ep_x_y[:, 0, :] = sheet_ep_x_y[:, y_size-2, :]
            sheet_ep_y_y[:, 0, :] = sheet_ep_y_y[:, y_size-2, :]
            sheet_ep_z_y[:, 0, :] = sheet_ep_z_y[:, y_size-2, :]
            sheet_sig_x_y[:, 0, :] = sheet_sig_x_y[:, y_size-2, :]
            sheet_sig_y_y[:, 0, :] = sheet_sig_y_y[:, y_size-2, :]
            sheet_sig_z_y[:, 0, :] = sheet_sig_z_y[:, y_size-2, :]

        if num_sheets_z>0:
            # z-boundary
            sheet_ep_x_z[:, :, 0] = sheet_ep_x_z[:, :, z_size-2]
            sheet_ep_y_z[:, :, 0] = sheet_ep_y_z[:, :, z_size-2]
            sheet_ep_z_z[:, :, 0] = sheet_ep_z_z[:, :, z_size-2]
            sheet_sig_x_z[:, :, 0] = sheet_sig_x_z[:, :, z_size-2]
            sheet_sig_y_z[:, :, 0] = sheet_sig_y_z[:, :, z_size-2]
            sheet_sig_z_z[:, :, 0] = sheet_sig_z_z[:, :, z_size-2]

            # y-boundary
            sheet_ep_x_z[:, 0, :] = sheet_ep_x_z[:, y_size-2, :]
            sheet_ep_y_z[:, 0, :] = sheet_ep_y_z[:, y_size-2, :]
            sheet_ep_z_z[:, 0, :] = sheet_ep_z_z[:, y_size-2, :]
            sheet_sig_x_z[:, 0, :] = sheet_sig_x_z[:, y_size-2, :]
            sheet_sig_y_z[:, 0, :] = sheet_sig_y_z[:, y_size-2, :]
            sheet_sig_z_z[:, 0, :] = sheet_sig_z_z[:, y_size-2, :]

    # For permittivity and electrical conductivity
    # X-component
    denom_x = 1.0 + sigma_x[:-1, :-1, :-1] * del_t / (2.0 * ep_0 * relative_ep_x[:-1, :-1, :-1])
    gax[:-1, :-1, :-1] = (1.0 - sigma_x[:-1, :-1, :-1] * del_t / (2.0 * ep_0 * relative_ep_x[:-1, :-1, :-1])) / denom_x
    gbx[:-1, :-1, :-1] = (del_t / (ep_0 * relative_ep_x[:-1, :-1, :-1])) / denom_x
    
    # Y-component
    denom_y = 1.0 + sigma_y[:-1, :-1, :-1] * del_t / (2.0 * ep_0 * relative_ep_y[:-1, :-1, :-1])
    gay[:-1, :-1, :-1] = (1.0 - sigma_y[:-1, :-1, :-1] * del_t / (2.0 * ep_0 * relative_ep_y[:-1, :-1, :-1])) / denom_y
    gby[:-1, :-1, :-1] = (del_t / (ep_0 * relative_ep_y[:-1, :-1, :-1])) / denom_y
    
    # Z-component
    denom_z = 1.0 + sigma_z[:-1, :-1, :-1] * del_t / (2.0 * ep_0 * relative_ep_z[:-1, :-1, :-1])
    gaz[:-1, :-1, :-1] = (1.0 - sigma_z[:-1, :-1, :-1] * del_t / (2.0 * ep_0 * relative_ep_z[:-1, :-1, :-1])) / denom_z
    gbz[:-1, :-1, :-1] = (del_t / (ep_0 * relative_ep_z[:-1, :-1, :-1])) / denom_z

    # Now lastly, we need to deal with sheet adjustments to ga,gb
    # Initialize sheet averages with default material values
    # Counters: 1 if a sheet is present in that direction, 0 otherwise

    sheet_ep_avg_x = relative_ep_x.clone()
    sheet_ep_avg_y = relative_ep_y.clone()
    sheet_ep_avg_z = relative_ep_z.clone()
    sheet_sig_avg_x = sigma_x.clone()
    sheet_sig_avg_y = sigma_y.clone()
    sheet_sig_avg_z = sigma_z.clone()

    # Counters: 1 if a sheet is present in that direction, 0 otherwise
    counter_x = th.zeros_like(relative_ep_x, device=device, dtype=th.bool)
    counter_y = th.zeros_like(relative_ep_y, device=device, dtype=th.bool)
    counter_z = th.zeros_like(relative_ep_z, device=device, dtype=th.bool)

    # --- x-normal sheets ---
    mask_y = (sheet_sig_y_x > 0) | (sheet_ep_y_x > 1)
    mask_z = (sheet_sig_z_x > 0) | (sheet_ep_z_x > 1)

    sheet_ep_avg_y[mask_y] += (sheet_thickness / del_x) * (sheet_ep_y_x[mask_y] - relative_ep_y[mask_y])
    sheet_sig_avg_y[mask_y] += (sheet_thickness / del_x) * (sheet_sig_y_x[mask_y] - sigma_y[mask_y])
    counter_y[mask_y] = True

    sheet_ep_avg_z[mask_z] += (sheet_thickness / del_x) * (sheet_ep_z_x[mask_z] - relative_ep_z[mask_z])
    sheet_sig_avg_z[mask_z] += (sheet_thickness / del_x) * (sheet_sig_z_x[mask_z] - sigma_z[mask_z])
    counter_z[mask_z] = True

    # --- y-normal sheets ---
    mask_x = (sheet_sig_x_y > 0) | (sheet_ep_x_y > 1)
    mask_z_y = (sheet_sig_z_y > 0) | (sheet_ep_z_y > 1)

    sheet_ep_avg_x[mask_x] += (sheet_thickness / del_y) * (sheet_ep_x_y[mask_x] - relative_ep_x[mask_x])
    sheet_sig_avg_x[mask_x] += (sheet_thickness / del_y) * (sheet_sig_x_y[mask_x] - sigma_x[mask_x])
    counter_x[mask_x] = True

    sheet_ep_avg_z[mask_z_y] += (sheet_thickness / del_y) * (sheet_ep_z_y[mask_z_y] - relative_ep_z[mask_z_y])
    sheet_sig_avg_z[mask_z_y] += (sheet_thickness / del_y) * (sheet_sig_z_y[mask_z_y] - sigma_z[mask_z_y])
    counter_z[mask_z_y] = True

    # --- z-normal sheets ---
    mask_x_z = (sheet_sig_x_z > 0) | (sheet_ep_x_z > 1)
    mask_y_z = (sheet_sig_y_z > 0) | (sheet_ep_y_z > 1)

    sheet_ep_avg_x[mask_x_z] += (sheet_thickness / del_z) * (sheet_ep_x_z[mask_x_z] - relative_ep_x[mask_x_z])
    sheet_sig_avg_x[mask_x_z] += (sheet_thickness / del_z) * (sheet_sig_x_z[mask_x_z] - sigma_x[mask_x_z])
    counter_x[mask_x_z] = True

    sheet_ep_avg_y[mask_y_z] += (sheet_thickness / del_z) * (sheet_ep_y_z[mask_y_z] - relative_ep_y[mask_y_z])
    sheet_sig_avg_y[mask_y_z] += (sheet_thickness / del_z) * (sheet_sig_y_z[mask_y_z] - sigma_y[mask_y_z])
    counter_y[mask_y_z] = True

    # x-direction
    mask = counter_x
    gax[mask] = (1.0 - sheet_sig_avg_x[mask]*del_t/(2.0*ep_0*sheet_ep_avg_x[mask])) / \
                (1.0 + sheet_sig_avg_x[mask]*del_t/(2.0*ep_0*sheet_ep_avg_x[mask]))
    gbx[mask] = (del_t/(ep_0*sheet_ep_avg_x[mask])) / \
                (1.0 + sheet_sig_avg_x[mask]*del_t/(2.0*ep_0*sheet_ep_avg_x[mask]))

    # y-direction
    mask = counter_y
    gay[mask] = (1.0 - sheet_sig_avg_y[mask]*del_t/(2.0*ep_0*sheet_ep_avg_y[mask])) / \
                (1.0 + sheet_sig_avg_y[mask]*del_t/(2.0*ep_0*sheet_ep_avg_y[mask]))
    gby[mask] = (del_t/(ep_0*sheet_ep_avg_y[mask])) / \
                (1.0 + sheet_sig_avg_y[mask]*del_t/(2.0*ep_0*sheet_ep_avg_y[mask]))

    # z-direction
    mask = counter_z
    gaz[mask] = (1.0 - sheet_sig_avg_z[mask]*del_t/(2.0*ep_0*sheet_ep_avg_z[mask])) / \
                (1.0 + sheet_sig_avg_z[mask]*del_t/(2.0*ep_0*sheet_ep_avg_z[mask]))
    gbz[mask] = (del_t/(ep_0*sheet_ep_avg_z[mask])) / \
                (1.0 + sheet_sig_avg_z[mask]*del_t/(2.0*ep_0*sheet_ep_avg_z[mask]))

    # masks are used in main fdtd for it to know if there is a sheet or not
    mask_x = ((sheet_sig_x_x > 0) | (sheet_ep_x_x > 1)).float()  # x-normal sheets
    mask_y = ((sheet_sig_y_y > 0) | (sheet_ep_y_y > 1)).float()  # y-normal sheets
    mask_z = ((sheet_sig_z_z > 0) | (sheet_ep_z_z > 1)).float()  # z-normal sheets

##########################################################
# Definitions for main FDTD loop - for calling individual sections to make code more readable
##########################################################
def update_H_fields(Hx, Hy, Hz, Ex, Ey, Ez, da, db, den_hx, den_hy, den_hz):
    """
    Update H field components using curl of E fields.
    
    Inputs (read): Ex, Ey, Ez, da, db, den_hx, den_hy, den_hz
    Outputs (modified in-place): Hx, Hy, Hz
    """
    Hx[:-1, :-1, :-1] = da * Hx[:-1, :-1, :-1] + db * (
        (Ey[:-1, :-1, 1:] - Ey[:-1, :-1, :-1]) * den_hz[None,None,:] +
        (Ez[:-1, :-1, :-1] - Ez[:-1, 1:, :-1]) * den_hy[None,:,None]
    )
    Hy[:-1, :-1, :-1] = da * Hy[:-1, :-1, :-1] + db * (
        (Ex[:-1, :-1, :-1] - Ex[:-1, :-1, 1:]) * den_hz[None,None,:] +
        (Ez[1:, :-1, :-1] - Ez[:-1, :-1, :-1]) * den_hx[:,None,None]
    )
    Hz[:-1, :-1, :-1] = da * Hz[:-1, :-1, :-1] + db * (
        (Ey[:-1, :-1, :-1] - Ey[1:, :-1, :-1]) * den_hx[:,None,None] +
        (Ex[:-1, 1:, :-1] - Ex[:-1, :-1, :-1]) * den_hy[None,:,None]
    )


def update_H_pml(Hx, Hy, Hz, Ex, Ey, Ez, db, del_x, del_y, del_z,
                 psi_Hyx_1, psi_Hyx_2, psi_Hzx_1, psi_Hzx_2,
                 psi_Hxy_1, psi_Hxy_2, psi_Hzy_1, psi_Hzy_2,
                 psi_Hxz_1, psi_Hxz_2, psi_Hyz_1, psi_Hyz_2,
                 bh_x_1, ch_x_1, bh_x_2, ch_x_2,
                 bh_y_1, ch_y_1, bh_y_2, ch_y_2,
                 bh_z_1, ch_z_1, bh_z_2, ch_z_2,
                 nxPML_1, nyPML_1, nzPML_1,
                 ii_indices_h, i_indices_h, i_plus_1_indices_h,
                 jj_indices_h, j_indices_h, j_plus_1_indices_h,
                 kk_indices_h, k_indices_h, k_plus_1_indices_h):
    """
    Apply PML boundary conditions to H fields.
    
    Inputs (read): Ex, Ey, Ez, db, del_x/y/z, bh_*, ch_*, PML parameters, indices
    Outputs (modified in-place): Hx, Hy, Hz, psi_* arrays
    """
    # X-direction PMLs
    # Hy update
    psi_Hyx_1[:,:-1,:-1] = bh_x_1[:,:-1,:-1] * psi_Hyx_1[:,:-1,:-1] + ch_x_1[:,:-1,:-1] * (Ez[1:(nxPML_1),:-1,:-1] - Ez[:(nxPML_1-1),:-1,:-1]) / del_x
    Hy[:(nxPML_1-1),:-1,:-1] = Hy[:(nxPML_1-1),:-1,:-1] + db * psi_Hyx_1[:,:-1,:-1]
    
    psi_Hyx_2[ii_indices_h,:-1,:-1] = (
        bh_x_2[ii_indices_h,:-1,:-1] * psi_Hyx_2[ii_indices_h,:-1,:-1] 
        + ch_x_2[ii_indices_h,:-1,:-1] * (Ez[i_plus_1_indices_h,:-1,:-1] - Ez[i_indices_h,:-1,:-1]) / del_x
    )
    Hy[i_indices_h,:-1,:-1] = Hy[i_indices_h,:-1,:-1] + db * psi_Hyx_2[ii_indices_h,:-1,:-1]
    
    # Hz update
    psi_Hzx_1[:,:,1:-1] = bh_x_1[:,:-1,1:-1] * psi_Hzx_1[:,:,1:-1] + ch_x_1[:,:-1,1:-1] * (Ey[:(nxPML_1-1),:-1,1:-1] - Ey[1:nxPML_1,:-1,1:-1]) / del_x
    Hz[:(nxPML_1-1),:-1,1:-1] = Hz[:(nxPML_1-1),:-1,1:-1] + db * psi_Hzx_1[:,:,1:-1]
    
    psi_Hzx_2[ii_indices_h,:,1:-1] = (
        bh_x_2[ii_indices_h,:-1,1:-1] * psi_Hzx_2[ii_indices_h,:,1:-1] 
        + ch_x_2[ii_indices_h,:-1,1:-1] * (Ey[i_indices_h,:-1,1:-1] - Ey[i_plus_1_indices_h,:-1,1:-1]) / del_x
    )
    Hz[i_indices_h,:-1,1:-1] = Hz[i_indices_h,:-1,1:-1] + db * psi_Hzx_2[ii_indices_h,:,1:-1]
    
    # Y-direction PMLs
    # Hx update
    psi_Hxy_1[:-1,:,:-1] = bh_y_1[:-1,:,:-1] * psi_Hxy_1[:-1,:,:-1] + ch_y_1[:-1,:,:-1] * (Ez[:-1,:(nyPML_1-1),:-1] - Ez[:-1,1:(nyPML_1),:-1]) / del_y
    Hx[:-1,:(nyPML_1-1),:-1] = Hx[:-1,:(nyPML_1-1),:-1] + db * psi_Hxy_1[:-1,:,:-1]
    
    psi_Hxy_2[:-1, jj_indices_h, :-1] = (
        bh_y_2[:-1, jj_indices_h, :-1] * psi_Hxy_2[:-1, jj_indices_h, :-1]
        + ch_y_2[:-1, jj_indices_h, :-1] * (Ez[:-1, j_indices_h, :-1] - Ez[:-1, j_plus_1_indices_h, :-1]) / del_y
    )
    Hx[:-1, j_indices_h, :-1] = Hx[:-1, j_indices_h, :-1] + db * psi_Hxy_2[:-1, jj_indices_h, :-1]
    
    # Hz update
    psi_Hzy_1[:,:,1:-1] = bh_y_1[:-1,:,1:-1] * psi_Hzy_1[:,:,1:-1] + ch_y_1[:-1,:,1:-1] * (Ex[:-1,1:nyPML_1,1:-1] - Ex[:-1,:(nyPML_1-1),1:-1]) / del_y
    Hz[:-1,:(nyPML_1-1),1:-1] = Hz[:-1,:(nyPML_1-1),1:-1] + db * psi_Hzy_1[:,:,1:-1]
    
    psi_Hzy_2[:,jj_indices_h,1:-1] = (
        bh_y_2[:-1,jj_indices_h,1:-1] * psi_Hzy_2[:,jj_indices_h,1:-1] 
        + ch_y_2[:-1,jj_indices_h,1:-1] * (Ex[:-1,j_plus_1_indices_h,1:-1] - Ex[:-1,j_indices_h,1:-1]) / del_y
    )
    Hz[:-1,j_indices_h,1:-1] = Hz[:-1,j_indices_h,1:-1] + db * psi_Hzy_2[:,jj_indices_h,1:-1]
    
    # Z-direction PMLs
    # Hx update
    psi_Hxz_1[:-1,:,:] = bh_z_1[:-1,:-1,:] * psi_Hxz_1[:-1,:,:] + ch_z_1[:-1,:-1,:] * (Ey[:-1,:-1,1:(nzPML_1)] - Ey[:-1,:-1,:(nzPML_1-1)]) / del_z
    Hx[:-1,:-1,:(nzPML_1-1)] = Hx[:-1,:-1,:(nzPML_1-1)] + db * psi_Hxz_1[:-1,:,:]
    
    psi_Hxz_2[:-1,:,kk_indices_h] = (
        bh_z_2[:-1,:-1,kk_indices_h] * psi_Hxz_2[:-1,:,kk_indices_h] 
        + ch_z_2[:-1,:-1,kk_indices_h] * (Ey[:-1,:-1,k_plus_1_indices_h] - Ey[:-1,:-1,k_indices_h]) / del_z
    )
    Hx[:-1,:-1,k_indices_h] = Hx[:-1,:-1,k_indices_h] + db * psi_Hxz_2[:-1,:,kk_indices_h]
    
    # Hy update
    psi_Hyz_1[:,:-1,:] = bh_z_1[:-1,:-1,:] * psi_Hyz_1[:,:-1,:] + ch_z_1[:-1,:-1,:] * (Ex[:-1,:-1,:(nzPML_1-1)] - Ex[:-1,:-1,1:(nzPML_1)]) / del_z
    Hy[:-1,:-1,:(nzPML_1-1)] = Hy[:-1,:-1,:(nzPML_1-1)] + db * psi_Hyz_1[:,:-1,:]
    
    psi_Hyz_2[:,:-1,kk_indices_h] = (
        bh_z_2[:-1,:-1,kk_indices_h] * psi_Hyz_2[:,:-1,kk_indices_h] 
        + ch_z_2[:-1,:-1,kk_indices_h] * (Ex[:-1,:-1,k_indices_h] - Ex[:-1,:-1,k_plus_1_indices_h]) / del_z
    )
    Hy[:-1,:-1,k_indices_h] = Hy[:-1,:-1,k_indices_h] + db * psi_Hyz_2[:,:-1,kk_indices_h]
    

def update_E_fields(Ex, Ey, Ez, Hx, Hy, Hz, gax, gay, gaz, gbx, gby, gbz, 
                    den_ex, den_ey, den_ez):
    """
    Update E field components using curl of H fields.
    
    Inputs (read): Hx, Hy, Hz, ga*, gb*, den_ex, den_ey, den_ez
    Outputs (modified in-place): Ex, Ey, Ez
    """
    Ex[:-1, 1:-1, 1:-1] = gax[:-1, 1:-1, 1:-1] * Ex[:-1, 1:-1, 1:-1] + gbx[:-1, 1:-1, 1:-1] * (
        (Hy[:-1, 1:-1, :-2] - Hy[:-1, 1:-1, 1:-1]) * den_ez[None,None,1:] +
        (Hz[:-1, 1:-1, 1:-1] - Hz[:-1, :-2, 1:-1]) * den_ey[None,1:,None]
    )
    Ey[1:-1, :-1, 1:-1] = gay[1:-1, :-1, 1:-1] * Ey[1:-1, :-1, 1:-1] + gby[1:-1, :-1, 1:-1] * (
        (Hx[1:-1, :-1, 1:-1] - Hx[1:-1, :-1, :-2]) * den_ez[None,None,1:] +
        (Hz[:-2, :-1, 1:-1] - Hz[1:-1, :-1, 1:-1]) * den_ex[1:,None,None]
    )
    Ez[1:-1, 1:-1, :-1] = gaz[1:-1, 1:-1, :-1] * Ez[1:-1, 1:-1, :-1] + gbz[1:-1, 1:-1, :-1] * (
        (Hx[1:-1, :-2, :-1] - Hx[1:-1, 1:-1, :-1]) * den_ey[None,1:,None] +
        (Hy[1:-1, 1:-1, :-1] - Hy[:-2, 1:-1, :-1]) * den_ex[1:,None,None]
    )


def update_E_pml(Ex, Ey, Ez, Hx, Hy, Hz, gbx, del_x, del_y, del_z,
                 psi_Eyx_1, psi_Eyx_2, psi_Ezx_1, psi_Ezx_2,
                 psi_Exy_1, psi_Exy_2, psi_Ezy_1, psi_Ezy_2,
                 psi_Exz_1, psi_Exz_2, psi_Eyz_1, psi_Eyz_2,
                 be_x_1, ce_x_1, be_x_2, ce_x_2,
                 be_y_1, ce_y_1, be_y_2, ce_y_2,
                 be_z_1, ce_z_1, be_z_2, ce_z_2,
                 nxPML_1, nyPML_1, nzPML_1,
                 ii_indices_e, i_indices_e, i_minus_1_indices_e,
                 jj_indices_e, j_indices_e, j_minus_1_indices_e,
                 kk_indices_e, k_indices_e, k_minus_1_indices_e):
    """
    Apply PML boundary conditions to E fields.
    
    Inputs (read): Hx, Hy, Hz, gbx, del_x/y/z, be_*, ce_*, PML parameters, indices
    Outputs (modified in-place): Ex, Ey, Ez, psi_* arrays
    """
    # X-direction PMLs
    # Ey update
    psi_Eyx_1[1:,:,1:-1] = be_x_1[1:,:-1,1:-1] * psi_Eyx_1[1:,:,1:-1] + ce_x_1[1:,:-1,1:-1] * (Hz[:(nxPML_1-1),:-1,1:-1] - Hz[1:(nxPML_1),:-1,1:-1]) / del_x
    Ey[1:(nxPML_1),:-1,1:-1] = Ey[1:(nxPML_1),:-1,1:-1] + gbx[1:(nxPML_1),:-1,1:-1] * psi_Eyx_1[1:,:,1:-1]
    
    psi_Eyx_2[ii_indices_e,:,1:-1] = (
        be_x_2[ii_indices_e,:-1,1:-1] * psi_Eyx_2[ii_indices_e,:,1:-1] 
        + ce_x_2[ii_indices_e,:-1,1:-1] * (Hz[i_minus_1_indices_e,:-1,1:-1] - Hz[i_indices_e,:-1,1:-1]) / del_x
    )
    Ey[i_indices_e,:-1,1:-1] = Ey[i_indices_e,:-1,1:-1] + gbx[i_indices_e,:-1,1:-1] * psi_Eyx_2[ii_indices_e,:,1:-1]
    
    # Ez update
    psi_Ezx_1[1:,1:-1,:-1] = be_x_1[1:,1:-1,:-1] * psi_Ezx_1[1:,1:-1,:-1] + ce_x_1[1:,1:-1,:-1] * (Hy[1:nxPML_1,1:-1,:-1] - Hy[:(nxPML_1-1),1:-1,:-1]) / del_x
    Ez[1:nxPML_1,1:-1,:-1] = Ez[1:nxPML_1,1:-1,:-1] + gbx[1:nxPML_1,1:-1,:-1] * psi_Ezx_1[1:,1:-1,:-1]
    
    psi_Ezx_2[ii_indices_e,1:-1,:-1] = (
        be_x_2[ii_indices_e,1:-1,:-1] * psi_Ezx_2[ii_indices_e,1:-1,:-1] 
        + ce_x_2[ii_indices_e,1:-1,:-1] * (Hy[i_indices_e,1:-1,:-1] - Hy[i_minus_1_indices_e,1:-1,:-1]) / del_x
    )
    Ez[i_indices_e,1:-1,:-1] = Ez[i_indices_e,1:-1,:-1] + gbx[i_indices_e,1:-1,:-1] * psi_Ezx_2[ii_indices_e,1:-1,:-1]
    
    # Y-direction PMLs
    # Ex update
    psi_Exy_1[:,1:,1:-1] = be_y_1[:-1,1:,1:-1] * psi_Exy_1[:,1:,1:-1] + ce_y_1[:-1,1:,1:-1] * (Hz[:-1,1:nyPML_1,1:-1] - Hz[:-1,:(nyPML_1-1),1:-1]) / del_y
    Ex[:-1,1:nyPML_1,1:-1] = Ex[:-1,1:nyPML_1,1:-1] + gbx[:-1,1:nyPML_1,1:-1] * psi_Exy_1[:,1:,1:-1]
    
    psi_Exy_2[:,jj_indices_e,1:-1] = (
        be_y_2[:-1,jj_indices_e,1:-1] * psi_Exy_2[:,jj_indices_e,1:-1] 
        + ce_y_2[:-1,jj_indices_e,1:-1] * (Hz[:-1,j_indices_e,1:-1] - Hz[:-1,j_minus_1_indices_e,1:-1]) / del_y
    )
    Ex[:-1,j_indices_e,1:-1] = Ex[:-1,j_indices_e,1:-1] + gbx[:-1,j_indices_e,1:-1] * psi_Exy_2[:,jj_indices_e,1:-1]
    
    # Ez update
    psi_Ezy_1[1:-1,1:,:-1] = be_y_1[1:-1,1:,:-1] * psi_Ezy_1[1:-1,1:,:-1] + ce_y_1[1:-1,1:,:-1] * (Hx[1:-1,:(nyPML_1-1),:-1] - Hx[1:-1,1:nyPML_1,:-1]) / del_y
    Ez[1:-1,1:nyPML_1,:-1] = Ez[1:-1,1:nyPML_1,:-1] + gbx[1:-1,1:nyPML_1,:-1] * psi_Ezy_1[1:-1,1:,:-1]
    
    psi_Ezy_2[1:-1,jj_indices_e,:-1] = (
        be_y_2[1:-1,jj_indices_e,:-1] * psi_Ezy_2[1:-1,jj_indices_e,:-1] 
        + ce_y_2[1:-1,jj_indices_e,:-1] * (Hx[1:-1,j_minus_1_indices_e,:-1] - Hx[1:-1,j_indices_e,:-1]) / del_y
    )
    Ez[1:-1,j_indices_e,:-1] = Ez[1:-1,j_indices_e,:-1] + gbx[1:-1,j_indices_e,:-1] * psi_Ezy_2[1:-1,jj_indices_e,:-1]
    
    # Z-direction PMLs
    # Ex update
    psi_Exz_1[:,1:-1,1:] = be_z_1[:-1,1:-1,1:] * psi_Exz_1[:,1:-1,1:] + ce_z_1[:-1,1:-1,1:] * (Hy[:-1,1:-1,:(nzPML_1-1)] - Hy[:-1,1:-1,1:nzPML_1]) / del_z
    Ex[:-1,1:-1,1:(nzPML_1)] = Ex[:-1,1:-1,1:(nzPML_1)] + gbx[:-1,1:-1,1:(nzPML_1)] * psi_Exz_1[:,1:-1,1:]
    
    psi_Exz_2[:,1:-1,kk_indices_e] = (
        be_z_2[:-1,1:-1,kk_indices_e] * psi_Exz_2[:,1:-1,kk_indices_e] 
        + ce_z_2[:-1,1:-1,kk_indices_e] * (Hy[:-1,1:-1,k_minus_1_indices_e] - Hy[:-1,1:-1,k_indices_e]) / del_z
    )
    Ex[:-1,1:-1,k_indices_e] = Ex[:-1,1:-1,k_indices_e] + gbx[:-1,1:-1,k_indices_e] * psi_Exz_2[:,1:-1,kk_indices_e]
    
    # Ey update
    psi_Eyz_1[1:-1,:,1:] = be_z_1[1:-1,:-1,1:] * psi_Eyz_1[1:-1,:,1:] + ce_z_1[1:-1,:-1,1:] * (Hx[1:-1,:-1,1:(nzPML_1)] - Hx[1:-1,:-1,:(nzPML_1-1)]) / del_z
    Ey[1:-1,:-1,1:(nzPML_1)] = Ey[1:-1,:-1,1:(nzPML_1)] + gbx[1:-1,:-1,1:(nzPML_1)] * psi_Eyz_1[1:-1,:,1:]
    
    psi_Eyz_2[1:-1,:,kk_indices_e] = (
        be_z_2[1:-1,:-1,kk_indices_e] * psi_Eyz_2[1:-1,:,kk_indices_e] 
        + ce_z_2[1:-1,:-1,kk_indices_e] * (Hx[1:-1,:-1,k_indices_e] - Hx[1:-1,:-1,k_minus_1_indices_e]) / del_z
    )
    Ey[1:-1,:-1,k_indices_e] = Ey[1:-1,:-1,k_indices_e] + gbx[1:-1,:-1,k_indices_e] * psi_Eyz_2[1:-1,:,kk_indices_e]


def compute_source_pulse(counter, del_t, pulse_type, t_spread, spread):
    """
    Compute the source pulse value at a given time step.
    
    Inputs: counter, del_t, pulse_type, t_spread, spread
    Returns: pulse value (scalar)
    """
    t = counter * del_t
    if pulse_type == 1:
        pulse = np.exp(-0.5 * ((t_spread - t) / spread) ** 2)
    elif pulse_type == 2:
        pulse = -1 * ((t_spread - t) / spread) * np.exp(0.5) * \
                np.exp(-0.5 * ((t_spread - t) / spread) ** 2)
    return pulse


def E_inc(weight, x, y, z):
    
    time_arg = (counter - 1.0) * del_t
    spatial_phase = (
        (th.sin(theta) * th.cos(phi) * (x - x_delay) * del_x) +
        (th.sin(theta) * th.sin(phi) * (y - y_delay) * del_y) +
        (th.cos(theta) * (z - z_delay) * del_z)
    ) / c
    term = (t_spread - time_arg + spatial_phase) / spread
    
    if pulse_type == 1:
        E_inc_val = weight * th.exp(-0.5 * term**2)
    elif pulse_type == 2:
        E_inc_val = -1.0 * weight * term * th.exp(th.tensor(0.5, device=device)) * th.exp(-0.5 * term**2)
    else:
        E_inc_val = th.zeros_like(x)
        
    return E_inc_val


def H_inc(weight, x, y, z):
    
    time_arg = (counter - 0.5) * del_t
    spatial_phase = (
        (th.sin(theta) * th.cos(phi) * (x - x_delay) * del_x) +
        (th.sin(theta) * th.sin(phi) * (y - y_delay) * del_y) +
        (th.cos(theta) * (z - z_delay) * del_z)
    ) / c
    term = (t_spread - time_arg + spatial_phase) / spread

    if pulse_type == 1:
        H_inc_val = weight * th.exp(-0.5 * term**2)
    elif pulse_type == 2:
        H_inc_val = -1.0 * weight * term * th.exp(th.tensor(0.5, device=device)) * th.exp(-0.5 * term**2)
    else:
        H_inc_val = th.zeros_like(x)

    return H_inc_val


def H_plane_waves(Hz, Hy, Hx, E_inc, WEx, WEy, WEz,
                                   xlow, xhigh, ylow, yhigh, zlow, zhigh,
                                   del_t, mu_0, del_x, del_y, del_z,
                                   xlow_wall, xhigh_wall, ylow_wall, 
                                   yhigh_wall, zlow_wall, zhigh_wall):
    # Pre-compute common factors
    coeff_x = del_t / (mu_0 * del_x)
    coeff_y = del_t / (mu_0 * del_y)
    coeff_z = del_t / (mu_0 * del_z)
    
    # ===== Y FACES =====
    # First Y face loop: Hz updates
    i_range = th.arange(xlow, xhigh, device=device)
    k_range = th.arange(zlow, zhigh + 1, device=device)
    i_grid, k_grid = th.meshgrid(i_range, k_range, indexing='ij')
    
    E_inc_vals = E_inc(WEx, i_grid + 0.5, ylow + 0.0, k_grid + 0.0)
    Hz[i_range, ylow-1, zlow:zhigh+1] -= coeff_y * E_inc_vals * ylow_wall
    
    E_inc_vals = E_inc(WEx, i_grid + 0.5, yhigh + 0.0, k_grid + 0.0)
    Hz[i_range, yhigh, zlow:zhigh+1] += coeff_y * E_inc_vals * yhigh_wall
    
    # Second Y face loop: Hx updates
    i_range = th.arange(xlow, xhigh + 1, device=device)
    k_range = th.arange(zlow, zhigh, device=device)
    i_grid, k_grid = th.meshgrid(i_range, k_range, indexing='ij')
    
    E_inc_vals = E_inc(WEz, i_grid + 0.0, ylow + 0.0, k_grid + 0.5)
    Hx[i_range, ylow-1, zlow:zhigh] += coeff_y * E_inc_vals * ylow_wall
    
    E_inc_vals = E_inc(WEz, i_grid + 0.0, yhigh + 0.0, k_grid + 0.5)
    Hx[i_range, yhigh, zlow:zhigh] -= coeff_y * E_inc_vals * yhigh_wall
    
    # ===== Z FACES =====
    # First Z face loop: Hy updates
    i_range = th.arange(xlow, xhigh, device=device)
    j_range = th.arange(ylow, yhigh + 1, device=device)
    i_grid, j_grid = th.meshgrid(i_range, j_range, indexing='ij')
    
    E_inc_vals = E_inc(WEx, i_grid + 0.5, j_grid + 0.0, zlow + 0.0)
    Hy[i_range, ylow:yhigh+1, zlow-1] += coeff_z * E_inc_vals * zlow_wall
    
    E_inc_vals = E_inc(WEx, i_grid + 0.5, j_grid + 0.0, zhigh + 0.0)
    Hy[i_range, ylow:yhigh+1, zhigh] -= coeff_z * E_inc_vals * zhigh_wall
    
    # Second Z face loop: Hx updates
    i_range = th.arange(xlow, xhigh + 1, device=device)
    j_range = th.arange(ylow, yhigh, device=device)
    i_grid, j_grid = th.meshgrid(i_range, j_range, indexing='ij')
    
    E_inc_vals = E_inc(WEy, i_grid + 0.0, j_grid + 0.5, zlow + 0.0)
    Hx[i_range, ylow:yhigh, zlow-1] -= coeff_z * E_inc_vals * zlow_wall
    
    E_inc_vals = E_inc(WEy, i_grid + 0.0, j_grid + 0.5, zhigh + 0.0)
    Hx[i_range, ylow:yhigh, zhigh] += coeff_z * E_inc_vals * zhigh_wall
    
    # ===== X FACES =====
    # First X face loop: Hz updates
    j_range = th.arange(ylow, yhigh, device=device)
    k_range = th.arange(zlow, zhigh + 1, device=device)
    j_grid, k_grid = th.meshgrid(j_range, k_range, indexing='ij')
    
    E_inc_vals = E_inc(WEy, xlow + 0.0, j_grid + 0.5, k_grid + 0.0)
    Hz[xlow-1, ylow:yhigh, zlow:zhigh+1] += coeff_x * E_inc_vals * xlow_wall
    
    E_inc_vals = E_inc(WEy, xhigh + 0.0, j_grid + 0.5, k_grid + 0.0)
    Hz[xhigh, ylow:yhigh, zlow:zhigh+1] -= coeff_x * E_inc_vals * xhigh_wall
    
    # Second X face loop: Hy updates
    j_range = th.arange(ylow, yhigh + 1, device=device)
    k_range = th.arange(zlow, zhigh, device=device)
    j_grid, k_grid = th.meshgrid(j_range, k_range, indexing='ij')
    
    E_inc_vals = E_inc(WEz, xlow + 0.0, j_grid + 0.0, k_grid + 0.5)
    Hy[xlow-1, ylow:yhigh+1, zlow:zhigh] -= coeff_x * E_inc_vals * xlow_wall
    
    E_inc_vals = E_inc(WEz, xhigh + 0.0, j_grid + 0.0, k_grid + 0.5)
    Hy[xhigh, ylow:yhigh+1, zlow:zhigh] += coeff_x * E_inc_vals * xhigh_wall


def E_plane_waves(Ez, Ey, Ex, H_inc, WHx, WHy, WHz,
                                    xlow, xhigh, ylow, yhigh, zlow, zhigh,
                                    del_t, ep_0, del_x, del_y, del_z,
                                    xlow_wall, xhigh_wall, ylow_wall, 
                                    yhigh_wall, zlow_wall, zhigh_wall):
    
    # Pre-compute common factors
    coeff_x = del_t / (ep_0 * del_x)
    coeff_y = del_t / (ep_0 * del_y)
    coeff_z = del_t / (ep_0 * del_z)
    
    # ===== Y FACES =====
    # First Y face loop: Ez updates
    i_range = th.arange(xlow, xhigh + 1, device=device)
    k_range = th.arange(zlow, zhigh, device=device)
    i_grid, k_grid = th.meshgrid(i_range, k_range, indexing='ij')
    
    H_inc_vals = H_inc(WHx, i_grid + 0.0, ylow - 0.5, k_grid + 0.5)
    Ez[i_range, ylow, zlow:zhigh] += coeff_y * H_inc_vals * ylow_wall
    
    H_inc_vals = H_inc(WHx, i_grid + 0.0, yhigh + 0.5, k_grid + 0.5)
    Ez[i_range, yhigh, zlow:zhigh] -= coeff_y * H_inc_vals * yhigh_wall
    
    # Second Y face loop: Ex updates
    i_range = th.arange(xlow, xhigh, device=device)
    k_range = th.arange(zlow, zhigh + 1, device=device)
    i_grid, k_grid = th.meshgrid(i_range, k_range, indexing='ij')
    
    H_inc_vals = H_inc(WHz, i_grid + 0.5, ylow - 0.5, k_grid + 0.0)
    Ex[i_range, ylow, zlow:zhigh+1] -= coeff_y * H_inc_vals * ylow_wall
    
    H_inc_vals = H_inc(WHz, i_grid + 0.5, yhigh + 0.5, k_grid + 0.0)
    Ex[i_range, yhigh, zlow:zhigh+1] += coeff_y * H_inc_vals * yhigh_wall
    
    # ===== Z FACES =====
    # First Z face loop: Ey updates
    i_range = th.arange(xlow, xhigh + 1, device=device)
    j_range = th.arange(ylow, yhigh, device=device)
    i_grid, j_grid = th.meshgrid(i_range, j_range, indexing='ij')
    
    H_inc_vals = H_inc(WHx, i_grid + 0.0, j_grid + 0.5, zlow - 0.5)
    Ey[i_range, ylow:yhigh, zlow] -= coeff_z * H_inc_vals * zlow_wall
    
    H_inc_vals = H_inc(WHx, i_grid + 0.0, j_grid + 0.5, zhigh + 0.5)
    Ey[i_range, ylow:yhigh, zhigh] += coeff_z * H_inc_vals * zhigh_wall
    
    # Second Z face loop: Ex updates
    i_range = th.arange(xlow, xhigh, device=device)
    j_range = th.arange(ylow, yhigh + 1, device=device)
    i_grid, j_grid = th.meshgrid(i_range, j_range, indexing='ij')
    
    H_inc_vals = H_inc(WHy, i_grid + 0.5, j_grid + 0.0, zlow - 0.5)
    Ex[i_range, ylow:yhigh+1, zlow] += coeff_z * H_inc_vals * zlow_wall
    
    H_inc_vals = H_inc(WHy, i_grid + 0.5, j_grid + 0.0, zhigh + 0.5)
    Ex[i_range, ylow:yhigh+1, zhigh] -= coeff_z * H_inc_vals * zhigh_wall
    
    # ===== X FACES =====
    # First X face loop: Ez updates
    j_range = th.arange(ylow, yhigh + 1, device=device)
    k_range = th.arange(zlow, zhigh, device=device)
    j_grid, k_grid = th.meshgrid(j_range, k_range, indexing='ij')
    
    H_inc_vals = H_inc(WHy, xlow - 0.5, j_grid + 0.0, k_grid + 0.5)
    Ez[xlow, ylow:yhigh+1, zlow:zhigh] -= coeff_x * H_inc_vals * xlow_wall
    
    H_inc_vals = H_inc(WHy, xhigh + 0.5, j_grid + 0.0, k_grid + 0.5)
    Ez[xhigh, ylow:yhigh+1, zlow:zhigh] += coeff_x * H_inc_vals * xhigh_wall
    
    # Second X face loop: Ey updates
    j_range = th.arange(ylow, yhigh, device=device)
    k_range = th.arange(zlow, zhigh + 1, device=device)
    j_grid, k_grid = th.meshgrid(j_range, k_range, indexing='ij')
    
    H_inc_vals = H_inc(WHz, xlow - 0.5, j_grid + 0.5, k_grid + 0.0)
    Ey[xlow, ylow:yhigh, zlow:zhigh+1] += coeff_x * H_inc_vals * xlow_wall
    
    H_inc_vals = H_inc(WHz, xhigh + 0.5, j_grid + 0.5, k_grid + 0.0)
    Ey[xhigh, ylow:yhigh, zlow:zhigh+1] -= coeff_x * H_inc_vals * xhigh_wall


def update_sheets_H(Ex,Ex_special,Ey,Ey_special,Ez,Ez_special,Hx,Hy,Hz,den_hx,den_hy,den_hz,sheet_thickness,del_x,del_y,del_z,mask_x,mask_y,mask_z):
    # -------------------------------
    # X-normal sheets update
    # -------------------------------
    if mask_x.any():
        # Hy update
        Hy[:-1, :-1, :-1] += db * (sheet_thickness / del_x) * mask_x[:-1, :-1, :-1] * (
            Ex[:-1, :-1, 1:] - Ex[:-1, :-1, :-1] + Ex_special[:-1, :-1, :-1] - Ex_special[:-1, :-1, 1:]
        ) * den_hz[None, None, :]

        # Hz update
        Hz[:-1, :-1, :-1] += db * (sheet_thickness / del_x) * mask_x[:-1, :-1, :-1] * (
            Ex[:-1, :-1, :-1] - Ex[:-1, 1:, :-1] + Ex_special[:-1, 1:, :-1] - Ex_special[:-1, :-1, :-1]
        ) * den_hy[None, :, None]

    # -------------------------------
    # Y-normal sheets update
    # -------------------------------
    if mask_y.any():
        # Hx update
        Hx[:-1, :-1, :-1] += db * (sheet_thickness / del_y) * mask_y[:-1, :-1, :-1] * (
            Ey[:-1, :-1, :-1] - Ey[:-1, :-1, 1:] + Ey_special[:-1, :-1, 1:] - Ey_special[:-1, :-1, :-1]
        ) * den_hz[None, None, :]

        # Hz update
        Hz[:-1, :-1, :-1] += db * (sheet_thickness / del_y) * mask_y[:-1, :-1, :-1] * (
            Ey[1:, :-1, :-1] - Ey[:-1, :-1, :-1] - Ey_special[1:, :-1, :-1] + Ey_special[:-1, :-1, :-1]
        ) * den_hx[:, None, None]

    # -------------------------------
    # Z-normal sheets update
    # -------------------------------
    if mask_z.any():
        # Hx update
        Hx[:-1, :-1, :-1] += db * (sheet_thickness / del_z) * mask_z[:-1, :-1, :-1] * (
            Ez[:-1, 1:, :-1] - Ez[:-1, :-1, :-1] - Ez_special[:-1, 1:, :-1] + Ez_special[:-1, :-1, :-1]
        ) * den_hy[None, :, None]

        # Hy update
        Hy[:-1, :-1, :-1] += db * (sheet_thickness / del_z) * mask_z[:-1, :-1, :-1] * (
            Ez[:-1, :-1, :-1] - Ez[1:, :-1, :-1] + Ez_special[1:, :-1, :-1] - Ez_special[:-1, :-1, :-1]
        ) * den_hx[:, None, None]


def update_sheets_E(Ex,Ex_special,Ey,Ey_special,Ez,Ez_special,Hx,Hy,Hz,den_ex,den_ey,den_ez,del_t,ep_0,mask_x,mask_y,mask_z):
    # -------------------------------
    # X-normal sheets: Ex_special
    # -------------------------------
    if mask_x.any():
        coef_x = 1.0 / (1.0 + sheet_sig_x_x[:-1, 1:-1, 1:-1] * del_t / (2.0 * ep_0 * sheet_ep_x_x[:-1, 1:-1, 1:-1]))
        Ex_special[:-1, 1:-1, 1:-1] = ((1.0 - sheet_sig_x_x[:-1, 1:-1, 1:-1] * del_t / (2.0 * ep_0 * sheet_ep_x_x[:-1, 1:-1, 1:-1])) *
                                    Ex_special[:-1, 1:-1, 1:-1]*coef_x + 
                                    (del_t / (ep_0 * sheet_ep_x_x[:-1, 1:-1, 1:-1])) * 
                                    ((Hz[:-1, 1:-1, 1:-1] - Hz[:-1, :-2, 1:-1]) * den_ey[None,1:,None] +
                                    (Hy[:-1, 1:-1, :-2] - Hy[:-1, 1:-1, 1:-1]) * den_ez[None,None,1:])
                                ) * mask_x[:-1, 1:-1, 1:-1]

    # -------------------------------
    # Y-normal sheets: Ey_special
    # -------------------------------
    if mask_y.any():
        coef_y = 1.0 / (1.0 + sheet_sig_y_y[1:-1, :-1, 1:-1]* del_t / (2.0 * ep_0 * sheet_ep_y_y[1:-1, :-1, 1:-1]))
        Ey_special[1:-1, :-1, 1:-1] = ((1.0 - sheet_sig_y_y[1:-1, :-1, 1:-1] * del_t / (2.0 * ep_0 * sheet_ep_y_y[1:-1, :-1, 1:-1])) *
                                    Ey_special[1:-1, :-1, 1:-1]*coef_y +
                                    (del_t / (ep_0 * sheet_ep_y_y[1:-1, :-1, 1:-1])) *
                                    ((Hz[:-2, :-1, 1:-1] - Hz[1:-1, :-1, 1:-1]) * den_ex[1:,None, None] +
                                    (Hx[1:-1, :-1, 1:-1] - Hx[1:-1, :-1, :-2]) * den_ez[None, None, 1:])
                                ) * mask_y[1:-1, :-1, 1:-1]

    # -------------------------------
    # Z-normal sheets: Ez_special
    # -------------------------------
    if mask_z.any():
        coef_z = 1.0 / (1.0 + sheet_sig_z_z[1:-1, 1:-1, :-1] * del_t / (2.0 * ep_0 * sheet_ep_z_z[1:-1, 1:-1, :-1]))
        Ez_special[1:-1, 1:-1, :-1] = ((1.0 - sheet_sig_z_z[1:-1, 1:-1, :-1] * del_t / (2.0 * ep_0 * sheet_ep_z_z[1:-1, 1:-1, :-1])) *
                                Ez_special[1:-1, 1:-1, :-1]*coef_z +
                                (del_t / (ep_0 * sheet_ep_z_z[1:-1, 1:-1, :-1])) *
                                ((Hy[1:-1, 1:-1, :-1] - Hy[:-2, 1:-1, :-1]) * den_ex[1:,None, None] +
                                (Hx[1:-1, :-2, :-1] - Hx[1:-1, 1:-1, :-1]) * den_ey[None, 1:,None])
                                ) * mask_z[1:-1, 1:-1, :-1]


##########################################################
# MAIN FDTD SOLVER
##########################################################
if compile_functions==True:
    print("\nCompiling functions with th.compile()...")
    update_H_fields = th.compile(update_H_fields)
    update_E_fields = th.compile(update_E_fields)
    update_H_pml = th.compile(update_H_pml)
    update_E_pml = th.compile(update_E_pml)
    # can do the same thing for any other definitions we throw in the loop 
    # but some can be quite slow to initialize here and for loop iteration counter=0
print("\nStarting FDTD solver...")
with GPUProfiler("FDTD Solver", device):

    # Start time loop
    for counter in range(time_steps):
        
        # Update H fields
        update_H_fields(Hx, Hy, Hz, Ex, Ey, Ez, da, db, den_hx, den_hy, den_hz)

        # update sheets for H fields if applicable
        update_sheets_H(Ex,Ex_special,Ey,Ey_special,Ez,Ez_special,Hx,Hy,Hz,den_hx,den_hy,den_hz,
                                sheet_thickness,del_x,del_y,del_z,mask_x,mask_y,mask_z)

        # Add plane wave excitation to H fields
        H_plane_waves(Hz, Hy, Hx, E_inc, WEx, WEy, WEz,
                                   xlow, xhigh, ylow, yhigh, zlow, zhigh,
                                   del_t, mu_0, del_x, del_y, del_z,
                                   xlow_wall, xhigh_wall, ylow_wall, 
                                   yhigh_wall, zlow_wall, zhigh_wall)
        
        # Apply PML to H fields
        update_H_pml(Hx, Hy, Hz, Ex, Ey, Ez, db, del_x, del_y, del_z,
                     psi_Hyx_1, psi_Hyx_2, psi_Hzx_1, psi_Hzx_2,
                     psi_Hxy_1, psi_Hxy_2, psi_Hzy_1, psi_Hzy_2,
                     psi_Hxz_1, psi_Hxz_2, psi_Hyz_1, psi_Hyz_2,
                     bh_x_1, ch_x_1, bh_x_2, ch_x_2,
                     bh_y_1, ch_y_1, bh_y_2, ch_y_2,
                     bh_z_1, ch_z_1, bh_z_2, ch_z_2,
                     nxPML_1, nyPML_1, nzPML_1,
                     ii_indices_h, i_indices_h, i_plus_1_indices_h,
                     jj_indices_h, j_indices_h, j_plus_1_indices_h,
                     kk_indices_h, k_indices_h, k_plus_1_indices_h)
        
        # Update E fields
        update_E_fields(Ex, Ey, Ez, Hx, Hy, Hz, gax, gay, gaz, gbx, gby, gbz,
                       den_ex, den_ey, den_ez)

        # update sheets for E fields if applicable
        update_sheets_E(Ex,Ex_special,Ey,Ey_special,Ez,Ez_special,Hx,Hy,Hz,den_ex,den_ey,den_ez,
                        del_t,ep_0,mask_x,mask_y,mask_z)

        # Add plane wave excitation to E fields
        E_plane_waves(Ez, Ey, Ex, H_inc, WHx, WHy, WHz,
                                    xlow, xhigh, ylow, yhigh, zlow, zhigh,
                                    del_t, ep_0, del_x, del_y, del_z,
                                    xlow_wall, xhigh_wall, ylow_wall, 
                                    yhigh_wall, zlow_wall, zhigh_wall)
        
        # Apply PML to E fields (using compiled version)
        update_E_pml(Ex, Ey, Ez, Hx, Hy, Hz, gbx, del_x, del_y, del_z,
                     psi_Eyx_1, psi_Eyx_2, psi_Ezx_1, psi_Ezx_2,
                     psi_Exy_1, psi_Exy_2, psi_Ezy_1, psi_Ezy_2,
                     psi_Exz_1, psi_Exz_2, psi_Eyz_1, psi_Eyz_2,
                     be_x_1, ce_x_1, be_x_2, ce_x_2,
                     be_y_1, ce_y_1, be_y_2, ce_y_2,
                     be_z_1, ce_z_1, be_z_2, ce_z_2,
                     nxPML_1, nyPML_1, nzPML_1,
                     ii_indices_e, i_indices_e, i_minus_1_indices_e,
                     jj_indices_e, j_indices_e, j_minus_1_indices_e,
                     kk_indices_e, k_indices_e, k_minus_1_indices_e)
        
        # Apply source
        # pulse = compute_source_pulse(counter, del_t, pulse_type, t_spread, spread)
        # pulse_tensor = th.tensor(pulse, dtype=dtype, device=device)
        # Ez[40, 40, 40] = pulse_tensor
        # input_data[counter] = pulse_tensor
        
        # Progress tracking
        if counter % 100 == 0:
            print(f"{counter} of {time_steps} time steps")

#this section will be heavily modified later on - just here for testing purposes now
##########################################################
# POST PROCESSING
##########################################################
with GPUProfiler("Post Processing", device):
    import numpy as np
    # Transfer only final results to CPU for plotting
    Ey_cpu = Ez.cpu().numpy()
    #input_cpu = input_data.cpu().numpy()
    
    # Example: Plot a slice
    #x=np.linspace(1,x_size,x_size)
    #y=np.linspace(1,y_size,y_size)
    #z=th.linspace(1,z_size,z_size)
    #X,Y=np.meshgrid(y,x)
    plt.pcolormesh(Ey_cpu[:,:,35], cmap='bwr', shading='gouraud')
    plt.colorbar()
    #plt.title('Ez field at z=40')
    plt.savefig('fdtd_result.png')
    # plt.close()

print("\n" + "="*60)
print("Simulation complete#")
print("="*60)

# Check GPU utilization
if device.type == 'cuda':
    print(f"\nFinal GPU memory usage: {th.cuda.memory_allocated() / 1e9:.2f} GB")
    print(f"Peak GPU memory usage: {th.cuda.max_memory_allocated() / 1e9:.2f} GB")