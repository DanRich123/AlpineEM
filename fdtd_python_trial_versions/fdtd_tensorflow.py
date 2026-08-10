import time as tp
from matplotlib import pyplot as plt
import tensorflow as tf
import math
import sys

import os
os.chdir('C:/Users/dari6475/Desktop/temp/prep_to_backup/fdtd_python')

##########################################################
# GPU SETUP AND VERIFICATION
##########################################################
def setup_gpu():
    """Setup and verify GPU availability"""
    gpus = tf.config.list_physical_devices('GPU')
    if gpus:
        try:
            # Enable memory growth to prevent TF from allocating all GPU memory at once
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
            
            print(f"✓ GPU detected: {gpus[0].name}")
            print(f"✓ TensorFlow version: {tf.__version__}")
            print(f"✓ Number of GPUs available: {len(gpus)}")
            
            return '/GPU:0'
        except RuntimeError as e:
            print(e)
            return '/CPU:0'
    else:
        print("⚠ GPU not available, using CPU (will be much slower)")
        return '/CPU:0'

device = setup_gpu()
dtype = tf.float32
jit_on=True

##########################################################
# PROFILING UTILITIES
##########################################################
class GPUProfiler:
    """Context manager for profiling GPU operations"""
    def __init__(self, name, device):
        self.name = name
        self.device = device
        
    def __enter__(self):
        # TensorFlow doesn't need explicit synchronization like CUDA
        self.start = tp.time()
        return self
        
    def __exit__(self, *args):
        # Force execution of any pending operations
        tf.test.experimental.sync_devices()
        self.elapsed = tp.time() - self.start
        print(f"{self.name}: {self.elapsed:.4f} seconds")
        if 'GPU' in self.device:
            # TensorFlow's memory info is less granular than PyTorch
            print(f"  GPU memory info: Use nvidia-smi for detailed stats")

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

    xlow=0
    ylow=0
    zlow=0
    xhigh=0
    yhigh=0
    zhigh=0

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
            with tf.device(device):
                materials_properties = tf.Variable(tf.zeros((num_total_materials, 3, 3)))
                sheet_properties = tf.Variable(tf.zeros((num_total_materials, 4, 3)))

        for ii in range(num_total_materials):
            material_type = f.readline()

            if material_type.strip() == 'volume':
                num_materials += 1

                type_ingored_for_now=f.readline()
                materials_properties[i, 0, 0].assign(int(f.readline()))

                vals = list(map(float, f.readline().split(',')))
                materials_properties[i, 1, 0].assign(vals[0])
                materials_properties[i, 1, 1].assign(vals[1])
                materials_properties[i, 1, 2].assign(vals[2])

                vals = list(map(float, f.readline().split(',')))
                materials_properties[i, 2, 0].assign(vals[0])
                materials_properties[i, 2, 1].assign(vals[1])
                materials_properties[i, 2, 2].assign(vals[2])

                i += 1

            elif material_type.strip() == 'sheet':
                num_sheet_materials += 1

                sheet_properties[j, 0, 0].assign(int(f.readline()))
                imped_value = float(f.readline())
                sheet_thickness = min(del_x, del_y, del_z) / 5000.0

                sheet_properties[j, 2, 0].assign(1.0)
                sheet_properties[j, 2, 1].assign(1.0)
                sheet_properties[j, 2, 2].assign(1.0)

                if imped_value == 0:
                    sheet_properties[j, 3, 0].assign(1E20)
                    sheet_properties[j, 3, 1].assign(1E20)
                    sheet_properties[j, 3, 2].assign(1E20)
                else:
                    val = 1.0 / (imped_value * sheet_thickness)
                    sheet_properties[j, 3, 0].assign(val)
                    sheet_properties[j, 3, 1].assign(val)
                    sheet_properties[j, 3, 2].assign(val)

                j += 1

        if num_materials > 0 and num_materials < num_total_materials:
            materials_properties = materials_properties[:num_materials]
        elif num_materials == 0 and materials_properties is not None:
            del materials_properties

        if num_sheet_materials > 0 and num_sheet_materials < num_total_materials:
            sheet_properties = sheet_properties[:num_sheet_materials]
        elif num_sheet_materials == 0 and sheet_properties is not None:
            del sheet_properties


        num_objects = int(f.readline())

        if num_objects > 0:
            object_type = [''] * num_objects
            with tf.device(device):
                spheres = tf.Variable(tf.zeros((num_objects, 3, 3)))
                blocks = tf.Variable(tf.zeros((num_objects, 3, 3)))
                cylinders = tf.Variable(tf.zeros((num_objects, 4, 3)))

        i = j = k = 0

        for ii in range(num_objects):
            object_type[ii] = f.readline()

            if object_type[ii].strip() == 'block':
                num_blocks += 1

                blocks[i, 0, 0].assign(int(f.readline()))

                vals = list(map(int, f.readline().split(',')))
                blocks[i, 1, 0].assign(vals[0])
                blocks[i, 1, 1].assign(vals[1])
                blocks[i, 1, 2].assign(vals[2])

                vals = list(map(int, f.readline().split(',')))
                blocks[i, 2, 0].assign(vals[0])
                blocks[i, 2, 1].assign(vals[1])
                blocks[i, 2, 2].assign(vals[2])

                i += 1

            elif object_type[ii].strip() == 'sphere':
                num_spheres += 1

                spheres[j, 0, 0].assign(int(f.readline()))

                vals = list(map(int, f.readline().split(',')))
                spheres[j, 1, 0].assign(vals[0])
                spheres[j, 1, 1].assign(vals[1])
                spheres[j, 1, 2].assign(vals[2])

                spheres[j, 2, 0].assign(int(f.readline()))
                j += 1

            elif object_type[ii].strip() == 'cylinder':
                num_cylinders += 1

                cylinders[k, 0, 0].assign(int(f.readline()))
                char_cylinders = f.readline().strip()

                if char_cylinders == 'x':
                    cylinders[k, 1, 0].assign(0)
                elif char_cylinders == 'y':
                    cylinders[k, 1, 0].assign(1)
                elif char_cylinders == 'z':
                    cylinders[k, 1, 0].assign(2)

                vals = list(map(int, f.readline().split(',')))
                cylinders[k, 2, 0].assign(vals[0])
                cylinders[k, 2, 1].assign(vals[1])
                cylinders[k, 2, 2].assign(vals[2])

                vals = list(map(int, f.readline().split(',')))
                cylinders[k, 3, 0].assign(vals[0])
                cylinders[k, 3, 1].assign(vals[1])

                k += 1


        num_sheets = int(f.readline())

        if num_sheets > 0:
            with tf.device(device):
                sheets_x = tf.Variable(tf.zeros((num_sheets, 4, 2)))
                sheets_y = tf.Variable(tf.zeros((num_sheets, 4, 2)))
                sheets_z = tf.Variable(tf.zeros((num_sheets, 4, 2)))

        i = j = k = 0

        for ii in range(num_sheets):
            sheet_position_index = f.readline().strip()

            if sheet_position_index == 'x':
                num_sheets_x += 1
                sheets_x[i, 0, 0].assign(int(f.readline()))
                sheets_x[i, 1, 0].assign(int(f.readline()))
                vals = list(map(int, f.readline().split(',')))
                sheets_x[i, 2, 0].assign(vals[0])
                sheets_x[i, 2, 1].assign(vals[1])
                vals = list(map(int, f.readline().split(',')))
                sheets_x[i, 3, 0].assign(vals[0])
                sheets_x[i, 3, 1].assign(vals[1])
                i += 1

            elif sheet_position_index == 'y':
                num_sheets_y += 1
                sheets_y[j, 0, 0].assign(int(f.readline()))
                sheets_y[j, 1, 0].assign(int(f.readline()))
                vals = list(map(int, f.readline().split(',')))
                sheets_y[j, 2, 0].assign(vals[0])
                sheets_y[j, 2, 1].assign(vals[1])
                vals = list(map(int, f.readline().split(',')))
                sheets_y[j, 3, 0].assign(vals[0])
                sheets_y[j, 3, 1].assign(vals[1])
                j += 1

            elif sheet_position_index == 'z':
                num_sheets_z += 1
                sheets_z[k, 0, 0].assign(int(f.readline()))
                sheets_z[k, 1, 0].assign(int(f.readline()))
                vals = list(map(int, f.readline().split(',')))
                sheets_z[k, 2, 0].assign(vals[0])
                sheets_z[k, 2, 1].assign(vals[1])
                vals = list(map(int, f.readline().split(',')))
                sheets_z[k, 3, 0].assign(vals[0])
                sheets_z[k, 3, 1].assign(vals[1])
                k += 1


        num_ports = int(f.readline())

        if num_ports > 0:
            with tf.device(device):
                ports = tf.Variable(tf.zeros((num_ports, 4, 4)))

        for i in range(num_ports):
            char_ports = f.readline().strip()

            if char_ports == 'x':
                ports[i, 0, 0].assign(0)
            elif char_ports == 'y':
                ports[i, 0, 0].assign(1)
            elif char_ports == 'z':
                ports[i, 0, 0].assign(2)

            imped_port = float(f.readline())
            speed_port = c

            ports[i, 1, 1].assign(imped_port / speed_port)
            ports[i, 1, 2].assign(1.0 / (imped_port * speed_port))

            vals = list(map(int, f.readline().split(',')))
            ports[i, 2, 0].assign(vals[0])
            ports[i, 2, 1].assign(vals[1])
            ports[i, 2, 2].assign(vals[2])

            vals = list(map(int, f.readline().split(',')))
            ports[i, 3, 0].assign(vals[0])
            ports[i, 3, 1].assign(vals[1])
            ports[i, 3, 2].assign(vals[2])

            if char_ports == 'x':
                ports[i, 3, 1].assign(ports[i, 3, 1] + 1)
                ports[i, 3, 2].assign(ports[i, 3, 2] + 1)
            elif char_ports == 'y':
                ports[i, 3, 0].assign(ports[i, 3, 0] + 1)
                ports[i, 3, 2].assign(ports[i, 3, 2] + 1)
            elif char_ports == 'z':
                ports[i, 3, 0].assign(ports[i, 3, 0] + 1)
                ports[i, 3, 1].assign(ports[i, 3, 1] + 1)

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
    with tf.device(device):
        theta = tf.constant(theta, dtype=dtype)
        phi = tf.constant(phi, dtype=dtype)
        pol = tf.constant(pol, dtype=dtype)

        theta_mirror = tf.constant(theta_mirror, dtype=dtype)
        phi_mirror = tf.constant(phi_mirror, dtype=dtype)

    # far field data uses coordinate system centered in the middle of the simulation space
    # currently only used in the PML version, not pbc version - pbc not present here anyway
    ic=int((x_size-1)/2.0) 
    jc=int((y_size-1)/2.0)
    kc=int((z_size-1)/2.0)

    # These are used in a number of places - mostly far fields
    if (type_sim==0):
        xlow=nxPML_1+buffer
        xhigh=x_size-nxPML_2-buffer
        ylow=1
        yhigh=y_size-1
        zlow=1
        zhigh=z_size-1
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
    WHx=plane_wave_amp*(tf.sin(pol)*tf.sin(phi)+tf.cos(pol)*tf.cos(theta)*tf.cos(phi))*(ep_0/mu_0)**(0.5)
    WHy=plane_wave_amp*(-1*tf.sin(pol)*tf.cos(phi)+tf.cos(pol)*tf.cos(theta)*tf.sin(phi))*(ep_0/mu_0)**0.5
    WHz=plane_wave_amp*(-1*tf.cos(pol)*tf.sin(theta))*(ep_0/mu_0)**0.5
    WEx=plane_wave_amp*(tf.cos(pol)*tf.sin(phi)-tf.sin(pol)*tf.cos(theta)*tf.cos(phi))
    WEy=plane_wave_amp*(-1*tf.cos(pol)*tf.cos(phi)-tf.sin(pol)*tf.cos(theta)*tf.sin(phi))
    WEz=plane_wave_amp*(tf.sin(pol)*tf.sin(theta))

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
    with tf.device(device):
        min_steps=2*t_spread/del_t+tf.sqrt(((tf.sin(theta)*tf.cos(phi)*(x_size)*del_x)**2 \
        +(tf.sin(theta)*tf.sin(phi)*(y_size)*del_y)**2+(tf.cos(theta)*(z_size)*del_z)**2))/(c*del_t)

##########################################################
# FDTD INITIALIZATION
##########################################################
with GPUProfiler("FDTD Initialization", device):

    # Initialize all arrays ON GPU with proper dtype
    with tf.device(device):

        # Main E,H Fields
        Ex = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        Ey = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        Ez = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        Hx = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        Hy = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        Hz = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))

        # Material cells
        relative_ep_x_cell =  tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        relative_ep_y_cell =  tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        relative_ep_z_cell =  tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sigma_x_cell =  tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sigma_y_cell =  tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sigma_z_cell =  tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))

        # Grid points
        relative_ep_x =  tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        relative_ep_y =  tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        relative_ep_z =  tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sigma_x =  tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sigma_y =  tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sigma_z =  tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))

        # Special field and sheet list arrays
        Ex_special = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        Ey_special = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        Ez_special = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))

        x_sheet_list = tf.Variable(tf.zeros((x_size - 1,), dtype=dtype))
        y_sheet_list = tf.Variable(tf.zeros((y_size - 1,), dtype=dtype))
        z_sheet_list = tf.Variable(tf.zeros((z_size - 1,), dtype=dtype))

        # Sheet materials (x-normal)
        sheet_ep_x_cell_x = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_ep_y_cell_x = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_ep_z_cell_x = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_x_cell_x = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_y_cell_x = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_z_cell_x = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))

        # Sheet grid points (x-normal)
        sheet_ep_x_x = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_ep_y_x = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_ep_z_x = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_x_x = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_y_x = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_z_x = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))

        # Sheet materials (y-normal)
        sheet_ep_x_cell_y = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_ep_y_cell_y = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_ep_z_cell_y = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_x_cell_y = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_y_cell_y = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_z_cell_y = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))

        # Sheet grid points y-normal)
        sheet_ep_x_y = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_ep_y_y = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_ep_z_y = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_x_y = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_y_y = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_z_y = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))

        # Sheet materials (z-normal)
        sheet_ep_x_cell_z = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_ep_y_cell_z = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_ep_z_cell_z = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_x_cell_z = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_y_cell_z = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_z_cell_z = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))

        # Sheet grid points (z-normal)
        sheet_ep_x_z = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_ep_y_z = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_ep_z_z = tf.Variable(tf.ones((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_x_z = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_y_z = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        sheet_sig_z_z = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))

        # For permeability (scalars for now, but will modify at some point)
        da = 1.0
        db = del_t / mu_0
        
        # Field aux for permittivity and electrical conductivity
        gax = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        gbx = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        gay = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        gby = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        gaz = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        gbz = tf.Variable(tf.zeros((x_size, y_size, z_size), dtype=dtype))
        
        # Initialize psi arrays
        psi_Ezx_1 = tf.Variable(tf.zeros((nxPML_1, y_size, z_size), dtype=dtype))
        psi_Ezx_2 = tf.Variable(tf.zeros((nxPML_2, y_size, z_size), dtype=dtype))
        psi_Hyx_1 = tf.Variable(tf.zeros((nxPML_1-1, y_size, z_size), dtype=dtype))
        psi_Hyx_2 = tf.Variable(tf.zeros((nxPML_2-1, y_size, z_size), dtype=dtype))
        psi_Ezy_1 = tf.Variable(tf.zeros((x_size, nyPML_1, z_size), dtype=dtype))
        psi_Ezy_2 = tf.Variable(tf.zeros((x_size, nyPML_2, z_size), dtype=dtype))
        psi_Hxy_1 = tf.Variable(tf.zeros((x_size, nyPML_1-1, z_size), dtype=dtype))
        psi_Hxy_2 = tf.Variable(tf.zeros((x_size, nyPML_2-1, z_size), dtype=dtype))
        psi_Hxz_1 = tf.Variable(tf.zeros((x_size, y_size-1, nzPML_1-1), dtype=dtype))
        psi_Hxz_2 = tf.Variable(tf.zeros((x_size, y_size-1, nzPML_2-1), dtype=dtype))
        psi_Hyz_1 = tf.Variable(tf.zeros((x_size-1, y_size, nzPML_1-1), dtype=dtype))
        psi_Hyz_2 = tf.Variable(tf.zeros((x_size-1, y_size, nzPML_2-1), dtype=dtype))
        psi_Exz_1 = tf.Variable(tf.zeros((x_size-1, y_size, nzPML_1), dtype=dtype))
        psi_Exz_2 = tf.Variable(tf.zeros((x_size-1, y_size, nzPML_2), dtype=dtype))
        psi_Eyz_1 = tf.Variable(tf.zeros((x_size, y_size-1, nzPML_1), dtype=dtype))
        psi_Eyz_2 = tf.Variable(tf.zeros((x_size, y_size-1, nzPML_2), dtype=dtype))
        psi_Hzx_1 = tf.Variable(tf.zeros((nxPML_1-1, y_size-1, z_size), dtype=dtype))
        psi_Eyx_1 = tf.Variable(tf.zeros((nxPML_1, y_size-1, z_size), dtype=dtype))
        psi_Hzx_2 = tf.Variable(tf.zeros((nxPML_2-1, y_size-1, z_size), dtype=dtype))
        psi_Eyx_2 = tf.Variable(tf.zeros((nxPML_2, y_size-1, z_size), dtype=dtype))
        psi_Hzy_1 = tf.Variable(tf.zeros((x_size-1, nyPML_1-1, z_size), dtype=dtype))
        psi_Exy_1 = tf.Variable(tf.zeros((x_size-1, nyPML_1, z_size), dtype=dtype))
        psi_Hzy_2 = tf.Variable(tf.zeros((x_size-1, nyPML_2-1, z_size), dtype=dtype))
        psi_Exy_2 = tf.Variable(tf.zeros((x_size-1, nyPML_2, z_size), dtype=dtype))
        
        # Initialize PML parameter arrays
        alphae_x_PML_1 = tf.Variable(tf.zeros(nxPML_1, dtype=dtype))
        sige_x_PML_1 = tf.Variable(tf.zeros(nxPML_1, dtype=dtype))
        kappae_x_PML_1 = tf.Variable(tf.zeros(nxPML_1, dtype=dtype))
        alphah_x_PML_1 = tf.Variable(tf.zeros(nxPML_1-1, dtype=dtype))
        sigh_x_PML_1 = tf.Variable(tf.zeros(nxPML_1-1, dtype=dtype))
        kappah_x_PML_1 = tf.Variable(tf.zeros(nxPML_1-1, dtype=dtype))
        
        alphae_x_PML_2 = tf.Variable(tf.zeros(nxPML_2, dtype=dtype))
        sige_x_PML_2 = tf.Variable(tf.zeros(nxPML_2, dtype=dtype))
        kappae_x_PML_2 = tf.Variable(tf.zeros(nxPML_2, dtype=dtype))
        alphah_x_PML_2 = tf.Variable(tf.zeros(nxPML_2-1, dtype=dtype))
        sigh_x_PML_2 = tf.Variable(tf.zeros(nxPML_2-1, dtype=dtype))
        kappah_x_PML_2 = tf.Variable(tf.zeros(nxPML_2-1, dtype=dtype))
        
        alphae_y_PML_1 = tf.Variable(tf.zeros(nyPML_1, dtype=dtype))
        sige_y_PML_1 = tf.Variable(tf.zeros(nyPML_1, dtype=dtype))
        kappae_y_PML_1 = tf.Variable(tf.zeros(nyPML_1, dtype=dtype))
        alphah_y_PML_1 = tf.Variable(tf.zeros(nyPML_1-1, dtype=dtype))
        sigh_y_PML_1 = tf.Variable(tf.zeros(nyPML_1-1, dtype=dtype))
        kappah_y_PML_1 = tf.Variable(tf.zeros(nyPML_1-1, dtype=dtype))
        
        alphae_y_PML_2 = tf.Variable(tf.zeros(nyPML_2, dtype=dtype))
        sige_y_PML_2 = tf.Variable(tf.zeros(nyPML_2, dtype=dtype))
        kappae_y_PML_2 = tf.Variable(tf.zeros(nyPML_2, dtype=dtype))
        alphah_y_PML_2 = tf.Variable(tf.zeros(nyPML_2-1, dtype=dtype))
        sigh_y_PML_2 = tf.Variable(tf.zeros(nyPML_2-1, dtype=dtype))
        kappah_y_PML_2 = tf.Variable(tf.zeros(nyPML_2-1, dtype=dtype))
        
        alphae_z_PML_1 = tf.Variable(tf.zeros(nzPML_1, dtype=dtype))
        sige_z_PML_1 = tf.Variable(tf.zeros(nzPML_1, dtype=dtype))
        kappae_z_PML_1 = tf.Variable(tf.zeros(nzPML_1, dtype=dtype))
        alphah_z_PML_1 = tf.Variable(tf.zeros(nzPML_1-1, dtype=dtype))
        sigh_z_PML_1 = tf.Variable(tf.zeros(nzPML_1-1, dtype=dtype))
        kappah_z_PML_1 = tf.Variable(tf.zeros(nzPML_1-1, dtype=dtype))
        
        alphae_z_PML_2 = tf.Variable(tf.zeros(nzPML_2, dtype=dtype))
        sige_z_PML_2 = tf.Variable(tf.zeros(nzPML_2, dtype=dtype))
        kappae_z_PML_2 = tf.Variable(tf.zeros(nzPML_2, dtype=dtype))
        alphah_z_PML_2 = tf.Variable(tf.zeros(nzPML_2-1, dtype=dtype))
        sigh_z_PML_2 = tf.Variable(tf.zeros(nzPML_2-1, dtype=dtype))
        kappah_z_PML_2 = tf.Variable(tf.zeros(nzPML_2-1, dtype=dtype))
        
        # Initialize b and c coefficient arrays
        be_x_1 = tf.Variable(tf.zeros((nxPML_1,y_size,z_size), dtype=dtype))
        ce_x_1 = tf.Variable(tf.zeros((nxPML_1,y_size,z_size), dtype=dtype))
        bh_x_1 = tf.Variable(tf.zeros((nxPML_1-1,y_size,z_size), dtype=dtype))
        ch_x_1 = tf.Variable(tf.zeros((nxPML_1-1,y_size,z_size), dtype=dtype))

        be_x_2 = tf.Variable(tf.zeros((nxPML_2,y_size,z_size), dtype=dtype))
        ce_x_2 = tf.Variable(tf.zeros((nxPML_2,y_size,z_size), dtype=dtype))
        bh_x_2 = tf.Variable(tf.zeros((nxPML_2-1,y_size,z_size), dtype=dtype))
        ch_x_2 = tf.Variable(tf.zeros((nxPML_2-1,y_size,z_size), dtype=dtype))

        be_y_1 = tf.Variable(tf.zeros((x_size,nyPML_1,z_size), dtype=dtype))
        ce_y_1 = tf.Variable(tf.zeros((x_size,nyPML_1,z_size), dtype=dtype))
        bh_y_1 = tf.Variable(tf.zeros((x_size,nyPML_1-1,z_size), dtype=dtype))
        ch_y_1 = tf.Variable(tf.zeros((x_size,nyPML_1-1,z_size), dtype=dtype))

        be_y_2 = tf.Variable(tf.zeros((x_size,nyPML_2,z_size), dtype=dtype))
        ce_y_2 = tf.Variable(tf.zeros((x_size,nyPML_2,z_size), dtype=dtype))
        bh_y_2 = tf.Variable(tf.zeros((x_size,nyPML_2-1,z_size), dtype=dtype))
        ch_y_2 = tf.Variable(tf.zeros((x_size,nyPML_2-1,z_size), dtype=dtype))

        be_z_1 = tf.Variable(tf.zeros((x_size,y_size,nzPML_1), dtype=dtype))
        ce_z_1 = tf.Variable(tf.zeros((x_size,y_size,nzPML_1), dtype=dtype))
        bh_z_1 = tf.Variable(tf.zeros((x_size,y_size,nzPML_1-1), dtype=dtype))
        ch_z_1 = tf.Variable(tf.zeros((x_size,y_size,nzPML_1-1), dtype=dtype))

        be_z_2 = tf.Variable(tf.zeros((x_size,y_size,nzPML_2), dtype=dtype))
        ce_z_2 = tf.Variable(tf.zeros((x_size,y_size,nzPML_2), dtype=dtype))
        bh_z_2 = tf.Variable(tf.zeros((x_size,y_size,nzPML_2-1), dtype=dtype))
        ch_z_2 = tf.Variable(tf.zeros((x_size,y_size,nzPML_2-1), dtype=dtype))

        # Denominators for update equations
        den_ex = tf.Variable(tf.fill((x_size-1,), 1/del_x))
        den_ey = tf.Variable(tf.fill((y_size-1,), 1/del_y))
        den_ez = tf.Variable(tf.fill((z_size-1,), 1/del_z))
        den_hx = tf.Variable(tf.fill((x_size-1,), 1/del_x))
        den_hy = tf.Variable(tf.fill((y_size-1,), 1/del_y))
        den_hz = tf.Variable(tf.fill((z_size-1,), 1/del_z))

        # Voltage and Current Aux Arrays if needed
        if num_ports>0:
            Voltage = tf.Variable(tf.zeros((num_ports, port_array_size), dtype=dtype))
            Current = tf.Variable(tf.zeros((num_ports, port_array_size), dtype=dtype))
            Voltage_out = tf.Variable(tf.zeros((num_ports, time_steps), dtype=dtype))
            Current_out = tf.Variable(tf.zeros((num_ports, time_steps), dtype=dtype))

        # Video Arrays (3D Fields) if needed
        if video_on=='yes':
            Ex_video = tf.Variable(tf.zeros((vid_size1, vid_size2, time_steps), dtype=dtype))
            Hx_video = tf.Variable(tf.zeros((vid_size1, vid_size2, time_steps), dtype=dtype))
            Ey_video = tf.Variable(tf.zeros((vid_size1, vid_size2, time_steps), dtype=dtype))
            Hy_video = tf.Variable(tf.zeros((vid_size1, vid_size2, time_steps), dtype=dtype))
            Ez_video = tf.Variable(tf.zeros((vid_size1, vid_size2, time_steps), dtype=dtype))
            Hz_video = tf.Variable(tf.zeros((vid_size1, vid_size2, time_steps), dtype=dtype))

        # Incident and Output Field Arrays (1D arrays) - various uses below
        incident = tf.Variable(tf.zeros((time_steps,), dtype=dtype))
        E_reflected = tf.Variable(tf.zeros((time_steps,), dtype=dtype))
        E_transmitted = tf.Variable(tf.zeros((time_steps,), dtype=dtype))

        # Far Field Post Processing Variables (2D Boundary Surfaces) if needed
        if num_far_field_angles > 0:
            # x-low/x-high faces (depend on y_size, z_size)
            My_xlow = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))
            Mz_xlow = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))
            Jy_xlow = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))
            Jz_xlow = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))
            My_xlow_oldt = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))
            Mz_xlow_oldt = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))
            Jy_xlow_oldt = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))
            Jz_xlow_oldt = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))

            My_xhigh = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))
            Mz_xhigh = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))
            Jy_xhigh = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))
            Jz_xhigh = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))
            My_xhigh_oldt = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))
            Mz_xhigh_oldt = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))
            Jy_xhigh_oldt = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))
            Jz_xhigh_oldt = tf.Variable(tf.zeros((y_size, z_size), dtype=dtype))

            # y-low/y-high faces (depend on x_size, z_size)
            Mx_ylow = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))
            Mz_ylow = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))
            Jx_ylow = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))
            Jz_ylow = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))
            Mx_ylow_oldt = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))
            Mz_ylow_oldt = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))
            Jx_ylow_oldt = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))
            Jz_ylow_oldt = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))

            Mx_yhigh = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))
            Mz_yhigh = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))
            Jx_yhigh = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))
            Jz_yhigh = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))
            Mx_yhigh_oldt = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))
            Mz_yhigh_oldt = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))
            Jx_yhigh_oldt = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))
            Jz_yhigh_oldt = tf.Variable(tf.zeros((x_size, z_size), dtype=dtype))

            # z-low/z-high faces (depend on x_size, y_size)
            Mx_zlow = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))
            My_zlow = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))
            Jx_zlow = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))
            Jy_zlow = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))
            Mx_zlow_oldt = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))
            My_zlow_oldt = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))
            Jx_zlow_oldt = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))
            Jy_zlow_oldt = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))

            Mx_zhigh = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))
            My_zhigh = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))
            Jx_zhigh = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))
            Jy_zhigh = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))
            Mx_zhigh_oldt = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))
            My_zhigh_oldt = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))
            Jx_zhigh_oldt = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))
            Jy_zhigh_oldt = tf.Variable(tf.zeros((x_size, y_size), dtype=dtype))

            # Far Field Summation and Output Variables (W, U, E)
            W_shape = (num_far_field_angles, len_far_field_arrays)
            E_out_shape = (num_far_field_angles, time_steps)

            Wx = tf.Variable(tf.zeros(W_shape, dtype=dtype))
            Wy = tf.Variable(tf.zeros(W_shape, dtype=dtype))
            Wz = tf.Variable(tf.zeros(W_shape, dtype=dtype))
            Ux = tf.Variable(tf.zeros(W_shape, dtype=dtype))
            Uy = tf.Variable(tf.zeros(W_shape, dtype=dtype))
            Uz = tf.Variable(tf.zeros(W_shape, dtype=dtype))
            W_theta = tf.Variable(tf.zeros(W_shape, dtype=dtype))
            W_phi = tf.Variable(tf.zeros(W_shape, dtype=dtype))
            U_theta = tf.Variable(tf.zeros(W_shape, dtype=dtype))
            U_phi = tf.Variable(tf.zeros(W_shape, dtype=dtype))

            # Combined outputs
            E_theta = tf.Variable(tf.zeros(W_shape, dtype=dtype))
            E_phi = tf.Variable(tf.zeros(W_shape, dtype=dtype))
            E_theta_out = tf.Variable(tf.zeros(E_out_shape, dtype=dtype))
            E_phi_out = tf.Variable(tf.zeros(E_out_shape, dtype=dtype))
            data_out_time = tf.Variable(tf.zeros((num_far_field_angles,), dtype=dtype))
        
        # Initialize arrays for saving data (on GPU for now)
        input_data = tf.Variable(tf.zeros(time_steps, dtype=dtype))
        
        # Pre-compute pulse time values (on GPU)
        time_values = tf.range(time_steps, dtype=dtype) * del_t

##########################################################
# PML SETUP
##########################################################
with GPUProfiler("PML Calculations", device):
    # Setup the parameters - probably a faster/better way to do this long term but copied and slightly modified from fortran
    # The PML are only 10 in size so it's pretty fast regardless of methods and cpu/gpu settings    
    for i in range(nxPML_1):
        sige_x_PML_1[i].assign(sig_x_max * ( (nxPML_1 - (i+1)) / (nxPML_1 - 1.0) )**m)
        alphae_x_PML_1[i].assign(alpha_x_max*(((i+1)-1)/(nxPML_1-1.0))**ma)
        kappae_x_PML_1[i].assign(1.0+(kappa_x_max-1.0)*((nxPML_1 - (i+1)) / (nxPML_1 - 1.0))**m)
        be_x_1[i,:,:].assign(tf.fill(be_x_1[i,:,:].shape, math.exp(-(sige_x_PML_1[i] / kappae_x_PML_1[i] + alphae_x_PML_1[i])*del_t/ep_0)))
        if ((sige_x_PML_1[i] == 0.0) and (alphae_x_PML_1[i] == 0.0) and (i == nxPML_1-1)):
            ce_x_1[i,:,:].assign(tf.zeros(ce_x_1[i,:,:].shape))
        else:
            ce_x_1[i,:,:].assign(tf.fill(ce_x_1[i,:,:].shape, sige_x_PML_1[i]*(be_x_1[i,0,0]-1.0)/(sige_x_PML_1[i]+kappae_x_PML_1[i]*alphae_x_PML_1[i]) / kappae_x_PML_1[i]))

    for i in range(nxPML_1-1):
        sigh_x_PML_1[i].assign(sig_x_max * ( (nxPML_1 - (i+1) - 0.5)/(nxPML_1-1.0))**m)
        alphah_x_PML_1[i].assign(alpha_x_max*(((i+1)-0.5)/(nxPML_1-1.0))**ma)
        kappah_x_PML_1[i].assign(1.0+(kappa_x_max-1.0)*((nxPML_1 - (i+1) - 0.5) / (nxPML_1 - 1.0))**m)
        bh_x_1[i,:,:].assign(tf.fill(bh_x_1[i,:,:].shape, math.exp(-(sigh_x_PML_1[i] / kappah_x_PML_1[i] + alphah_x_PML_1[i])*del_t/ep_0)))
        ch_x_1[i,:,:].assign(tf.fill(ch_x_1[i,:,:].shape, sigh_x_PML_1[i]*(bh_x_1[i,0,0]-1.0)/(sigh_x_PML_1[i]+kappah_x_PML_1[i]*alphah_x_PML_1[i])/kappah_x_PML_1[i]))

    for i in range(nxPML_2):
        sige_x_PML_2[i].assign(sig_x_max * ( (nxPML_2 - (i+1)) / (nxPML_2 - 1.0) )**m)
        alphae_x_PML_2[i].assign(alpha_x_max*(((i+1)-1)/(nxPML_2-1.0))**ma)
        kappae_x_PML_2[i].assign(1.0+(kappa_x_max-1.0)*((nxPML_2 - (i+1)) / (nxPML_2 - 1.0))**m)
        be_x_2[i,:,:].assign(tf.fill(be_x_2[i,:,:].shape, math.exp(-(sige_x_PML_2[i] / kappae_x_PML_2[i] +alphae_x_PML_2[i])*del_t/ep_0)))
        if ((sige_x_PML_2[i] == 0.0) and (alphae_x_PML_2[i] == 0.0) and (i == nxPML_2-1)):
            ce_x_2[i,:,:].assign(tf.zeros(ce_x_2[i,:,:].shape))
        else:
            ce_x_2[i,:,:].assign(tf.fill(ce_x_2[i,:,:].shape, sige_x_PML_2[i]*(be_x_2[i,0,0]-1.0)/(sige_x_PML_2[i]+kappae_x_PML_2[i]*alphae_x_PML_2[i]) / kappae_x_PML_2[i]))

    for i in range(nxPML_2-1):
        sigh_x_PML_2[i].assign(sig_x_max * ( (nxPML_2 - (i+1) - 0.5)/(nxPML_2-1.0))**m)
        alphah_x_PML_2[i].assign(alpha_x_max*(((i+1)-0.5)/(nxPML_2-1.0))**ma)
        kappah_x_PML_2[i].assign(1.0+(kappa_x_max-1.0)*((nxPML_2 - (i+1) - 0.5) / (nxPML_2 - 1.0))**m)
        bh_x_2[i,:,:].assign(tf.fill(bh_x_2[i,:,:].shape, math.exp(-(sigh_x_PML_2[i] / kappah_x_PML_2[i] +alphah_x_PML_2[i])*del_t/ep_0)))
        ch_x_2[i,:,:].assign(tf.fill(ch_x_2[i,:,:].shape, sigh_x_PML_2[i]*(bh_x_2[i,0,0]-1.0)/(sigh_x_PML_2[i]+kappah_x_PML_2[i]*alphah_x_PML_2[i])/kappah_x_PML_2[i]))

    for j in range(nyPML_1):
        sige_y_PML_1[j].assign(sig_y_max * ( (nyPML_1 - (j+1)) / (nyPML_1 - 1.0) )**m)
        alphae_y_PML_1[j].assign(alpha_y_max*(((j+1)-1)/(nyPML_1-1.0))**ma)
        kappae_y_PML_1[j].assign(1.0+(kappa_y_max-1.0)*((nyPML_1 - (j+1)) / (nyPML_1 - 1.0))**m)
        be_y_1[:,j,:].assign(tf.fill(be_y_1[:,j,:].shape, math.exp(-(sige_y_PML_1[j] / kappae_y_PML_1[j] +alphae_y_PML_1[j])*del_t/ep_0)))
        if ((sige_y_PML_1[j] == 0.0) and (alphae_y_PML_1[j] == 0.0) and (j == nyPML_1-1)):
            ce_y_1[:,j,:].assign(tf.zeros(ce_y_1[:,j,:].shape))
        else:
            ce_y_1[:,j,:].assign(tf.fill(ce_y_1[:,j,:].shape, sige_y_PML_1[j]*(be_y_1[0,j,0]-1.0)/(sige_y_PML_1[j]+kappae_y_PML_1[j]*alphae_y_PML_1[j]) / kappae_y_PML_1[j]))

    for j in range(nyPML_1-1):
        sigh_y_PML_1[j].assign(sig_y_max * ( (nyPML_1 - (j+1) - 0.5)/(nyPML_1-1.0))**m)
        alphah_y_PML_1[j].assign(alpha_y_max*(((j+1)-0.5)/(nyPML_1-1.0))**ma)
        kappah_y_PML_1[j].assign(1.0+(kappa_y_max-1.0)*((nyPML_1 - (j+1) - 0.5) / (nyPML_1 - 1.0))**m)
        bh_y_1[:,j,:].assign(tf.fill(bh_y_1[:,j,:].shape, math.exp(-(sigh_y_PML_1[j] / kappah_y_PML_1[j] +alphah_y_PML_1[j])*del_t/ep_0)))
        ch_y_1[:,j,:].assign(tf.fill(ch_y_1[:,j,:].shape, sigh_y_PML_1[j]*(bh_y_1[0,j,0]-1.0)/(sigh_y_PML_1[j]+kappah_y_PML_1[j]*alphah_y_PML_1[j])/ kappah_y_PML_1[j]))

    for j in range(nyPML_2):
        sige_y_PML_2[j].assign(sig_y_max * ( (nyPML_2 - (j+1)) / (nyPML_2 - 1.0) )**m)
        alphae_y_PML_2[j].assign(alpha_y_max*(((j+1)-1)/(nyPML_2-1.0))**ma)
        kappae_y_PML_2[j].assign(1.0+(kappa_y_max-1.0)*((nyPML_2 - (j+1)) / (nyPML_2 - 1.0))**m)
        be_y_2[:,j,:].assign(tf.fill(be_y_2[:,j,:].shape, math.exp(-(sige_y_PML_2[j] / kappae_y_PML_2[j] +alphae_y_PML_2[j])*del_t/ep_0)))
        if ((sige_y_PML_2[j] == 0.0) and(alphae_y_PML_2[j] == 0.0) and(j == nyPML_2-1)):
            ce_y_2[:,j,:].assign(tf.zeros(ce_y_2[:,j,:].shape))
        else:
            ce_y_2[:,j,:].assign(tf.fill(ce_y_2[:,j,:].shape, sige_y_PML_2[j]*(be_y_2[0,j,0]-1.0)/(sige_y_PML_2[j]+kappae_y_PML_2[j]*alphae_y_PML_2[j])/kappae_y_PML_2[j]))

    for j in range(nyPML_2-1):
        sigh_y_PML_2[j].assign(sig_y_max * ( (nyPML_2 - (j+1) - 0.5)/(nyPML_2-1.0))**m)
        alphah_y_PML_2[j].assign(alpha_y_max*(((j+1)-0.5)/(nyPML_2-1.0))**ma)
        kappah_y_PML_2[j].assign(1.0+(kappa_y_max-1.0)*((nyPML_2 - (j+1) - 0.5) / (nyPML_2 - 1.0))**m)
        bh_y_2[:,j,:].assign(tf.fill(bh_y_2[:,j,:].shape, math.exp(-(sigh_y_PML_2[j] / kappah_y_PML_2[j] +alphah_y_PML_2[j])*del_t/ep_0)))
        ch_y_2[:,j,:].assign(tf.fill(ch_y_2[:,j,:].shape, sigh_y_PML_2[j]*(bh_y_2[0,j,0]-1.0)/(sigh_y_PML_2[j]+kappah_y_PML_2[j]*alphah_y_PML_2[j])/kappah_y_PML_2[j]))

    for k in range(nzPML_1):
        sige_z_PML_1[k].assign(sig_z_max * ( (nzPML_1 - (k+1)) / (nzPML_1 - 1.0) )**m)
        alphae_z_PML_1[k].assign(alpha_z_max*(((k+1)-1)/(nzPML_1-1.0))**ma)
        kappae_z_PML_1[k].assign(1.0+(kappa_z_max-1.0)*((nzPML_1 - (k+1)) / (nzPML_1 - 1.0))**m)
        be_z_1[:,:,k].assign(tf.fill(be_z_1[:,:,k].shape, math.exp(-(sige_z_PML_1[k] / kappae_z_PML_1[k] +alphae_z_PML_1[k])*del_t/ep_0)))
        if ((sige_z_PML_1[k] == 0.0) and (alphae_z_PML_1[k] == 0.0) and(k == nzPML_1-1)):
            ce_z_1[:,:,k].assign(tf.zeros(ce_z_1[:,:,k].shape))
        else:
            ce_z_1[:,:,k].assign(tf.fill(ce_z_1[:,:,k].shape, sige_z_PML_1[k]*(be_z_1[0,0,k]-1.0)/(sige_z_PML_1[k]+kappae_z_PML_1[k]*alphae_z_PML_1[k]) / kappae_z_PML_1[k]))

    for k in range(nzPML_1-1):
        sigh_z_PML_1[k].assign(sig_z_max * ( (nzPML_1 - (k+1) - 0.5)/(nzPML_1-1.0))**m)
        alphah_z_PML_1[k].assign(alpha_z_max*(((k+1)-0.5)/(nzPML_1-1.0))**ma)
        kappah_z_PML_1[k].assign(1.0+(kappa_z_max-1.0)*((nzPML_1 - (k+1) - 0.5) / (nzPML_1 - 1.0))**m)
        bh_z_1[:,:,k].assign(tf.fill(bh_z_1[:,:,k].shape, math.exp(-(sigh_z_PML_1[k] / kappah_z_PML_1[k] +alphah_z_PML_1[k])*del_t/ep_0)))
        ch_z_1[:,:,k].assign(tf.fill(ch_z_1[:,:,k].shape, sigh_z_PML_1[k]*(bh_z_1[0,0,k]-1.0)/(sigh_z_PML_1[k]+kappah_z_PML_1[k]*alphah_z_PML_1[k])/ kappah_z_PML_1[k]))

    for k in range(nzPML_2):
        sige_z_PML_2[k].assign(sig_z_max * ( (nzPML_2 - (k+1)) / (nzPML_2 - 1.0) )**m)
        alphae_z_PML_2[k].assign(alpha_z_max*(((k+1)-1)/(nzPML_2-1.0))**ma)
        kappae_z_PML_2[k].assign(1.0+(kappa_z_max-1.0)*((nzPML_2 - (k+1)) / (nzPML_2 - 1.0))**m)
        be_z_2[:,:,k].assign(tf.fill(be_z_2[:,:,k].shape, math.exp(-(sige_z_PML_2[k] / kappae_z_PML_2[k] +alphae_z_PML_2[k])*del_t/ep_0)))
        if ((sige_z_PML_2[k] == 0.0) and (alphae_z_PML_2[k] == 0.0) and(k == nzPML_2-1)):
            ce_z_2[:,:,k].assign(tf.zeros(ce_z_2[:,:,k].shape))
        else:
            ce_z_2[:,:,k].assign(tf.fill(ce_z_2[:,:,k].shape, sige_z_PML_2[k]*(be_z_2[0,0,k]-1.0)/(sige_z_PML_2[k]+kappae_z_PML_2[k]*alphae_z_PML_2[k])/ kappae_z_PML_2[k]))

    for k in range(nzPML_2-1):
        sigh_z_PML_2[k].assign(sig_z_max * ( (nzPML_2 - (k+1) - 0.5)/(nzPML_2-1.0))**m)
        alphah_z_PML_2[k].assign(alpha_z_max*(((k+1)-0.5)/(nzPML_2-1.0))**ma)
        kappah_z_PML_2[k].assign(1.0+(kappa_z_max-1.0)*((nzPML_2 - (k+1) - 0.5) / (nzPML_2 - 1.0))**m)
        bh_z_2[:,:,k].assign(tf.fill(bh_z_2[:,:,k].shape, math.exp(-(sigh_z_PML_2[k] / kappah_z_PML_2[k] + alphah_z_PML_2[k])*del_t/ep_0)))
        ch_z_2[:,:,k].assign(tf.fill(ch_z_2[:,:,k].shape, sigh_z_PML_2[k]*(bh_z_2[0,0,k]-1.0)/(sigh_z_PML_2[k]+kappah_z_PML_2[k]*alphah_z_PML_2[k]) / kappah_z_PML_2[k]))

##########################################################
# Denominator SETUP
##########################################################
with GPUProfiler("Denominator Calculations", device):
    
    den_hx[0 : nxPML_1 - 1].assign(1.0 / (kappah_x_PML_1[:nxPML_1 - 1] * del_x))
    den_hx[x_size - nxPML_2 : x_size - 1].assign(1.0 / (tf.reverse(kappah_x_PML_2[:nxPML_2 - 1], axis=[0]) * del_x))

    if type_sim==1:
        den_hy[0 : nyPML_1 - 1].assign(1.0 / (kappah_y_PML_1[:nyPML_1 - 1] * del_y))
        den_hy[y_size - nyPML_2 : y_size - 1].assign(1.0 / (tf.reverse(kappah_y_PML_2[:nyPML_2 - 1], axis=[0]) * del_y))

        den_hz[0 : nzPML_1 - 1].assign(1.0 / (kappah_z_PML_1[:nzPML_1 - 1] * del_z))
        den_hz[z_size - nzPML_2 : z_size - 1].assign(1.0 / (tf.reverse(kappah_z_PML_2[:nzPML_2 - 1], axis=[0]) * del_z))

    den_ex[0 : nxPML_1].assign(1.0 / (kappae_x_PML_1[:nxPML_1] * del_x))
    den_ex[x_size - nxPML_2 : x_size - 1].assign(1.0 / (tf.reverse(kappae_x_PML_2[1:nxPML_2], axis=[0]) * del_x))

    if type_sim==1:
        den_ey[0 : nyPML_1].assign(1.0 / (kappae_y_PML_1[:nyPML_1] * del_y))
        den_ey[y_size - nyPML_2 : y_size - 1].assign(1.0 / (tf.reverse(kappae_y_PML_2[1:nyPML_2], axis=[0]) * del_y))

        den_ez[0 : nzPML_1].assign(1.0 / (kappae_z_PML_1[:nzPML_1] * del_z))
        den_ez[z_size - nzPML_2 : z_size - 1].assign(1.0 / (tf.reverse(kappae_z_PML_2[1:nzPML_2], axis=[0]) * del_z))

##########################################################
# Geometry SETUP
##########################################################
with GPUProfiler("Geometry Calculations", device):
    """
    This section is long and messy
    Key items we keep when done with this section are:
    ga,gb variables (x,y,z)
    masks for main fdtd loop for sheets
    ep,sig matrices for sheets (xx,yy,zz) for special subcell update in E,H
    """

    # First we need to put all blocks,spheres,cylinders, and sheets into material cells
    ep_sig_tensors = [
        relative_ep_x_cell, relative_ep_y_cell, relative_ep_z_cell,
        sigma_x_cell, sigma_y_cell, sigma_z_cell
    ]

    # Grid indices (using reshape instead of view)
    gi = tf.range(relative_ep_x_cell.shape[0])[:, tf.newaxis, tf.newaxis]
    gj = tf.range(relative_ep_x_cell.shape[1])[tf.newaxis, :, tf.newaxis]
    gk = tf.range(relative_ep_x_cell.shape[2])[tf.newaxis, tf.newaxis, :]

    ii, jj, kk = 0, 0, 0
    for counter in range(num_objects):
        obj_type = object_type[counter].numpy().decode('utf-8').strip() if hasattr(object_type[counter], 'numpy') else object_type[counter].strip()
        
        # --- BLOCK ---
        if obj_type == 'block':
            # m_id logic: find index where material property matches
            m_id = tf.where(materials_properties[:, 0, 0] == blocks[ii, 0, 0])[0][0]
            
            i_s, i_e = int(blocks[ii, 1, 0]), int(blocks[ii, 1, 0] + blocks[ii, 2, 0])
            j_s, j_e = int(blocks[ii, 1, 1] + pbc_shift), int(blocks[ii, 1, 1] + blocks[ii, 2, 1] + pbc_shift)
            k_s, k_e = int(blocks[ii, 1, 2] + pbc_shift), int(blocks[ii, 1, 2] + blocks[ii, 2, 2] + pbc_shift)
            
            vals = tf.concat([materials_properties[m_id, 1, :], materials_properties[m_id, 2, :]], axis=0)
            
            for idx, tensor in enumerate(ep_sig_tensors):
                # In TF, we use assign for in-place updates on Variables
                update_shape = [i_e - i_s, j_e - j_s, k_e - k_s]
                tensor[i_s:i_e, j_s:j_e, k_s:k_e].assign(tf.fill(update_shape, vals[idx]))
            
            ii += 1

        # --- SPHERE ---
        elif obj_type == 'sphere':
            m_id = tf.where(materials_properties[:, 0, 0] == spheres[jj, 0, 0])[0][0]
            rad = int(spheres[jj, 2, 0])
            ci, cj, ck = int(spheres[jj, 1, 0]), int(spheres[jj, 1, 1]), int(spheres[jj, 1, 2])
            
            s_i = slice(ci - rad, ci + rad)
            s_j = slice(cj - rad + pbc_shift, cj + rad + pbc_shift)
            s_k = slice(ck - rad + pbc_shift, ck + rad + pbc_shift)
            
            dist_sq = (gi[s_i] - ci)**2 + (gj[:, s_j, :] - cj)**2 + (gk[:, :, s_k] - ck)**2
            mask = dist_sq <= rad**2
            
            vals = tf.concat([materials_properties[m_id, 1, :], materials_properties[m_id, 2, :]], axis=0)
            
            for idx, tensor in enumerate(ep_sig_tensors):
                sub_volume = tensor[s_i, s_j, s_k]
                # Replace sub-volume with masked values
                updated_sub = tf.where(mask, vals[idx], sub_volume)
                tensor[s_i, s_j, s_k].assign(updated_sub)
                
            jj += 1

        # --- CYLINDER ---
        elif obj_type == 'cylinder':
            m_id = tf.where(materials_properties[:, 0, 0] == cylinders[kk, 0, 0])[0][0]
            axis, length, rad = int(cylinders[kk, 1, 0]), int(cylinders[kk, 3, 0]), int(cylinders[kk, 3, 1])
            bi, bj, bk = cylinders[kk, 2, 0], cylinders[kk, 2, 1], cylinders[kk, 2, 2]

            D, H, W = relative_ep_x_cell.shape

            if axis == 0: # X-Axis
                s_i = slice(max(0, int(bi)), min(D, int(bi + length)))
                s_j = slice(max(0, int(bj - rad + pbc_shift)), min(H, int(bj + rad + pbc_shift)))
                s_k = slice(max(0, int(bk - rad + pbc_shift)), min(W, int(bk + rad + pbc_shift)))
                mask = (gj[:, s_j, :] - bj)**2 + (gk[:, :, s_k] - bk)**2 <= rad**2
            elif axis == 1: # Y-Axis
                s_i = slice(max(0, int(bi - rad)), min(D, int(bi + rad)))
                s_j = slice(max(0, int(bj + pbc_shift)), min(H, int(bj + length + pbc_shift)))
                s_k = slice(max(0, int(bk - rad + pbc_shift)), min(W, int(bk + rad + pbc_shift)))
                mask = (gi[s_i] - bi)**2 + (gk[:, :, s_k] - bk)**2 <= rad**2
            else: # Z-Axis
                s_i = slice(max(0, int(bi - rad)), min(D, int(bi + rad)))
                s_j = slice(max(0, int(bj - rad + pbc_shift)), min(H, int(bj + rad + pbc_shift)))
                s_k = slice(max(0, int(bk + pbc_shift)), min(W, int(bk + length + pbc_shift)))
                mask = (gi[s_i] - bi)**2 + (gj[:, s_j, :] - bj)**2 <= rad**2

            vals = tf.concat([materials_properties[m_id, 1, :], materials_properties[m_id, 2, :]], axis=0)
            
            for idx, tensor in enumerate(ep_sig_tensors):
                sub_volume = tensor[s_i, s_j, s_k]
                # broadcast_to ensures the 2D mask matches the 3D sub_volume
                current_mask = tf.broadcast_to(mask, sub_volume.shape)
                updated_sub = tf.where(current_mask, vals[idx], sub_volume)
                tensor[s_i, s_j, s_k].assign(updated_sub)
                
            kk += 1

    # Helper to find material IDs for sheets
    def get_material_ids(sheets, properties, dtype):
        # Broadcast comparison
        matches = tf.equal(sheets[:, 0:1, 0], properties[:, 0, 0])
        # Cast boolean to the requested dtype before argmax
        return tf.argmax(tf.cast(matches, dtype=dtype), axis=1)

    # Pre-calculate material IDs for all three orientations
    if num_sheets_x > 0:
        m_ids_x = get_material_ids(sheets_x, sheet_properties, dtype)
    if num_sheets_y > 0:
        m_ids_y = get_material_ids(sheets_y, sheet_properties, dtype)
    if num_sheets_z > 0:
        m_ids_z = get_material_ids(sheets_z, sheet_properties, dtype)

    # --- X-Oriented Sheets ---
    for ii in range(num_sheets_x):
        i = int(sheets_x[ii, 1, 0])
        j_s, j_e = int(sheets_x[ii, 2, 0] + pbc_shift), int(sheets_x[ii, 2, 0] + sheets_x[ii, 3, 0] + pbc_shift)
        k_s, k_e = int(sheets_x[ii, 2, 1] + pbc_shift), int(sheets_x[ii, 2, 1] + sheets_x[ii, 3, 1] + pbc_shift)
        
        mid = m_ids_x[ii]
        # In TF, we must match the slice shape during assignment
        shape = (j_e - j_s, k_e - k_s)
        
        sheet_ep_x_cell_x[i, j_s:j_e, k_s:k_e].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 2, 0], dtype)))
        sheet_ep_y_cell_x[i, j_s:j_e, k_s:k_e].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 2, 1], dtype)))
        sheet_ep_z_cell_x[i, j_s:j_e, k_s:k_e].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 2, 2], dtype)))
        sheet_sig_x_cell_x[i, j_s:j_e, k_s:k_e].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 3, 0], dtype)))
        sheet_sig_y_cell_x[i, j_s:j_e, k_s:k_e].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 3, 1], dtype)))
        sheet_sig_z_cell_x[i, j_s:j_e, k_s:k_e].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 3, 2], dtype)))

    # --- Y-Oriented Sheets ---
    for ii in range(num_sheets_y):
        j = int(sheets_y[ii, 1, 0] + pbc_shift)
        i_s, i_e = int(sheets_y[ii, 2, 0]), int(sheets_y[ii, 2, 0] + sheets_y[ii, 3, 0])
        k_s, k_e = int(sheets_y[ii, 2, 1] + pbc_shift), int(sheets_y[ii, 2, 1] + sheets_y[ii, 3, 1] + pbc_shift)
        
        mid = m_ids_y[ii]
        shape = (i_e - i_s, k_e - k_s)
        
        sheet_ep_x_cell_y[i_s:i_e, j, k_s:k_e].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 2, 0], dtype)))
        sheet_ep_y_cell_y[i_s:i_e, j, k_s:k_e].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 2, 1], dtype)))
        sheet_ep_z_cell_y[i_s:i_e, j, k_s:k_e].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 2, 2], dtype)))
        sheet_sig_x_cell_y[i_s:i_e, j, k_s:k_e].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 3, 0], dtype)))
        sheet_sig_y_cell_y[i_s:i_e, j, k_s:k_e].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 3, 1], dtype)))
        sheet_sig_z_cell_y[i_s:i_e, j, k_s:k_e].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 3, 2], dtype)))

    # --- Z-Oriented Sheets ---
    for ii in range(num_sheets_z):
        k = int(sheets_z[ii, 1, 0] + pbc_shift)
        i_s, i_e = int(sheets_z[ii, 2, 0]), int(sheets_z[ii, 2, 0] + sheets_z[ii, 3, 0])
        j_s, j_e = int(sheets_z[ii, 2, 1] + pbc_shift), int(sheets_z[ii, 2, 1] + sheets_z[ii, 3, 1] + pbc_shift)
        
        mid = m_ids_z[ii]
        shape = (i_e - i_s, j_e - j_s)
        
        sheet_ep_x_cell_z[i_s:i_e, j_s:j_e, k].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 2, 0], dtype)))
        sheet_ep_y_cell_z[i_s:i_e, j_s:j_e, k].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 2, 1], dtype)))
        sheet_ep_z_cell_z[i_s:i_e, j_s:j_e, k].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 2, 2], dtype)))
        sheet_sig_x_cell_z[i_s:i_e, j_s:j_e, k].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 3, 0], dtype)))
        sheet_sig_y_cell_z[i_s:i_e, j_s:j_e, k].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 3, 1], dtype)))
        sheet_sig_z_cell_z[i_s:i_e, j_s:j_e, k].assign(tf.fill(shape, tf.cast(sheet_properties[mid, 3, 2], dtype)))

    # now if pbc we need to account for this in material cells before moving to yee cell creation
    if type_sim == 0:
        # z-direction PBC for i,j
        relative_ep_x_cell[:, :, 0].assign(relative_ep_x_cell[:, :, z_size-2])
        relative_ep_y_cell[:, :, 0].assign(relative_ep_y_cell[:, :, z_size-2])
        relative_ep_z_cell[:, :, 0].assign(relative_ep_z_cell[:, :, z_size-2])
        sigma_x_cell[:, :, 0].assign(sigma_x_cell[:, :, z_size-2])
        sigma_y_cell[:, :, 0].assign(sigma_y_cell[:, :, z_size-2])
        sigma_z_cell[:, :, 0].assign(sigma_z_cell[:, :, z_size-2])

        sheet_ep_x_cell_x[:, :, 0].assign(sheet_ep_x_cell_x[:, :, z_size-2])
        sheet_ep_y_cell_x[:, :, 0].assign(sheet_ep_y_cell_x[:, :, z_size-2])
        sheet_ep_z_cell_x[:, :, 0].assign(sheet_ep_z_cell_x[:, :, z_size-2])
        sheet_sig_x_cell_x[:, :, 0].assign(sheet_sig_x_cell_x[:, :, z_size-2])
        sheet_sig_y_cell_x[:, :, 0].assign(sheet_sig_y_cell_x[:, :, z_size-2])
        sheet_sig_z_cell_x[:, :, 0].assign(sheet_sig_z_cell_x[:, :, z_size-2])

        sheet_ep_x_cell_y[:, :, 0].assign(sheet_ep_x_cell_y[:, :, z_size-2])
        sheet_ep_y_cell_y[:, :, 0].assign(sheet_ep_y_cell_y[:, :, z_size-2])
        sheet_ep_z_cell_y[:, :, 0].assign(sheet_ep_z_cell_y[:, :, z_size-2])
        sheet_sig_x_cell_y[:, :, 0].assign(sheet_sig_x_cell_y[:, :, z_size-2])
        sheet_sig_y_cell_y[:, :, 0].assign(sheet_sig_y_cell_y[:, :, z_size-2])
        sheet_sig_z_cell_y[:, :, 0].assign(sheet_sig_z_cell_y[:, :, z_size-2])

        sheet_ep_x_cell_z[:, :, 0].assign(sheet_ep_x_cell_z[:, :, z_size-2])
        sheet_ep_y_cell_z[:, :, 0].assign(sheet_ep_y_cell_z[:, :, z_size-2])
        sheet_ep_z_cell_z[:, :, 0].assign(sheet_ep_z_cell_z[:, :, z_size-2])
        sheet_sig_x_cell_z[:, :, 0].assign(sheet_sig_x_cell_z[:, :, z_size-2])
        sheet_sig_y_cell_z[:, :, 0].assign(sheet_sig_y_cell_z[:, :, z_size-2])
        sheet_sig_z_cell_z[:, :, 0].assign(sheet_sig_z_cell_z[:, :, z_size-2])

        # y-direction PBC for i,k
        relative_ep_x_cell[:, 0, :].assign(relative_ep_x_cell[:, y_size-2, :])
        relative_ep_y_cell[:, 0, :].assign(relative_ep_y_cell[:, y_size-2, :])
        relative_ep_z_cell[:, 0, :].assign(relative_ep_z_cell[:, y_size-2, :])
        sigma_x_cell[:, 0, :].assign(sigma_x_cell[:, y_size-2, :])
        sigma_y_cell[:, 0, :].assign(sigma_y_cell[:, y_size-2, :])
        sigma_z_cell[:, 0, :].assign(sigma_z_cell[:, y_size-2, :])

        sheet_ep_x_cell_x[:, 0, :].assign(sheet_ep_x_cell_x[:, y_size-2, :])
        sheet_ep_y_cell_x[:, 0, :].assign(sheet_ep_y_cell_x[:, y_size-2, :])
        sheet_ep_z_cell_x[:, 0, :].assign(sheet_ep_z_cell_x[:, y_size-2, :])
        sheet_sig_x_cell_x[:, 0, :].assign(sheet_sig_x_cell_x[:, y_size-2, :])
        sheet_sig_y_cell_x[:, 0, :].assign(sheet_sig_y_cell_x[:, y_size-2, :])
        sheet_sig_z_cell_x[:, 0, :].assign(sheet_sig_z_cell_x[:, y_size-2, :])

        sheet_ep_x_cell_y[:, 0, :].assign(sheet_ep_x_cell_y[:, y_size-2, :])
        sheet_ep_y_cell_y[:, 0, :].assign(sheet_ep_y_cell_y[:, y_size-2, :])
        sheet_ep_z_cell_y[:, 0, :].assign(sheet_ep_z_cell_y[:, y_size-2, :])
        sheet_sig_x_cell_y[:, 0, :].assign(sheet_sig_x_cell_y[:, y_size-2, :])
        sheet_sig_y_cell_y[:, 0, :].assign(sheet_sig_y_cell_y[:, y_size-2, :])
        sheet_sig_z_cell_y[:, 0, :].assign(sheet_sig_z_cell_y[:, y_size-2, :])

        sheet_ep_x_cell_z[:, 0, :].assign(sheet_ep_x_cell_z[:, y_size-2, :])
        sheet_ep_y_cell_z[:, 0, :].assign(sheet_ep_y_cell_z[:, y_size-2, :])
        sheet_ep_z_cell_z[:, 0, :].assign(sheet_ep_z_cell_z[:, y_size-2, :])
        sheet_sig_x_cell_z[:, 0, :].assign(sheet_sig_x_cell_z[:, y_size-2, :])
        sheet_sig_y_cell_z[:, 0, :].assign(sheet_sig_y_cell_z[:, y_size-2, :])
        sheet_sig_z_cell_z[:, 0, :].assign(sheet_sig_z_cell_z[:, y_size-2, :])

    # now material cells to yee cell creation
    # Bulk materials
    relative_ep_x[1:-1, 1:-1, 1:-1].assign(
        (
            relative_ep_x_cell[1:-1, 1:-1, 1:-1] +
            relative_ep_x_cell[1:-1, 0:-2, 1:-1] +
            relative_ep_x_cell[1:-1, 1:-1, 0:-2] +
            relative_ep_x_cell[1:-1, 0:-2, 0:-2]
        ) / 4.0
    )

    relative_ep_y[1:-1, 1:-1, 1:-1].assign(
        (
            relative_ep_y_cell[1:-1, 1:-1, 1:-1] +
            relative_ep_y_cell[0:-2, 1:-1, 1:-1] +
            relative_ep_y_cell[1:-1, 1:-1, 0:-2] +
            relative_ep_y_cell[0:-2, 1:-1, 0:-2]
        ) / 4.0
    )

    relative_ep_z[1:-1, 1:-1, 1:-1].assign(
        (
            relative_ep_z_cell[1:-1, 1:-1, 1:-1] +
            relative_ep_z_cell[0:-2, 1:-1, 1:-1] +
            relative_ep_z_cell[1:-1, 0:-2, 1:-1] +
            relative_ep_z_cell[0:-2, 0:-2, 1:-1]
        ) / 4.0
    )

    sigma_x[1:-1, 1:-1, 1:-1].assign(
        (
            sigma_x_cell[1:-1, 1:-1, 1:-1] +
            sigma_x_cell[1:-1, 0:-2, 1:-1] +
            sigma_x_cell[1:-1, 1:-1, 0:-2] +
            sigma_x_cell[1:-1, 0:-2, 0:-2]
        ) / 4.0
    )

    sigma_y[1:-1, 1:-1, 1:-1].assign(
        (
            sigma_y_cell[1:-1, 1:-1, 1:-1] +
            sigma_y_cell[0:-2, 1:-1, 1:-1] +
            sigma_y_cell[1:-1, 1:-1, 0:-2] +
            sigma_y_cell[0:-2, 1:-1, 0:-2]
        ) / 4.0
    )

    sigma_z[1:-1, 1:-1, 1:-1].assign(
        (
            sigma_z_cell[1:-1, 1:-1, 1:-1] +
            sigma_z_cell[0:-2, 1:-1, 1:-1] +
            sigma_z_cell[1:-1, 0:-2, 1:-1] +
            sigma_z_cell[0:-2, 0:-2, 1:-1]
        ) / 4.0
    )

    # These next 3 if statements just save time for material to yee if nothing it's 1 and 0s anyway in all places
    # Now put sheets into yee cells from material cells
    if num_sheets_x>0:
        # x-normal sheets
        sheet_ep_x_x[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_ep_x_cell_x[1:-1, 1:-1, 1:-1] +
                sheet_ep_x_cell_x[1:-1, 0:-2, 1:-1] +
                sheet_ep_x_cell_x[1:-1, 1:-1, 0:-2] +
                sheet_ep_x_cell_x[1:-1, 0:-2, 0:-2]
            ) / 4.0
        )

        sheet_ep_y_x[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_ep_y_cell_x[1:-1, 1:-1, 1:-1] +
                sheet_ep_y_cell_x[1:-1, 1:-1, 0:-2] +
                sheet_ep_y_cell_x[1:-1, 1:-1, 1:-1] +
                sheet_ep_y_cell_x[1:-1, 1:-1, 0:-2]
            ) / 4.0
        )

        sheet_ep_z_x[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_ep_z_cell_x[1:-1, 1:-1, 1:-1] +
                sheet_ep_z_cell_x[1:-1, 0:-2, 1:-1] +
                sheet_ep_z_cell_x[1:-1, 1:-1, 1:-1] +
                sheet_ep_z_cell_x[1:-1, 0:-2, 1:-1]
            ) / 4.0
        )

        sheet_sig_x_x[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_sig_x_cell_x[1:-1, 1:-1, 1:-1] +
                sheet_sig_x_cell_x[1:-1, 0:-2, 1:-1] +
                sheet_sig_x_cell_x[1:-1, 1:-1, 0:-2] +
                sheet_sig_x_cell_x[1:-1, 0:-2, 0:-2]
            ) / 4.0
        )

        sheet_sig_y_x[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_sig_y_cell_x[1:-1, 1:-1, 1:-1] +
                sheet_sig_y_cell_x[1:-1, 1:-1, 0:-2] +
                sheet_sig_y_cell_x[1:-1, 1:-1, 1:-1] +
                sheet_sig_y_cell_x[1:-1, 1:-1, 0:-2]
            ) / 4.0
        )

        sheet_sig_z_x[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_sig_z_cell_x[1:-1, 1:-1, 1:-1] +
                sheet_sig_z_cell_x[1:-1, 0:-2, 1:-1] +
                sheet_sig_z_cell_x[1:-1, 1:-1, 1:-1] +
                sheet_sig_z_cell_x[1:-1, 0:-2, 1:-1]
            ) / 4.0
        )

    if num_sheets_y>0:
        # y-normal sheets
        sheet_ep_x_y[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_ep_x_cell_y[1:-1, 1:-1, 1:-1] +
                sheet_ep_x_cell_y[1:-1, 1:-1, 0:-2] +
                sheet_ep_x_cell_y[1:-1, 1:-1, 1:-1] +
                sheet_ep_x_cell_y[1:-1, 1:-1, 0:-2]
            ) / 4.0
        )

        sheet_ep_y_y[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_ep_y_cell_y[1:-1, 1:-1, 1:-1] +
                sheet_ep_y_cell_y[1:-1, 1:-1, 0:-2] +
                sheet_ep_y_cell_y[0:-2, 1:-1, 1:-1] +
                sheet_ep_y_cell_y[0:-2, 1:-1, 0:-2]
            ) / 4.0
        )

        sheet_ep_z_y[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_ep_z_cell_y[1:-1, 1:-1, 1:-1] +
                sheet_ep_z_cell_y[0:-2, 1:-1, 1:-1] +
                sheet_ep_z_cell_y[1:-1, 1:-1, 1:-1] +
                sheet_ep_z_cell_y[0:-2, 1:-1, 1:-1]
            ) / 4.0
        )

        sheet_sig_x_y[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_sig_x_cell_y[1:-1, 1:-1, 1:-1] +
                sheet_sig_x_cell_y[1:-1, 1:-1, 0:-2] +
                sheet_sig_x_cell_y[1:-1, 1:-1, 1:-1] +
                sheet_sig_x_cell_y[1:-1, 1:-1, 0:-2]
            ) / 4.0
        )

        sheet_sig_y_y[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_sig_y_cell_y[1:-1, 1:-1, 1:-1] +
                sheet_sig_y_cell_y[1:-1, 1:-1, 0:-2] +
                sheet_sig_y_cell_y[0:-2, 1:-1, 1:-1] +
                sheet_sig_y_cell_y[0:-2, 1:-1, 0:-2]
            ) / 4.0
        )

        sheet_sig_z_y[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_sig_z_cell_y[1:-1, 1:-1, 1:-1] +
                sheet_sig_z_cell_y[0:-2, 1:-1, 1:-1] +
                sheet_sig_z_cell_y[1:-1, 1:-1, 1:-1] +
                sheet_sig_z_cell_y[0:-2, 1:-1, 1:-1]
            ) / 4.0
        )

    if num_sheets_z>0:
        # z-normal sheets
        sheet_ep_x_z[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_ep_x_cell_z[1:-1, 1:-1, 1:-1] +
                sheet_ep_x_cell_z[1:-1, 0:-2, 1:-1] +
                sheet_ep_x_cell_z[1:-1, 1:-1, 1:-1] +
                sheet_ep_x_cell_z[1:-1, 0:-2, 1:-1]
            ) / 4.0
        )

        sheet_ep_y_z[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_ep_y_cell_z[1:-1, 1:-1, 1:-1] +
                sheet_ep_y_cell_z[0:-2, 1:-1, 1:-1] +
                sheet_ep_y_cell_z[1:-1, 1:-1, 1:-1] +
                sheet_ep_y_cell_z[0:-2, 1:-1, 1:-1]
            ) / 4.0
        )

        sheet_ep_z_z[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_ep_z_cell_z[1:-1, 1:-1, 1:-1] +
                sheet_ep_z_cell_z[0:-2, 1:-1, 1:-1] +
                sheet_ep_z_cell_z[1:-1, 0:-2, 1:-1] +
                sheet_ep_z_cell_z[0:-2, 0:-2, 1:-1]
            ) / 4.0
        )

        sheet_sig_x_z[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_sig_x_cell_z[1:-1, 1:-1, 1:-1] +
                sheet_sig_x_cell_z[1:-1, 0:-2, 1:-1] +
                sheet_sig_x_cell_z[1:-1, 1:-1, 1:-1] +
                sheet_sig_x_cell_z[1:-1, 0:-2, 1:-1]
            ) / 4.0
        )

        sheet_sig_y_z[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_sig_y_cell_z[1:-1, 1:-1, 1:-1] +
                sheet_sig_y_cell_z[0:-2, 1:-1, 1:-1] +
                sheet_sig_y_cell_z[1:-1, 1:-1, 1:-1] +
                sheet_sig_y_cell_z[0:-2, 1:-1, 1:-1]
            ) / 4.0
        )

        sheet_sig_z_z[1:-1, 1:-1, 1:-1].assign(
            (
                sheet_sig_z_cell_z[1:-1, 1:-1, 1:-1] +
                sheet_sig_z_cell_z[0:-2, 1:-1, 1:-1] +
                sheet_sig_z_cell_z[1:-1, 0:-2, 1:-1] +
                sheet_sig_z_cell_z[0:-2, 0:-2, 1:-1]
            ) / 4.0
        )

    # now if pbc exists we need to apply the same condition like we did for material cells
    if type_sim == 0:
        # z-boundary
        relative_ep_x[:, :, 0].assign(relative_ep_x[:, :, z_size-2])
        relative_ep_y[:, :, 0].assign(relative_ep_y[:, :, z_size-2])
        relative_ep_z[:, :, 0].assign(relative_ep_z[:, :, z_size-2])
        sigma_x[:, :, 0].assign(sigma_x[:, :, z_size-2])
        sigma_y[:, :, 0].assign(sigma_y[:, :, z_size-2])
        sigma_z[:, :, 0].assign(sigma_z[:, :, z_size-2])

        # y-boundary
        relative_ep_x[:, 0, :].assign(relative_ep_x[:, y_size-2, :])
        relative_ep_y[:, 0, :].assign(relative_ep_y[:, y_size-2, :])
        relative_ep_z[:, 0, :].assign(relative_ep_z[:, y_size-2, :])
        sigma_x[:, 0, :].assign(sigma_x[:, y_size-2, :])
        sigma_y[:, 0, :].assign(sigma_y[:, y_size-2, :])
        sigma_z[:, 0, :].assign(sigma_z[:, y_size-2, :])

        if num_sheets_x > 0:
            sheet_ep_x_x[:, :, 0].assign(sheet_ep_x_x[:, :, z_size-2])
            sheet_ep_y_x[:, :, 0].assign(sheet_ep_y_x[:, :, z_size-2])
            sheet_ep_z_x[:, :, 0].assign(sheet_ep_z_x[:, :, z_size-2])
            sheet_sig_x_x[:, :, 0].assign(sheet_sig_x_x[:, :, z_size-2])
            sheet_sig_y_x[:, :, 0].assign(sheet_sig_y_x[:, :, z_size-2])
            sheet_sig_z_x[:, :, 0].assign(sheet_sig_z_x[:, :, z_size-2])

            sheet_ep_x_x[:, 0, :].assign(sheet_ep_x_x[:, y_size-2, :])
            sheet_ep_y_x[:, 0, :].assign(sheet_ep_y_x[:, y_size-2, :])
            sheet_ep_z_x[:, 0, :].assign(sheet_ep_z_x[:, y_size-2, :])
            sheet_sig_x_x[:, 0, :].assign(sheet_sig_x_x[:, y_size-2, :])
            sheet_sig_y_x[:, 0, :].assign(sheet_sig_y_x[:, y_size-2, :])
            sheet_sig_z_x[:, 0, :].assign(sheet_sig_z_x[:, y_size-2, :])

        if num_sheets_y > 0:
            sheet_ep_x_y[:, :, 0].assign(sheet_ep_x_y[:, :, z_size-2])
            sheet_ep_y_y[:, :, 0].assign(sheet_ep_y_y[:, :, z_size-2])
            sheet_ep_z_y[:, :, 0].assign(sheet_ep_z_y[:, :, z_size-2])
            sheet_sig_x_y[:, :, 0].assign(sheet_sig_x_y[:, :, z_size-2])
            sheet_sig_y_y[:, :, 0].assign(sheet_sig_y_y[:, :, z_size-2])
            sheet_sig_z_y[:, :, 0].assign(sheet_sig_z_y[:, :, z_size-2])

            sheet_ep_x_y[:, 0, :].assign(sheet_ep_x_y[:, y_size-2, :])
            sheet_ep_y_y[:, 0, :].assign(sheet_ep_y_y[:, y_size-2, :])
            sheet_ep_z_y[:, 0, :].assign(sheet_ep_z_y[:, y_size-2, :])
            sheet_sig_x_y[:, 0, :].assign(sheet_sig_x_y[:, y_size-2, :])
            sheet_sig_y_y[:, 0, :].assign(sheet_sig_y_y[:, y_size-2, :])
            sheet_sig_z_y[:, 0, :].assign(sheet_sig_z_y[:, y_size-2, :])

        if num_sheets_z > 0:
            sheet_ep_x_z[:, :, 0].assign(sheet_ep_x_z[:, :, z_size-2])
            sheet_ep_y_z[:, :, 0].assign(sheet_ep_y_z[:, :, z_size-2])
            sheet_ep_z_z[:, :, 0].assign(sheet_ep_z_z[:, :, z_size-2])
            sheet_sig_x_z[:, :, 0].assign(sheet_sig_x_z[:, :, z_size-2])
            sheet_sig_y_z[:, :, 0].assign(sheet_sig_y_z[:, :, z_size-2])
            sheet_sig_z_z[:, :, 0].assign(sheet_sig_z_z[:, :, z_size-2])

            sheet_ep_x_z[:, 0, :].assign(sheet_ep_x_z[:, y_size-2, :])
            sheet_ep_y_z[:, 0, :].assign(sheet_ep_y_z[:, y_size-2, :])
            sheet_ep_z_z[:, 0, :].assign(sheet_ep_z_z[:, y_size-2, :])
            sheet_sig_x_z[:, 0, :].assign(sheet_sig_x_z[:, y_size-2, :])
            sheet_sig_y_z[:, 0, :].assign(sheet_sig_y_z[:, y_size-2, :])
            sheet_sig_z_z[:, 0, :].assign(sheet_sig_z_z[:, y_size-2, :])

    # X-component
    denom_x = 1.0 + sigma_x[:-1, :-1, :-1] * del_t / (
        2.0 * ep_0 * relative_ep_x[:-1, :-1, :-1]
    )

    gax[:-1, :-1, :-1].assign(
        (1.0 - sigma_x[:-1, :-1, :-1] * del_t / (
            2.0 * ep_0 * relative_ep_x[:-1, :-1, :-1]
        )) / denom_x
    )

    gbx[:-1, :-1, :-1].assign(
        (del_t / (ep_0 * relative_ep_x[:-1, :-1, :-1])) / denom_x
    )

    # Y-component
    denom_y = 1.0 + sigma_y[:-1, :-1, :-1] * del_t / (
        2.0 * ep_0 * relative_ep_y[:-1, :-1, :-1]
    )

    gay[:-1, :-1, :-1].assign(
        (1.0 - sigma_y[:-1, :-1, :-1] * del_t / (
            2.0 * ep_0 * relative_ep_y[:-1, :-1, :-1]
        )) / denom_y
    )

    gby[:-1, :-1, :-1].assign(
        (del_t / (ep_0 * relative_ep_y[:-1, :-1, :-1])) / denom_y
    )

    # Z-component
    denom_z = 1.0 + sigma_z[:-1, :-1, :-1] * del_t / (
        2.0 * ep_0 * relative_ep_z[:-1, :-1, :-1]
    )

    gaz[:-1, :-1, :-1].assign(
        (1.0 - sigma_z[:-1, :-1, :-1] * del_t / (
            2.0 * ep_0 * relative_ep_z[:-1, :-1, :-1]
        )) / denom_z
    )

    gbz[:-1, :-1, :-1].assign(
        (del_t / (ep_0 * relative_ep_z[:-1, :-1, :-1])) / denom_z
    )

    # Now lastly, we need to deal with sheet adjustments to ga,gb
    # Initialize sheet averages with default material values

    # Create Variables instead of tensors
    sheet_ep_avg_x = tf.Variable(relative_ep_x)
    sheet_ep_avg_y = tf.Variable(relative_ep_y)
    sheet_ep_avg_z = tf.Variable(relative_ep_z)
    sheet_sig_avg_x = tf.Variable(sigma_x)
    sheet_sig_avg_y = tf.Variable(sigma_y)
    sheet_sig_avg_z = tf.Variable(sigma_z)

    # Counters
    counter_x = tf.Variable(tf.zeros_like(relative_ep_x, dtype=tf.bool))
    counter_y = tf.Variable(tf.zeros_like(relative_ep_y, dtype=tf.bool))
    counter_z = tf.Variable(tf.zeros_like(relative_ep_z, dtype=tf.bool))

    # --- x-normal sheets ---
    mask_y = tf.logical_or(sheet_sig_y_x > 0, sheet_ep_y_x > 1)
    mask_z = tf.logical_or(sheet_sig_z_x > 0, sheet_ep_z_x > 1)

    sheet_ep_avg_y.assign(
        tf.where(
            mask_y,
            sheet_ep_avg_y + (sheet_thickness / del_x) *
            (sheet_ep_y_x - relative_ep_y),
            sheet_ep_avg_y
        )
    )

    sheet_sig_avg_y.assign(
        tf.where(
            mask_y,
            sheet_sig_avg_y + (sheet_thickness / del_x) *
            (sheet_sig_y_x - sigma_y),
            sheet_sig_avg_y
        )
    )

    counter_y.assign(tf.logical_or(counter_y, mask_y))

    sheet_ep_avg_z.assign(
        tf.where(
            mask_z,
            sheet_ep_avg_z + (sheet_thickness / del_x) *
            (sheet_ep_z_x - relative_ep_z),
            sheet_ep_avg_z
        )
    )

    sheet_sig_avg_z.assign(
        tf.where(
            mask_z,
            sheet_sig_avg_z + (sheet_thickness / del_x) *
            (sheet_sig_z_x - sigma_z),
            sheet_sig_avg_z
        )
    )

    counter_z.assign(tf.logical_or(counter_z, mask_z))

    # --- y-normal sheets ---
    mask_x = tf.logical_or(sheet_sig_x_y > 0, sheet_ep_x_y > 1)
    mask_z_y = tf.logical_or(sheet_sig_z_y > 0, sheet_ep_z_y > 1)

    sheet_ep_avg_x.assign(
        tf.where(
            mask_x,
            sheet_ep_avg_x + (sheet_thickness / del_y) *
            (sheet_ep_x_y - relative_ep_x),
            sheet_ep_avg_x
        )
    )

    sheet_sig_avg_x.assign(
        tf.where(
            mask_x,
            sheet_sig_avg_x + (sheet_thickness / del_y) *
            (sheet_sig_x_y - sigma_x),
            sheet_sig_avg_x
        )
    )

    counter_x.assign(tf.logical_or(counter_x, mask_x))

    sheet_ep_avg_z.assign(
        tf.where(
            mask_z_y,
            sheet_ep_avg_z + (sheet_thickness / del_y) *
            (sheet_ep_z_y - relative_ep_z),
            sheet_ep_avg_z
        )
    )

    sheet_sig_avg_z.assign(
        tf.where(
            mask_z_y,
            sheet_sig_avg_z + (sheet_thickness / del_y) *
            (sheet_sig_z_y - sigma_z),
            sheet_sig_avg_z
        )
    )

    counter_z.assign(tf.logical_or(counter_z, mask_z_y))

    # --- z-normal sheets ---
    mask_x_z = tf.logical_or(sheet_sig_x_z > 0, sheet_ep_x_z > 1)
    mask_y_z = tf.logical_or(sheet_sig_y_z > 0, sheet_ep_y_z > 1)

    sheet_ep_avg_x.assign(
        tf.where(
            mask_x_z,
            sheet_ep_avg_x + (sheet_thickness / del_z) *
            (sheet_ep_x_z - relative_ep_x),
            sheet_ep_avg_x
        )
    )

    sheet_sig_avg_x.assign(
        tf.where(
            mask_x_z,
            sheet_sig_avg_x + (sheet_thickness / del_z) *
            (sheet_sig_x_z - sigma_x),
            sheet_sig_avg_x
        )
    )

    counter_x.assign(tf.logical_or(counter_x, mask_x_z))

    sheet_ep_avg_y.assign(
        tf.where(
            mask_y_z,
            sheet_ep_avg_y + (sheet_thickness / del_z) *
            (sheet_ep_y_z - relative_ep_y),
            sheet_ep_avg_y
        )
    )

    sheet_sig_avg_y.assign(
        tf.where(
            mask_y_z,
            sheet_sig_avg_y + (sheet_thickness / del_z) *
            (sheet_sig_y_z - sigma_y),
            sheet_sig_avg_y
        )
    )

    counter_y.assign(tf.logical_or(counter_y, mask_y_z))

    # x-direction
    mask = counter_x
    gax.assign(
        tf.where(
            mask,
            (1.0 - sheet_sig_avg_x * del_t / (2.0 * ep_0 * sheet_ep_avg_x)) /
            (1.0 + sheet_sig_avg_x * del_t / (2.0 * ep_0 * sheet_ep_avg_x)),
            gax
        )
    )

    gbx.assign(
        tf.where(
            mask,
            (del_t / (ep_0 * sheet_ep_avg_x)) /
            (1.0 + sheet_sig_avg_x * del_t / (2.0 * ep_0 * sheet_ep_avg_x)),
            gbx
        )
    )

    # y-direction
    mask = counter_y
    gay.assign(
        tf.where(
            mask,
            (1.0 - sheet_sig_avg_y * del_t / (2.0 * ep_0 * sheet_ep_avg_y)) /
            (1.0 + sheet_sig_avg_y * del_t / (2.0 * ep_0 * sheet_ep_avg_y)),
            gay
        )
    )

    gby.assign(
        tf.where(
            mask,
            (del_t / (ep_0 * sheet_ep_avg_y)) /
            (1.0 + sheet_sig_avg_y * del_t / (2.0 * ep_0 * sheet_ep_avg_y)),
            gby
        )
    )

    # z-direction
    mask = counter_z
    gaz.assign(
        tf.where(
            mask,
            (1.0 - sheet_sig_avg_z * del_t / (2.0 * ep_0 * sheet_ep_avg_z)) /
            (1.0 + sheet_sig_avg_z * del_t / (2.0 * ep_0 * sheet_ep_avg_z)),
            gaz
        )
    )

    gbz.assign(
        tf.where(
            mask,
            (del_t / (ep_0 * sheet_ep_avg_z)) /
            (1.0 + sheet_sig_avg_z * del_t / (2.0 * ep_0 * sheet_ep_avg_z)),
            gbz
        )
    )

    # masks used in main FDTD
    mask_x = tf.cast(tf.logical_or(sheet_sig_x_x > 0, sheet_ep_x_x > 1), dtype=dtype)
    mask_y = tf.cast(tf.logical_or(sheet_sig_y_y > 0, sheet_ep_y_y > 1), dtype=dtype)
    mask_z = tf.cast(tf.logical_or(sheet_sig_z_z > 0, sheet_ep_z_z > 1), dtype=dtype)

##########################################################
# Definitions for main FDTD loop - for calling individual sections to make code more readable
##########################################################
@tf.function(jit_compile=jit_on)
def update_h_fields(Hx, Hy, Hz, Ex, Ey, Ez, da, db, den_hx, den_hy, den_hz):

    # Update Hx
    Hx[:-1, :-1, :-1].assign(da * Hx[:-1, :-1, :-1] + db * (
        (Ey[:-1, :-1, 1:] - Ey[:-1, :-1, :-1]) * den_hz[None,None,:] +
        (Ez[:-1, :-1, :-1] - Ez[:-1, 1:, :-1]) * den_hy[None,:,None]
    ))

    # Update Hy
    Hy[:-1, :-1, :-1].assign(da * Hy[:-1, :-1, :-1] + db * (
        (Ex[:-1, :-1, :-1] - Ex[:-1, :-1, 1:]) * den_hz[None,None,:] +
        (Ez[1:, :-1, :-1] - Ez[:-1, :-1, :-1]) * den_hx[:,None,None]
    ))

    # Update Hz
    Hz[:-1, :-1, :-1].assign(da * Hz[:-1, :-1, :-1] + db * (
        (Ey[:-1, :-1, :-1] - Ey[1:, :-1, :-1]) * den_hx[:,None,None] +
        (Ex[:-1, 1:, :-1] - Ex[:-1, :-1, :-1]) * den_hy[None,:,None]
    ))

@tf.function(jit_compile=jit_on)
def update_h_pml(Hx, Hy, Hz, Ex, Ey, Ez, db, del_x, del_y, del_z,
                 psi_Hyx_1, psi_Hyx_2, psi_Hzx_1, psi_Hzx_2,
                 psi_Hxy_1, psi_Hxy_2, psi_Hzy_1, psi_Hzy_2,
                 psi_Hxz_1, psi_Hxz_2, psi_Hyz_1, psi_Hyz_2,
                 bh_x_1, ch_x_1, bh_x_2, ch_x_2,
                 bh_y_1, ch_y_1, bh_y_2, ch_y_2,
                 bh_z_1, ch_z_1, bh_z_2, ch_z_2,
                 nxPML_1, nyPML_1, nzPML_1, nxPML_2, nyPML_2, nzPML_2,
                 x_size, y_size, z_size):

    # =================================================================
    # Hx UPDATES
    # =================================================================
    if type_sim!=0:
        # --- bottom y for Hx ---
        psi_Hxy_1[:-1, :, :-1].assign(bh_y_1[:-1, :, :-1] * psi_Hxy_1[:-1, :, :-1] + ch_y_1[:-1, :, :-1] * (Ez[:-1, :nyPML_1-1, :-1] - Ez[:-1, 1:nyPML_1, :-1]) / del_y)
        Hx[:-1, :nyPML_1-1, :-1].assign(Hx[:-1, :nyPML_1-1, :-1] + db * psi_Hxy_1[:-1, :, :-1])

        # --- top y for Hx ---
        new_psi_Hxy_2 = (
            bh_y_2[:-1, :nyPML_2-1, :-1] * psi_Hxy_2[:-1, :nyPML_2-1, :-1] + 
            ch_y_2[:-1, :nyPML_2-1, :-1] * (tf.reverse(Ez[:-1, y_size-nyPML_2 : y_size-1, :-1], [1]) - 
                                            tf.reverse(Ez[:-1, y_size-nyPML_2+1 : y_size, :-1], [1])) / del_y
        )
        psi_Hxy_2[:-1, :nyPML_2-1, :-1].assign(new_psi_Hxy_2)
        Hx[:-1, y_size-nyPML_2 : y_size-1, :-1].assign(Hx[:-1, y_size-nyPML_2 : y_size-1, :-1] + db * tf.reverse(new_psi_Hxy_2, axis=[1]))

        # --- bottom z for Hx ---
        psi_Hxz_1[:-1, :, :].assign(bh_z_1[:-1, :-1, :] * psi_Hxz_1[:-1, :, :] + ch_z_1[:-1, :-1, :] * (Ey[:-1, :-1, 1:nzPML_1] - Ey[:-1, :-1, :nzPML_1-1]) / del_z)
        Hx[:-1, :-1, :nzPML_1-1].assign(Hx[:-1, :-1, :nzPML_1-1] + db * psi_Hxz_1[:-1, :, :])

        # --- top z for Hx ---
        new_psi_Hxz_2 = (
            bh_z_2[:-1, :-1, :nzPML_2-1] * psi_Hxz_2[:-1, :, :nzPML_2-1] + 
            ch_z_2[:-1, :-1, :nzPML_2-1] * (tf.reverse(Ey[:-1, :-1, z_size-nzPML_2+1 : z_size], [2]) - 
                                            tf.reverse(Ey[:-1, :-1, z_size-nzPML_2 : z_size-1], [2])) / del_z
        )
        psi_Hxz_2[:-1, :, :nzPML_2-1].assign(new_psi_Hxz_2)
        Hx[:-1, :-1, z_size-nzPML_2 : z_size-1].assign(Hx[:-1, :-1, z_size-nzPML_2 : z_size-1] + db * tf.reverse(new_psi_Hxz_2, axis=[2]))

    # =================================================================
    # Hy UPDATES 
    # =================================================================

    # --- bottom x for Hy ---
    psi_Hyx_1[:, :-1, :-1].assign(bh_x_1[:, :-1, :-1] * psi_Hyx_1[:, :-1, :-1] + ch_x_1[:, :-1, :-1] * (Ez[1:nxPML_1, :-1, :-1] - Ez[:nxPML_1-1, :-1, :-1]) / del_x)
    Hy[:nxPML_1-1, :-1, :-1].assign(Hy[:nxPML_1-1, :-1, :-1] + db * psi_Hyx_1[:, :-1, :-1])

    # --- top x for Hy ---
    new_psi_Hyx_2 = (
        bh_x_2[:nxPML_2-1, :-1, :-1] * psi_Hyx_2[:nxPML_2-1, :-1, :-1] + 
        ch_x_2[:nxPML_2-1, :-1, :-1] * (tf.reverse(Ez[x_size-nxPML_2+1 : x_size, :-1, :-1], [0]) - 
                                        tf.reverse(Ez[x_size-nxPML_2 : x_size-1, :-1, :-1], [0])) / del_x
    )
    psi_Hyx_2[:nxPML_2-1, :-1, :-1].assign(new_psi_Hyx_2)
    Hy[x_size-nxPML_2 : x_size-1, :-1, :-1].assign(Hy[x_size-nxPML_2 : x_size-1, :-1, :-1] + db * tf.reverse(new_psi_Hyx_2, axis=[0]))

    if type_sim!=0:

        # --- bottom z for Hy ---
        psi_Hyz_1[:, :-1, :].assign(bh_z_1[:-1, :-1, :] * psi_Hyz_1[:, :-1, :] + ch_z_1[:-1, :-1, :] * (Ex[:-1, :-1, :nzPML_1-1] - Ex[:-1, :-1, 1:nzPML_1]) / del_z)
        Hy[:-1, :-1, :nzPML_1-1].assign(Hy[:-1, :-1, :nzPML_1-1] + db * psi_Hyz_1[:, :-1, :])

        # --- top z for Hy ---
        new_psi_Hyz_2 = (
            bh_z_2[:-1, :-1, :nzPML_2-1] * psi_Hyz_2[:, :-1, :nzPML_2-1] + 
            ch_z_2[:-1, :-1, :nzPML_2-1] * (tf.reverse(Ex[:-1, :-1, z_size-nzPML_2 : z_size-1], [2]) - 
                                            tf.reverse(Ex[:-1, :-1, z_size-nzPML_2+1 : z_size], [2])) / del_z
        )
        psi_Hyz_2[:, :-1, :nzPML_2-1].assign(new_psi_Hyz_2)
        Hy[:-1, :-1, z_size-nzPML_2 : z_size-1].assign(Hy[:-1, :-1, z_size-nzPML_2 : z_size-1] + db * tf.reverse(new_psi_Hyz_2, axis=[2]))

    # =================================================================
    # Hz UPDATES 
    # =================================================================

    # --- bottom x for Hz ---
    psi_Hzx_1[:, :, 1:-1].assign(bh_x_1[:, :-1, 1:-1] * psi_Hzx_1[:, :, 1:-1] + ch_x_1[:, :-1, 1:-1] * (Ey[:nxPML_1-1, :-1, 1:-1] - Ey[1:nxPML_1, :-1, 1:-1]) / del_x)
    Hz[:nxPML_1-1, :-1, 1:-1].assign(Hz[:nxPML_1-1, :-1, 1:-1] + db * psi_Hzx_1[:, :, 1:-1])

    # --- top x for Hz ---
    new_psi_Hzx_2 = (
        bh_x_2[:nxPML_2-1, :-1, 1:-1] * psi_Hzx_2[:nxPML_2-1, :, 1:-1] + 
        ch_x_2[:nxPML_2-1, :-1, 1:-1] * (tf.reverse(Ey[x_size-nxPML_2 : x_size-1, :-1, 1:-1], [0]) - 
                                        tf.reverse(Ey[x_size-nxPML_2+1 : x_size, :-1, 1:-1], [0])) / del_x
    )
    psi_Hzx_2[:nxPML_2-1, :, 1:-1].assign(new_psi_Hzx_2)
    Hz[x_size-nxPML_2 : x_size-1, :-1, 1:-1].assign(Hz[x_size-nxPML_2 : x_size-1, :-1, 1:-1] + db * tf.reverse(new_psi_Hzx_2, axis=[0]))

    if type_sim!=0:

        # --- bottom y for Hz ---
        psi_Hzy_1[:, :, 1:-1].assign(bh_y_1[:-1, :, 1:-1] * psi_Hzy_1[:, :, 1:-1] + ch_y_1[:-1, :, 1:-1] * (Ex[:-1, 1:nyPML_1, 1:-1] - Ex[:-1, :nyPML_1-1, 1:-1]) / del_y)
        Hz[:-1, :nyPML_1-1, 1:-1].assign(Hz[:-1, :nyPML_1-1, 1:-1] + db * psi_Hzy_1[:, :, 1:-1])

        # --- top y for Hz ---
        new_psi_Hzy_2 = (
            bh_y_2[:-1, :nyPML_2-1, 1:-1] * psi_Hzy_2[:, :nyPML_2-1, 1:-1] + 
            ch_y_2[:-1, :nyPML_2-1, 1:-1] * (tf.reverse(Ex[:-1, y_size-nyPML_2+1 : y_size, 1:-1], [1]) - 
                                            tf.reverse(Ex[:-1, y_size-nyPML_2 : y_size-1, 1:-1], [1])) / del_y
        )
        psi_Hzy_2[:, :nyPML_2-1, 1:-1].assign(new_psi_Hzy_2)
        Hz[:-1, y_size-nyPML_2 : y_size-1, 1:-1].assign(Hz[:-1, y_size-nyPML_2 : y_size-1, 1:-1] + db * tf.reverse(new_psi_Hzy_2, axis=[1]))

@tf.function(jit_compile=jit_on)
def update_e_fields(Ex, Ey, Ez, Hx, Hy, Hz, gax, gay, gaz, gbx, gby, gbz, 
                    den_ex, den_ey, den_ez):

    # Update Ex
    Ex[:-1, 1:-1, 1:-1].assign(gax[:-1, 1:-1, 1:-1] * Ex[:-1, 1:-1, 1:-1] + gbx[:-1, 1:-1, 1:-1] * (
        (Hy[:-1, 1:-1, :-2] - Hy[:-1, 1:-1, 1:-1]) * den_ez[None,None,1:] +
        (Hz[:-1, 1:-1, 1:-1] - Hz[:-1, :-2, 1:-1]) * den_ey[None,1:,None]
    ))
    
    # Update Ey
    Ey[1:-1, :-1, 1:-1].assign(gay[1:-1, :-1, 1:-1] * Ey[1:-1, :-1, 1:-1] + gby[1:-1, :-1, 1:-1] * (
        (Hx[1:-1, :-1, 1:-1] - Hx[1:-1, :-1, :-2]) * den_ez[None,None,1:] +
        (Hz[:-2, :-1, 1:-1] - Hz[1:-1, :-1, 1:-1]) * den_ex[1:,None,None]
    ))
    
    # Update Ez
    Ez[1:-1, 1:-1, :-1].assign(gaz[1:-1, 1:-1, :-1] * Ez[1:-1, 1:-1, :-1] + gbz[1:-1, 1:-1, :-1] * (
        (Hx[1:-1, :-2, :-1] - Hx[1:-1, 1:-1, :-1]) * den_ey[None,1:,None] +
        (Hy[1:-1, 1:-1, :-1] - Hy[:-2, 1:-1, :-1]) * den_ex[1:,None,None]
    ))

@tf.function(jit_compile=jit_on)
def update_e_pml(Ex, Ey, Ez, Hx, Hy, Hz, gbx, gby, gbz, del_x, del_y, del_z,
                 psi_Eyx_1, psi_Eyx_2, psi_Ezx_1, psi_Ezx_2,
                 psi_Exy_1, psi_Exy_2, psi_Ezy_1, psi_Ezy_2,
                 psi_Exz_1, psi_Exz_2, psi_Eyz_1, psi_Eyz_2,
                 be_x_1, ce_x_1, be_x_2, ce_x_2,
                 be_y_1, ce_y_1, be_y_2, ce_y_2,
                 be_z_1, ce_z_1, be_z_2, ce_z_2,
                 nxPML_1, nyPML_1, nzPML_1, nxPML_2, nyPML_2, nzPML_2, 
                 x_size, y_size, z_size):
    
    # =================================================================
    # Ex UPDATES
    # =================================================================

    if type_sim!=0:
        # --- bottom y for Ex ---
        psi_Exy_1[:, 1:, 1:-1].assign(be_y_1[:-1, 1:, 1:-1] * psi_Exy_1[:, 1:, 1:-1] + ce_y_1[:-1, 1:, 1:-1] * (Hz[:-1, 1:nyPML_1, 1:-1] - Hz[:-1, :nyPML_1-1, 1:-1]) / del_y)
        Ex[:-1, 1:nyPML_1, 1:-1].assign(Ex[:-1, 1:nyPML_1, 1:-1] + gbx[:-1, 1:nyPML_1, 1:-1] * psi_Exy_1[:, 1:, 1:-1])

        # --- top y for Ex ---
        new_psi_Exy_2 = (
            be_y_2[:-1, 1:nyPML_2, 1:-1] * psi_Exy_2[:, 1:nyPML_2, 1:-1] + 
            ce_y_2[:-1, 1:nyPML_2, 1:-1] * (tf.reverse(Hz[:-1, y_size-nyPML_2 : y_size-1, 1:-1], [1]) - 
                                            tf.reverse(Hz[:-1, y_size-nyPML_2-1 : y_size-2, 1:-1], [1])) / del_y
        )
        psi_Exy_2[:, 1:nyPML_2, 1:-1].assign(new_psi_Exy_2)
        Ex[:-1, y_size-nyPML_2 : y_size-1, 1:-1].assign(Ex[:-1, y_size-nyPML_2 : y_size-1, 1:-1] + gbx[:-1, y_size-nyPML_2 : y_size-1, 1:-1] * tf.reverse(new_psi_Exy_2, [1]))

        # --- bottom z for Ex ---
        psi_Exz_1[:, 1:-1, 1:].assign(be_z_1[:-1, 1:-1, 1:] * psi_Exz_1[:, 1:-1, 1:] + ce_z_1[:-1, 1:-1, 1:] * (Hy[:-1, 1:-1, :nzPML_1-1] - Hy[:-1, 1:-1, 1:nzPML_1]) / del_z)
        Ex[:-1, 1:-1, 1:nzPML_1].assign(Ex[:-1, 1:-1, 1:nzPML_1] + gbx[:-1, 1:-1, 1:nzPML_1] * psi_Exz_1[:, 1:-1, 1:])

        # --- top z for Ex ---
        new_psi_Exz_2 = (
            be_z_2[:-1, 1:-1, 1:nzPML_2] * psi_Exz_2[:, 1:-1, 1:nzPML_2] + 
            ce_z_2[:-1, 1:-1, 1:nzPML_2] * (tf.reverse(Hy[:-1, 1:-1, z_size-nzPML_2-1 : z_size-2], [2]) - 
                                            tf.reverse(Hy[:-1, 1:-1, z_size-nzPML_2 : z_size-1], [2])) / del_z
        )
        psi_Exz_2[:, 1:-1, 1:nzPML_2].assign(new_psi_Exz_2)
        Ex[:-1, 1:-1, z_size-nzPML_2 : z_size-1].assign(Ex[:-1, 1:-1, z_size-nzPML_2 : z_size-1] + gbx[:-1, 1:-1, z_size-nzPML_2 : z_size-1] * tf.reverse(new_psi_Exz_2, [2]))

    # =================================================================
    # Ey UPDATES
    # =================================================================

    # --- bottom x for Ey ---
    psi_Eyx_1[1:, :, 1:-1].assign(be_x_1[1:, :-1, 1:-1] * psi_Eyx_1[1:, :, 1:-1] + ce_x_1[1:, :-1, 1:-1] * (Hz[:nxPML_1-1, :-1, 1:-1] - Hz[1:nxPML_1, :-1, 1:-1]) / del_x)
    Ey[1:nxPML_1, :-1, 1:-1].assign(Ey[1:nxPML_1, :-1, 1:-1] + gby[1:nxPML_1, :-1, 1:-1] * psi_Eyx_1[1:, :, 1:-1])

    # --- top x for Ey ---
    new_psi_Eyx_2 = (
        be_x_2[1:nxPML_2, :-1, 1:-1] * psi_Eyx_2[1:nxPML_2, :, 1:-1] + 
        ce_x_2[1:nxPML_2, :-1, 1:-1] * (tf.reverse(Hz[x_size-nxPML_2-1 : x_size-2, :-1, 1:-1], [0]) - 
                                        tf.reverse(Hz[x_size-nxPML_2 : x_size-1, :-1, 1:-1], [0])) / del_x
    )
    psi_Eyx_2[1:nxPML_2, :, 1:-1].assign(new_psi_Eyx_2)
    Ey[x_size-nxPML_2 : x_size-1, :-1, 1:-1].assign(Ey[x_size-nxPML_2 : x_size-1, :-1, 1:-1] + gby[x_size-nxPML_2 : x_size-1, :-1, 1:-1] * tf.reverse(new_psi_Eyx_2, [0]))

    if type_sim!=0:

        # --- bottom y for Ey ---
        psi_Eyz_1[1:-1, :, 1:].assign(be_z_1[1:-1, :-1, 1:] * psi_Eyz_1[1:-1, :, 1:] + ce_z_1[1:-1, :-1, 1:] * (Hx[1:-1, :-1, 1:nzPML_1] - Hx[1:-1, :-1, :nzPML_1-1]) / del_z)
        Ey[1:-1, :-1, 1:nzPML_1].assign(Ey[1:-1, :-1, 1:nzPML_1] + gby[1:-1, :-1, 1:nzPML_1] * psi_Eyz_1[1:-1, :, 1:])

        # --- top y for Ey ---
        new_psi_Eyz_2 = (
            be_z_2[1:-1, :-1, 1:nzPML_2] * psi_Eyz_2[1:-1, :, 1:nzPML_2] + 
            ce_z_2[1:-1, :-1, 1:nzPML_2] * (tf.reverse(Hx[1:-1, :-1, z_size-nzPML_2 : z_size-1], [2]) - 
                                            tf.reverse(Hx[1:-1, :-1, z_size-nzPML_2-1 : z_size-2], [2])) / del_z
        )
        psi_Eyz_2[1:-1, :, 1:nzPML_2].assign(new_psi_Eyz_2)
        Ey[1:-1, :-1, z_size-nzPML_2 : z_size-1].assign(Ey[1:-1, :-1, z_size-nzPML_2 : z_size-1] + gby[1:-1, :-1, z_size-nzPML_2 : z_size-1] * tf.reverse(new_psi_Eyz_2, [2]))

    # =================================================================
    # Ez UPDATES
    # =================================================================

    # --- bottom x for Ez ---
    psi_Ezx_1[1:, 1:-1, :-1].assign(be_x_1[1:, 1:-1, :-1] * psi_Ezx_1[1:, 1:-1, :-1] + ce_x_1[1:, 1:-1, :-1] * (Hy[1:nxPML_1, 1:-1, :-1] - Hy[:nxPML_1-1, 1:-1, :-1]) / del_x)
    Ez[1:nxPML_1, 1:-1, :-1].assign(Ez[1:nxPML_1, 1:-1, :-1] + gbz[1:nxPML_1, 1:-1, :-1] * psi_Ezx_1[1:, 1:-1, :-1])

    # --- top x for Ez ---
    new_psi_Ezx_2 = (
        be_x_2[1:nxPML_2, 1:-1, :-1] * psi_Ezx_2[1:nxPML_2, 1:-1, :-1] + 
        ce_x_2[1:nxPML_2, 1:-1, :-1] * (tf.reverse(Hy[x_size-nxPML_2 : x_size-1, 1:-1, :-1], [0]) - 
                                        tf.reverse(Hy[x_size-nxPML_2-1 : x_size-2, 1:-1, :-1], [0])) / del_x
    )
    psi_Ezx_2[1:nxPML_2, 1:-1, :-1].assign(new_psi_Ezx_2)
    Ez[x_size-nxPML_2 : x_size-1, 1:-1, :-1].assign(Ez[x_size-nxPML_2 : x_size-1, 1:-1, :-1] + gbz[x_size-nxPML_2 : x_size-1, 1:-1, :-1] * tf.reverse(new_psi_Ezx_2, [0]))

    if type_sim!=0:

        # --- bottom y for Ez ---
        psi_Ezy_1[1:-1, 1:, :-1].assign(be_y_1[1:-1, 1:, :-1] * psi_Ezy_1[1:-1, 1:, :-1] + ce_y_1[1:-1, 1:, :-1] * (Hx[1:-1, :nyPML_1-1, :-1] - Hx[1:-1, 1:nyPML_1, :-1]) / del_y)
        Ez[1:-1, 1:nyPML_1, :-1].assign(Ez[1:-1, 1:nyPML_1, :-1] + gbz[1:-1, 1:nyPML_1, :-1] * psi_Ezy_1[1:-1, 1:, :-1])

        # --- top y for Ez ---
        new_psi_Ezy_2 = (
            be_y_2[1:-1, 1:nyPML_2, :-1] * psi_Ezy_2[1:-1, 1:nyPML_2, :-1] + 
            ce_y_2[1:-1, 1:nyPML_2, :-1] * (tf.reverse(Hx[1:-1, y_size-nyPML_2-1 : y_size-2, :-1], [1]) - 
                                            tf.reverse(Hx[1:-1, y_size-nyPML_2 : y_size-1, :-1], [1])) / del_y
        )
        psi_Ezy_2[1:-1, 1:nyPML_2, :-1].assign(new_psi_Ezy_2)
        Ez[1:-1, y_size-nyPML_2 : y_size-1, :-1].assign(Ez[1:-1, y_size-nyPML_2 : y_size-1, :-1] + gbz[1:-1, y_size-nyPML_2 : y_size-1, :-1] * tf.reverse(new_psi_Ezy_2, [1]))

#don't compile, this is for testing only
def compute_source_pulse(counter, del_t, pulse_type, t_spread, spread):
    """
    Compute the source pulse value at a given time step.
    
    Inputs: counter, del_t, pulse_type, t_spread, spread
    Returns: pulse value (scalar)
    """
    import numpy as np
    t = counter * del_t
    if pulse_type == 1:
        pulse = np.exp(-0.5 * ((t_spread - t) / spread) ** 2)
    elif pulse_type == 2:
        pulse = -1 * ((t_spread - t) / spread) * np.exp(0.5) * \
                np.exp(-0.5 * ((t_spread - t) / spread) ** 2)
    Ez[40, 40, 40].assign(pulse)
    input_data[counter].assign(pulse)

@tf.function(jit_compile=jit_on)
def E_inc(weight, x, y, z,counter_val):
    
    time_arg = (counter_val - 1.0) * del_t
    spatial_phase = (
        (tf.sin(theta) * tf.cos(phi) * (x - x_delay) * del_x) +
        (tf.sin(theta) * tf.sin(phi) * (y - y_delay) * del_y) +
        (tf.cos(theta) * (z - z_delay) * del_z)
    ) / c
    term = (t_spread - time_arg + spatial_phase) / spread
    
    if pulse_type == 1:
        E_inc_val = weight * tf.exp(-0.5 * term**2)
    elif pulse_type == 2:
        E_inc_val = -1.0 * weight * term * tf.exp(tf.constant(0.5, dtype=x.dtype)) * tf.exp(-0.5 * term**2)
    else:
        E_inc_val = tf.zeros_like(x)
        
    return E_inc_val

@tf.function(jit_compile=jit_on)
def H_inc(weight, x, y, z,counter_val):
    
    time_arg = (counter_val - 0.5) * del_t
    spatial_phase = (
        (tf.sin(theta) * tf.cos(phi) * (x - x_delay) * del_x) +
        (tf.sin(theta) * tf.sin(phi) * (y - y_delay) * del_y) +
        (tf.cos(theta) * (z - z_delay) * del_z)
    ) / c
    term = (t_spread - time_arg + spatial_phase) / spread

    if pulse_type == 1:
        H_inc_val = weight * tf.exp(-0.5 * term**2)
    elif pulse_type == 2:
        H_inc_val = -1.0 * weight * term * tf.exp(tf.constant(0.5, dtype=x.dtype)) * tf.exp(-0.5 * term**2)
    else:
        H_inc_val = tf.zeros_like(x)

    return H_inc_val

@tf.function(jit_compile=jit_on)
def h_plane_waves(Hz, Hy, Hx, E_inc, WEx, WEy, WEz,
                  xlow, xhigh, ylow, yhigh, zlow, zhigh,
                  del_t, mu_0, del_x, del_y, del_z,
                  xlow_wall, xhigh_wall, ylow_wall, 
                  yhigh_wall, zlow_wall, zhigh_wall,counter_val):

    # Pre-compute common factors
    coeff_x = del_t / (mu_0 * del_x)
    coeff_y = del_t / (mu_0 * del_y)
    coeff_z = del_t / (mu_0 * del_z)
    
    if type_sim==1:
        # ===== Y FACES =====
        # First Y face loop: Hz updates
        i_range = tf.range(xlow, xhigh, dtype=tf.int32)
        k_range = tf.range(zlow, zhigh + 1, dtype=tf.int32)
        i_grid, k_grid = tf.meshgrid(i_range, k_range, indexing='ij')
        
        E_inc_vals = E_inc(WEx, tf.cast(i_grid, tf.float32) + 0.5, 
                        tf.cast(ylow, tf.float32) + 0.0, 
                        tf.cast(k_grid, tf.float32) + 0.0, counter_val)
        updates = -(coeff_y * E_inc_vals * ylow_wall)
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.fill([tf.size(i_grid)], ylow-1),
                        tf.reshape(k_grid, [-1])], axis=1)
        Hz.assign(tf.tensor_scatter_nd_add(Hz, indices, tf.reshape(updates, [-1])))
        
        E_inc_vals = E_inc(WEx, tf.cast(i_grid, tf.float32) + 0.5, 
                        tf.cast(yhigh, tf.float32) + 0.0, 
                        tf.cast(k_grid, tf.float32) + 0.0, counter_val)
        updates = coeff_y * E_inc_vals * yhigh_wall
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.fill([tf.size(i_grid)], yhigh),
                        tf.reshape(k_grid, [-1])], axis=1)
        Hz.assign(tf.tensor_scatter_nd_add(Hz, indices, tf.reshape(updates, [-1])))
        
        # Second Y face loop: Hx updates
        i_range = tf.range(xlow, xhigh + 1, dtype=tf.int32)
        k_range = tf.range(zlow, zhigh, dtype=tf.int32)
        i_grid, k_grid = tf.meshgrid(i_range, k_range, indexing='ij')
        
        E_inc_vals = E_inc(WEz, tf.cast(i_grid, tf.float32) + 0.0, 
                        tf.cast(ylow, tf.float32) + 0.0, 
                        tf.cast(k_grid, tf.float32) + 0.5, counter_val)
        updates = coeff_y * E_inc_vals * ylow_wall
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.fill([tf.size(i_grid)], ylow-1),
                        tf.reshape(k_grid, [-1])], axis=1)
        Hx.assign(tf.tensor_scatter_nd_add(Hx, indices, tf.reshape(updates, [-1])))
        
        E_inc_vals = E_inc(WEz, tf.cast(i_grid, tf.float32) + 0.0, 
                        tf.cast(yhigh, tf.float32) + 0.0, 
                        tf.cast(k_grid, tf.float32) + 0.5, counter_val)
        updates = -(coeff_y * E_inc_vals * yhigh_wall)
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.fill([tf.size(i_grid)], yhigh),
                        tf.reshape(k_grid, [-1])], axis=1)
        Hx.assign(tf.tensor_scatter_nd_add(Hx, indices, tf.reshape(updates, [-1])))
        
        # ===== Z FACES =====
        # First Z face loop: Hy updates
        i_range = tf.range(xlow, xhigh, dtype=tf.int32)
        j_range = tf.range(ylow, yhigh + 1, dtype=tf.int32)
        i_grid, j_grid = tf.meshgrid(i_range, j_range, indexing='ij')
        
        E_inc_vals = E_inc(WEx, tf.cast(i_grid, tf.float32) + 0.5, 
                        tf.cast(j_grid, tf.float32) + 0.0, 
                        tf.cast(zlow, tf.float32) + 0.0, counter_val)
        updates = coeff_z * E_inc_vals * zlow_wall
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.reshape(j_grid, [-1]),
                        tf.fill([tf.size(i_grid)], zlow-1)], axis=1)
        Hy.assign(tf.tensor_scatter_nd_add(Hy, indices, tf.reshape(updates, [-1])))
        
        E_inc_vals = E_inc(WEx, tf.cast(i_grid, tf.float32) + 0.5, 
                        tf.cast(j_grid, tf.float32) + 0.0, 
                        tf.cast(zhigh, tf.float32) + 0.0, counter_val)
        updates = -(coeff_z * E_inc_vals * zhigh_wall)
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.reshape(j_grid, [-1]),
                        tf.fill([tf.size(i_grid)], zhigh)], axis=1)
        Hy.assign(tf.tensor_scatter_nd_add(Hy, indices, tf.reshape(updates, [-1])))
        
        # Second Z face loop: Hx updates
        i_range = tf.range(xlow, xhigh + 1, dtype=tf.int32)
        j_range = tf.range(ylow, yhigh, dtype=tf.int32)
        i_grid, j_grid = tf.meshgrid(i_range, j_range, indexing='ij')
        
        E_inc_vals = E_inc(WEy, tf.cast(i_grid, tf.float32) + 0.0, 
                        tf.cast(j_grid, tf.float32) + 0.5, 
                        tf.cast(zlow, tf.float32) + 0.0, counter_val)
        updates = -(coeff_z * E_inc_vals * zlow_wall)
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.reshape(j_grid, [-1]),
                        tf.fill([tf.size(i_grid)], zlow-1)], axis=1)
        Hx.assign(tf.tensor_scatter_nd_add(Hx, indices, tf.reshape(updates, [-1])))
        
        E_inc_vals = E_inc(WEy, tf.cast(i_grid, tf.float32) + 0.0, 
                        tf.cast(j_grid, tf.float32) + 0.5, 
                        tf.cast(zhigh, tf.float32) + 0.0, counter_val)
        updates = coeff_z * E_inc_vals * zhigh_wall
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.reshape(j_grid, [-1]),
                        tf.fill([tf.size(i_grid)], zhigh)], axis=1)
        Hx.assign(tf.tensor_scatter_nd_add(Hx, indices, tf.reshape(updates, [-1])))
    
    # ===== X FACES =====
    # First X face loop: Hz updates
    j_range = tf.range(ylow, yhigh, dtype=tf.int32)
    k_range = tf.range(zlow, zhigh + 1, dtype=tf.int32)
    j_grid, k_grid = tf.meshgrid(j_range, k_range, indexing='ij')
    
    E_inc_vals = E_inc(WEy, tf.cast(xlow, tf.float32) + 0.0, 
                       tf.cast(j_grid, tf.float32) + 0.5, 
                       tf.cast(k_grid, tf.float32) + 0.0, counter_val)
    updates = coeff_x * E_inc_vals * xlow_wall
    indices = tf.stack([tf.fill([tf.size(j_grid)], xlow-1),
                       tf.reshape(j_grid, [-1]),
                       tf.reshape(k_grid, [-1])], axis=1)
    Hz.assign(tf.tensor_scatter_nd_add(Hz, indices, tf.reshape(updates, [-1])))
    
    E_inc_vals = E_inc(WEy, tf.cast(xhigh, tf.float32) + 0.0, 
                       tf.cast(j_grid, tf.float32) + 0.5, 
                       tf.cast(k_grid, tf.float32) + 0.0, counter_val)
    updates = -(coeff_x * E_inc_vals * xhigh_wall)
    indices = tf.stack([tf.fill([tf.size(j_grid)], xhigh),
                       tf.reshape(j_grid, [-1]),
                       tf.reshape(k_grid, [-1])], axis=1)
    Hz.assign(tf.tensor_scatter_nd_add(Hz, indices, tf.reshape(updates, [-1])))
    
    # Second X face loop: Hy updates
    j_range = tf.range(ylow, yhigh + 1, dtype=tf.int32)
    k_range = tf.range(zlow, zhigh, dtype=tf.int32)
    j_grid, k_grid = tf.meshgrid(j_range, k_range, indexing='ij')
    
    E_inc_vals = E_inc(WEz, tf.cast(xlow, tf.float32) + 0.0, 
                       tf.cast(j_grid, tf.float32) + 0.0, 
                       tf.cast(k_grid, tf.float32) + 0.5, counter_val)
    updates = -(coeff_x * E_inc_vals * xlow_wall)
    indices = tf.stack([tf.fill([tf.size(j_grid)], xlow-1),
                       tf.reshape(j_grid, [-1]),
                       tf.reshape(k_grid, [-1])], axis=1)
    Hy.assign(tf.tensor_scatter_nd_add(Hy, indices, tf.reshape(updates, [-1])))
    
    E_inc_vals = E_inc(WEz, tf.cast(xhigh, tf.float32) + 0.0, 
                       tf.cast(j_grid, tf.float32) + 0.0, 
                       tf.cast(k_grid, tf.float32) + 0.5, counter_val)
    updates = coeff_x * E_inc_vals * xhigh_wall
    indices = tf.stack([tf.fill([tf.size(j_grid)], xhigh),
                       tf.reshape(j_grid, [-1]),
                       tf.reshape(k_grid, [-1])], axis=1)
    Hy.assign(tf.tensor_scatter_nd_add(Hy, indices, tf.reshape(updates, [-1])))

@tf.function(jit_compile=jit_on)
def e_plane_waves(Ez, Ey, Ex, H_inc, WHx, WHy, WHz,
                  xlow, xhigh, ylow, yhigh, zlow, zhigh,
                  del_t, ep_0, del_x, del_y, del_z,
                  xlow_wall, xhigh_wall, ylow_wall, 
                  yhigh_wall, zlow_wall, zhigh_wall,counter_val):
    
    # Pre-compute common factors
    coeff_x = del_t / (ep_0 * del_x)
    coeff_y = del_t / (ep_0 * del_y)
    coeff_z = del_t / (ep_0 * del_z)
    
    if type_sim==1:
        # ===== Y FACES =====
        # First Y face loop: Ez updates
        i_range = tf.range(xlow, xhigh + 1, dtype=tf.int32)
        k_range = tf.range(zlow, zhigh, dtype=tf.int32)
        i_grid, k_grid = tf.meshgrid(i_range, k_range, indexing='ij')
        
        H_inc_vals = H_inc(WHx, tf.cast(i_grid, tf.float32) + 0.0, 
                        tf.cast(ylow, tf.float32) - 0.5, 
                        tf.cast(k_grid, tf.float32) + 0.5, counter_val)
        updates = coeff_y * H_inc_vals * ylow_wall
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.fill([tf.size(i_grid)], ylow),
                        tf.reshape(k_grid, [-1])], axis=1)
        Ez.assign(tf.tensor_scatter_nd_add(Ez, indices, tf.reshape(updates, [-1])))
        
        H_inc_vals = H_inc(WHx, tf.cast(i_grid, tf.float32) + 0.0, 
                        tf.cast(yhigh, tf.float32) + 0.5, 
                        tf.cast(k_grid, tf.float32) + 0.5, counter_val)
        updates = -(coeff_y * H_inc_vals * yhigh_wall)
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.fill([tf.size(i_grid)], yhigh),
                        tf.reshape(k_grid, [-1])], axis=1)
        Ez.assign(tf.tensor_scatter_nd_add(Ez, indices, tf.reshape(updates, [-1])))
        
        # Second Y face loop: Ex updates
        i_range = tf.range(xlow, xhigh, dtype=tf.int32)
        k_range = tf.range(zlow, zhigh + 1, dtype=tf.int32)
        i_grid, k_grid = tf.meshgrid(i_range, k_range, indexing='ij')
        
        H_inc_vals = H_inc(WHz, tf.cast(i_grid, tf.float32) + 0.5, 
                        tf.cast(ylow, tf.float32) - 0.5, 
                        tf.cast(k_grid, tf.float32) + 0.0, counter_val)
        updates = -(coeff_y * H_inc_vals * ylow_wall)
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.fill([tf.size(i_grid)], ylow),
                        tf.reshape(k_grid, [-1])], axis=1)
        Ex.assign(tf.tensor_scatter_nd_add(Ex, indices, tf.reshape(updates, [-1])))
        
        H_inc_vals = H_inc(WHz, tf.cast(i_grid, tf.float32) + 0.5, 
                        tf.cast(yhigh, tf.float32) + 0.5, 
                        tf.cast(k_grid, tf.float32) + 0.0, counter_val)
        updates = coeff_y * H_inc_vals * yhigh_wall
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.fill([tf.size(i_grid)], yhigh),
                        tf.reshape(k_grid, [-1])], axis=1)
        Ex.assign(tf.tensor_scatter_nd_add(Ex, indices, tf.reshape(updates, [-1])))
        
        # ===== Z FACES =====
        # First Z face loop: Ey updates
        i_range = tf.range(xlow, xhigh + 1, dtype=tf.int32)
        j_range = tf.range(ylow, yhigh, dtype=tf.int32)
        i_grid, j_grid = tf.meshgrid(i_range, j_range, indexing='ij')
        
        H_inc_vals = H_inc(WHx, tf.cast(i_grid, tf.float32) + 0.0, 
                        tf.cast(j_grid, tf.float32) + 0.5, 
                        tf.cast(zlow, tf.float32) - 0.5, counter_val)
        updates = -(coeff_z * H_inc_vals * zlow_wall)
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.reshape(j_grid, [-1]),
                        tf.fill([tf.size(i_grid)], zlow)], axis=1)
        Ey.assign(tf.tensor_scatter_nd_add(Ey, indices, tf.reshape(updates, [-1])))
        
        H_inc_vals = H_inc(WHx, tf.cast(i_grid, tf.float32) + 0.0, 
                        tf.cast(j_grid, tf.float32) + 0.5, 
                        tf.cast(zhigh, tf.float32) + 0.5, counter_val)
        updates = coeff_z * H_inc_vals * zhigh_wall
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.reshape(j_grid, [-1]),
                        tf.fill([tf.size(i_grid)], zhigh)], axis=1)
        Ey.assign(tf.tensor_scatter_nd_add(Ey, indices, tf.reshape(updates, [-1])))
        
        # Second Z face loop: Ex updates
        i_range = tf.range(xlow, xhigh, dtype=tf.int32)
        j_range = tf.range(ylow, yhigh + 1, dtype=tf.int32)
        i_grid, j_grid = tf.meshgrid(i_range, j_range, indexing='ij')
        
        H_inc_vals = H_inc(WHy, tf.cast(i_grid, tf.float32) + 0.5, 
                        tf.cast(j_grid, tf.float32) + 0.0, 
                        tf.cast(zlow, tf.float32) - 0.5, counter_val)
        updates = coeff_z * H_inc_vals * zlow_wall
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.reshape(j_grid, [-1]),
                        tf.fill([tf.size(i_grid)], zlow)], axis=1)
        Ex.assign(tf.tensor_scatter_nd_add(Ex, indices, tf.reshape(updates, [-1])))
        
        H_inc_vals = H_inc(WHy, tf.cast(i_grid, tf.float32) + 0.5, 
                        tf.cast(j_grid, tf.float32) + 0.0, 
                        tf.cast(zhigh, tf.float32) + 0.5, counter_val)
        updates = -(coeff_z * H_inc_vals * zhigh_wall)
        indices = tf.stack([tf.reshape(i_grid, [-1]), 
                        tf.reshape(j_grid, [-1]),
                        tf.fill([tf.size(i_grid)], zhigh)], axis=1)
        Ex.assign(tf.tensor_scatter_nd_add(Ex, indices, tf.reshape(updates, [-1])))
        
    # ===== X FACES =====
    # First X face loop: Ez updates
    j_range = tf.range(ylow, yhigh + 1, dtype=tf.int32)
    k_range = tf.range(zlow, zhigh, dtype=tf.int32)
    j_grid, k_grid = tf.meshgrid(j_range, k_range, indexing='ij')
    
    H_inc_vals = H_inc(WHy, tf.cast(xlow, tf.float32) - 0.5, 
                       tf.cast(j_grid, tf.float32) + 0.0, 
                       tf.cast(k_grid, tf.float32) + 0.5, counter_val)
    updates = -(coeff_x * H_inc_vals * xlow_wall)
    indices = tf.stack([tf.fill([tf.size(j_grid)], xlow),
                       tf.reshape(j_grid, [-1]),
                       tf.reshape(k_grid, [-1])], axis=1)
    Ez.assign(tf.tensor_scatter_nd_add(Ez, indices, tf.reshape(updates, [-1])))
    
    H_inc_vals = H_inc(WHy, tf.cast(xhigh, tf.float32) + 0.5, 
                       tf.cast(j_grid, tf.float32) + 0.0, 
                       tf.cast(k_grid, tf.float32) + 0.5, counter_val)
    updates = coeff_x * H_inc_vals * xhigh_wall
    indices = tf.stack([tf.fill([tf.size(j_grid)], xhigh),
                       tf.reshape(j_grid, [-1]),
                       tf.reshape(k_grid, [-1])], axis=1)
    Ez.assign(tf.tensor_scatter_nd_add(Ez, indices, tf.reshape(updates, [-1])))
    
    # Second X face loop: Ey updates
    j_range = tf.range(ylow, yhigh, dtype=tf.int32)
    k_range = tf.range(zlow, zhigh + 1, dtype=tf.int32)
    j_grid, k_grid = tf.meshgrid(j_range, k_range, indexing='ij')
    
    H_inc_vals = H_inc(WHz, tf.cast(xlow, tf.float32) - 0.5, 
                       tf.cast(j_grid, tf.float32) + 0.5, 
                       tf.cast(k_grid, tf.float32) + 0.0, counter_val)
    updates = coeff_x * H_inc_vals * xlow_wall
    indices = tf.stack([tf.fill([tf.size(j_grid)], xlow),
                       tf.reshape(j_grid, [-1]),
                       tf.reshape(k_grid, [-1])], axis=1)
    Ey.assign(tf.tensor_scatter_nd_add(Ey, indices, tf.reshape(updates, [-1])))
    
    H_inc_vals = H_inc(WHz, tf.cast(xhigh, tf.float32) + 0.5, 
                       tf.cast(j_grid, tf.float32) + 0.5, 
                       tf.cast(k_grid, tf.float32) + 0.0, counter_val)
    updates = -(coeff_x * H_inc_vals * xhigh_wall)
    indices = tf.stack([tf.fill([tf.size(j_grid)], xhigh),
                       tf.reshape(j_grid, [-1]),
                       tf.reshape(k_grid, [-1])], axis=1)
    Ey.assign(tf.tensor_scatter_nd_add(Ey, indices, tf.reshape(updates, [-1])))

@tf.function(jit_compile=jit_on)
def update_sheets_h(Ex, Ex_special,
                    Ey, Ey_special,
                    Ez, Ez_special,
                    Hx, Hy, Hz,
                    den_hx, den_hy, den_hz,
                    sheet_thickness,
                    del_x, del_y, del_z,
                    mask_x, mask_y, mask_z,
                    db):

    # -------------------------------
    # X-normal sheets update
    # -------------------------------
    if tf.reduce_any(tf.cast(mask_x, tf.bool)):
        # Hy update
        delta_Hy = (
            db * (sheet_thickness / del_x) *
            mask_x[:-1, :-1, :-1] *
            (
                Ex[:-1, :-1, 1:] - Ex[:-1, :-1, :-1] +
                Ex_special[:-1, :-1, :-1] - Ex_special[:-1, :-1, 1:]
            ) *
            den_hz[None, None, :]
        )
        Hy.assign_add(tf.pad(delta_Hy, [[0, 1], [0, 1], [0, 1]]))

        # Hz update
        delta_Hz = (
            db * (sheet_thickness / del_x) *
            mask_x[:-1, :-1, :-1] *
            (
                Ex[:-1, :-1, :-1] - Ex[:-1, 1:, :-1] +
                Ex_special[:-1, 1:, :-1] - Ex_special[:-1, :-1, :-1]
            ) *
            den_hy[None, :, None]
        )
        Hz.assign_add(tf.pad(delta_Hz, [[0, 1], [0, 1], [0, 1]]))

    # -------------------------------
    # Y-normal sheets update
    # -------------------------------
    if tf.reduce_any(tf.cast(mask_y, tf.bool)):
        # Hx update
        delta_Hx = (
            db * (sheet_thickness / del_y) *
            mask_y[:-1, :-1, :-1] *
            (
                Ey[:-1, :-1, :-1] - Ey[:-1, :-1, 1:] +
                Ey_special[:-1, :-1, 1:] - Ey_special[:-1, :-1, :-1]
            ) *
            den_hz[None, None, :]
        )
        Hx.assign_add(tf.pad(delta_Hx, [[0, 1], [0, 1], [0, 1]]))

        # Hz update
        delta_Hz_y = (
            db * (sheet_thickness / del_y) *
            mask_y[:-1, :-1, :-1] *
            (
                Ey[1:, :-1, :-1] - Ey[:-1, :-1, :-1] -
                Ey_special[1:, :-1, :-1] + Ey_special[:-1, :-1, :-1]
            ) *
            den_hx[:, None, None]
        )
        Hz.assign_add(tf.pad(delta_Hz_y, [[0, 1], [0, 1], [0, 1]]))

    # -------------------------------
    # Z-normal sheets update
    # -------------------------------
    if tf.reduce_any(tf.cast(mask_z, tf.bool)):
        # Hx update
        delta_Hx_z = (
            db * (sheet_thickness / del_z) *
            mask_z[:-1, :-1, :-1] *
            (
                Ez[:-1, 1:, :-1] - Ez[:-1, :-1, :-1] -
                Ez_special[:-1, 1:, :-1] + Ez_special[:-1, :-1, :-1]
            ) *
            den_hy[None, :, None]
        )
        Hx.assign_add(tf.pad(delta_Hx_z, [[0, 1], [0, 1], [0, 1]]))

        # Hy update
        delta_Hy_z = (
            db * (sheet_thickness / del_z) *
            mask_z[:-1, :-1, :-1] *
            (
                Ez[:-1, :-1, :-1] - Ez[1:, :-1, :-1] + # Fixed index logic here
                Ez_special[1:, :-1, :-1] - Ez_special[:-1, :-1, :-1]
            ) *
            den_hx[:, None, None]
        )
        Hy.assign_add(tf.pad(delta_Hy_z, [[0, 1], [0, 1], [0, 1]]))

@tf.function(jit_compile=jit_on)
def update_sheets_e(Ex_special,
                    Ey_special,
                    Ez_special,
                    Hx, Hy, Hz,
                    den_ex, den_ey, den_ez,
                    del_t, ep_0,
                    mask_x, mask_y, mask_z,
                    sheet_sig_x_x, sheet_ep_x_x,
                    sheet_sig_y_y, sheet_ep_y_y,
                    sheet_sig_z_z, sheet_ep_z_z):

    # -------------------------------
    # X-normal sheets: Ex_special
    # -------------------------------
    if tf.reduce_any(tf.cast(mask_x, tf.bool)):
        sig = sheet_sig_x_x[:-1, 1:-1, 1:-1]
        ep_r = sheet_ep_x_x[:-1, 1:-1, 1:-1]
        
        coef_loss_minus = (1.0 - sig * del_t / (2.0 * ep_0 * ep_r))
        coef_loss_plus = (1.0 + sig * del_t / (2.0 * ep_0 * ep_r))
        
        update_x = (
            (coef_loss_minus * Ex_special[:-1, 1:-1, 1:-1] +
            (del_t / (ep_0 * ep_r)) * (
                (Hz[:-1, 1:-1, 1:-1] - Hz[:-1, :-2, 1:-1]) * den_ey[None, 1:, None] +
                (Hy[:-1, 1:-1, :-2] - Hy[:-1, 1:-1, 1:-1]) * den_ez[None, None, 1:]
            )) / coef_loss_plus
        ) * mask_x[:-1, 1:-1, 1:-1]
        
        Ex_special.assign(tf.pad(update_x, [[0, 1], [1, 1], [1, 1]]))

    # -------------------------------
    # Y-normal sheets: Ey_special
    # -------------------------------
    if tf.reduce_any(tf.cast(mask_y, tf.bool)):
        sig = sheet_sig_y_y[1:-1, :-1, 1:-1]
        ep_r = sheet_ep_y_y[1:-1, :-1, 1:-1]
        
        coef_loss_minus = (1.0 - sig * del_t / (2.0 * ep_0 * ep_r))
        coef_loss_plus = (1.0 + sig * del_t / (2.0 * ep_0 * ep_r))

        update_y = (
            (coef_loss_minus * Ey_special[1:-1, :-1, 1:-1] +
            (del_t / (ep_0 * ep_r)) * (
                (Hz[:-2, :-1, 1:-1] - Hz[1:-1, :-1, 1:-1]) * den_ex[1:, None, None] +
                (Hx[1:-1, :-1, 1:-1] - Hx[1:-1, :-1, :-2]) * den_ez[None, None, 1:]
            )) / coef_loss_plus
        ) * mask_y[1:-1, :-1, 1:-1]

        Ey_special.assign(tf.pad(update_y, [[1, 1], [0, 1], [1, 1]]))

    # -------------------------------
    # Z-normal sheets: Ez_special
    # -------------------------------
    if tf.reduce_any(tf.cast(mask_z, tf.bool)):
        sig = sheet_sig_z_z[1:-1, 1:-1, :-1]
        ep_r = sheet_ep_z_z[1:-1, 1:-1, :-1]
        
        coef_loss_minus = (1.0 - sig * del_t / (2.0 * ep_0 * ep_r))
        coef_loss_plus = (1.0 + sig * del_t / (2.0 * ep_0 * ep_r))

        update_z = (
            (coef_loss_minus * Ez_special[1:-1, 1:-1, :-1] +
            (del_t / (ep_0 * ep_r)) * (
                (Hy[1:-1, 1:-1, :-1] - Hy[:-2, 1:-1, :-1]) * den_ex[1:, None, None] +
                (Hx[1:-1, :-2, :-1] - Hx[1:-1, 1:-1, :-1]) * den_ey[None, 1:, None]
            )) / coef_loss_plus
        ) * mask_z[1:-1, 1:-1, :-1]

        Ez_special.assign(tf.pad(update_z, [[1, 1], [1, 1], [0, 1]]))

@tf.function(jit_compile=jit_on)
def h_pbc(Hx,Hy,Hz):
    Hz[:,-2,:].assign(Hz[:,0,:])
    Hx[:,-2,:].assign(Hx[:,0,:])

    Hy[:,:,-2].assign(Hy[:,:,0])
    Hx[:,:,-2].assign(Hx[:,:,0])

@tf.function(jit_compile=jit_on)
def e_pbc(Ex,Ey,Ez):
    Ez[:,0,:].assign(Ez[:,-2,:])
    Ex[:,0,:].assign(Ex[:,-2,:])
    
    Ey[:,:,0].assign(Ey[:,:,-2])
    Ex[:,:,0].assign(Ex[:,:,-2])

##########################################################
# MAIN FDTD SOLVER
##########################################################
print("\nStarting FDTD solver...")
with GPUProfiler("FDTD Solver", device):
    step = tf.Variable(0.0, dtype=tf.float32)
    # Start time loop
    for counter in range(time_steps):
        step.assign(counter)
        
        # Update H fields (using compiled version)
        update_h_fields(Hx, Hy, Hz, Ex, Ey, Ez, da, db, den_hx, den_hy, den_hz)

        # update sheets for H fields if applicable
        update_sheets_h(Ex, Ex_special,
            Ey, Ey_special,
            Ez, Ez_special,
            Hx, Hy, Hz,
            den_hx, den_hy, den_hz,
            sheet_thickness,
            del_x, del_y, del_z,
            mask_x, mask_y, mask_z,
            db)

        if type_sim==1:
            # Add H field contributions to plane waves
            h_plane_waves(Hz, Hy, Hx, E_inc, WEx, WEy, WEz,
                        xlow, xhigh, ylow, yhigh, zlow, zhigh,
                        del_t, mu_0, del_x, del_y, del_z,
                        xlow_wall, xhigh_wall, ylow_wall, 
                        yhigh_wall, zlow_wall, zhigh_wall,step)
        if type_sim==0:
            # Add H field contributions to plane waves
            h_plane_waves(Hz, Hy, Hx, E_inc, 0, 0, 0,
                        xlow, xhigh, ylow, yhigh, zlow, zhigh+1,
                        del_t, mu_0, del_x, del_y, del_z,
                        xlow_wall, xhigh_wall, ylow_wall, 
                        yhigh_wall, zlow_wall, zhigh_wall,step) 
            h_plane_waves(Hz, Hy, Hx, E_inc, 0, 0, 0,
                        xlow, xhigh, ylow, yhigh+1, zlow, zhigh,
                        del_t, mu_0, del_x, del_y, del_z,
                        xlow_wall, xhigh_wall, ylow_wall, 
                        yhigh_wall, zlow_wall, zhigh_wall,step)  
 
        # Apply PML to H fields (using compiled version)
        update_h_pml(Hx, Hy, Hz, Ex, Ey, Ez, db, del_x, del_y, del_z,
                 psi_Hyx_1, psi_Hyx_2, psi_Hzx_1, psi_Hzx_2,
                 psi_Hxy_1, psi_Hxy_2, psi_Hzy_1, psi_Hzy_2,
                 psi_Hxz_1, psi_Hxz_2, psi_Hyz_1, psi_Hyz_2,
                 bh_x_1, ch_x_1, bh_x_2, ch_x_2,
                 bh_y_1, ch_y_1, bh_y_2, ch_y_2,
                 bh_z_1, ch_z_1, bh_z_2, ch_z_2,
                 nxPML_1, nyPML_1, nzPML_1, nxPML_2, nyPML_2, nzPML_2,
                 x_size, y_size, z_size)
        
        #apply pbc condition (if all pml this is still fine)
        h_pbc(Hx,Hy,Hz)
        
        # Update E fields (using compiled version)
        update_e_fields(Ex, Ey, Ez, Hx, Hy, Hz, gax, gay, gaz, gbx, gby, gbz,
                       den_ex, den_ey, den_ez)

        # update sheets for E fields if applicable
        update_sheets_e(Ex_special, Ey_special, Ez_special,
                Hx, Hy, Hz,
                den_ex, den_ey, den_ez,
                del_t, ep_0,
                mask_x, mask_y, mask_z,
                sheet_sig_x_x, sheet_ep_x_x,
                sheet_sig_y_y, sheet_ep_y_y,
                sheet_sig_z_z, sheet_ep_z_z
            )

        if type_sim==1:
            # Add E field plane wave contributions
            e_plane_waves(Ez, Ey, Ex, H_inc, WHx, WHy, WHz,
                        xlow, xhigh, ylow, yhigh, zlow, zhigh,
                        del_t, ep_0, del_x, del_y, del_z,
                        xlow_wall, xhigh_wall, ylow_wall, 
                        yhigh_wall, zlow_wall, zhigh_wall,step)
        if type_sim==0:
            e_plane_waves(Ez, Ey, Ex, H_inc, 0, WHy, 0,
                        xlow, xhigh, ylow, yhigh-1, zlow, zhigh-1,
                        del_t, ep_0, del_x, del_y, del_z,
                        xlow_wall, xhigh_wall, ylow_wall, 
                        yhigh_wall, zlow_wall, zhigh_wall,step)
        
        # Apply PML to E fields (using compiled version)
        update_e_pml(Ex, Ey, Ez, Hx, Hy, Hz, gbx, gby, gbz, del_x, del_y, del_z,
                 psi_Eyx_1, psi_Eyx_2, psi_Ezx_1, psi_Ezx_2,
                 psi_Exy_1, psi_Exy_2, psi_Ezy_1, psi_Ezy_2,
                 psi_Exz_1, psi_Exz_2, psi_Eyz_1, psi_Eyz_2,
                 be_x_1, ce_x_1, be_x_2, ce_x_2,
                 be_y_1, ce_y_1, be_y_2, ce_y_2,
                 be_z_1, ce_z_1, be_z_2, ce_z_2,
                 nxPML_1, nyPML_1, nzPML_1, nxPML_2, nyPML_2, nzPML_2, 
                 x_size, y_size, z_size)
        
        #apply pbc condition (if all pml this is still fine)
        e_pbc(Hx,Hy,Hz)
        
        # Apply source as a test for PML only
        #compute_source_pulse(counter, del_t, pulse_type, t_spread, spread)

        
        # Progress tracking
        if counter % 100 == 0:
            print(f"{counter} of {time_steps} time steps")

##########################################################
# POST PROCESSING
##########################################################
with GPUProfiler("Post Processing", device):
    import numpy as np
    # Transfer only final results to CPU for plotting
    E_test = Ez.numpy()
    #input_cpu = input_data.cpu().numpy()
    
    # Example: Plot a slice
    x=np.linspace(1,x_size,x_size)
    y=np.linspace(1,y_size,y_size)
    z=np.linspace(1,z_size,z_size)
    X,Y=np.meshgrid(y,x)
    #use this one for pulse checking for pml
    #plt.pcolormesh(X,Y,10*np.log10(np.abs(E_test[:,:,40])), cmap='bwr',vmin=-90, vmax=-30, shading='gouraud')
    #use this one for plane wave checking
    plt.pcolormesh(E_test[:,:,15], cmap='bwr', shading='gouraud')
    plt.colorbar()
    #plt.title('Ez field at z=40')
    plt.savefig('fdtd_result00.png')
    # plt.close()

print("\n" + "="*60)
print("Simulation complete#")
print("="*60)

# Check GPU utilization
if 'GPU' in device:
    # TensorFlow doesn't have as direct memory reporting as PyTorch
    # Use tf.config.experimental.get_memory_info for similar functionality
    try:
        memory_info = tf.config.experimental.get_memory_info(device)
        print(f"\nCurrent GPU memory usage: {memory_info['current'] / 1e9:.2f} GB")
        print(f"Peak GPU memory usage: {memory_info['peak'] / 1e9:.2f} GB")
    except:
        # Fallback if get_memory_info is not available
        print(f"\nGPU memory info not available. Use nvidia-smi for detailed memory usage.")
        print("Run: !nvidia-smi")