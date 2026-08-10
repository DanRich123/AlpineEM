import numpy as np
import os
from pyevtk.hl import gridToVTK
import shutil

# ----------------------------
# User setup
# ----------------------------
#inputs file name for read-in
inputs_filename = 'inputs.txt'
#output directory name for VTR files
output_dir = "unit cell design"
#main VTM output file name
filename = 'main geometry' # main VTM file
#label if using spice or kmax versions
spice = False
kmax = True

# ----------------------------
# Setup directory
# ----------------------------
# move to the current directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
# Check if the output directory exists, then delete it and its contents if it does
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
# Create the fresh directory
os.makedirs(output_dir)
# check if the output vtm exists, if so remove it too before I make the new one below
vtm_file = os.path.join(script_dir, "{}.vtm".format(filename))
if os.path.isfile(vtm_file):
    os.remove(vtm_file)

#----------------------------
# Read-in fdtd inputs file
#----------------------------
# most of this is unused as the actual geometry is drawn using fortran outputs
# we need to interate through the entire file to get the info we need:
# filename to know what binary is named
# sizes (nx,ny,nz)
# pbc vs pml for type simulation
# ports that we do need
# IGP yes or no
with open(inputs_filename, "r") as file:
    # Read in necessary info here
    nx,ny,nz = map(int, file.readline().split(','))
    step_size_x, step_size_y, step_size_z = map(float, file.readline().split(','))
    time_steps = int(file.readline())
    time_reduction_factor = float(file.readline())
    f_center = float(file.readline())
    pulse_type = int(file.readline())
    b_xlow,b_xhigh,b_ylow,b_yhigh,b_zlow,b_zhigh= map(float, file.readline().split(','))

    igp_on=file.readline().strip()
    if igp_on=='IGP':
        igp_value = file.readline().strip()
        igp_location=int(file.readline())

    if kmax==False:
        excitation_type=file.readline().strip()
        if excitation_type=='antenna':
            antenna_amp, antenna_time_delay = map(float, file.readline().split(','))
            port_excited = int(file.readline())
        if excitation_type=='plane wave':
            plane_wave_amp, plane_wave_time_delay = map(float, file.readline().split(','))
            theta, phi, pol = map(float, file.readline().split(','))

    if kmax==True:
        excitation_type=file.readline().strip()
        if excitation_type=='plane wave':
            plane_wave_amp, plane_wave_time_delay = map(float, file.readline().split(','))
            TE_or_TM=file.readline()
            plus_minus=file.readline()
        if excitation_type=='antenna':
            antenna_amp, antenna_time_delay = map(float, file.readline().split(','))
            port_excited = int(file.readline())
        k1,k2=map(float, file.readline().split(','))

    # Read and convert far-field angles
    num_far_field_angles = int(file.readline())
    for _ in range(num_far_field_angles):
        angle_theta, angle_phi = map(float, file.readline().split(','))

    # Read video settings
    video_on = file.readline().strip()
    if video_on=='no':
        video_on=0
    if video_on=='yes':
        video_on=1
    if video_on == 1:
        slice_value = file.readline().strip()
        slice_location = int(file.readline())

    filename_fdtd = file.readline().strip()

    # Read bulk and sheet material properties
    num_total_materials = int(file.readline()) 
    for _ in range(num_total_materials):
        name=file.readline().strip()
        if name=='volume':
            bulk_mat_type=file.readline().strip()
            materials_id = int(file.readline())
            properties_1 = list(map(float, file.readline().split(',')))
            properties_2 = list(map(float, file.readline().split(',')))
            if bulk_mat_type=='plasma':
                num_poles=int(file.readline())
                for i in range(num_poles):
                    plasma_properties1=list(map(float, file.readline().split(',')))
                    plasma_properties2=list(map(float, file.readline().split(',')))
        if name=='sheet':
            sheet_materials_id = int(file.readline())
            props_1 = float(file.readline())
            #props_2 = list(map(float, file.readline().split(',')))
            #props_3 = list(map(float, file.readline().split(',')))

    binary_fdtd_bulkop_filename_yesno = file.readline().strip()
    if (binary_fdtd_bulkop_filename_yesno=='yes'):
        binary_fdtd_bulkop_filename=file.readline().strip()

    # Read bulk blocks, spheres, and cylinders - cylinders to come
    num_total_objects=int(file.readline())
    object_type=[]
    for _ in range(num_total_objects):
        name=file.readline().strip()
        object_type.append(name)
        if name=='block':
            block_id = int(file.readline())
            block_props_1 = list(map(int, file.readline().split(',')))
            block_props_2 = list(map(int, file.readline().split(',')))
        if name=='sphere':
            sphere_id = int(file.readline())
            sphere_props_1 = list(map(int, file.readline().split(',')))
            sphere_props_2 = int(file.readline())
        if name=='cylinder':
            cylinder_id=int(file.readline())
            cyl_ori=file.readline().strip()
            cylinder_props_2=list(map(int, file.readline().split(',')))
            cylinder_props_3=list(map(int, file.readline().split(',')))

    # Read sheets
    num_total_sheets=int(file.readline())
    for _ in range(num_total_sheets):
        name=file.readline().strip()
        sheet_id = int(file.readline())
        prop_1 = int(file.readline())
        prop_2 = list(map(int, file.readline().split(',')))
        prop_3 = list(map(int, file.readline().split(',')))

    # put ports into correct formatting for stuff below - combine normal and spice together
    raw_ports = []

    # Read ports
    total_num_ports=int(file.readline())
    for _ in range(total_num_ports):
        type_port=file.readline().strip()
        name=file.readline().strip()
        port_props_1 = list(map(float, file.readline().split(',')))
        port_props_2 = list(map(int, file.readline().split(',')))
        port_props_3 = list(map(int, file.readline().split(',')))
        if type_port=='basic':
            xcc=0
            ycc=0
            zcc=0
            if name=='x':
                if port_props_3[1]>0:
                    ycc=0.5
                if port_props_3[2]>0:
                    zcc=0.5
            if name=='y':
                if port_props_3[0]>0:
                    xcc=0.5
                if port_props_3[2]>0:
                    zcc=0.5
            if name=='z':
                if port_props_3[0]>0:
                    xcc=0.5
                if port_props_3[1]>0:
                    ycc=0.5
            offsets = {
            'x': (1, 1.5-ycc, 1.5-zcc),
            'y': (1.5-xcc, 1, 1.5-zcc),
            'z': (1.5-xcc, 1.5-ycc, 1)}
            current_offset = offsets.get(name, (1.5, 1.5, 1.5))
            props_2_adjusted = tuple(val - off for val, off in zip(port_props_2[0:3], current_offset))
            new_entry = props_2_adjusted + tuple(port_props_3[0:3]) + (name,)
            raw_ports.append(new_entry)
        if type_port=='gridded':
            #this code is a litte odd - i did it to keep using prior lumped port structure formatting and shifted things here accordingly
            if name=='-x' or name=='+x':
                port_props_3.append(0)
                port_props_3[2]=port_props_3[1]-1
                port_props_3[1]=port_props_3[0]-1
                port_props_3[0]=0
                if name=='+x':
                    offsets = {
                    'x': (1.5, 0.5, 0.5)}
                if name=='-x':
                    offsets = {
                    'x': (2.5, 0.5, 0.5)}
            if name=='-y' or name=='+y':
                port_props_3.append(0)
                port_props_3[2]=port_props_3[1]-1
                port_props_3[1]=0
                port_props_3[0]=port_props_3[0]-1
                if name=='+y':
                    offsets = {
                    'y': (0.5, 1.5, 0.5)}
                if name=='-y':
                    offsets = {
                    'y': (0.5, 2.5, 0.5)}
            if name=='-z' or name=='+z':
                port_props_3.append(0)
                port_props_3[1]=port_props_3[1]-1
                port_props_3[0]=port_props_3[0]-1
                if name=='+z':
                    offsets = {
                    'z': (0.5, 0.5, 1.5)}
                if name=='-z':
                    offsets = {
                    'z': (0.5, 0.5, 2.5)}
            current_offset = offsets.get(name[1:], (0.0, 0.0, 0.0))
            props_2_adjusted = tuple(val - off for val, off in zip(port_props_2[0:3], current_offset))
            new_entry = props_2_adjusted + tuple(port_props_3[0:3]) + (name[1:],)
            raw_ports.append(new_entry)
            grid_file_name=file.readline().strip()

    num_spice_ports = int(file.readline())
    if spice:
        netlist_name = file.readline().strip()
        for _ in range(num_spice_ports):
            type_port=file.readline().strip()
            name=file.readline().strip()
            port_props_2 = list(map(int, file.readline().split(',')))
            port_props_3 = list(map(int, file.readline().split(',')))
            if type_port=='basic':
                bias_voltage = float(file.readline())
                v_name = file.readline().strip()
                I_name = file.readline().strip()
                Cap_name = file.readline().strip()
                xcc=0
                ycc=0
                zcc=0
                if name=='x':
                    if port_props_3[1]>0:
                        ycc=0.5
                    if port_props_3[2]>0:
                        zcc=0.5
                if name=='y':
                    if port_props_3[0]>0:
                        xcc=0.5
                    if port_props_3[2]>0:
                        zcc=0.5
                if name=='z':
                    if port_props_3[0]>0:
                        xcc=0.5
                    if port_props_3[1]>0:
                        ycc=0.5
                offsets = {
                'x': (1, 1.5-ycc, 1.5-zcc),
                'y': (1.5-xcc, 1, 1.5-zcc),
                'z': (1.5-xcc, 1.5-ycc, 1)}
                current_offset = offsets.get(name, (1.5, 1.5, 1.5))
                props_2_adjusted = tuple(val - off for val, off in zip(port_props_2[0:3], current_offset))
                new_entry = props_2_adjusted + tuple(port_props_3[0:3]) + (name,)
                raw_ports.append(new_entry)
            if type_port=='gridded':
                grid_file_name=file.readline().strip()
                bias_voltage = float(file.readline())
                v_name = file.readline().strip()
                I_name = file.readline().strip()
                #this code is a litte odd - i did it to keep using prior lumped port structure formatting and shifted things here accordingly
                if name=='-x' or name=='+x':
                    port_props_3.append(0)
                    port_props_3[2]=port_props_3[1]-1
                    port_props_3[1]=port_props_3[0]-1
                    port_props_3[0]=0
                    if name=='+x':
                        offsets = {
                        'x': (1.5, 0.5, 0.5)}
                    if name=='-x':
                        offsets = {
                        'x': (2.5, 0.5, 0.5)}
                if name=='-y' or name=='+y':
                    port_props_3.append(0)
                    port_props_3[2]=port_props_3[1]-1
                    port_props_3[1]=0
                    port_props_3[0]=port_props_3[0]-1
                    if name=='+y':
                        offsets = {
                        'y': (0.5, 1.5, 0.5)}
                    if name=='-y':
                        offsets = {
                        'y': (0.5, 2.5, 0.5)}
                if name=='-z' or name=='+z':
                    port_props_3.append(0)
                    port_props_3[1]=port_props_3[1]-1
                    port_props_3[0]=port_props_3[0]-1
                    if name=='+z':
                        offsets = {
                        'z': (0.5, 0.5, 1.5)}
                    if name=='-z':
                        offsets = {
                        'z': (0.5, 0.5, 2.5)}
                current_offset = offsets.get(name[1:], (0.0, 0.0, 0.0))
                props_2_adjusted = tuple(val - off for val, off in zip(port_props_2[0:3], current_offset))
                new_entry = props_2_adjusted + tuple(port_props_3[0:3]) + (name[1:],)
                raw_ports.append(new_entry)

# ----------------------------
# Load Fortran binaries and create other needed geometry arrays
# ----------------------------
data = np.fromfile("{}".format(filename_fdtd[:-4]+'_geometry.bin'), dtype=np.float32)
data = data.reshape((4, nx, ny, nz), order='F')   # Fortran order
volume, sheets_x, sheets_y, sheets_z = [np.ascontiguousarray(a) for a in data]
#replace all zeros with base number for blank spaces (unused) that will be lower number to filter in paraview
base_number=-6.0
igp_val=-5.0
bound_val=-4.0
pml_val=-3.0
port_ind_val=-2.0
port_val=-1.0
#so we need user to select >=0 for any numbers for sheets and volume materials
volume[volume == 0] = base_number
sheets_x[sheets_x == 0] = base_number
sheets_y[sheets_y == 0] = base_number
sheets_z[sheets_z == 0] = base_number
#create bounding box
bounding=np.zeros((nx,ny,nz))
#don't need to fill 5's first since all filled in
bounding[:,:,:]=bound_val
#now PML
pml_box=np.zeros((nx,ny,nz))
pml_box[:,:,:]=base_number
if (b_xlow==0 and b_xhigh==0):
    pml_box[0:10,:,:]=pml_val
    pml_box[nx-10:,:,:]=pml_val
if (b_ylow==0 and b_yhigh==0):
    pml_box[:,0:10,:]=pml_val
    pml_box[:,ny-10:,:]=pml_val
if (b_zlow==0 and b_zhigh==0):
    pml_box[:,:,0:10]=pml_val
    pml_box[:,:,nz-10:]=pml_val
if igp_on=='IGP':
    if igp_value=='+x':
        pml_box[igp_location:,:,:]=base_number
    if igp_value=='+y':
        pml_box[:,igp_location:,:]=base_number
    if igp_value=='+z':
        pml_box[:,:,igp_location:]=base_number
    if igp_value=='-x':
        pml_box[:igp_location,:,:]=base_number
    if igp_value=='-y':
        pml_box[:,:igp_location,:]=base_number
    if igp_value=='-z':
        pml_box[:,:,:igp_location]=base_number

# ----------------------------
# Create the grid for all volumes, sheets, and ports (in that order)
# ----------------------------
x = np.linspace(0.5, nx + 0.5, nx + 1)
y = np.linspace(0.5, ny + 0.5, ny + 1)
z = np.linspace(0.5, nz + 0.5, nz + 1)

gridToVTK(
    os.path.join(output_dir, "bounding box"),
    x, y, z,
    cellData={"MaterialID": bounding}
)

gridToVTK(
    os.path.join(output_dir, "pml box"),
    x, y, z,
    cellData={"MaterialID": pml_box}
)

gridToVTK(
    os.path.join(output_dir, "volume"),
    x, y, z,
    cellData={"MaterialID": volume}
)

def create_thin_grid(nx, ny, nz, sheet_data, axis='x', thickness=0.01):
    # 1. Define base edge coordinates for non-sheet axes (centered)
    coords = {
        'x': np.linspace(0.5, nx + 0.5, nx + 1),
        'y': np.linspace(0.5, ny + 0.5, ny + 1),
        'z': np.linspace(0.5, nz + 0.5, nz + 1)
    }
    
    # 2. Define the sheet axis with "stuttered" coordinates
    # For axis 'x', sheets sit at 0.5, 1.5, 2.5...
    # We create a tiny cell of 'thickness' centered on those interface values.
    n_primary = {"x": nx, "y": ny, "z": nz}[axis]
    half_thick = thickness / 2.0
    
    thin_coords = []
    for i in range(n_primary):
        interface_pos = i + 0.5
        thin_coords.append(interface_pos - half_thick)
        thin_coords.append(interface_pos + half_thick)
    
    # Add a final coordinate to close the last gap if needed, 
    # though usually sheets correspond 1:1 with cell indices.
    thin_coords.append(n_primary + 0.5)
    coords[axis] = np.array(thin_coords)
    
    # 3. Create the data array
    # VTK requirement: Data shape must be (len(coords)-1) for each axis
    shape = (len(coords['x']) - 1, len(coords['y']) - 1, len(coords['z']) - 1)
    expanded_data = np.full(shape, np.nan)
    
    # 4. Inject data into the thin slots
    # The sheets now occupy the even indices (0, 2, 4...) along the primary axis
    if axis == 'x':
        expanded_data[0::2, :, :] = sheet_data
    elif axis == 'y':
        expanded_data[:, 0::2, :] = sheet_data
    elif axis == 'z':
        expanded_data[:, :, 0::2] = sheet_data
        
    return coords['x'], coords['y'], coords['z'], expanded_data

# --- Apply to X Sheets ---
xc, yc, zc, data_x = create_thin_grid(nx, ny, nz, sheets_x, axis='x')
gridToVTK(os.path.join(output_dir, "sheets_x"), xc, yc, zc, cellData={"MaterialID": data_x})

# --- Apply to Y Sheets ---
xc, yc, zc, data_y = create_thin_grid(nx, ny, nz, sheets_y, axis='y')
gridToVTK(os.path.join(output_dir, "sheets_y"), xc, yc, zc, cellData={"MaterialID": data_y})

# --- Apply to Z Sheets ---
xc, yc, zc, data_z = create_thin_grid(nx, ny, nz, sheets_z, axis='z')
gridToVTK(os.path.join(output_dir, "sheets_z"), xc, yc, zc, cellData={"MaterialID": data_z})

if igp_on=='IGP':
    igp_grid=np.zeros((nx,ny,nz))
    igp_grid[:,:,:]=base_number
    if igp_value=='+x' or igp_value=='-x':
        igp_grid[igp_location-1,:,:]=igp_val
        xc, yc, zc, data_x = create_thin_grid(nx, ny, nz, igp_grid, axis='x')
        gridToVTK(os.path.join(output_dir, "IGP"), xc, yc, zc, cellData={"MaterialID": data_x})
    if igp_value=='+y' or igp_value=='-y':
        igp_grid[:,igp_location-1,:]=igp_val
        xc, yc, zc, data_y = create_thin_grid(nx, ny, nz, igp_grid, axis='y')
        gridToVTK(os.path.join(output_dir, "IGP"), xc, yc, zc, cellData={"MaterialID": data_y})
    if igp_value=='+z' or igp_value=='-z':
        igp_grid[:,:,igp_location-1]=igp_val
        xc, yc, zc, data_z = create_thin_grid(nx, ny, nz, igp_grid, axis='x')
        gridToVTK(os.path.join(output_dir, "IGP"), xc, yc, zc, cellData={"MaterialID": data_z})

def create_port_grid(port_data, start_idx, nx_p, ny_p, nz_p, thickness=0.01):
    dims = {'x': nx_p, 'y': ny_p, 'z': nz_p}
    starts = {'x': start_idx[0], 'y': start_idx[1], 'z': start_idx[2]}
    coords = {}
    
    for axis in ['x', 'y', 'z']:
        n = dims[axis]
        s = starts[axis]
        
        if n == 0:
            center = float(s + 1.0)
            coords[axis] = np.array([
                center - thickness/2.0, 
                center + thickness/2.0
            ], dtype=np.float32)
        else:
            coords[axis] = np.linspace(s + 0.5, s + n + 0.5, n + 1, dtype=np.float32)

    return coords['x'], coords['y'], coords['z']

def create_direction_indicator(sx, sy, sz, nx_p, ny_p, nz_p, direction, thickness=0.02):
    def axis_coords(start, size):
        if size == 0:
            center = float(start + 1.0)
        else:
            center = start + 0.5 + size / 2.0
        return center

    cx = axis_coords(sx, nx_p)
    cy = axis_coords(sy, ny_p)
    cz = axis_coords(sz, nz_p)

    #def span(start, size):
    #   if size == 0:
    #        center = float(start + 1.0)
    #        return np.array([center - 1.0, center + 1.0], dtype=np.float32)
    #    else:
    #        return np.array([start + 0.5 - 1.0, start + size + 0.5 + 1.0], dtype=np.float32)
        
    def span(start, size):
        padding = 0.5 
        if size == 0:
            center = float(start + 1.0)
            return np.array([center - padding, center + padding], dtype=np.float32)
        else:
            # This spans from the start face to the end face, plus 0.5 on each side
            return np.array([start + 0.5 - padding, start + size + 0.5 + padding], dtype=np.float32)

    sx_span = span(sx, nx_p)
    sy_span = span(sy, ny_p)
    sz_span = span(sz, nz_p)

    if direction == 'x':
        px = np.array([cx - thickness/2.0, cx + thickness/2.0], dtype=np.float32)
        py = sy_span
        pz = sz_span
    elif direction == 'y':
        px = sx_span
        py = np.array([cy - thickness/2.0, cy + thickness/2.0], dtype=np.float32)
        pz = sz_span
    elif direction == 'z':
        px = sx_span
        py = sy_span
        pz = np.array([cz - thickness/2.0, cz + thickness/2.0], dtype=np.float32)

    data = np.full((1, 1, 1), port_ind_val, dtype=np.float32)
    return px, py, pz, data

port_vtr_names = []

for i, (sx, sy, sz, nx_p, ny_p, nz_p, direction) in enumerate(raw_ports):
    # --- Port itself (unchanged) ---
    data_shape = (max(1, nx_p), max(1, ny_p), max(1, nz_p))
    p_data = np.full(data_shape, port_val, dtype=np.float32)
    px, py, pz = create_port_grid(p_data, (sx, sy, sz), nx_p, ny_p, nz_p)
    vtr_name = f"port_{i}"
    gridToVTK(os.path.join(output_dir, vtr_name), px, py, pz,
              cellData={"MaterialID": p_data})
    port_vtr_names.append(f"{vtr_name}.vtr")

    # --- Direction indicator sheet (new) ---
    dx, dy, dz, d_data = create_direction_indicator(sx, sy, sz, nx_p, ny_p, nz_p, direction)
    dir_name = f"port_{i}_dir"
    gridToVTK(os.path.join(output_dir, dir_name), dx, dy, dz,
              cellData={"MaterialID": d_data})
    port_vtr_names.append(f"{dir_name}.vtr")

# ----------------------------
# Create the .vtm Wrapper
# ----------------------------
vtm_file = os.path.join(script_dir, "{}.vtm".format(filename))

all_files = [
    "bounding box.vtr", "pml box.vtr", "volume.vtr", 
    "sheets_x.vtr", "sheets_y.vtr", "sheets_z.vtr"
]
all_files.extend(port_vtr_names)

if igp_on=='IGP':
    all_files = [
    "bounding box.vtr", "pml box.vtr", "volume.vtr", 
    "sheets_x.vtr", "sheets_y.vtr", "sheets_z.vtr", "IGP.vtr"
    ]
    all_files.extend(port_vtr_names)

with open(vtm_file, "w") as f:
    f.write('<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">\n')
    f.write('  <vtkMultiBlockDataSet>\n')
    
    for i, fname in enumerate(all_files):
        f.write(f'    <DataSet index="{i}" file="{output_dir}/{fname}"/>\n')
        
    f.write('  </vtkMultiBlockDataSet>\n')
    f.write('</VTKFile>')
