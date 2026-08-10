import numpy as np

#np.int32 will match the standard precision in Fortran
#np.int64 will match the standard (8) precision in Fortran
#note though that real(8) and double precision are not technically the same thing
#same rules for floats which we need here

#first setup the array
x_size=65
y_size=30
z_size=30

#create the array
#zeros will be ignored by fdtd - specific to this optional read in section
data=np.zeros((x_size,y_size,z_size), dtype=np.float32)

#draw custom geometry here
#recall that the last number in the slice is not used
#this format will match the way we upload in normal fdtd py script
start_x=35
start_y=1
start_z=1
step_x=8
step_y=30
step_z=30
material_ID_number=5
data[start_x-1:start_x+step_x-1,start_y-1:start_y+step_y-1,start_z-1:start_z+step_z-1]=material_ID_number

#testing zone
#print(data[29:31,29:31,29:31])

#save
data.flatten(order='F').tofile('optional_geom_bulk.bin')