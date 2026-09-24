import numpy as np

#np.int32 will match the standard precision in Fortran
#np.int64 will match the standard (8) precision in Fortran
#note though that real(8) and double precision are not technically the same thing
#same rules for floats which we need here

#first setup the array
x_size=600
y_size=600
z_size=350

#create the array
#zeros will be ignored by fdtd - specific to this optional read in section
data=np.zeros((x_size,y_size,z_size), dtype=np.float32)

#draw custom geometry here
#recall that the last number in the slice is not used
#this format will match the way we upload in normal fdtd py script

TEM_line_base=np.load('waveguide_geometry.npy')
TEM_line_base[TEM_line_base == 1] = 6

material_ID_number=6
data[50:50+458,50:50+480,50:50+242]=TEM_line_base

#testing zone
#print(data[29:31,29:31,29:31])

#save
data.flatten(order='F').tofile('optional_geom_bulk.bin')