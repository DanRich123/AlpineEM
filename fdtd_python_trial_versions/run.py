import os
import sys
import numpy as np
#os.chdir('C:/Users/dari6475/Desktop/temp')

### General Notes ####################################
#this doesn't have to be done in python.
#good rule of thumb is 20 cells from PML for any geometric item - pml is 10 cells thick.
#pulse shape rules for freq range that is trustable is as follows:
#rule is for gaussian, 0 to 2.15*f_parameter. similar for diff gaussian but it's 2.76*f_parameter and 0.06*f_parameter limits.
#need min 10 grid points (not cells) per wavlength in incident direction but 20 recommended - consider dielectrics and oblique angles
#incident angle phi is any angle 0<=phi<360.
#incident angle theta is any angle 0<=theta<=180.
#incident polarization is any angle 0<=pol<360. 90 is -theta E-polarized and 0 is -phi E-polarized.
#good rule of thumb for number of time steps is 1000-4000 but you can examine the output time series to get a better idea.
#use .format command to easily input variables and enable more dynamic control (loops) over what you write into the text file.
#######################################################

#generate needed text file. Creates one if it doesn't exist and replaces if it does.
f=open('{}'.format(sys.argv[1]), 'w')

f.write('100,100,100\n') #xsize,ysize,zsize of number of cells
f.write('0.5E-3,0.5E-3,0.5E-3\n') #cell size in meters - recommended to use the same size for each one.
f.write('300\n') #time steps
f.write('5E9\n') #frequency parameter of output data - see above notes
f.write('2\n') #pulse type - 1 is a normalized gaussian pulse and 2 is a normalized differentiated gaussian pulse
f.write('pml\n') #'uc' or 'pml' - 'uc' is a unit cell version for infinite media (PML +-x and PBC y,z) and 'pml' is PML on all boundaries for free space imitation

f.write('plane wave\n') #excitation - 'plane wave' or 'antenna'

#if plane wave is chosen:
f.write('1,0\n') #amplitude and time delay (V/m and seconds)
#if unit cell use theta,phi,pol of (90,0,xx) for plane wave
f.write('90,0,270\n') #incident angles and polarization in degrees (theta,phi,pol)
#if antenna is chosen:
#f.write('0.001,0\n') #amplitude and time delay (V and seconds)
#f.write('4\n') #internal port number to excite (ports numbered by submission order below - 1,2,3,...,n)

f.write('2\n') #num far field angles we want returned data for - both pols returned
#if >0 then add this line each time
f.write('90,0\n') #ang1 (theta,phi)
f.write('90,180\n') #ang2 (theta,phi)
#many angles with ease:
#angle_array=np.linspace(0,180,11)
#for i in range(len(angle_array)):
    #f.write('{},0\n'.format(angle_array[i]))

f.write('none\n') #'none' or 'IGP' for infinite ground plane condition if using all PML version with infinite metal plane
#if 'IGP' then these lines:
#f.write('-x\n') #-x,+x,-y,+y,-z,+z indicating which wall is PEC
#f.write('40\n') #height location of the plane. Note: you must draw a metal block at this location as well.

f.write('no\n') #'yes' or 'no' to video 
#if 1, then add lines of code below
#f.write('z\n') #plane to make the cut
#f.write('15\n') #location/height of the plane to make the cut

f.write('data.dat\n') #name the output file

f.write('4\n') #number of materials
#do these lines for each bulk volume material
f.write('volume\n')
f.write('5\n') #material idenfitication number
f.write('2.2,2.2,2.2\n') #relative epsilon
f.write('0,0,0\n') #sigma
#do these lines for each bulk volume material
f.write('volume\n')
f.write('6\n') #material idenfitication number
f.write('1,1,1\n') #relative epsilon
f.write('1E8,1E8,1E8\n') #sigma
#do these lines for each sheet material
f.write('sheet\n')
f.write('4\n') #identifier for sheet material
f.write('0\n') #surface impedance
#do these lines for each sheet material
f.write('sheet\n')
f.write('7\n') #identifier for sheet material
f.write('50\n') #surface impedance

f.write('2\n') #number of objects
#do these lines for each block
f.write('block\n') #state if 'block' or 'sphere' or 'cylinder'
f.write('5\n') #material identification number
f.write('35,31,31\n') #where it starts, 1 is first location
f.write('8,10,10\n') #cells long in x, y, z directions
#do these lines for each block
f.write('block\n') #state if 'block' or 'sphere' or 'cylinder'
f.write('6\n') #material identification number
f.write('41,31,31\n') #where it starts, 1 is first location
f.write('2,10,10\n') #cells long in x, y, z directions
#do these lines for each sphere
#f.write('sphere\n') #state if 'block' or 'sphere' or 'cylinder'
#f.write('5\n') #material identification number
#f.write('85,85,85\n') #location for center of sphere - cube on the right side of middle.
#f.write('50\n') #sphere radius in cells
#do these lines for each cylinder
#f.write('cylinder\n') #state if 'block' or 'sphere' or 'cylinder'
#f.write('5\n') #material identification number
#f.write('x\n') #major axis of the cylinder - 0,1,2 for x,y,z
#f.write('10,15,15\n') #starting location for height, center for circle portion in other two directions (x,y) or (y,z) or (x,z)
#f.write('10,5\n') #distance in major axis, radius of circle portion

f.write('8\n') #number of sheets
#do these lines for each sheet
f.write('x\n') #x,y,z normal to the sheet surface
f.write('4\n') #idenfitication number for sheet material
f.write('35\n') #sheet location/height in normal direction
f.write('35,35\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('3,10\n') #cells along other two directions (x,y) or (y,z) or (x,z)
#do these lines for each sheet
f.write('x\n') #x,y,z normal to the sheet surface
f.write('4\n') #idenfitication number for sheet material
f.write('35\n') #sheet location/height in normal direction
f.write('35,46\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('3,10\n') #cells along other two directions (x,y) or (y,z) or (x,z)
#do these lines for each sheet
f.write('x\n') #x,y,z normal to the sheet surface
f.write('7\n') #idenfitication number for sheet material
f.write('35\n') #sheet location/height in normal direction
f.write('50,35\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('2,10\n') #cells along other two directions (x,y) or (y,z) or (x,z)
#do these lines for each sheet
f.write('x\n') #x,y,z normal to the sheet surface
f.write('7\n') #idenfitication number for sheet material
f.write('35\n') #sheet location/height in normal direction
f.write('50,46\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('2,10\n') #cells along other two directions (x,y) or (y,z) or (x,z)
#do these lines for each sheet
f.write('x\n') #x,y,z normal to the sheet surface
f.write('4\n') #idenfitication number for sheet material
f.write('36\n') #sheet location/height in normal direction
f.write('35,35\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('10,3\n') #cells along other two directions (x,y) or (y,z) or (x,z)
#do these lines for each sheet
f.write('x\n') #x,y,z normal to the sheet surface
f.write('4\n') #idenfitication number for sheet material
f.write('36\n') #sheet location/height in normal direction
f.write('46,35\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('10,3\n') #cells along other two directions (x,y) or (y,z) or (x,z)
#do these lines for each sheet
f.write('z\n') #x,y,z normal to the sheet surface
f.write('7\n') #idenfitication number for sheet material
f.write('45\n') #sheet location/height in normal direction
f.write('65,45\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('2,3\n') #cells along other two directions (x,y) or (y,z) or (x,z)
#do these lines for each sheet
f.write('z\n') #x,y,z normal to the sheet surface
f.write('7\n') #idenfitication number for sheet material
f.write('45\n') #sheet location/height in normal direction
f.write('38,35\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('3,3\n') #cells along other two directions (x,y) or (y,z) or (x,z)

f.write('0\n') #number of internal feed ports - spice ports are seperate
"""
#do these lines for each feed/port
f.write('z\n') #x,y,z for E-field direction
f.write('50\n') #port impedance (must be real)
f.write('35,5,15\n') #starting location of the port
f.write('0,0,1\n') #port size in all directions
#do these lines for each feed/port
f.write('z\n') #x,y,z for E-field direction
f.write('50\n') #port impedance (must be real)
f.write('35,20,15\n') #starting location of the port
f.write('0,0,1\n') #port size in all directions
#do these lines for each feed/port
f.write('y\n') #x,y,z for E-field direction
f.write('50\n') #port impedance (must be real)
f.write('36,15,5\n') #starting location of the port
f.write('0,1,0\n') #port size in all directions
#do these lines for each feed/port
f.write('x\n') #x,y,z for E-field direction
f.write('50\n') #port impedance (must be real)
f.write('37,15,15\n') #starting location of the port
f.write('1,0,0\n') #port size in all directions
"""

f.close()