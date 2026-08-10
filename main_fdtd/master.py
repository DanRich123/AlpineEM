import os
import time as tp
import numpy as np
#os.chdir('C:/Users/dari6475/Desktop/temp')

#this is the master py script that demonstrates how to use all features with ample notes

### General Notes ####################################
#this script generates a text file and then runs the executable from the precompiled fdtd solver.
#this doesn't have to be done in python.
#there are 4 versions - normal or kmax - w/ or w/ out spice utlilized
#good rule of thumb is 20 cells from PML for any geometric item - pml is 10 cells thick.
#need min 10 grid points (not cells) per wavlength but 20 recommended - consider dielectrics and oblique angles.
#good rule of thumb for number of time steps is 1000-4000 but you can examine the output time series to get a better idea.
#use .format command or sys.argv to easily input variables and enable more dynamic control (loops) over what you write to text file.
#gridded ports require a non-dispersive medium at the port location, and they assume whatever mode is at the cross section is non-dispersive as well.

#if using the normal version#
#pulse shape rules for freq range that is trustable is as follows: (-20dB field, -40dB in power)
#rule is for gaussian, 0 to 2.146*f_parameter. similar for diff gaussian but it's 3.087*f_parameter and 0.183*f_parameter limits.
#incident angle phi is any angle 0<=phi<360.
#incident angle theta is any angle 0<=theta<=180.
#incident polarization is any angle 0<=pol<360. 90 is -theta E-polarized and 0 is -phi E-polarized.

#if using the kmax version#
#frequency range that is trustable is chosen based on ky,kz values and the frequency parameter.
#range is roughly 0.241*f+c/(2pi)*sqrt(k1^2+k2^2) up to 1.759*f+c/(2pi)*sqrt(k1^2+k2^2) for -20dB field, -40dB in power.
#avoid using k values near 2pi/d - where d is the length in meters of the direction of k vector - also true for multiples of 2pi/d.
#this is because at integer multiples of 2pi/d there is a known numerical resonance behavior at the center of a Brillouin zone that can cause convergence issues.

#if using spice version#
#make sure you install ngspice appropriately - there are guides with this document.
#if using kmax+spice - as noted below - you need to use duplicate ports in the fdtd grid with different spice identities
#this is because kmax uses complex time domain data and spice doesn't support transient sims w/ complex data.
#ex: 1 port desired; 2 same location w/ different nomenclature ports entered here, 2 identical w/ different nomeclature circuits entered into netlist, 1 port produced from post processor
#note the two corresponding ports must be submitted one after another in sequential order.
#spice ports of any type can't have a dispersive medium at the port location in the fdtd grid.

#there are instructions present elsewhere for using openMP or ACC.
#nothing new required beyond using the right solver name here.
#note acc is not recommended in general for kmax because of elongated arrays and complex data types
#######################################################

#select the solver and name the inputs txt file to create - both used at bottom of script to submit via command line for convenience
solver='FDTD'
inputs_name='inputs.txt'

#generate needed text file. Creates one if it doesn't exist and replaces if it does.
f=open('{}'.format(inputs_name), 'w')

f.write('65,30,30\n') #xsize,ysize,zsize of number of cells
f.write('0.5E-3,0.5E-3,0.5E-3\n') #cell size in meters - recommended to use the same size for each one.
f.write('4000\n') #time steps
f.write('1\n') #time step reduction factor - CFL*0.98/factor
f.write('5E9\n') #frequency parameter of output data - see above notes
f.write('2\n') #pulse type - normal: 1 Gauss, 2 Diff Gauss; kmax: any value, currently unused

f.write('0,0,1,1,1,1\n') #boundary types 0-pml, 1-pbc (-x,+x,-y,+y,-z,+z)
f.write('none\n') #'none' or 'IGP' for infinite ground plane condition - not needed for a unit cell (2 pbc pairs)
#if 'IGP' then these lines:
#f.write('-x\n') #-x,+x,-y,+y,-z,+z indicating which wall is PEC
#f.write('40\n') #cube location of the plane. Note: you must draw a metal block at this location as well.

#if using normal fdtd solver:
f.write('plane wave\n') #excitation - 'plane wave' or 'antenna'
#if 'plane wave':
f.write('1,0\n') #amplitude and time delay (V/m and seconds)
f.write('90,0,270\n') #incident angles and polarization in degrees (theta,phi,pol)
#if 'antenna':
#f.write('0.001,0\n') #amplitude and time delay (V and seconds)
#f.write('4\n') #internal port number to excite (ports numbered by submission order below - 1,2,3,...,n)

#if using kmax fdtd solver:
#f.write('plane wave\n') #excitation - 'plane wave' or 'antenna'
#if 'plane wave':
#f.write('1,0\n') #amplitude and time shift (V/m and seconds)
#f.write('TM\n') #orientation of TE or TM
#f.write('+\n') #propogation direction
#if 'antenna':
#f.write('0.001,0\n') #amplitude and time shift (V and seconds)
#f.write('1\n') #internal port number to excite (ports numbered by submission order below - 1,2,3,...,n)
#for either, excitation specify k:
#f.write('0,0\n') #k1 and k2 in rad/m (ky,kz),(kx,kz), or (kx,ky) - used in both excitation schemes for boundary phase shift
##speciate note###
#if k1,k2=0,0 then TM or TE selects '1' directed transverse field component for plane wave case 
#if k1,k2=0,0 for antenna case then output TE will correspond to E1,H2 and output TM will coresopnd to E2,H1 (+- depending on prop direction)
#variable1 = float(sys.argv[1]) #set to retrieve the k values submitted via command line for ease in looping
#variable2 = float(sys.argv[2])

f.write('2\n') #num far field angles we want returned data for - both pols returned
#far field phase center is cubic center at the center of the grid - round down if sizes are odd. I.e. 65/2 yields 32.
#if IGP is used, phase center will be centered at the ground plane cube
#if >0 then add this line each time
f.write('90,0\n') #ang1 (theta,phi)
f.write('90,180\n') #ang2 (theta,phi)
#many angles with ease:
#angle_array=np.linspace(0,180,11)
#for i in range(len(angle_array)):
    #f.write('{},0\n'.format(angle_array[i]))

f.write('no\n') #'yes' or 'no' to export fields for videos or otherwise
#if yes, then add lines of code below
#f.write('y\n') #plane to make the cut
#f.write('16\n') #location/height of the plane to make the cut - cube centered cut

f.write('data.dat\n') #name the output file, this will also influence the name of the geometry output file name. 

f.write('4\n') #number of materials, these ID numbers must be greater than zero
#do these lines for each bulk volume material
f.write('volume\n')
f.write('simple\n') #type of material
f.write('5\n') #material idenfitication number
f.write('2.2,2.2,2.2\n') #relative epsilon
f.write('0,0,0\n') #sigma
#do these lines for each bulk volume material
f.write('volume\n')
f.write('simple\n') #type of material
f.write('6\n') #material idenfitication number
f.write('1,1,1\n') #relative epsilon
f.write('1E8,1E8,1E8\n') #sigma
#do these lines for each sheet material
f.write('sheet\n')
f.write('4\n') #identifier for sheet material
f.write('0\n') #surface impedance ohms/sq
#do these lines for each sheet material
f.write('sheet\n')
f.write('7\n') #identifier for sheet material
f.write('50\n') #surface impedance ohms/sq
#do these lines for each bulk volume material
#f.write('volume\n')
#f.write('plasma\n') #type of material
#f.write('9\n') #material idenfitication number
#f.write('3,3,3\n') #relative epsilon
#f.write('0.1,0.1,0.1\n') #sigma
#f.write('3\n') #number of poles if not 'simple'
#f.write('30E9,30E9,30E9\n') #plasma angular frequency term
#f.write('10E9,10E9,10E9\n') #plasma collisional angular frequency gamma
#f.write('40E9,40E9,40E9\n') #plasma angular frequency term
#f.write('5E9,5E9,5E9\n') #plasma collisional angular frequency gamma
#f.write('35E9,35E9,35E9\n') #plasma angular frequency term
#f.write('20E9,20E9,20E9\n') #plasma collisional angular frequency gamma

#optional geometry section
f.write('no\n') # 'yes' or 'no' to use an optional bulk geometry file
#f.write('optional_geom_bulk.bin\n') # if 'yes', what is the name of the binary file

f.write('2\n') #number of objects
#do these lines for each block
f.write('block\n') #state if 'block' or 'sphere' or 'cylinder'
f.write('5\n') #material identification number
f.write('35,1,1\n') #where it starts, 1 is first location
f.write('8,30,30\n') #cells long in x, y, z directions
#do these lines for each block
f.write('block\n') #state if 'block' or 'sphere' or 'cylinder'
f.write('6\n') #material identification number
f.write('41,1,1\n') #where it starts, 1 is first location
f.write('2,30,30\n') #cells long in x, y, z directions
#do these lines for each sphere
#f.write('sphere\n') #state if 'block' or 'sphere' or 'cylinder'
#f.write('5\n') #material identification number
#f.write('85,85,85\n') #location for center of sphere - cube on the right side of middle.
#f.write('50\n') #sphere radius in cells
#do these lines for each cylinder
#f.write('cylinder\n') #state if 'block' or 'sphere' or 'cylinder'
#f.write('5\n') #material identification number
#f.write('x\n') #major axis of the cylinder - x,y,z
#f.write('30,15,15\n') #starting location (x,y,z) - center has cube on the right (forward) side of middle. 
#f.write('5,5\n') #distance in major axis, radius of circle portion. distances go out radially and in height from starting location.

f.write('8\n') #number of sheets
#do these lines for each sheet
f.write('x\n') #x,y,z normal to the sheet surface
f.write('4\n') #idenfitication number for sheet material
f.write('35\n') #sheet location/height in normal direction - lower side of this cube
f.write('5,5\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('3,10\n') #cells along other two directions (x,y) or (y,z) or (x,z)
#do these lines for each sheet
f.write('x\n') #x,y,z normal to the sheet surface
f.write('4\n') #idenfitication number for sheet material
f.write('35\n') #sheet location/height in normal direction - lower side of this cube
f.write('5,16\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('3,10\n') #cells along other two directions (x,y) or (y,z) or (x,z)
#do these lines for each sheet
f.write('x\n') #x,y,z normal to the sheet surface
f.write('7\n') #idenfitication number for sheet material
f.write('35\n') #sheet location/height in normal direction - lower side of this cube
f.write('20,5\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('2,10\n') #cells along other two directions (x,y) or (y,z) or (x,z)
#do these lines for each sheet
f.write('x\n') #x,y,z normal to the sheet surface
f.write('7\n') #idenfitication number for sheet material
f.write('35\n') #sheet location/height in normal direction - lower side of this cube
f.write('20,16\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('2,10\n') #cells along other two directions (x,y) or (y,z) or (x,z)
#do these lines for each sheet
f.write('x\n') #x,y,z normal to the sheet surface
f.write('4\n') #idenfitication number for sheet material
f.write('36\n') #sheet location/height in normal direction - lower side of this cube
f.write('5,5\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('10,3\n') #cells along other two directions (x,y) or (y,z) or (x,z)
#do these lines for each sheet
f.write('x\n') #x,y,z normal to the sheet surface
f.write('4\n') #idenfitication number for sheet material
f.write('36\n') #sheet location/height in normal direction
f.write('16,5\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('10,3\n') #cells along other two directions (x,y) or (y,z) or (x,z)
#do these lines for each sheet
f.write('z\n') #x,y,z normal to the sheet surface
f.write('7\n') #idenfitication number for sheet material
f.write('15\n') #sheet location/height in normal direction
f.write('35,15\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('2,3\n') #cells along other two directions (x,y) or (y,z) or (x,z)
#do these lines for each sheet
f.write('z\n') #x,y,z normal to the sheet surface
f.write('7\n') #idenfitication number for sheet material
f.write('15\n') #sheet location/height in normal direction
f.write('38,15\n') #1,1 is first location, where it starts in other two directions (x,y) or (y,z) or (x,z)
f.write('3,3\n') #cells along other two directions (x,y) or (y,z) or (x,z)

f.write('4\n') #number of internal lumped feed ports - spice lumped ports are seperate
#do these lines for each feed/port
f.write('basic\n') #"basic" or "gridded" - basic is an omnidirectional injector, gridded is 1 way mode injector.
f.write('z\n') #x,y,z for E-field direction
f.write('50\n') #port impedance (must be real) 
f.write('35,5,15\n') #starting location of the port (bottom corner of the cube)
f.write('0,3,1\n') #port size in all directions
#do these lines for each feed/port
f.write('basic\n') #"basic" or "gridded" - basic is an omnidirectional injector, gridded is 1 way mode injector.
f.write('z\n') #x,y,z for E-field direction
f.write('50\n') #port impedance (must be real) 
f.write('35,20,15\n') #starting location of the port (bottom corner of the cube)
f.write('0,2,1\n') #port size in all directions
#do these lines for each feed/port
f.write('basic\n') #"basic" or "gridded" - basic is an omnidirectional injector, gridded is 1 way mode injector.
f.write('y\n') #x,y,z for E-field direction
f.write('50\n') #port impedance (must be real) 
f.write('36,15,5\n') #starting location of the port (bottom corner of the cube)
f.write('0,1,3\n') #port size in all directions
#do these lines for each feed/port
f.write('basic\n') #"basic" or "gridded" - basic is an omnidirectional injector, gridded is 1 way mode injector.
f.write('x\n') #x,y,z for E-field direction
f.write('50\n') #port impedance (must be real) 
f.write('37,15,15\n') #starting location of the port (bottom corner of the cube)
f.write('1,3,0\n') #port size in all directions
#do these lines for each feed/port
#f.write('gridded\n') #"basic" or "gridded" - basic is an omnidirectional injector, gridded is 1 way mode injector.
#f.write('-x\n') #-x,+x,-y,+y,-z,+z for wave propogation direction
#f.write('48\n') #port impedance (must be real)
#f.write('47,25,25\n') #bottom left cube of the mask
#f.write('50,50\n') #size of grid file (in cells)
#f.write('gridded_feed.bin\n') #mode data - can be any mode type or geometry but only valid if simulated frequencies support this mode or if we invoke linearity and ignore some.

####REQUIRES SPICE COMPILATION IF USED BEYOND ZERO PORTS#######
f.write('0\n') #number of spice lumped ports connected to fdtd grid - 1 netlist for the entire circuit(s)
#special note - if using kmax+spice refer to notes above for special usage
#do this line if spice ports>0
#f.write('fdtd_netlist.cir\n') #name of netlist file in this folder
#do these lines for each spice port in the netlist
#f.write('basic\n') #"basic" or "gridded" - basic is an omnidirectional injector, gridded is 1 way mode injector.
#f.write('z\n') #x,y,z for E-field direction
#f.write('35,20,15\n') #starting location of the port (bottom corner of the cube)
#f.write('0,2,1\n') #port size in all directions
#f.write('0\n') #DC bias voltage if any - 0 if none
#f.write('out2\n') #name of the port to link to fdtd
#f.write('I2\n') #name of the current source associated with the port linked to fdtd
#f.write('Cp2\n') #name of the capacitance linked to this fdtd cell - will update on fdtd side to be correct value
#do these lines for each spice port in the netlist
#f.write('gridded\n') #"basic" or "gridded" - basic is an omnidirectional injector, gridded is 1 way mode injector.
#f.write('-x\n') #-x,+x,-y,+y,-z,+z for wave propogation direction
#f.write('47,25,25\n') #bottom left cube of the mask
#f.write('50,50\n') #size of grid file (in cells)
#f.write('gridded_feed.bin\n') #mode data - can be any mode type or geometry but only valid if simulated frequencies support this mode or if we invoke linearity and ignore some.
#f.write('0\n') #DC bias voltage if any - 0 if none
#f.write('out2\n') #name of the port to link to fdtd
#f.write('I2\n') #name of the current source associated with the port linked to fdtd
#do these lines for each spice port in the netlist
#f.write('gridded\n') #"basic" or "gridded" - basic is an omnidirectional injector, gridded is 1 way mode injector.
#f.write('-x\n') #-x,+x,-y,+y,-z,+z for wave propogation direction
#f.write('47,25,25\n') #bottom left cube of the mask
#f.write('50,50\n') #size of grid file (in cells)
#f.write('gridded_feed.bin\n') #mode data - can be any mode type or geometry but only valid if simulated frequencies support this mode or if we invoke linearity and ignore some.
#f.write('0\n') #DC bias voltage if any - 0 if none
#f.write('out3\n') #name of the port to link to fdtd
#f.write('I3\n') #name of the current source associated with the port linked to fdtd

f.close()

#section to write the netlist from here with all PWL's input ready to go as dummy variables
#make sure to define all variables here to align with fdtd naming scheme
#recall kmax+spice requires duplicate ports with different naming and idnetities for real+imag splitting
#several examples below:

#normal example w/ basic ports:
#times=4000 #get from above
#factor=1 #get from above
#base_time_step=9.4365835E-13 #get from a clear case or similar
#f=open('fdtd_netlist.cir','w')
#f.write('*Validation of two 50 ohm ports with FDTD\n')
#f.write('I1 0 out1 PWL(0 0)\n')
#f.write('I2 0 out2 PWL(0 0)\n')
#f.write('.tran 1E-13 {}\n'.format(times*base_time_step/factor))
#f.write('Cp1 out1 0 1\n') #Cp1 and Cp2 values are overwritten by FDTD solver - they are grid capacitances
#f.write('Cp2 out2 0 1\n')
#f.write('Rterm1 out1 0 50\n')
#f.write('Rterm2 out2 0 50\n')
#f.write('.end\n')
#f.close()

#normal example w/ basic ports:
#times=4000 #get from above
#factor=1 #get from above
#base_time_step=9.4365835E-13 #get from a clear case or similar
#f=open('fdtd_netlist.cir','w')
#f.write("*SPICE two port circuit for MACOM MAVR-011020 varactors\n")
#f.write("I1 0 out_port1 PWL(0 0)\n")
#f.write("I2 0 out_port2 PWL(0 0)\n")
#f.write(".tran 1e-13 {}\n".format(times*base_time_step/factor))
#f.write("*The following capacitance and its parameters are set by the FDTD cell size\n")
#f.write(".param e0=8.85e-12\n")
#f.write(".param er=1\n")
#f.write(".param A=0.3048e-3*0.3048e-3\n")
#f.write(".param dl=0.3048e-3\n")
#f.write(".param Cshunt=er*e0*A/dl\n")
#f.write("C_fdtd1 out_port1 0 {Cshunt}\n")
#f.write("C_fdtd2 out_port2 0 {Cshunt}\n")
#f.write("*Two separate bias voltages\n")
#f.write('Vbias1 2_port1 0 SIN({} 0 7E9)\n'.format(bias_1)) # (SPICE_Net_Name CIR_LOC1 CIR_LOC2 SIN(BIAS AMP FREQ TIME_DELAY)
#f.write('Vbias2 2_port2 0 SIN({} 0 7E9)\n'.format(bias_2)) # (SPICE_Net_Name CIR_LOC1 CIR_LOC2 SIN(BIAS AMP FREQ TIME_DELAY)
#f.write("*Large RF chokes to DC bias the varactors\n")
#f.write("Lbias1 out_port1 2_port1 50u\n")
#f.write("Lbias2 out_port2 2_port2 50u\n")
#f.write("*Place the varactors\n")
#f.write("X1 0 out_port1 mavr011020\n")
#f.write("X2 0 out_port2 mavr011020\n")
#f.write("*Confirmed in measurement that MAVR011020 has negligible series inductance\n")
#f.write("*Output voltage across each varactor written to spice_intermed.txt\n")
#f.write(".endc\n")
#f.write(".subckt mavr011020 A K\n")
#f.write("Rs A 1 13.2\n")
#f.write("*Use parasitic capacitance value from mavr000120; same package\n")
#f.write("Cp A K 0.1394e-12\n")
#f.write("D1 1 K Diode_Model\n")
#f.write("*If encounter problem \"trouble with x1:diode_model-instance d.x1.d1,\" ... set M=1.0001, not M=1\n")
#f.write(".model Diode_Model D(IS=1e-14 RS=0 N=2 TT=3e-09 CJO=0.233p VJ=2.4 M=1.0001 EG=1.42 XTITUN=3 KF=0 AF=1 FC=0.5 BV=20 IBV=1e-05 IKF=0)\n")
#f.write(".ends\n")
#f.write(".end\n")
#f.close()

#normal example w/ gridded ports:
#times=4000 #get from above
#factor=1 #get from above
#base_time_step=9.4365835E-13 #get from a clear case or similar
#f=open('fdtd_netlist.cir','w')
#f.write('*Fully Stabilized Ideal Non-Foster Co-Simulation Netlist\n\n')
#f.write('.options method=gear reltol=1e-4 chgtol=1e-15\n\n')
#f.write('I2 0 out2 PWL(0 0)\n\n')
#f.write('C_interface out2 match_node 1pF\n\n')
#f.write('R_lp match_node v_filter 1\n')
#f.write('C_lp v_filter 0 1pF\n\n')
#f.write('B_neg_cap 0 match_node I=-0.2p * ddt(v(v_filter))\n')
#f.write('R_dc_anchor match_node 0 1meg\n\n')
#f.write('R2 match_node 0 50\n\n')
#f.write('.tran 1E-13 {} UIC\n'.format(times * base_time_step / factor))
#f.write('.end\n')
#f.close()

#kmax example for gridded ports:
#times=4000 #get from above
#factor=1 #get from above
#base_time_step=9.4365835E-13 #get from a clear case or similar
#f=open('fdtd_netlist.cir','w')
#f.write('*Validation of coax operation with FDTD\n')
#f.write('I2 0 out2 PWL(0 0)\n')
#f.write('I3 0 out3 PWL(0 0)\n')
#f.write('.tran 1E-13 {}\n'.format(times*base_time_step/factor))
#f.write('R2 out2 0 36.68\n')
#f.write('R3 out3 0 36.68\n')
#f.write('.end\n')
#f.close()

#kmax example w/ basic ports
#times=4000 #get from above
#factor=1 #get from above
#base_time_step=9.4365835E-13 #get from a clear case or similar
#f=open('fdtd_netlist.cir','w')
#f.write('*Validation of two 50 ohm ports with FDTD\n')
#f.write('I1 0 out1 PWL(0 0)\n')
#f.write('I2 0 out2 PWL(0 0)\n')
#f.write('I12 0 out12 PWL(0 0)\n')
#f.write('I22 0 out22 PWL(0 0)\n')
#f.write('.tran 1E-13 {}\n'.format(times*base_time_step/factor))
#f.write('Cp1 out1 0 1\n') #Cp1 and Cp2 values are overwritten by FDTD solver - they are grid capacitances
#f.write('Cp2 out2 0 1\n')
#f.write('Cp12 out12 0 1\n') #Cp12 and Cp22 values are overwritten by FDTD solver - they are grid capacitances
#f.write('Cp22 out22 0 1\n')
#f.write('Rterm1 out1 0 50\n')
#f.write('Rterm2 out2 0 50\n')
#f.write('Rterm12 out12 0 50\n')
#f.write('Rterm22 out22 0 50\n')
#f.write('.end\n')
#f.close()

##runs fdtd and clocks the time to run in minutes
start_time=tp.time()
os.system('./{} {}'.format(solver,inputs_name))
end_time=tp.time()
print('Total simulation time is ',(end_time-start_time)/60, ' minutes') #in minutes time to run everything