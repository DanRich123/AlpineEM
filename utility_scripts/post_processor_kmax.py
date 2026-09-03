import numpy as np
import os
from matplotlib import pyplot as plt
import time as tp
start_time=tp.time()

#SEE LICENSE FILE FOR LICENSE INFORMATION

#This is for the kmax + optional spice
#see other version for normal + optional spice

####################################################################
####SETUP INFORMATION HERE##########################################
####ONLY THESE LINES NEED CHANGING TO RUN###########################
#os.chdir('C:/Users/foldername')
simulation_type='plane wave' #'antenna' or 'plane wave' type of simulation
#input 3 file names - not all will necessarily be used
data_filename='data-k.dat'
clear_filename='clear-k.dat'
metal_filename='metal-k.dat'
#padding factor (integer) for increasing output # of points - time domain interpolation via padding at the end of the time sequence
#it can cause ripples in output data if padding>0 and it's EM fields are not converged in the time domain
padding=0
#decide if metal should be used - effects phase centering - useful for measurement comparison
use_metal=False
#output file names - not all will ncessarily be used
S_parameter_output_file_name='S_parameters_k.csv'
Scattering_output_file_name='Scattering_Far_Field_k.csv'
Antenna_gain_output_file_name='Realized_antenna_gain_k.csv'
#select dB rule - larger values get you more data from a single simulation
#20 is standard and trustable
#30 becomes more qualitative at wider freq and larger angles
#40 is very qualitative at wider freq and larger angles
dB_rule=20
#NOTE, spice port impedance is not known a priori, the post processor defaults to 50 ohms and makes a note in the csv output files that it needs to be modified directly by the user if needed.
#NOTE, if using spice port excitation, this needs to be uploaded here so that this script knows to use that incident wave instead of the .dat file's zero'd out wave.
spice_excitation=False
#if spice_excitation is True specify these two items, otherwise not needed and unused.
spice_port_number=1 # tell the script which of the spice ports it is (1,2,3,...n) for the submission order of the spice ports into the master.py - recall the pair submission is removed here (i.e port 3&4 in fdtd is port 2 here)
if spice_excitation==True:
    spice_incident_wave=np.load('incident.npy') # load in the incident wave, at the right fdtd time steps (1-xxx), as a numpy array
####################################################################
####################################################################

########IMPORT INFORMATION FROM .DAT FILES##########################
time_step=np.loadtxt(data_filename,usecols=3,skiprows=4,max_rows=1)
num_time_steps=int(np.loadtxt(data_filename,usecols=5,skiprows=5,max_rows=1))
pulse_freq=np.loadtxt(data_filename,usecols=3, skiprows=7,max_rows=1)
num_angles=int(np.loadtxt(data_filename, usecols=6, skiprows=13, max_rows=1))
num_ports=int(np.loadtxt(data_filename, usecols=5, skiprows=14, max_rows=1))
num_spice_ports=int(np.loadtxt(data_filename, usecols=5, skiprows=15, max_rows=1))
wave_port_num=4 #two TE and two TM regardless of excitation type
sim_type_x=int(np.loadtxt(data_filename,usecols=2,skiprows=16,max_rows=1))
sim_type_y=int(np.loadtxt(data_filename,usecols=3,skiprows=16,max_rows=1))
sim_type_z=int(np.loadtxt(data_filename,usecols=4,skiprows=16,max_rows=1))
k_direction=int(np.loadtxt(data_filename,usecols=5,skiprows=16,max_rows=1))
size_x=int(np.loadtxt(data_filename,usecols=3,skiprows=0,max_rows=1))
size_y=int(np.loadtxt(data_filename,usecols=3,skiprows=1,max_rows=1))
size_z=int(np.loadtxt(data_filename,usecols=3,skiprows=2,max_rows=1))
step_size_y=np.loadtxt(data_filename,usecols=4, skiprows=3,max_rows=1)
step_size_z=np.loadtxt(data_filename,usecols=5, skiprows=3,max_rows=1)
step_size_x=np.loadtxt(data_filename,usecols=3,skiprows=3,max_rows=1)
if (sim_type_y+sim_type_z==2):
    area=size_z*step_size_z*size_y*step_size_y
if (sim_type_x+sim_type_z==2):
    area=size_z*step_size_z*size_x*step_size_x
if (sim_type_x+sim_type_y==2):
    area=size_y*step_size_y*size_x*step_size_x
imped_free=376.730313
k1=np.loadtxt(data_filename,usecols=0,skiprows=18,max_rows=1)
k2=np.loadtxt(data_filename,usecols=1,skiprows=18,max_rows=1)
mode_type=np.loadtxt(data_filename,usecols=10,skiprows=19,max_rows=1)
if simulation_type=='antenna':
    port_number_str = np.loadtxt(data_filename,usecols=5, skiprows=20,max_rows=1)
    port_number_int = int(port_number_str)
    if spice_excitation==True:
        port_number_int = num_ports + spice_port_number
maxrows=num_time_steps
######################################################################

######First import all data into arrays and setup any needed arrays####
inc_data=np.loadtxt(data_filename, skiprows=21, usecols=0, max_rows=maxrows)
inc_data=inc_data.reshape(-1)

inc_data_imag=np.loadtxt(data_filename, skiprows=21, usecols=1, max_rows=maxrows)
inc_data_imag=inc_data_imag.reshape(-1)
inc_data=inc_data+1j*inc_data_imag

if spice_excitation==True:
    inc_data=spice_incident_wave

#plt.figure()
#plt.plot(np.real(inc_data))
#plt.plot(np.imag(inc_data))
#plt.savefig('inc-kmax')

wave_port_time_array=np.zeros((3,wave_port_num,num_time_steps), dtype=complex)
far_field_time_array=np.zeros((2,int(num_angles*2),num_time_steps), dtype=complex)
v_ports_time_array=np.zeros((int(num_ports+num_spice_ports), num_time_steps), dtype=complex)
impedance_array=[]
far_field_angles=[]

#create skip_start for skipping everything up to this
skip_start=21+maxrows
#plt.figure()
for i in range(wave_port_num):
    skip_start=skip_start+1
    temp=np.loadtxt(data_filename,skiprows=skip_start, usecols=0, max_rows=maxrows)
    temp2=temp.reshape(-1)
    #plt.plot(temp2)
    temp_imag=np.loadtxt(data_filename,skiprows=skip_start, usecols=1, max_rows=maxrows)
    temp2_imag=temp_imag.reshape(-1)
    wave_port_time_array[0,i]=temp2+temp2_imag*1j
    if simulation_type=='plane wave':
        temp=np.loadtxt(clear_filename,skiprows=skip_start, usecols=0, max_rows=maxrows)
        temp2=temp.reshape(-1)
        temp_imag=np.loadtxt(clear_filename,skiprows=skip_start, usecols=1, max_rows=maxrows)
        temp2_imag=temp_imag.reshape(-1)
        #plt.plot(temp2)
        wave_port_time_array[1,i]=temp2+temp2_imag*1j
        if use_metal==True:
            temp=np.loadtxt(metal_filename,skiprows=skip_start, usecols=0, max_rows=maxrows)
            temp2=temp.reshape(-1)
            temp_imag=np.loadtxt(metal_filename,skiprows=skip_start, usecols=1, max_rows=maxrows)
            temp2_imag=temp_imag.reshape(-1)
            wave_port_time_array[2,i]=temp2+temp2_imag*1j
    skip_start+=maxrows

#plt.savefig('waves')

#far field E fields w/ 1 added for first far field header
for i in range(int(num_angles*2)):
    skip_start=skip_start+1
    far1=np.loadtxt(data_filename,usecols=0, skiprows=skip_start,max_rows=1)
    far2=np.loadtxt(data_filename,usecols=1, skiprows=skip_start,max_rows=1)
    if i % 2 == 0:
        far_field_angles.append([far1,far2])
    skip_start=skip_start+1
    temp=np.loadtxt(data_filename,skiprows=skip_start, usecols=0, max_rows=maxrows)
    temp2=temp.reshape(-1)
    temp_imag=np.loadtxt(data_filename,skiprows=skip_start, usecols=1, max_rows=maxrows)
    temp2_imag=temp_imag.reshape(-1)
    far_field_time_array[0,i]=temp2+temp2_imag*1j
    if simulation_type=='plane wave':
        temp=np.loadtxt(clear_filename,skiprows=skip_start, usecols=0, max_rows=maxrows)
        temp2=temp.reshape(-1)
        temp_imag=np.loadtxt(clear_filename,skiprows=skip_start, usecols=1, max_rows=maxrows)
        temp2_imag=temp_imag.reshape(-1)
        far_field_time_array[1,i]=temp2+temp2_imag*1j
    skip_start+=maxrows

#clear for far field phase correction - self contained - dont need to pull from both imag and real - should be the same
skip_start_clear=skip_start+1
time_corrections_clear=[]
for i in range(int(num_angles)):
    if simulation_type=='plane wave':
        time_corrections_clear.append(np.loadtxt(clear_filename,skiprows=skip_start_clear,max_rows=1))
        skip_start_clear+=1
#now for the incident field for phase centering - only applies to plane waves and we only need it from the clear case
if (num_angles>0):
    if simulation_type=='plane wave':
        Ex_real=np.loadtxt(clear_filename,skiprows=skip_start_clear+1,usecols=0,max_rows=num_time_steps)
        Ex_real=Ex_real.reshape(-1)
        Ex_imag=np.loadtxt(clear_filename,skiprows=skip_start_clear+1,usecols=1,max_rows=num_time_steps)
        Ex_imag=Ex_imag.reshape(-1)
        Ex=Ex_real+1j*Ex_imag
        skip_start_clear=skip_start_clear+num_time_steps
        Ey_real=np.loadtxt(clear_filename,skiprows=skip_start_clear+1,usecols=0,max_rows=num_time_steps)
        Ey_real=Ey_real.reshape(-1)
        Ey_imag=np.loadtxt(clear_filename,skiprows=skip_start_clear+1,usecols=1,max_rows=num_time_steps)
        Ey_imag=Ey_imag.reshape(-1)
        Ey=Ey_real+1j*Ey_imag
        skip_start_clear=skip_start_clear+num_time_steps
        Ez_real=np.loadtxt(clear_filename,skiprows=skip_start_clear+1,usecols=0,max_rows=num_time_steps)
        Ez_real=Ez_real.reshape(-1)
        Ez_imag=np.loadtxt(clear_filename,skiprows=skip_start_clear+1,usecols=1,max_rows=num_time_steps)
        Ez_imag=Ez_imag.reshape(-1)
        Ez=Ez_real+1j*Ez_imag

#voltages - w/ 1 added to skip first header
#plt.figure()
for i in range(int(num_ports+num_spice_ports)):
    skip_start=skip_start+1
    impedance_array.append(np.loadtxt(data_filename,skiprows=skip_start, max_rows=1))
    skip_start=skip_start+1
    temp=np.loadtxt(data_filename,skiprows=skip_start, usecols=0, max_rows=maxrows)
    temp2=temp.reshape(-1)
    #plt.plot(temp2)
    temp_imag=np.loadtxt(data_filename,skiprows=skip_start, usecols=1, max_rows=maxrows)
    temp2_imag=temp_imag.reshape(-1)
    v_ports_time_array[i]=temp2+temp2_imag*1j
    skip_start+=maxrows
#plt.savefig('volts')

#data for far field phase correction - self contained - imag and real should be the same
skip_start=skip_start+1
time_corrections=[]
for i in range(int(num_angles)):
    time_corrections.append(np.loadtxt(data_filename,skiprows=skip_start,max_rows=1))
    skip_start+=1
#######################################################################

######Pad all and concatenate as needed################################
pad=np.linspace(0+0*1j,0+0*1j,len(inc_data)*padding)
inc_data=np.concatenate((inc_data,pad))

wave_port_time2_array=np.zeros((3,wave_port_num,int((padding+1)*num_time_steps)), dtype=complex)
far_field_time2_array=np.zeros((2,int(num_angles*2),int((padding+1)*num_time_steps)),dtype=complex)
v_ports_time2_array=np.zeros((int(num_ports+num_spice_ports),int((padding+1)*num_time_steps)),dtype=complex)

#first waveports
for i in range(wave_port_num):
    for j in range(3):
        wave_port_time2_array[j][i]=np.concatenate((wave_port_time_array[j][i],pad))

#far field E fields
for i in range(int(num_angles*2)):
    for j in range(2):
        far_field_time2_array[j][i]=np.concatenate((far_field_time_array[j][i],pad))

#voltages
for i in range(int(num_ports+num_spice_ports)):
    v_ports_time2_array[i]=np.concatenate((v_ports_time_array[i],pad))

#clear case incident for phase centering if it exists:
if (num_angles>0):
    if simulation_type=='plane wave':
        Ex=np.concatenate((Ex,pad))
        Ey=np.concatenate((Ey,pad))
        Ez=np.concatenate((Ez,pad))
#######################################################################

#####FFT all data sets#################################################
inc_data_f=np.fft.fft(inc_data)
freq=np.fft.fftfreq(len(inc_data),time_step)

wave_port_freq_array=np.zeros((3,wave_port_num,int(num_time_steps*(padding+1))), dtype=complex)
far_field_freq_array=np.zeros((2,int(num_angles*2), int((padding+1)*num_time_steps)), dtype=complex)
v_ports_freq_array=np.zeros((int(num_ports+num_spice_ports), int((padding+1)*num_time_steps)), dtype=complex)

#first waveports
for i in range(wave_port_num):
    for j in range(3):
        wave_port_freq_array[j][i]=np.fft.fft(wave_port_time2_array[j][i])

#far field E fields
for i in range(int(num_angles*2)):
    for j in range(2):
        far_field_freq_array[j][i]=np.fft.fft(far_field_time2_array[j][i])
    
#voltages
for i in range(int(num_ports+num_spice_ports)):
    v_ports_freq_array[i]=np.fft.fft(v_ports_time2_array[i])
#######################################################################
       
####Calculations section###############################################
##first S parameters and then far field quantities of interest
s_parameters_array=[]
scatter_array=[]
k = 2 * np.pi * freq / 299792458
cos_ang = np.ones_like(k, dtype=float)
if k1 != 0 or k2 != 0:
    arg = np.maximum(k**2 - k1**2 - k2**2, 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        cos_ang = np.where(k != 0, np.sqrt(arg) / k, 1)
cos_ang = np.clip(cos_ang, 1e-12, 1.0)  # Prevent zero or negative values, I won't use the invalid regions anyway later on

if simulation_type=='plane wave':
    #first waveports
    #notes# 
    #recieved-clear for that port type then divide by clear - metal for the inc port type for the respective refl or trans
    #if metal is not used, it's set to zero so no harm in keeping as is, just wanted time and memory
    if mode_type==0: #TE inc
        #TE and then TM refl
        s_parameters_array.append(-1*(wave_port_freq_array[0][0]-wave_port_freq_array[1][0])/(wave_port_freq_array[2][0]-wave_port_freq_array[1][0]))
        s_parameters_array.append(-1*(wave_port_freq_array[0][1]-wave_port_freq_array[1][1])/(wave_port_freq_array[2][0]-wave_port_freq_array[1][0])/cos_ang)
        #TE amd then TM trans
        s_parameters_array.append(-1*(wave_port_freq_array[0][2]-wave_port_freq_array[2][2])/(wave_port_freq_array[1][2]-wave_port_freq_array[2][2]))
        s_parameters_array.append(-1*(wave_port_freq_array[0][3]-wave_port_freq_array[2][3])/(wave_port_freq_array[1][2]-wave_port_freq_array[2][2])/cos_ang)
        #need an incident wave equivlanet since r(k) and we don't want to save entire array from fortran
        #this uses clear of inc pol at the plane of measurement for de-embedding
        inc_place_holder=wave_port_freq_array[1][0]
    if mode_type==1: #TM inc
        #TE and then TM refl
        s_parameters_array.append(-1*(wave_port_freq_array[0][0]-wave_port_freq_array[1][0])/(wave_port_freq_array[2][1]-wave_port_freq_array[1][1])*cos_ang)
        s_parameters_array.append(-1*(wave_port_freq_array[0][1]-wave_port_freq_array[1][1])/(wave_port_freq_array[2][1]-wave_port_freq_array[1][1]))
        #TE amd then TM trans
        s_parameters_array.append(-1*(wave_port_freq_array[0][2]-wave_port_freq_array[2][2])/(wave_port_freq_array[1][3]-wave_port_freq_array[2][3])*cos_ang)
        s_parameters_array.append(-1*(wave_port_freq_array[0][3]-wave_port_freq_array[2][3])/(wave_port_freq_array[1][3]-wave_port_freq_array[2][3]))
        #need an incident wave equivlanet since r(k) and we don't want to save entire array from fortran
        #this uses clear of inc pol at the plane of measurement for de-embedding
        inc_place_holder=wave_port_freq_array[1][1]
    #then lumped ports
    for i in range(num_ports+num_spice_ports):
        port_pwave=v_ports_freq_array[i]/np.sqrt(np.real(impedance_array[i]))
        if mode_type==0: #then it's a TE port
            waveport_pwave=inc_place_holder*np.sqrt(area/(np.real(imped_free)))
        if mode_type==1: #then it's a TM port
            waveport_pwave=inc_place_holder*np.sqrt(area/(np.real(imped_free)*cos_ang**2))
        s_parameters_array.append(port_pwave/waveport_pwave)

    #now scattered far fields
    #first we determine the incident field specifically from the clear case phase centering 
    #this is a lot of work but allows us to avoid manual phase corrections.
    #can always revert back to inc_data_f if we need the right magnitude for sure but don't care about phase centering
    import numpy.lib.scimath as scimath
    k3=scimath.sqrt(k**2 - k1**2 - k2**2)
    k12=np.sqrt(k1**2+k2**2)
    for i in range(len(k)):
        if k[i]==0:
            k[i]=1E-6
    if (sim_type_y+sim_type_z==2):
        if mode_type==0: #then TE so E fields are added to H fields where all H trans fields are always positive when k is positive, H is always negative when k is negative
            if (k1==0 and k2==0):
                coeff_ex=0
                coeff_ey=2*k_direction-1
                coeff_ez=0
            else:
                coeff_ex=0
                coeff_ey=k2/k12*(2*k_direction-1)
                coeff_ez=-k1/k12*(2*k_direction-1)
        if mode_type==1: #then TM so H fields are added to E fields where all E trans fields are always positive when k is positive, E is always negative when k is negative
            if (k1==0 and k2==0):
                coeff_ex=0
                coeff_ey=0
                coeff_ez=2*k_direction-1
            else:
                coeff_ex=k12/k*(2*k_direction-1)*-1.0
                coeff_ey=k3*k1/(k*k12)
                coeff_ez=k3*k2/(k*k12)
    if (sim_type_x+sim_type_z==2):
        if mode_type==0: #then TE so E fields are added to H fields where all H trans fields are always positive when k is positive, H is always negative when k is negative
            if (k1==0 and k2==0):
                coeff_ex=(2*k_direction-1)*-1
                coeff_ey=0
                coeff_ez=0
            else:
                coeff_ex=-k2/k12*(2*k_direction-1)
                coeff_ey=0
                coeff_ez=k1/k12*(2*k_direction-1)
        if mode_type==1: #then TM so H fields are added to E fields where all E trans fields are always positive when k is positive, E is always negative when k is negative
            if (k1==0 and k2==0):
                coeff_ex=0
                coeff_ey=0
                coeff_ez=2*k_direction-1
            else:
                coeff_ex=k3*k1/(k*k12)
                coeff_ey=k12/k*(2*k_direction-1)*-1.0
                coeff_ez=k3*k2/(k*k12)
    if (sim_type_x+sim_type_y==2):
        if mode_type==0: #then TE so E fields are added to H fields where all H trans fields are always positive when k is positive, H is always negative when k is negative
            if (k1==0 and k2==0):
                coeff_ex=0
                coeff_ey=(2*k_direction-1)*-1
                coeff_ez=0
            else:
                coeff_ex=k2/k12*(2*k_direction-1)
                coeff_ey=-k1/k12*(2*k_direction-1)
                coeff_ez=0
        if mode_type==1: #then TM so H fields are added to E fields where all E trans fields are always positive when k is positive, E is always negative when k is negative
            if (k1==0 and k2==0):
                coeff_ex=2*k_direction-1
                coeff_ey=0
                coeff_ez=0
            else:
                coeff_ex=k3*k1/(k*k12)
                coeff_ey=k3*k2/(k*k12)
                coeff_ez=k12/k*(2*k_direction-1)*-1.0

    inc_ff_pc_time = Ex*coeff_ex+Ey*coeff_ey+Ez*coeff_ez
    inc_ff_pc_freq = np.fft.fft(inc_ff_pc_time)

    timeplace=0
    for i in range(int(2*num_angles)):
        if i > 0 and i % 2 == 0:
            timeplace+=1
        normal_phase=np.exp(1j*2*np.pi*freq*(time_corrections[timeplace])-1j*np.pi/2)
        clear_phase=np.exp(1j*2*np.pi*freq*(time_corrections_clear[timeplace])-1j*np.pi/2)
        scatter_array.append((far_field_freq_array[0][i]*normal_phase-far_field_freq_array[1][i]*clear_phase)/(inc_ff_pc_freq))

if simulation_type=='antenna':
    #first waveports
    #note inc wave here is exact so no issues with r(k) like in plane wave case
    for i in range(wave_port_num):
        port_pwave=inc_data_f/np.sqrt(np.real(impedance_array[port_number_int-1]))
        if i==0 or i==2: #then it's a TE port
            waveport_pwave=wave_port_freq_array[0,i]*np.sqrt(area/(np.real(imped_free)))
        if i==1 or i==3: #then it's a TM port
            waveport_pwave=wave_port_freq_array[0,i]*np.sqrt(area/(np.real(imped_free)*cos_ang**2))
        s_parameters_array.append(waveport_pwave/port_pwave) 
    for i in range(num_ports+num_spice_ports):
        if i==(port_number_int-1):
            s_parameters_array.append((v_ports_freq_array[i]-inc_data_f)/inc_data_f)
            #temp=(v_ports_freq_array[i]-inc_data_f)/inc_data_f
            #imped_out=impededance_array[i]*(1+temp)/(1-temp)
        if i!=(port_number_int-1):
            s_parameters_array.append((v_ports_freq_array[i]/np.sqrt(np.real(impedance_array[i])))/(inc_data_f/np.sqrt(np.real(impedance_array[port_number_int-1]))))
    #now gain from scattered field
    for i in range(int(2*num_angles)):
        scatter_array.append((far_field_freq_array[0][i]))

#######################################################################

#####Output the arrays to csv files####################################
#first establish bounds for freq range that is valid and send data out
max_amp=np.max(np.abs(inc_data_f))
# Calculate the dB down boundary.
# A drop of 30 dB means the power is reduced by a factor of 10^(30/10) = 1000
dB_val=10**(dB_rule/10)
db_down_boundary = max_amp / dB_val
# Find the indeces where the y data is above the dB down boundary
index = np.where((np.abs(inc_data_f) >= db_down_boundary) & (freq >= 0.0))
# Find the minimum and maximum indeces
l = np.min(index)
h = np.max(index)

#start with S parameter data, if any exist
if (num_ports+num_spice_ports+wave_port_num)>0:
    out_put_data=[]
    out_put_data.append(freq[l:h+1]/1E9)
    for i in range(wave_port_num+num_ports+num_spice_ports):
        out_put_data.append(20*np.log10(np.abs(s_parameters_array[i][l:h+1])))
        out_put_data.append(np.angle(s_parameters_array[i][l:h+1]))

    header_parts = ["Frequency (GHz)"]
    if simulation_type=='plane wave':
        header_parts.append(f"Wave Refl TE (dB)")
        header_parts.append(f"Wave Refl TE (rad)")
        header_parts.append(f"Wave Refl TM (dB)")
        header_parts.append(f"Wave Refl TM (rad)")
        header_parts.append(f"Wave Trans TE (dB)")
        header_parts.append(f"Wave Trans TE (rad)")
        header_parts.append(f"Wave Trans TM (dB)")
        header_parts.append(f"Wave Trans TM (rad)")
    if simulation_type=='antenna':
        if spice_excitation==False:
            header_parts.append(f"S_TEwaveport1:lumpedport{port_number_int} (dB)")
            header_parts.append(f"S_TEwaveport1:lumpedport{port_number_int} (rad)")
            header_parts.append(f"S_TMwaveport1:lumpedport{port_number_int} (dB)")
            header_parts.append(f"S_TMwaveport1:lumpedport{port_number_int} (rad)")
            header_parts.append(f"S_TEwaveport2:lumpedport{port_number_int} (dB)")
            header_parts.append(f"S_TEwaveport2:lumpedport{port_number_int} (rad)")
            header_parts.append(f"S_TMwaveport2:lumpedport{port_number_int} (dB)")
            header_parts.append(f"S_TMwaveport2:lumpedport{port_number_int} (rad)")
        if spice_excitation==True:
            header_parts.append(f"S_TEwaveport1:S_spiceport{spice_port_number} (dB) - needs renormalizing Z_spice")
            header_parts.append(f"S_TEwaveport1:S_spiceport{spice_port_number} (rad) - needs renormalizing Z_spice")
            header_parts.append(f"S_TMwaveport1:S_spiceport{spice_port_number} (dB) - needs renormalizing Z_spice")
            header_parts.append(f"S_TMwaveport1:S_spiceport{spice_port_number} (rad) - needs renormalizing Z_spice")
            header_parts.append(f"S_TEwaveport2:S_spiceport{spice_port_number} (dB) - needs renormalizing Z_spice")
            header_parts.append(f"S_TEwaveport2:S_spiceport{spice_port_number} (rad) - needs renormalizing Z_spice")
            header_parts.append(f"S_TMwaveport2:S_spiceport{spice_port_number} (dB) - needs renormalizing Z_spice")
            header_parts.append(f"S_TMwaveport2:S_spiceport{spice_port_number} (rad) - needs renormalizing Z_spice")
    for i in range(num_ports):
        if simulation_type=='plane wave':
            s_parameter_name=f"S_lumpedport{i+1}:inc_waveport"
            header_parts.append(f"{s_parameter_name} (dB)")
            header_parts.append(f"{s_parameter_name} (rad)")
        if simulation_type=='antenna':
            if spice_excitation==False:
                s_parameter_name = f"S{i+1}:{port_number_int}"
                header_parts.append(f"{s_parameter_name} (dB)")
                header_parts.append(f"{s_parameter_name} (rad)")
            if spice_excitation==True:
                s_parameter_name = f"S_lumpedport{i+1}:spiceport{spice_port_number}"
                header_parts.append(f"{s_parameter_name} (dB) - needs renormalizing Z_spice")
                header_parts.append(f"{s_parameter_name} (rad) - needs renormalizing Z_spice")
    for i in range(num_spice_ports):
        if simulation_type=='plane wave':
            s_parameter_name=f"S_spiceport{i+1}:inc_waveport"
            header_parts.append(f"{s_parameter_name} (dB) - needs renormalizing Z_spice")
            header_parts.append(f"{s_parameter_name} (rad) - needs renormalizing Z_spice")
        if simulation_type=='antenna':
            if spice_excitation==False:
                s_parameter_name = f"S_spiceport{i+1}:lumpedport{port_number_int}"
                header_parts.append(f"{s_parameter_name} (dB) - needs renormalizing Z_spice")
                header_parts.append(f"{s_parameter_name} (rad) - needs renormalizing Z_spice")
            if spice_excitation==True:
                s_parameter_name = f"S_spiceport{i+1}:spiceport{spice_port_number}"
                header_parts.append(f"{s_parameter_name} (dB) - needs renormalizing Z_spice (unless #n:#n)")
                header_parts.append(f"{s_parameter_name} (rad) - needs renormalizing Z_spice (unless #n:#n)")

    
    header_text = ", ".join(header_parts)
    # Save the arrays with a header, without comment character
    out_put_data=np.transpose(out_put_data)
    np.savetxt(S_parameter_output_file_name, out_put_data, delimiter=",", header=header_text, comments="")

#now do far field parameters of interest
if num_angles>0:
    out_put_data=[]
    out_put_data.append(freq[l:h+1]/1E9)
    header_parts = ["Frequency (GHz)"]
    if simulation_type=='plane wave':
        for i in range(int(num_angles*2)):
                out_put_data.append(20*np.log10(np.abs(scatter_array[i][l:h+1]))+10*np.log10(4*np.pi))
                out_put_data.append(np.angle(scatter_array[i][l:h+1]))
        for i in range(int(num_angles)):
            header_parts.append(f"Theta Pol: Theta={far_field_angles[i][0]} & Phi={far_field_angles[i][1]} (dBsm)")
            header_parts.append(f"Theta Pol: Theta={far_field_angles[i][0]} & Phi={far_field_angles[i][1]} (rad of rEfar/Einc - phase centered at ff origin)") 
            header_parts.append(f"Phi Pol: Theta={far_field_angles[i][0]} & Phi={far_field_angles[i][1]} (dBsm)")
            header_parts.append(f"Phi Pol: Theta={far_field_angles[i][0]} & Phi={far_field_angles[i][1]} (rad of rEfar/Einc - phase centered at ff origin)") 
        header_text = ", ".join(header_parts)
        # Save the arrays with a header, without comment character
        out_put_data=np.transpose(out_put_data)
        np.savetxt(Scattering_output_file_name, out_put_data, delimiter=",", header=header_text, comments="")

    if simulation_type=='antenna':
        P_ant=(np.abs(inc_data_f[l:h+1])**2/(2*np.real(impedance_array[port_number_int-1])))
        timeplace=0
        for i in range(int(num_angles*2)):
            out_put_data.append(10*np.log10(4*np.pi*np.abs(scatter_array[i][l:h+1])**2/(2*imped_free*P_ant)))
            if i > 0 and i % 2 == 0:
                timeplace+=1
            normal_phase=np.exp(1j*2*np.pi*freq[l:h+1]*(time_corrections[timeplace])-1j*np.pi/2)
            out_put_data.append(np.angle(scatter_array[i][l:h+1]*normal_phase/inc_data_f[l:h+1]))
        for i in range(int(num_angles)):
            if spice_excitation==False:
                header_parts.append(f"Theta Pol: Theta={far_field_angles[i][0]} & Phi={far_field_angles[i][1]} (dB)")
                header_parts.append(f"Theta Pol: Theta={far_field_angles[i][0]} & Phi={far_field_angles[i][1]} (rad of rEfar/Vinc - eEfar phase centered at ff origin & Vinc phase centered on fdtd grid port location)") 
                header_parts.append(f"Phi Pol: Theta={far_field_angles[i][0]} & Phi={far_field_angles[i][1]} (dB)")
                header_parts.append(f"Phi Pol: Theta={far_field_angles[i][0]} & Phi={far_field_angles[i][1]} (rad of rEfar/Vinc- eEfar phase centered at ff origin & Vinc phase centered on fdtd grid port location)") 
            if spice_excitation==True:
                header_parts.append(f"Theta Pol: Theta={far_field_angles[i][0]} & Phi={far_field_angles[i][1]} (dB) - needs renormalizing Z_spice")
                header_parts.append(f"Theta Pol: Theta={far_field_angles[i][0]} & Phi={far_field_angles[i][1]} (rad of rEfar/Vinc - eEfar phase centered at ff origin & Vinc phase centered on fdtd grid port location) - needs renormalizing Z_spice") 
                header_parts.append(f"Phi Pol: Theta={far_field_angles[i][0]} & Phi={far_field_angles[i][1]} (dB) - needs renormalizing Z_spice")
                header_parts.append(f"Phi Pol: Theta={far_field_angles[i][0]} & Phi={far_field_angles[i][1]} (rad of rEfar/Vinc- eEfar phase centered at ff origin & Vinc phase centered on fdtd grid port location) - needs renormalizing Z_spice") 
        header_text = ", ".join(header_parts)
        out_put_data=np.transpose(out_put_data)
        # Save the arrays with a header, without comment character
        np.savetxt(Antenna_gain_output_file_name, out_put_data, delimiter=",", header=header_text, comments="")
    
########################################################################################
end_time=tp.time()
print('Time to post process data is {} seconds'.format(end_time-start_time))