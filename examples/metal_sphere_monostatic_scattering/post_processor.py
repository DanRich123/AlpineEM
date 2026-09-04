import numpy as np
import os
from matplotlib import pyplot as plt
import time as tp
start_time=tp.time()

#SEE LICENSE FILE FOR LICENSE INFORMATION

#This is for the normal + optional spice
#see other version for kmax + optional spice

####################################################################
####SETUP INFORMATION HERE##########################################
####ONLY THESE LINES NEED CHANGING TO RUN###########################
#os.chdir('C:/Users/foldername')
simulation_type='plane wave' #'antenna' or 'plane wave' type of simulation
#input 4 file names - not all will necessarily be used depending on several factors
data_filename='data.dat'
clear_filename='clear.dat' # clear but can/should include IGP if present in data.dat case
metal_filename='metal.dat'
IGP_full_clear_filename='IGP_full_clear.dat' # full clear, no IGP even if present in data.dat - used for far field inc amptlitude and phase centering
#padding factor (integer) for increasing output # of points - time domain interpolation via padding at the end of the time sequence
#it can cause ripples in output data if padding>0 and it's EM fields are not converged in the time domain
padding=10
#decide if metal should be used (only for unit cells, otherwise ignored) - effects phase centering - useful for measurement comparison and higher fidelity accuracy of amplitudes
use_metal=False
#IGP full clear needed for far field cal if IGP used (only for plane waves, otherwise ignored) - used for incident wave amplitude and phase centering for far fields
use_IGP_full_clear=False
#output file names - not all will ncessarily be used
S_parameter_output_file_name='S_parameters.csv'
Scattering_output_file_name='Scattering_Far_Field.csv'
Antenna_gain_output_file_name='Realized_antenna_gain.csv'
#NOTE, spice port impedance is not known a priori, the post processor defaults to 50 ohms and makes a note in the csv output files that it needs to be modified directly by the user if needed.
#NOTE, if using spice port excitation, this needs to be uploaded here so that this script knows to use that incident wave instead of the .dat file's zero'd out wave.
spice_excitation=False
#if spice_excitation is True specify these two items, otherwise not needed and unused.
spice_port_number=1 # tell the script which of the spice ports it is (1,2,3,...n) for the submission order of the spice ports into the master.py
if spice_excitation==True:
    spice_incident_wave=np.load('incident.npy') # load in the incident wave, at the right fdtd time steps (1-xxx), as a numpy array
####################################################################
####################################################################

########IMPORT INFORMATION FROM .DAT FILES##########################
wave_port_num=2 #default, can get overwritten
time_step=np.loadtxt(data_filename,usecols=3,skiprows=4,max_rows=1)
num_time_steps=int(np.loadtxt(data_filename,usecols=5,skiprows=5,max_rows=1))
pulse_freq=np.loadtxt(data_filename,usecols=3, skiprows=7,max_rows=1)
pulse_freq=pulse_freq/1E9
pulse_type=np.loadtxt(data_filename, usecols=3, skiprows=8,max_rows=1)
if pulse_type==1:
    pulse_bounds=[0.0,2.146]
if pulse_type==2:
    pulse_bounds=[0.183,3.087]
num_angles=int(np.loadtxt(data_filename, usecols=6, skiprows=13, max_rows=1))
num_ports=int(np.loadtxt(data_filename, usecols=5, skiprows=14, max_rows=1))
num_spice_ports=int(np.loadtxt(data_filename, usecols=5, skiprows=15, max_rows=1))
sim_type_x=int(np.loadtxt(data_filename,usecols=2,skiprows=16,max_rows=1))
sim_type_y=int(np.loadtxt(data_filename,usecols=3,skiprows=16,max_rows=1))
sim_type_z=int(np.loadtxt(data_filename,usecols=4,skiprows=16,max_rows=1))
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
if (sim_type_x+sim_type_y+sim_type_z<2):
    wave_port_num=0
    area=1 #output using this is then per area - effective aperature area
imped_free=376.730313
maxrows=num_time_steps
if simulation_type=='antenna':
    wave_port_num=0
    port_number_str = np.loadtxt(data_filename,usecols=5, skiprows=17,max_rows=1)
    port_number_int = int(port_number_str)
    if spice_excitation==True:
        port_number_int = num_ports + spice_port_number
#######################################################################

######First import all data into arrays and setup any needed arrays####
inc_data=np.loadtxt(data_filename, skiprows=18,max_rows=maxrows)
inc_data=inc_data.reshape(-1)

if spice_excitation==True:
    inc_data=spice_incident_wave

wave_port_time_array=np.zeros((3,wave_port_num,num_time_steps))
far_field_time_array=np.zeros((2,int(num_angles*2),num_time_steps))
v_ports_time_array=np.zeros((int(num_ports+num_spice_ports), num_time_steps))
impedance_array=[]
far_field_angles=[]

#create skip_start for skipping everything up to this
skip_start=18+maxrows

for i in range(wave_port_num):
    skip_start=skip_start+1
    temp=np.loadtxt(data_filename,skiprows=skip_start, max_rows=maxrows)
    wave_port_time_array[0,i]=temp.reshape(-1)
    if simulation_type=='plane wave':
        temp_clear=np.loadtxt(clear_filename,skiprows=skip_start, max_rows=maxrows)
        wave_port_time_array[1,i]=temp_clear.reshape(-1)
        if use_metal==True:
            temp_metal=np.loadtxt(metal_filename,skiprows=skip_start, max_rows=maxrows)
            wave_port_time_array[2,i]=temp_metal.reshape(-1)
    skip_start+=maxrows

#far field E fields w/ 1 added for first far field header
for i in range(int(num_angles*2)):
    skip_start=skip_start+1
    far1=np.loadtxt(data_filename,usecols=0, skiprows=skip_start,max_rows=1)
    far2=np.loadtxt(data_filename,usecols=1, skiprows=skip_start,max_rows=1)
    if i % 2 == 0:
        far_field_angles.append([far1,far2])
    skip_start=skip_start+1
    temp=np.loadtxt(data_filename,skiprows=skip_start, max_rows=maxrows)
    far_field_time_array[0,i]=temp.reshape(-1)
    if simulation_type=='plane wave':
        temp=np.loadtxt(clear_filename,skiprows=skip_start, max_rows=maxrows)
        far_field_time_array[1,i]=temp.reshape(-1)
    skip_start+=maxrows

#clear for far field phase correction - self contained
skip_start_clear=skip_start+1
time_corrections_clear=[]
for i in range(int(num_angles)):
    if simulation_type=='plane wave':
        time_corrections_clear.append(np.loadtxt(clear_filename,skiprows=skip_start_clear,max_rows=1))
        skip_start_clear+=1    
#now for the incident field for phase centering - only applies to plane waves and we only need it from the clear case
if (num_angles>0):
    if simulation_type=='plane wave':
        if (use_IGP_full_clear==False):
            inc_ff_pc=np.loadtxt(clear_filename,skiprows=skip_start_clear+1,max_rows=num_time_steps)
            inc_ff_pc=inc_ff_pc.reshape(-1)
        if (use_IGP_full_clear==True):
            inc_ff_pc=np.loadtxt(IGP_full_clear_filename,skiprows=skip_start_clear+1,max_rows=num_time_steps)
            inc_ff_pc=inc_ff_pc.reshape(-1)

#voltages - w/ 1 added to skip first header
for i in range(int(num_ports+num_spice_ports)):
    skip_start=skip_start+1
    impedance_array.append(np.loadtxt(data_filename,skiprows=skip_start,max_rows=1))
    skip_start=skip_start+1
    temp=np.loadtxt(data_filename,skiprows=skip_start, max_rows=maxrows)
    v_ports_time_array[i]=temp.reshape(-1)
    skip_start+=maxrows

#data for far field phase correction - self contained
skip_start=skip_start+1
time_corrections=[]
for i in range(int(num_angles)):
    time_corrections.append(np.loadtxt(data_filename,skiprows=skip_start,max_rows=1))
    skip_start+=1
#######################################################################

######Pad all and concatenate as needed################################
pad=np.linspace(0,0,len(inc_data)*padding)
inc_data=np.concatenate((inc_data,pad))

wave_port_time2_array=np.zeros((3,wave_port_num,int((padding+1)*num_time_steps)))
far_field_time2_array=np.zeros((2,int(num_angles*2),int((padding+1)*num_time_steps)))
v_ports_time2_array=np.zeros((int(num_ports+num_spice_ports),int((padding+1)*num_time_steps)))

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
        inc_ff_pc=np.concatenate((inc_ff_pc,pad))
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

#clear case incident for phase centering if it exists:
if (num_angles>0):
    if simulation_type=='plane wave':
        inc_ff_pc_freq=np.fft.fft(inc_ff_pc)
#######################################################################
       
####Calculations section###############################################
##first S parameters and then far field quantities of interest
s_parameters_array=[]
scatter_array=[]

if simulation_type=='plane wave':
    #first waveports
    for i in range(wave_port_num):
        #first is reflection and second is transmission
        #if metal is not used - then the array was defaulted to zero so all is fine - it becomes as if it's unused
        if i==0:
            s_parameters_array.append(-1*(wave_port_freq_array[0][i]-wave_port_freq_array[1][i])/(wave_port_freq_array[2][i]-wave_port_freq_array[1][i]))
        if i>0:
            s_parameters_array.append(-1*(wave_port_freq_array[0][i]-wave_port_freq_array[2][i])/(wave_port_freq_array[1][i]-wave_port_freq_array[2][i]))
    #then lumped ports w/ the right phase centering
    for i in range(num_ports+num_spice_ports):
        port_pwave=v_ports_freq_array[i]/np.sqrt(np.real(impedance_array[i]))
        if (wave_port_num!=0):
            inc_place_holder=wave_port_freq_array[1][0]
            waveport_pwave=inc_place_holder*np.sqrt(area/(np.real(imped_free)))
        else:
            waveport_pwave=inc_data_f*np.sqrt(area/(np.real(imped_free))) #area set to 1 for this calculation only
        s_parameters_array.append(port_pwave/waveport_pwave)
    
    #now scattered far fields
    timeplace=0
    for i in range(int(2*num_angles)):
        if i > 0 and i % 2 == 0:
            timeplace+=1
        normal_phase=np.exp(1j*2*np.pi*freq*(time_corrections[timeplace])-1j*np.pi/2)
        clear_phase=np.exp(1j*2*np.pi*freq*(time_corrections_clear[timeplace])-1j*np.pi/2)
        scatter_array.append((far_field_freq_array[0][i]*normal_phase-far_field_freq_array[1][i]*clear_phase)/(inc_ff_pc_freq))

if simulation_type=='antenna':
    #first waveports
    for i in range(wave_port_num):
        port_pwave=inc_data_f/np.sqrt(np.real(impedance_array[port_number_int-1]))
        waveport_pwave=wave_port_freq_array[0,i]*np.sqrt(area/(np.real(imped_free)))
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
low=pulse_bounds[0]*pulse_freq*1E9
high=pulse_bounds[1]*pulse_freq*1E9
record_freq=[]
for i in range(len(freq)):
    if freq[i]>=low and freq[i]<=high:
        record_freq.append(i)
l=np.min(record_freq)
h=np.max(record_freq)

#start with S parameter data, if any exist
if (num_ports+num_spice_ports+wave_port_num)>0:
    out_put_data=[]
    out_put_data.append(freq[l:h+1]/1E9)
    for i in range(wave_port_num+num_ports+num_spice_ports):
        out_put_data.append(20*np.log10(np.abs(s_parameters_array[i][l:h+1])))
        out_put_data.append(np.angle(s_parameters_array[i][l:h+1]))

    header_parts = ["Frequency (GHz)"]
    if wave_port_num>0:
        if simulation_type=='plane wave':
            header_parts.append(f"Wave Refl (dB)")
            header_parts.append(f"Wave Refl (rad)")
            header_parts.append(f"Wave Trans (dB)")
            header_parts.append(f"Wave Trans (rad)")
        if simulation_type=='antenna':
            if spice_excitation==False:
                header_parts.append(f"S_waveport1:lumpedport{port_number_int} (dB)")
                header_parts.append(f"S_waveport1:lumpedport{port_number_int} (rad)")
                header_parts.append(f"S_waveport2:lumpedport{port_number_int} (dB)")
                header_parts.append(f"S_waveport2:lumpedport{port_number_int} (rad)")
            if spice_excitation==True:
                header_parts.append(f"S_waveport1:S_spiceport{spice_port_number} (dB) - needs renormalizing Z_spice")
                header_parts.append(f"S_waveport1:S_spiceport{spice_port_number} (rad) - needs renormalizing Z_spice")
                header_parts.append(f"S_waveport2:S_spiceport{spice_port_number} (dB) - needs renormalizing Z_spice")
                header_parts.append(f"S_waveport2:S_spiceport{spice_port_number} (rad) - needs renormalizing Z_spice")
    for i in range(num_ports):
        if simulation_type=='plane wave':
            if (wave_port_num!=0):
                s_parameter_name=f"S_lumpedport{i+1}:inc_waveport"
            else:
                s_parameter_name=f"A_e for port {i+1} - Vport/Einc phase center info: Vport fdtd grid (centered) port location & Einc location"
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
            if (wave_port_num!=0):
                s_parameter_name=f"S_spiceport{i+1}:inc_waveport"
            else:
                s_parameter_name=f"A_e for spice port {i+1} - Vport/Einc phase center info: Vport fdtd grid (centered) port location & Einc location"
            header_parts.append(f"{s_parameter_name} (dB) - needs renormalizaing Z_spice")
            header_parts.append(f"{s_parameter_name} (rad) - needs renormalizaing Z_spice")
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