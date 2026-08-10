##this is a trial adjoint optimization script
#I will import an S parameter, then use that to determine the adjoint source
#Then I will get the fields from both simulations to determine the gradient using a known formula for this type of problem
import numpy as np
from matplotlib import pyplot as plt
import sys
import os

#load 'truth' simulation - will uses this to be our goal
S=np.loadtxt('S_parameters_truth.csv', delimiter=',', skiprows=1)
S=np.transpose(S)
freq=S[0]
S11_db=S[1]
S11_ang=S[2]
S11_ang = np.unwrap(S11_ang)

#save plot to verify it looks as expected
plt.figure()
plt.plot(freq,S11_db)
plt.savefig('truth_s11')

#switch to real and imaginary parts for convenience
S11=10**(S11_db/20)*np.exp(1j*S11_ang)

#get the S parameter interpolated to only frequency of interest for this example
f_0=5E9
S11_real = np.interp(f_0/1E9, freq, np.real(S11))
S11_imag = np.interp(f_0/1E9, freq, np.imag(S11))
S11_at_f0 = S11_real + 1j * S11_imag
#print(S11_at_f0)

#how many times to iterate
number_of_interations=30
#norms for learning rates
ep_norm=8.85418782E-12 * (2*np.pi*f_0)**2 * (0.5E-3)**2
sigma_norm=(2*np.pi*f_0) * (0.5E-3)**2
#initial guess - isotropic permittivity and conductivity
#doing 2,5 in both layers to check
epsilon=[3,1.5] #relative permittivity
sigma=[0.5,1] #S/m
#value for gradient decent convergence - learning rates - Newton's method circumvents this arbitrary requirement
alpha_ep=[0.2/ep_norm,0.2/ep_norm]
alpha_sigma=[0.2/sigma_norm,0.2/sigma_norm]
#I want history when done
epsilon_history=[]
sigma_history=[]
J_history=[]
F_history=[]

for i in range(number_of_interations):
    #record epsilon and sigma for this iteration
    epsilon_history.append([epsilon[0],epsilon[1]])
    sigma_history.append([sigma[0],sigma[1]])
    #run first fdtd run and post process to create S parameter
    os.system('python fdtd.py {} {} {} {} {} {}'.format(1,0,epsilon[0],sigma[0],epsilon[1],sigma[1]))
    os.system('python post_processor.py')
    #then import S parameter and put in real+im format at f_0
    S_run1=np.loadtxt('S_parameters.csv', delimiter=',', skiprows=1)
    S_run1=np.transpose(S_run1)
    S11_db_run1=S_run1[1]
    S11_ang_run1=S_run1[2]
    S11_ang_run1 = np.unwrap(S11_ang_run1)
    S11_run1=10**(S11_db_run1/20)*np.exp(1j*S11_ang_run1)
    S11_run1 = np.interp(f_0/1E9,freq,np.real(S11_run1)) + 1j*np.interp(f_0/1E9,freq,np.imag(S11_run1))
    #then import the fields for later - will overwrite with data next fdtd run to save memory outside of python
    file_path = 'data_video_Ez.bin'
    #Read the binary file into a NumPy array
    field_run1 = np.fromfile(file_path, dtype=np.float32)
    #need remove first and last elements - not sure why - based on the way binaries are saved initially
    field_run1=field_run1[1:]
    field_run1=field_run1[:-1]
    #needs to be backwards from the way it was made - can do this in a more automatic way but this is fine for now
    field_run1=field_run1.reshape(4000,30,65)
    data_fft_run1=np.fft.fft(field_run1,axis=0)
    #then determine the adjoint source via scoring
    J=(np.abs(S11_at_f0-S11_run1))**2
    J_history.append(J)
    #adjoint source relative to the incident wave amplitude - so I can just scale my original signal
    Adj_source=-2*np.conj(S11_at_f0 - S11_run1)
    Adj_source_A=np.abs(Adj_source)
    Adj_source_phase=np.angle(Adj_source)
    # Ensure the phase is in the range [0, 2pi]
    if Adj_source_phase < 0:
        Adj_source_phase += 2 * np.pi
    time_delay_from_phase = Adj_source_phase / (2 * np.pi * f_0)
    #print(Adj_source_A)
    #print(np.abs(time_delay_from_phase))
    #then run fdtd again with amplitude and phase change to form adjoint source
    os.system('python fdtd.py {} {} {} {} {} {}'.format(Adj_source_A,time_delay_from_phase,epsilon[0],sigma[0],epsilon[1],sigma[1]))
    """
    #if we want the s parameters for the adjoint source, we actually need to run a metal and clear case too
    #would need to restructure with names to avoid overwriting things
    #we don't need them though
    os.system('python post_processor.py')
    #then import S parameter and get real+img at f_0
    S_run2=np.loadtxt('S_parameters.csv', delimiter=',', skiprows=1)
    S_run2=np.transpose(S_run2)
    S11_db_run2=S_run2[1]
    S11_ang_run2=S_run2[2]
    S11_ang_run2 = np.unwrap(S11_ang_run2)
    S11_run2=10**(S11_db_run2/20)*np.exp(1j*S11_ang_run2)
    S11_run2 = np.interp(f_0/1E9,freq,np.real(S11_run2)) + 1j*np.interp(f_0/1E9,freq,np.imag(S11_run2))   
    """ 
    #then import fields again
    #Read the binary file into a NumPy array
    field_run2 = np.fromfile(file_path, dtype=np.float32)
    #need remove first and last elements - not sure why - based on the way binaries are saved initially
    field_run2=field_run2[1:]
    field_run2=field_run2[:-1]
    #needs to be backwards from the way it was made - can do this in a more automatic way but this is fine for now
    field_run2=field_run2.reshape(4000,30,65)
    data_fft_run2=np.fft.fft(field_run2,axis=0)
    #use both fields to get the gradient
    #we are interested in all 2 materials so the math is a little more simple
    n_time = field_run1.shape[0]
    dt = 9.4365835E-13
    freq_vids = np.fft.fftfreq(n_time, dt) / 1E9  # no fftshift
    half = n_time // 2  # round down
    freq_vids_half = freq_vids[:half]  # 0 to just below Nyquist

    #mat 1
    slice_fwd_1 = data_fft_run1[:, 15, 34:39]
    slice_adj_1 = data_fft_run2[:, 15, 34:39]
    overlap_freq_1 = np.sum(slice_fwd_1 * slice_adj_1, axis=1)
    total_overlap_f0_1 = np.interp(f_0/1E9, freq_vids_half, np.real(overlap_freq_1[:half])) + 1j*np.interp(f_0/1E9, freq_vids_half, np.imag(overlap_freq_1[:half]))
    F1 = np.real(total_overlap_f0_1) * (0.5E-3)**3 * 8.85418782E-12 * (2*np.pi*f_0)**2
    F2 = -1.0 * (2*np.pi*f_0) * np.imag(total_overlap_f0_1) * (0.5E-3)**3
    #mat 2
    slice_fwd_2 = data_fft_run1[:, 15, 39:44]
    slice_adj_2 = data_fft_run2[:, 15, 39:44]
    overlap_freq_2 = np.sum(slice_fwd_2 * slice_adj_2, axis=1)
    total_overlap_f0_2 = np.interp(f_0/1E9, freq_vids_half, np.real(overlap_freq_2[:half])) + 1j*np.interp(f_0/1E9, freq_vids_half, np.imag(overlap_freq_2[:half]))
    F3 = np.real(total_overlap_f0_2) * (0.5E-3)**3 * 8.85418782E-12 * (2*np.pi*f_0)**2 
    F4 = -1.0 * (2*np.pi*f_0) * np.imag(total_overlap_f0_2) * (0.5E-3)**3

    #all
    F_history.append([F1,F2,F3,F4])
    #mat 1
    epsilon[0]=epsilon[0]-alpha_ep[0]*F1
    sigma[0]=sigma[0]-alpha_sigma[0]*F2
    #mat 2
    epsilon[1]=epsilon[1]-alpha_ep[1]*F3
    sigma[1]=sigma[1]-alpha_sigma[1]*F4

    #set new epsilon and protect if too low or too high - same for sigma
    for k in range(len(epsilon)):
        if epsilon[k] < 1:
            epsilon[k]=1
        if epsilon[k] > 10:
            epsilon[k]=10
        if sigma[k] < 0:
            sigma[k]=0
        if sigma[k] > 10:
            sigma[k]=10

    #moving plots within the iterations so they will update at each iterations so I can track more easily.
    x=np.linspace(1,i+1,i+1)
    plt.figure()
    plt.plot(x,np.abs(J_history), '-o')
    plt.title('Square Error Vs. Iteration')
    plt.xlabel('Iteration Number')
    plt.ylabel('Square Error (Max of 1.0 possible)')
    plt.yscale('log')
    plt.grid()
    plt.tight_layout()
    plt.savefig('J_history')
    plt.close()

    plt.figure()
    plt.plot(x,np.transpose(F_history)[0]*alpha_ep[0], '-o', label=r'$\epsilon$ #1')
    plt.plot(x,np.transpose(F_history)[1]*alpha_sigma[0], '-o', label=r'$\sigma$ #1')
    plt.plot(x,np.transpose(F_history)[2]*alpha_ep[1], '-o', label=r'$\epsilon$ #2')
    plt.plot(x,np.transpose(F_history)[3]*alpha_sigma[1], '-o', label=r'$\sigma$ #2')
    plt.title('Gradient Vs. Iteration')
    plt.xlabel('Iteration Number')
    plt.ylabel(r'Gradient * learning rate $\alpha$')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('F_history')
    plt.close()

    plt.figure()
    plt.plot(x,np.transpose(epsilon_history)[0], '-o', label=r'$\epsilon$ #1')
    plt.plot(x,np.transpose(epsilon_history)[1], '-o', label=r'$\epsilon$ #2')
    plt.title('Permittivity Vs. Iteration')
    #x_0=np.linspace(1,1,number_of_interations)
    #y=5.0*x_0
    #plt.plot(x,y, '--', label='truth/goal')
    plt.xlabel('Iteration Number')
    plt.ylabel('Permittivity (relative)')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('ep_history')
    plt.close()

    plt.figure()
    plt.plot(x,np.transpose(sigma_history)[0], '-o', label=r'$\sigma$ #1')
    plt.plot(x,np.transpose(sigma_history)[1], '-o', label=r'$\sigma$ #2')
    plt.title('Conductivity Vs. Iteration')
    #x_0=np.linspace(1,1,number_of_interations)
    #y=5.0*x_0
    #plt.plot(x,y, '--', label='truth/goal')
    plt.xlabel('Iteration Number')
    plt.ylabel('Conductivity (S/m)')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('sigma_history')
    plt.close()

print('done')

