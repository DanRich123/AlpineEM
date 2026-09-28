import numpy as np

time_step = 1.8873167E-12 
num_steps = 15000
time = np.linspace(1,15000,15000)
time = time*time_step
source = 0.5*0.001 * np.sin(2*np.pi*10E9*time) * np.exp(-np.log(2) * ((time - 150E-12)/50E-12)**2)
np.save('incident', source)