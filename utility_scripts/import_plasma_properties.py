import os
import time as tp
import numpy as np

#this script just shows how to import plasma properties from import_tecplot_and_interpolate.py
#it also shows how you would submit the materials with a few additional ones that might be relevant

#import some processed data
materials_id=np.load('material_properties.npy')
#num plasmas includes the radome number too so it's 1 extra
num_plasmas=len(materials_id[:])-1
#each index is a set of material properties for all species
#these lines can be used if the loss is the same for all species
#materials_id_omega_p_sum=np.linspace(0,0,num_plasmas)
#for i in range(num_plasmas):
    #materials_id_omega_p_sum[i]=np.sqrt(materials_id[i,0]**2+materials_id[i,2]**2+materials_id[i,4]**2+materials_id[i,6]**2+materials_id[i,8]**2+materials_id[i,10]**2)
#print(materials_id_omega_p_sum)

f.write('{}\n'.format(num_plasmas+4)) #number of materials
#do these lines for each bulk volume material
for i in range(num_plasmas):
    f.write('volume\n')
    f.write('plasma\n') #type of material
    f.write('{}\n'.format(i+1)) #material idenfitication number
    f.write('1,1,1\n') #relative epsilon
    f.write('0.0,0.0,0.0\n') #sigma
    f.write('6\n') #number of poles if not 'simple'
    f.write('{},{},{}\n'.format(materials_id[i,0],materials_id[i,0],materials_id[i,0])) #plasma angular frequency term
    f.write('{},{},{}\n'.format(materials_id[i,1],materials_id[i,1],materials_id[i,1])) #plasma angular collision frequency term
    f.write('{},{},{}\n'.format(materials_id[i,2],materials_id[i,2],materials_id[i,2])) #plasma angular frequency term
    f.write('{},{},{}\n'.format(materials_id[i,3],materials_id[i,3],materials_id[i,3])) #plasma angular collision frequency term
    f.write('{},{},{}\n'.format(materials_id[i,4],materials_id[i,4],materials_id[i,4])) #plasma angular frequency term
    f.write('{},{},{}\n'.format(materials_id[i,5],materials_id[i,5],materials_id[i,5])) #plasma angular collision frequency term
    f.write('{},{},{}\n'.format(materials_id[i,6],materials_id[i,6],materials_id[i,6])) #plasma angular frequency term
    f.write('{},{},{}\n'.format(materials_id[i,7],materials_id[i,7],materials_id[i,7])) #plasma angular collision frequency term
    f.write('{},{},{}\n'.format(materials_id[i,8],materials_id[i,8],materials_id[i,8])) #plasma angular frequency term
    f.write('{},{},{}\n'.format(materials_id[i,9],materials_id[i,9],materials_id[i,9])) #plasma angular collision frequency term
    f.write('{},{},{}\n'.format(materials_id[i,10],materials_id[i,10],materials_id[i,10])) #plasma angular frequency term
    f.write('{},{},{}\n'.format(materials_id[i,11],materials_id[i,11],materials_id[i,11])) #plasma angular collision frequency term
#radome
f.write('volume\n')
f.write('simple\n')
f.write('{}\n'.format(num_plasmas+1)) #material idenfitication number
f.write('3,3,3\n') #relative epsilon
f.write('0.001,0.001,0.001\n') #sigma
#metal
f.write('volume\n')
f.write('simple\n')
f.write('{}\n'.format(num_plasmas+2)) #material idenfitication number
f.write('1,1,1\n') #relative epsilon
f.write('1E8,1E8,1E8\n') #sigma
#do these lines for each sheet material
f.write('sheet\n')
f.write('{}\n'.format(num_plasmas+3)) #material idenfitication number for sheets
f.write('0.01\n') #surface impedance
#air
f.write('volume\n')
f.write('simple\n')
f.write('{}\n'.format(num_plasmas+4)) #material idenfitication number
f.write('1,1,1\n') #relative epsilon
f.write('0,0,0\n') #sigma
