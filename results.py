# -*- coding: utf-8 -*-

import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#file_output_TRRM = "./Outputs/cask_lattice_RR_P3_adjoint.json" 
file_output_TRRM = "./Outputs/cask_lattice_LSP2_adjoint.json"    

# returns JSON object as python dictionary
with open(file_output_TRRM) as f:
    dictionary = json.load(f)

# get values
tally   = dictionary.get('fluxMap')
resp    = dictionary.get('Response')
res     = np.array(tally.get('fluxMap'))
flux    = np.sum(res[:,:,:,:,0],3)
flux_z  = flux[17,:,:]
flux_y  = flux[:,12,:]

print (np.shape(res))
print (np.shape(flux))

# convert 8x array into a WW array
weights1  = np.zeros([17,24,64])
weights2  = np.zeros([17,12,64])
weights3  = np.zeros([17,12,32])

for n in range(0,17):
    weights1[n,:,:] = flux[2*n,:,:] + flux[2*n+1,:,:] 

for n in range(0,12):
    weights2[:,n,:] =  weights1[:,2*n,:] + weights1[:,2*n+1,:] 
    
for n in range(0,32):
    weights3[:,:,n] = weights2[:,:,2*n] + weights2[:,:,2*n+1] 

WWmap  = resp/weights3
#WWmap = WWmap/np.max(WWmap)

WW    = np.reshape(WWmap, np.size(WWmap), order='C')
####################################

file_name = "WW_trrm5.txt"

f = open(file_name, 'w')
f.write("map { type multiMap; maps (map_x map_y map_z); map_x { type spaceMap; axis x; grid lin; min -320.0; max 320.0; N 32; } map_y { type spaceMap; axis y; grid lin; min -120.0; max 120.0; N 12; }  map_z { type spaceMap; axis z; grid lin; min -170.0; max 170.0; N 17; } } \n")
f.write("constSurvival 2.0; \n")

f.write("wLower ( ")
for x in range(np.size(WW)):
    f.write(str(WW[x]*0.5) + " ")    
f.write(" ); \n")

f.write("wUpper ( ")
for x in range(np.size(WW)):
    f.write(str(WW[x]*2) + " ")    
f.write(" ); \n")

f.close()

#######################################

# plot results
print (sum(sum(WWmap[:,6,:])))
print (np.shape(WWmap[:,6,:]))
plt.figure()
plt.pcolor(WWmap[:,6,:], cmap='Spectral', norm=LogNorm(vmin=1.0e-15, vmax=1.0))#))
cbar = plt.colorbar()
plt.xlabel('X')
plt.ylabel('Z')
plt.title('P0')
cbar.set_label('CADIS weight windows [-]', rotation=270, labelpad=25)
#plt.savefig("wwCADISP3.pdf", format="pdf")
plt.show()

