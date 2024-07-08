# -*- coding: utf-8 -*-
# line above has to be the first line of the code when executing in Linux shell
import os
import json
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
# -------- Creating main input file
# ----------------------------------------------------------------------------------------------------------------------

def makeWW(n):

    file_output_TRRM = "./caskRR5_" + str(n) + ".json"    
    
    # returns JSON object as python dictionary
    with open(file_output_TRRM) as f:
        dictionary = json.load(f)
    
    # get values
    tally   = dictionary.get('fluxMap')
    res     = np.array(tally.get('fluxMap'))
    flux    = np.sum(res[:,:,:,:,0],3)
    resp    = dictionary.get('Response')

    # construct weight windows from output
    WWmap  = resp/flux
    WW     = np.reshape(WWmap, np.size(WWmap), order='C')
    
    # write WW file
    file_name = "ww_P1" + str(n) + ".txt" #change for scattering order.
    
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


# ----------------------------------------------------------------------------------------------------------------------
# -------- Creating main input file
# ----------------------------------------------------------------------------------------------------------------------
def create_main_input_file (sim_nr, input_name):       # create new main input file
    file_input_template = "./" + input_name  				      # load the template 
    file_input = input_name + str(sim_nr)

    f = open(file_input_template, 'r')
    template = f.readlines()
    f.close()
    
    for x in range(len(template)):
        if "XXX" in template[x]:
            template[x] = template[x].replace("XXX", str(sim_nr))

    f = open(file_input, 'w')
    for x in range(len(template)):
        f.write(str(template[x]))
    f.close()

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# run scone
# ----------------------------------------------------------------------------------------------------------------------

def run_sbatch(input_name):

    os.system("sbatch /home/rn438/SCONE/InputFiles/TRRM/hpc_script_TRRM /home/rn438/SCONE/InputFiles/TRRM/" + input_name)

#-------------------------------------------------------------------------------
#-------------------------------      MAIN      --------------------------------
#-------------------------------------------------------------------------------

path = "/home/rn438/SCONE/InputFiles/TRRM"   # MAIN FOLDER PATH
os.chdir(path)

# input file name
input_name = "cask_P1"

for sim_nr in range(1,21):  

    #makeWW(sim_nr)
    create_main_input_file(sim_nr, input_name)
    run_sbatch(input_name + str(sim_nr))



