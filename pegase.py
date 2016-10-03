###################################
#
# Python wrapper for PÉGASE (http://www2.iap.fr/pegase/pegasehr/index.html)
#
###################################
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as p
import subprocess as s
import glob as g
import os

prefix_run="test"
name_SSPs_input="SSPs_input"
s.call(["rm",name_SSPs_input])
SSPs_input=open(name_SSPs_input,"w")
SSPs_input.write("2\n0.08\n120\nA\ny\n3\n"+prefix_run) # IMF - lower mass - higher mass - SN model - Stellar wind - Spectral library
SSPs_input.close()
s.call("./SSPs_HR <"+name_SSPs_input,shell=True)
#os.system("./SSPs_HR < "+name_SSPs_input)
tracks=g.glob("*Z*.dat")
Zs=sorted([t[12:-4] for t in tracks])
scenario_name="scenario"
s.call(["rm",scenario_name])
name_scenario_input="scenario_input"
s.call(["rm",name_scenario_input])
scenario_input=open(name_scenario_input,"w")
scenario_input.write(scenario_name+"\n")
scenario_input.write(prefix_run+"_SSPs.dat\n")
scenario_input.write("0.05\n3\n")
props="n0y0nn0" # no infall - single burst SFH - evolution of metallicity - no substellar objects - no winds/nebular emission - no extinction
spectra_list=[]
for i,z in enumerate(Zs):
    s.call(["rm","spectra_"+str(i)+".fits"])
    scenario_input.write("spectra_"+str(i)+".fits\n")
    spectra_list.append("spectra_"+str(i)+".fits")
    scenario_input.write(str(z)+"\n")
    dump=[scenario_input.write(pr+"\n") for pr in props]
scenario_input.write("end")
scenario_input.close()
print scenario_name
s.call("./scenarios_HR <"+name_scenario_input,shell=True)
name_spectra_input="spectra_input"
s.call(["rm",name_spectra_input])
spectra_input=open(name_spectra_input,"w")
spectra_input.write(scenario_name)
spectra_input.close()
s.call("./spectra_HR <"+name_spectra_input,shell=True)
# REDSHIFT SPECTRA
name_color_input="color_input"
for i,sp in enumerate(spectra_list):
    s.call(["rm",name_color_input])
    s.call(["rm","color_"+str(i)+".dat"])
    color_input=open(name_color_input,"w")
    color_input.write(sp+"\n")
    color_input.write("color_"+str(i)+".dat")
    color_input.close()
    s.call("./colors_HR <"+name_color_input,shell=True)

