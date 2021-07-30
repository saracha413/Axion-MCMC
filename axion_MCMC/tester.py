
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import utilities
import multiprocessing as mp
from multiprocessing import Pool
from multiprocessing import Process
import time


params = utilities.read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/Axion-MCMC/axion_MCMC/')

#ls, Cls, Dls = utilities.get_power(input_pars)

#starting chain
JSD_current = 0.1
p_current = [params['log10_axion_ac'], params['log10_fraction_axion_ac'], params['omega_cdm'], params['H0']] #whatever params you're varying in MCMC
stdDevs = [float(params['log10_axion_ac'])*0.05, float(params['log10_fraction_axion_ac'])*0.05, float(params['omega_cdm'])*0.05, float(params['H0'])*0.05] #standard deviation for params
stdDevs = [str(x) for x in stdDevs]

outFile = 'vary_ac_fEDE_wCDM_H0'+str(time.time())+'.txt'
   
with open(outFile, 'a') as fileObject:
    line = np.append(p_current, JSD_current)
    line = [float(x) for x in line]
    np.savetxt(fileObject,np.transpose(line),delimiter=',',newline = ' ')
    fileObject.write('\n')
