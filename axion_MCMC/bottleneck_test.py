import numpy as np
from utilities import *
import os, os.path, copy
import signal
from classy import Class
from classy import CosmoComputationError

import mcmc as MC

from tqdm import *


from utilities import *
import matplotlib

import multiprocessing as mp
from multiprocessing import Process, TimeoutError

import numpy.random as rd
import scipy.stats as st
import copy, time, os
from datetime import datetime as dt


import matplotlib.pyplot as plt

import functools



# check for line_profiler or memory_profiler in the local scope, both
# are injected by their respective tools or they're absent
# if these tools aren't being used (in which case we need to substite
# a dummy @profile decorator)
##from High Performance Python
if 'line_profiler' not in dir() and 'profile' not in dir():
    def profile(func):
        return func

#instructions for user-defined exceptions from here https://www.programiz.com/python-programming/user-defined-exception
class Error(Exception):
    """Base class for other exceptions"""
    pass

class ParamValueError(Error):
    """Raised when strange parameter values prevent AxiCLASS from convering in reasonable time"""
    pass


#borrowed from https://stackoverflow.com/questions/492519/timeout-on-a-function-call
def handler(signum, frame):
    #raise ParamValueError
    pass




#============================================
# main MCMC function
#============================================

#@profile
#def mcmc(params, num_steps):
def mcmc(num_burn_in, l_min, l_max, n_axion, log10_axion_ac, log10_fraction_axion_ac, omega_cdm, H0, num_steps):
    """
          params = {'num_burn_in': l_min':   'l_max':  'n_axion':
                    log10_axion_ac':   'log10_fraction_axion_ac':  'omega_cdm':  'H0': }
    """

    model_pars = read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/Axion-MCMC/axion_MCMC/')
    model_pars['n_axion'] = n_axion

    log10_axion_ac = float(model_pars['log10_axion_ac'])
    log10_fraction_axion_ac = float(model_pars['log10_fraction_axion_ac'])
    omega_cdm = float(model_pars['omega_cdm'])
    H0 = float(model_pars['H0'])


    log10_axion_ac = np.random.normal(log10_axion_ac, abs(log10_axion_ac*0.1))
    log10_fraction_axion_ac = np.random.normal(log10_fraction_axion_ac, abs(log10_fraction_axion_ac*0.1))
    omega_cdm = np.random.normal(omega_cdm, abs(omega_cdm*0.05))
    H0 = np.random.normal(H0, abs(H0*0.1))

    start = time.time()


    directory = os.getcwd()
    fileName = os.path.join(directory,'planck/planck_tt_spectrum_2018.txt')
    _, Dl_data, _, _ = np.loadtxt(fileName, unpack = True)
    Dl_data = Dl_data[l_min-2:l_max-1] #truncation checked to match l_max and l_min


    #initiatiate
    l_model, Cl_model, Dl_model = get_power(model_pars, l_min, l_max)

    #outFile = name+'.txt'

    JSD_current = JSD(Dl_model, Dl_data)
    p_current = [model_pars['log10_axion_ac'], model_pars['log10_fraction_axion_ac'], model_pars['omega_cdm'], model_pars['H0']] #whatever params you're varying in MCMC
    stdDevs = [float(model_pars['log10_axion_ac'])*0.05, float(model_pars['log10_fraction_axion_ac'])*0.05, float(model_pars['omega_cdm'])*0.05, float(model_pars['H0'])*0.05] #standard deviation for params
    stdDevs = [abs(x) for x in stdDevs]
    stdDevs = [str(x) for x in stdDevs] #Python have trouble with type of data, recasting float to str

    sampling_result = []
    line = np.append(p_current, JSD_current)
    line = [float(x) for x in line] #Python have trouble with type of data, recasting str to float
    sampling_result.append(line)

    text = 'Chain # '+str(os.getpid())

    pos = os.getpid()

    for t in range(num_steps):
        signal.signal(signal.SIGALRM, handler)
        signal.alarm(30)

        try:

            p_propose = np.random.normal(p_current, stdDevs)
            ##reset params array

            ####TO-DO: write this fxn to use whatever variable param you want. hard-coded for now
            model_pars['log10_axion_ac'] = p_propose[0]
            model_pars['log10_fraction_axion_ac'] = p_propose[1]
            model_pars['omega_cdm'] = p_propose[2]
            model_pars['H0'] = p_propose[3]

            l_propose, Cl_propose, Dl_propose = get_power(model_pars, l_min, l_max)
            JSD_propose = JSD(Dl_propose, Dl_data)
            x = JSD_propose/JSD_current

            if x < 1+np.random.uniform():
                p_current = p_propose
                JSD_current = JSD_propose

        except (ParamValueError, CosmoComputationError):
            pass

        #moving outside try block so **always** writes to file even if times out
        if t > num_burn_in:
            line = np.append(p_current, JSD_current)
            line = [float(x) for x in line] #Python have trouble with type of data, recasting str to float
            sampling_result.append(line)

        if t%10 == 0:
            print("Process ", pos, ' has completed ', t ,' steps.')


    end = time.time()

    return sampling_result





 #========================================================
 # Main.py contents here
 #========================================================
if __name__ == '__main__':

    l_min = 90
    l_max = 2000
    num_steps = 100
    num_chains = 3
    num_burn_in = 0
    saveFile = False
    runFromFile = False


    n_axion = 3

    model_pars = read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/Axion-MCMC/axion_MCMC/')
    model_pars['n_axion'] = n_axion

    params = {'num_burn_in': num_burn_in,
              'l_min':  l_min, 'l_max': l_max, 'model_pars':model_pars }
    #get parameters to start chains from

    if runFromFile:
        #get full params from file
        par0, par1, par2, par3, par4 = np.loadtxt(inFileName, unpack = True)
        #find location of best Djs and start there
        idx = np.argmin(par4)

        log10_axion_ac_IN = par0[idx]
        log10_fraction_axion_ac_IN = par1[idx]
        omega_cdm_IN = par2[idx]
        H0_IN = par3[idx]
    else:
        log10_axion_ac_IN = float(params['model_pars']['log10_axion_ac'])
        log10_fraction_axion_ac_IN = float(params['model_pars']['log10_fraction_axion_ac'])
        omega_cdm_IN = float(params['model_pars']['omega_cdm'])
        H0_IN = float(params['model_pars']['H0'])
        
        
    ######THIS CELL RUNS THE MCMC##########

    #freeze_support()

    #print('Making ', num_steps, ' samples per chain for ',  num_chains, ' chains. Burn-in rate: ', num_burn_in/num_steps)
    pool = mp.Pool(processes=num_chains)
    print('Starting')
    n_trials_per_process = [num_steps] * num_chains
    start = time.time()
    #with mp.Pool(processes = num_chains, initializer=tqdm.set_lock, initargs=(tqdm.get_lock(),)) as pool:
    #pool.map(functools.partial(mcmc, num_burn_in, l_min, l_max, n_axion, log10_axion_ac_IN, log10_fraction_axion_ac_IN, omega_cdm_IN, H0_IN), n_trials_per_process)
    #total_sampling_result = mcmc(num_burn_in, l_min, l_max, n_axion, log10_axion_ac_IN, log10_fraction_axion_ac_IN, omega_cdm_IN, H0_IN, 10)
    total_sampling_result = pool.map(MC.mcmc, n_trials_per_process)
    end = time.time()
    pool.close()
    print('total exec time: ', end-start)







