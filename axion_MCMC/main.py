#adapted from embarassingly paralell MCMC using Pool here:
# https://linuxtut.com/en/02d66723d03d7f977bb8/


from mcmc import *
from utilities import *
import matplotlib
import signal
from multiprocessing import Process, TimeoutError
from classy import Class

import numpy as np
import numpy.random as rd
import scipy.stats as st
import copy, time, os
from datetime import datetime as dt

from multiprocessing import Pool, freeze_support

import matplotlib.pyplot as plt

import functools


if __name__ == '__main__':

    l_min = 90
    l_max = 2000
    num_steps = 1
    num_chains = 1
    num_burn_in = 0
    saveFile = False
    fileName = 'Sept-30_tests.txt'
    runFromFile = False
    inFileName = 'Sept-21_runs.txt'


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

    freeze_support()

    print('Making ', num_steps, ' samples per chain for ',  num_chains, ' chains. Burn-in rate: ', num_burn_in/num_steps)
    pool = Pool(processes=num_chains)
    n_trials_per_process = [num_steps] * num_chains
    print('pool is open')
    start = time.time()
    total_sampling_result = pool.map(mcmc, n_trials_per_process)
    #total_sampling_result = pool.map(functools.partial(mcmc, num_burn_in, l_min, l_max, n_axion, log10_axion_ac_IN, log10_fraction_axion_ac_IN, omega_cdm_IN, H0_IN), n_trials_per_process)
    end = time.time()
    pool.close()
    print('total exec time: ', end-start)



    log10_axion_ac = np.zeros((num_chains, len(total_sampling_result[0])))
    log10_fraction_axion_ac = np.zeros((num_chains, len(total_sampling_result[0])))
    omega_cdm = np.zeros((num_chains, len(total_sampling_result[0])))
    H0 = np.zeros((num_chains, len(total_sampling_result[0])))
    Djs = np.zeros((num_chains, len(total_sampling_result[0])))

    for i in range(num_chains):
        log10_axion_ac[i] = [col[0] for col in total_sampling_result[i]] #extract column
        log10_fraction_axion_ac[i] = [col[1] for col in total_sampling_result[i]]
        omega_cdm[i] = [col[2] for col in total_sampling_result[i]]
        H0[i] = [col[3] for col in total_sampling_result[i]]
        Djs[i] = [col[4] for col in total_sampling_result[i]]


        
        

    if saveFile:
        #create total array to combine all the chain data
        big_arr = np.zeros((5, len(log10_axion_ac[0])*num_chains))

        big_arr[0] = np.concatenate(([log10_axion_ac[i] for i in range(num_chains)]))
        big_arr[1] = np.concatenate(([log10_fraction_axion_ac[i] for i in range(num_chains)]))
        big_arr[2] = np.concatenate(([omega_cdm[i] for i in range(num_chains)]))
        big_arr[3] = np.concatenate(([H0[i] for i in range(num_chains)]))
        big_arr[4] = np.concatenate(([Djs[i] for i in range(num_chains)]))
        #save run to file
        with open(fileName, 'a') as fileObject:
            np.savetxt(fileObject, np.transpose(big_arr))
                
                
    steps = np.arange(0,len(total_sampling_result[0]))
    for i in range(num_chains):
        plt.plot(steps, log10_axion_ac[i])
    plt.show()


