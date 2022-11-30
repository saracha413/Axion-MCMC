from utilities import *
#from mcmc import *

import numpy as np
import multiprocessing as mp
from classy import Class
from classy import CosmoComputationError, CosmoSevereError

import os, os.path, copy

from interruptingcow import timeout
import tqdm


#NOTE IF THIS WORKS YOU CAN REMOVE PARS_LIST FROM ARGUMENTS!!!!
def change_starting_place(pars_list, bounds):

    new_pars_list = np.zeros(len(pars_list))
    for i in range(len(pars_list)):
        new_pars_list[i] = np.random.uniform(bounds[i][0], bounds[i][1])

    return new_pars_list
    
    
#calculate Gelman-Rubin statistic to test for chain convergence
#take an array of arrays of samples for each chain
def grubin(samples, burn_in_steps, num_chains, num_steps):
    L = num_steps - burn_in_steps #np.zeros(num_chains)
    chain_mean = np.zeros(num_chains) #mean w/in chain
    

    for i in range(num_chains):
        if type(samples[i]) == np.float64:
            print('Chain number ', str(i), ' has just one sample. More steps need to be run for convergence.')
            chain_mean[i] = samples[i]
            
        elif len(samples[i]) == 0:
            print('Chain number ', str(i), ' is empty. More steps need to be run for convergence.')
            chain_mean[i] = float("NaN")
        else:
            #remove burn-in steps
            #samples[i] = samples[burn_in_steps:len(samples[i])-1] #BURN-IN STEPS ALREADY REMOVED IN MCMC
            #get mean within chain
            #chain_mean[i] = 1/L[i] * np.sum(samples[i])
            chain_mean[i] = 1/L * np.sum(samples[i])
            
            
    grand_mean = 1/num_chains * np.sum(chain_mean) #mean of means
    
    s = np.zeros(num_chains)
    for i in range(num_chains):
        s[i] = 1/(L-1) * np.sum((samples[i]-chain_mean[i])**2) #variance within chains
    B = L/(num_chains-1) * np.sum((chain_mean[i] - grand_mean)**2) #variance between chains
    W = 1/num_chains * np.sum(s)**2 #weighted variance within chains
    R = ((L-1)/L * W + 1/L * B)/W #Gelman-Rubin statistic
    #print(L, W, B, R)
    
    return R

def func(num_it):

    num_burn_in = 0
    l_min = 90
    l_max = 2000
    n_axion = 3


    #set min and max for sampled params
    #note this kind of setting a prior :
    param_ranges =  {'log10_axion_ac':[-4.6, -3], 'log10_fraction_axion_ac': [-2, -0.82], 'omega_cdm': [0.095, 0.15], 'H0': [65, 80]}
    


    #specify input param file
    model_pars = read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/Axion-MCMC/axion_MCMC/')
    model_pars['n_axion'] = n_axion

    sampling_result = []



    p_current = [model_pars['log10_axion_ac'], model_pars['log10_fraction_axion_ac'], model_pars['omega_cdm'], model_pars['H0']] #whatever params you're varying in MCMC
    stdDevs = [str(abs(float(x)*0.1)) for x in p_current]

    fEDE_stdDev = 10**float(p_current[1])*0.1
    #print('std is ', fEDE_stdDev)

    #assign new starting places for each chain
    #note since p_current is an array not a dict, have to convert valeus in param_ranges to a list as well
    p_current = change_starting_place(p_current, list(param_ranges.values()))
    


    model_pars['log10_axion_ac'] = p_current[0]
    model_pars['log10_fraction_axion_ac'] = p_current[1]
    model_pars['omega_cdm'] = p_current[2]
    model_pars['H0'] = p_current[3]


    directory = os.getcwd()
    fileName = os.path.join(directory,'planck/planck_tt_spectrum_2018.txt')
    _, Dl_data, _, _ = np.loadtxt(fileName, unpack = True)
    Dl_data = Dl_data[l_min-2:l_max-1] #truncation checked to match l_max and l_min


    _, _, Dl_init = get_power(model_pars, l_min, l_max)
    JSD_current = JSD(Dl_init, Dl_data)
    
    #start = time.time()

    i = 0
    num_accept = 0 #number of steps accepted

    #for i in range(num_it):
    while i < num_it:

        #need to take step that's flat in fEDE not log10(fEDE), so treat separately from rest of params
        fEDE_current = 10**float(p_current[1])
        fEDE_propose = np.random.normal(fEDE_current, fEDE_stdDev)


        #print('fEDE_propose is ', fEDE_propose)


        p_propose = np.random.normal(p_current, stdDevs)
        p_propose[1] = np.log10(fEDE_propose) #overwrite with fEDE chosen from gaussian distrib
        #print('log10fEDE is ', p_propose[1])


        ####TO-DO: write this fxn to use whatever variable param you want. hard-coded for now
        model_pars['log10_axion_ac'] = p_propose[0]
        model_pars['log10_fraction_axion_ac'] = p_propose[1]
        model_pars['omega_cdm'] = p_propose[2]
        model_pars['H0'] = p_propose[3]



        #NOTE: you can organize this more pythonically by changing param_range to be a dict and using the key from model_pars
        allParamsWithinRange = True
        for param in param_ranges:
            if model_pars[param] < param_ranges[param][0] or model_pars[param] > param_ranges[param][1]  or np.isnan(model_pars[param]):
                allParamsWithinRange = False


        if allParamsWithinRange:

            try:
                with timeout(60, exception=RuntimeError):
                    _, _, Dl_propose = get_power(model_pars, l_min, l_max)

                    JSD_propose = JSD(Dl_propose, Dl_data)
                    x = JSD_propose/JSD_current
                    
                    #Metropolis-Hastings acceptance
                    
                    if x < 0.5+np.random.uniform():
                        p_current = p_propose
                        JSD_current = JSD_propose
                        num_accept = num_accept + 1


            except RuntimeError:
                pass

            except CosmoComputationError:
                print('Computation error at step number ', i)
                print('Param values for log10_axion_ac, log10_fraction_axion_ac, omega_cdm, and H0 are', p_propose[0], p_propose[1], p_propose[2], p_propose[3])

            except CosmoSevereError:
                print('Params causing Classy to crash. log10_axion_ac is ', model_pars['log10_axion_ac'], ' while log10_fraction_axion_ac is ', model_pars['log10_fraction_axion_ac'])
                print('fEDE is ', fEDE_propose)



            line = np.append(p_current, JSD_current)
            sampling_result.append(line)

            i=i+1

    frac_accepted = num_accept/num_it

    print('Acceptance fraction is ', frac_accepted)
    return sampling_result


if __name__ == '__main__':
    num_chains = 12
    num_steps = 1000
    saveFile = True
    #if saveFile:
    fileName = 'Nov-16.txt'



    print('starting code')


    pool = mp.Pool(processes=24)
    print('The pool is open! Come swim!')

    n_trials_per_process = [num_steps] * num_chains


    total_sampling_result = pool.map(func, n_trials_per_process)
    #total_sampling_result = []
    #for _ in tqdm.tqdm(pool.imap_unordered(func, n_trials_per_process), total=num_steps*num_chains):
    #    total_sampling_result.append(_)



    pool.close()
    print('The pool is now closed.')


    log10_axion_ac = np.zeros((num_chains, len(total_sampling_result[0])))
    log10_fraction_axion_ac = np.zeros((num_chains, len(total_sampling_result[0])))
    omega_cdm = np.zeros((num_chains, len(total_sampling_result[0])))
    H0 = np.zeros((num_chains, len(total_sampling_result[0])))
    Djs = np.zeros((num_chains, len(total_sampling_result[0])))
    
    num_params = 4+1 #number of parameters you're varying plus Djs
    #collection of all the chains (for G-R calculation)
    super_chain_array = np.zeros((num_chains, num_params, num_steps))
    
    for i in range(num_chains):

        log10_axion_ac[i] = [col[0] for col in total_sampling_result[i]] #extract column
        log10_fraction_axion_ac[i] = [col[1] for col in total_sampling_result[i]]
        omega_cdm[i] = [col[2] for col in total_sampling_result[i]]
        H0[i] = [col[3] for col in total_sampling_result[i]]
        Djs[i] = [col[4] for col in total_sampling_result[i]]
        chain_arr = np.array([log10_axion_ac[i], log10_fraction_axion_ac[i], omega_cdm[i], H0[i], Djs[i]])
        super_chain_array[i] = chain_arr
        chainFileName = fileName.split('.tx')[0] + '#'+str(i)+'.txt'
        np.savetxt(chainFileName, np.transpose(chain_arr))
        
        
    
    burn_in_steps = 0

    gr_log10_axion_ac = grubin(super_chain_array[:,0,:], burn_in_steps, num_chains, num_steps)
    gr_log10_fraction_ac = grubin(super_chain_array[:,1,:], burn_in_steps, num_chains, num_steps)
    gr_omega_cdm = grubin(super_chain_array[:,2,:], burn_in_steps, num_chains, num_steps)
    gr_H0 = grubin(super_chain_array[:,3,:], burn_in_steps, num_chains, num_steps)
    
    print('Gelman-Rubin statistic for log10_axion_ac is ',gr_log10_axion_ac)
    print('Gelman-Rubin statistic for log10_fraction_ac is ', gr_log10_fraction_ac)
    print('Gelman-Rubin statistic for omega_cdm is ', gr_omega_cdm)
    print('Gelman-Rubin statistic for H0 is ', gr_H0)
