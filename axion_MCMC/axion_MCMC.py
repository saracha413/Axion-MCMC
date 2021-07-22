
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import utilities
import multiprocessing as mp
from multiprocessing import Pool
from multiprocessing import Process
import time

#Bubble sorting algortihm from https://realpython.com/sorting-algorithms-python/#the-bubble-sort-algorithm-in-python

def reorganize(log_ac_vals, JSD_vals):
    new_log_ac_vals, new_JSD_vals = log_ac_vals, JSD_vals
    n = len(log_ac_vals)
    
    for i in range(n):
        already_sorted = True
        
        for j in range(n-i-1):
            if new_log_ac_vals[j] > new_log_ac_vals[j+1]:
                new_log_ac_vals[j], new_log_ac_vals[j+1] = new_log_ac_vals[j+1], new_log_ac_vals[j]
                new_JSD_vals[j], new_JSD_vals[j+1] = new_JSD_vals[j+1], new_JSD_vals[j]
                
                already_sorted = False
                
        if already_sorted:
            break
            
    return new_log_ac_vals, new_JSD_vals
    
    
    
    
    
    
    




if __name__ == '__main__':

    input_pars = utilities.read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/Axion-MCMC/axion_MCMC/')




    #doing this to reduce memory being hogged
    #og_ac = float(input_pars['log10_a_c'])
    #og_frac_ac = float(input_pars['fraction_fld_ac'])
    og_log10_axion_ac = float(input_pars['log10_axion_ac'])
    og_log10_fraction_axion_ac = float(input_pars['log10_fraction_axion_ac'])
    og_omega_cdm = float(input_pars['omega_cdm'])
    og_H0 = float(input_pars['H0'])

    pars = input_pars
    pars['n_axion'] = 2
    
    num_walkers = 10

    start = time.time()
    
    processlist = []

    for i in range(num_walkers):

        process = Process(target=utilities.MCMC_run, args=(pars,))
        processlist.append(process)
        process.start()

    for process in processlist:
        process.join()

    print('All done! Completed in ', time.time()-start, ' seconds.')
    
    
    
    
#    pool = mp.Pool(mp.cpu_count()+2)
#
#    num_walkers = 10
#
#
#    with Pool(processes = 10) as pool:
#        for _ in range(num_walkers):
#
#            #for each _, create a "walker"...
#            #... at some point in parameter space chosen by the random fxn
#
#
#            #allow 10% standard dev for now
#            pars['log10_axion_ac'] = np.random.normal(og_log10_axion_ac, abs(og_log10_axion_ac*0.1))
#            pars['log10_fraction_axion_ac'] = np.random.normal(og_log10_fraction_axion_ac, abs(og_log10_fraction_axion_ac*0.1))
#            pars['omega_cdm'] = np.random.normal(og_omega_cdm, abs(og_omega_cdm*0.1))
#            pars['H0'] = np.random.normal(og_H0, abs(og_H0*0.1))
#
#
#            pool.apply_async(utilities.MCMC_run, (pars,))
#
#            print('working?')
#
#
#    ####WHY IS THIS HERE????
#    #pool.apply_async(utilities.MCMC_run, args=(input_pars,))
#        print('exited')
#        pool.close()
#        pool.join()

#    print('All done! Completed in ', time.time()-start, ' seconds.')


#log_ac_n2, JSD_n2 = np.loadtxt('vary_ac_naxion=2.txt', unpack=True)
#log_ac_n3, JSD_n3 = np.loadtxt('vary_ac_naxion=3.txt', unpack=True)
#log_ac_n10, JSD_n10 = np.loadtxt('vary_ac_naxion=10.txt', unpack=True)




#MCMC_run(input_pars, 200, 'log_ac_n=2')

#log_ac_vals, JSD_vals = np.loadtxt('test_log_ac.txt', unpack=True)


#plt.scatter(log_ac_vals, JSD_vals)
#plt.xlabel('log10_ac')
#plt.ylabel('JSD')
