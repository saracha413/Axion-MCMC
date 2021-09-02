
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import utilities_OOP
import multiprocessing as mp
from multiprocessing import Pool
from multiprocessing import Process
import time

from datetime import date


from os import listdir
from os.path import isfile, join

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
        s[i] = 1/(L-1) * np.sum(samples[i]-chain_mean[i]) #variance within chains
    B = L/(num_chains-1) * np.sum(chain_mean[i] - grand_mean) #variance between chains
    
    W = 1/num_chains * np.sum(s)**2 #weighted variance within chains
    R = ((L-1)/L * W + 1/L * B)/W #Gelman-Rubin statistic
    
    return R




if __name__ == '__main__':


    run_MCMC = True
    run_GR_test = True
    file_combine = True


    input_pars = utilities_OOP.read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/Axion-MCMC/axion_MCMC/')


    

    #doing this to reduce memory being hogged
    #og_ac = float(input_pars['log10_a_c'])
    #og_frac_ac = float(input_pars['fraction_fld_ac'])
    og_log10_axion_ac = float(input_pars['log10_axion_ac'])
    og_log10_fraction_axion_ac = float(input_pars['log10_fraction_axion_ac'])
    og_omega_cdm = float(input_pars['omega_cdm'])
    og_H0 = float(input_pars['H0'])
    og_djs = 10

    pars = input_pars
    pars['n_axion'] = 2
    
    num_walkers = 10

    start = time.time()
    

    burn_in_steps = 30
    
    gr_test = 50
    counter = 1
    #TEMPORARILY CHANGING THIS
    #chains are converging too fast, want to get more data
    steps_per_it = 2000 #50 #number of steps in MCMC before checking G-R test

    while gr_test > 1.2:
    
        #list of processes and list of chains
        #serve similar purposes, just need to refer to different types of objects
        processlist = []
        chain_list = []
        
        #lists to hold parameter values for G-R test
        log10_axion_ac_list = []
        log10_fraction_ac_list = []
        omega_cdm_list = []
        H0_list = []
        Djs_list = []
        
        
        #initialize and run MCMC chains
        for i in range(num_walkers):
            if run_MCMC:
                name=str(i)
                if counter > 1:
                    og_log10_axion_ac = best_vals[i][0]
                    og_log10_fraction_axion_ac = best_vals[i][1]
                    og_omega_cdm = best_vals[i][2]
                    og_H0 = best_vals[i][3]
                    og_djs = best_vals[i][4]
        
        
                #allow 10% standard dev for now
                pars['log10_axion_ac'] = np.random.normal(og_log10_axion_ac, abs(og_log10_axion_ac*0.1))
                pars['log10_fraction_axion_ac'] = np.random.normal(og_log10_fraction_axion_ac, abs(og_log10_fraction_axion_ac*0.1))
                pars['omega_cdm'] = np.random.normal(og_omega_cdm, abs(og_omega_cdm*0.1))
                pars['H0'] = np.random.normal(og_H0, abs(og_H0*0.1))
        
                chain_list.append(utilities_OOP.Chain(str(i), pars, steps_per_it))
                process = Process(target=chain_list[i].MCMC_run, args=(chain_list[i],))
                processlist.append(process)
                process.start()
            #join MCMC chains to get results
            for process in processlist:
                process.join()

            #best values of the 4 params + djs for each chain
            best_vals = np.zeros((6, num_walkers))
        
        if run_GR_test:
            #get data for Gelman-Rubin test
            #TO-DO: is it faster to read this data from a file or have the Chain objects save their entries so I can access them here?
            for n in range(num_walkers):
                path = '/Users/saravannah/Axion-MCMC/axion_MCMC/'
                fileName ='vary_ac_fEDE_wCDM_H0-chain#'+str(n)+'.txt'#'vary_ac_fEDE_wCDM_H0-chain#'+chain_list[n-1].name+'.txt'
                #unpack file of MCMC chains and assign them to the correct varriables
                log10_axion_ac, log10_fraction_ac, omega_cdm, H0, Djs =  np.loadtxt(path+fileName, unpack=True)
                log10_axion_ac_list.append(log10_axion_ac)
                log10_fraction_ac_list.append(log10_fraction_ac)
                omega_cdm_list.append(omega_cdm)
                H0_list.append(H0)
                Djs_list.append(Djs)
                
                best_vals[n] = np.array([log10_axion_ac[len(log10_axion_ac)-1], log10_fraction_ac[len(log10_fraction_ac)-1], omega_cdm[len(omega_cdm) - 1], H0[len(H0) - 1], Djs[len(Djs)-1]])
            
            num_steps = steps_per_it * counter
            #perform G-R test on all params
            gr_log10_axion_ax = grubin(log10_axion_ac_list, burn_in_steps, num_walkers, num_steps)
            gr_log10_fraction_ac = grubin(log10_fraction_ac_list, burn_in_steps, num_walkers, num_steps)
            gr_omega_cdm = grubin(omega_cdm_list, burn_in_steps, num_walkers, num_steps)
            gr_omega_cdm = grubin(H0_list, burn_in_steps, num_walkers, num_steps)
            print('Gelman-Rubin statistic for log10_axion_ac is ',gr_log10_axion_ax )
            print('Gelman-Rubin statistic for log10_fraction_ac is ', gr_log10_fraction_ac)
            print('Gelman-Rubin statistic for omega_CDM is ', gr_omega_cdm)
            print('Gelman-Rubin statistic for H0 is ',gr_omega_cdm )

            #set gr_test to max param value to check if all params satisy GR test
            gr_test = max(gr_log10_axion_ax, gr_log10_fraction_ac, gr_omega_cdm, gr_omega_cdm)
            
            #increment counter to run another steps_per_it samples
            counter = counter+1
        
    print('All done! Completed in ', time.time()-start, ' seconds.')
    
    
    if file_combine:
        print('Combining files...')
        
        path = '/Users/saravannah/Axion-MCMC/axion_MCMC/'
        files = [f for f in listdir(path) if isfile(join(path, f))]
        
        today = date.today()

        newFileName = path+'vary_ac_fEDE_wCDM-'+str(today)+'.txt'
        
        
        for fileName in files:
            with open(newFileName, 'a') as fileObject:
                #there's got to be a more efficient way to do this...
                if 'vary_ac_fEDE' in fileName and fileName != 'vary_ac_fEDE_wCDM.txt':
                    log10_axion_ac, log10_fraction_axion_ac, omega_cdm, H0, Djs = np.loadtxt(path+fileName, unpack=True)
                    for i in range(len(log10_axion_ac)):
                        line = [log10_axion_ac[i], log10_fraction_axion_ac[i], omega_cdm[i], H0[i], Djs[i]]
                        line = [float(x) for x in line] #Python have trouble with type of data, recasting str to float
                        np.savetxt(fileObject, np.transpose(line),delimiter=',',newline = ' ')
                        fileObject.write('\n')
            fileObject.close()
        

