
import numpy as np
from utilities_NEW import *
import os, os.path




####TO-DO: add truncation scheme

def initiate(params, l_min, l_max):
    
    directory = os.getcwd()
    fileName = os.path.join(directory,'planck/planck_tt_spectrum_2018.txt')
    l_data, Dl_data, Dl_data_err_lo, Dl_data_err_hi = np.loadtxt(fileName, unpack = True)
    
    ####TO-DO: add a_c to params (check axiCLASS)
    #log10_axion_ac = -3.7
    
    l_model, Cl_model, Dl_model = get_power(params, l_min, l_max)

    
    ##TO-DO: remove indexing when you add truncation scheme
    return Dl_model, Dl_data[l_min-1:l_max] #, l_max



#============================================
# main MCMC function
#============================================


def mcmc(params):
    """
          params = {'num_steps':  'num_burn_in':  'name': 'l_min':   'l_max':  'n_axion':
                    log10_axion_ac':   'log10_fraction_axion_ac':  'omega_cdm':  'H0': }
    """

    num_steps = params['num_steps']
    num_burn_in = params['num_burn_in']
    name = params['name']
    l_min = params['l_min']
    l_max = params['l_max']
    model_pars = params['model_pars']

    #model_pars = read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/Axion-MCMC/axion_MCMC/')
    #model_pars['n_axion'] = params['n_axion']
    #model_pars['log10_axion_ac'] = params['log10_axion_ac']
    #model_pars['log10_fraction_axion_ac'] = params['log10_fraction_axion_ac']
    #model_pars['omega_cdm'] = params['omega_cdm']
    #model_pars['H0'] = params['H0']


    Dl_model, Dl_data = initiate(model_pars, l_min, l_max)

    outFile = name+'.txt'

    #starting chain
    print('Dl model is ', Dl_model)
    print('Dl data is ', Dl_data)
    JSD_current = JSD(Dl_model, Dl_data)
    p_current = [model_pars['log10_axion_ac'], model_pars['log10_fraction_axion_ac'], model_pars['omega_cdm'], model_pars['H0']] #whatever params you're varying in MCMC
    stdDevs = [float(model_pars['log10_axion_ac'])*0.05, float(model_pars['log10_fraction_axion_ac'])*0.05, float(model_pars['omega_cdm'])*0.05, float(model_pars['H0'])*0.05] #standard deviation for params
    stdDevs = [abs(x) for x in stdDevs]
    stdDevs = [str(x) for x in stdDevs] #Python have trouble with type of data, recasting float to str

    for t in range(num_steps):

        p_propose = np.random.normal(p_current, stdDevs)
        ##reset params array

        ####TO-DO: write this fxn to use whatever variable param you want. hard-coded for now
        model_pars['log10_axion_ac'] = p_propose[0]
        model_pars['log10_fraction_axion_ac'] = p_propose[1]
        model_pars['omega_cdm'] = p_propose[2]
        model_pars['H0'] = p_propose[3]

        l_propose, Cl_propose, Dl_propose = get_power(params, l_min, l_max)
        JSD_propose = JSD(Dl_propose, Dl_data)
        x = JSD_propose/JSD_current

        if x < 1+np.random.uniform():
            p_current = p_propose
            JSD_current = JSD_propose
        if t > burn_in_steps:
            with open(outFile, 'a') as fileObject:
                line = np.append(p_current, JSD_current)
                line = [float(x) for x in line] #Python have trouble with type of data, recasting str to float
                np.savetxt(fileObject,np.transpose(line),delimiter=',',newline = ' ')
                fileObject.write('\n')
    fileObject.close()





