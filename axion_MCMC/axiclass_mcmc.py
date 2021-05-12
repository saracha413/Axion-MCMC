###Sara Vannah


# import necessary modules
#%matplotlib inline
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from scipy.stats import norm
import math

import os, os.path



#------------------------------------------------------------------------------------------------------
#                                    Helper functions
#------------------------------------------------------------------------------------------------------

##This is from AxiCLASS
def is_number(s):
# ---------------------------------- This func checks whether a thing is a number. Found online
    try:
        float(s)
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

    return False
    
#function for ignoring comments in ini files
def ignore_comment(line):
    if '#' in line:
        #save all elements up to the #
        line = line[:line.find('#')]
    if '*' in line:
        line = ''

    return line
    
    
##This is from AxiCLASS

def read_ini_file(inifile, loc = ''):

# Function to read ini file and save it in a dictionary that can be passed to classy
# Takes the required argument inifile = filename with extension
# Takes the optional argument loc = location of your ini file, ending in a '/'
# Returns dictionary of everything contained in your ini file
#
    inivals = {}

    with open(loc + inifile) as f: # opening initialisation file as f
        content = f.readlines() # reading the initialisation file and turning it into a list

    q = {} # initialise q as an empty dictionary
    for s in content: # iterates over lines in .ini file and saves information, turning numbers into floats from strings
        #SV --- added this skip over commented sections
        s = ignore_comment(s)
        if s != '':
            if is_number(s[s.find('=')+2:]):
                q[s[:s.find(' =')]] = float(s[s.find('=')+2:])
            else:
                q[s[:s.find(' =')]] = s[s.find('=')+2:-1]

    q.pop('')
    return q # inivals dict has dict of initial values at key given by 'original'

#------------------------------------------------------------------------------------------------------
#                               Information theory
#------------------------------------------------------------------------------------------------------

#calculate modal fraction for a power spectrum
def modal(dls):
    modal = dls/np.sum(dls)
    
    return modal
    
#calculate Jensen-Shannon divergence
def JSD(mod_Dl, dat_Dl):
    p, q = modal(mod_Dl), modal(dat_Dl)
    r = 1/2 * (p+q)
    
    return 1/2 * np.nansum(p*np.log(p/r)) + 1/2 * np.nansum(q*np.log(q/r))


#------------------------------------------------------------------------------------------------------
#                               Cosmology
#------------------------------------------------------------------------------------------------------

###TO-DO: remove hardcoding here
#input_pars = read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/axiclass/')


#get power spectrum from CLASS for a given parameter set
#currently set to TT only
#modified from https://github.com/lesgourg/class_public/wiki/Python-wrapper
def get_power(params):
    
    l_max = 2000
    
    #create an instance of CLASS wrapper w/correct params
    cosmo = Class()
    cosmo.set(params)
    #cosmo.set({'output':'tCl,pCl,lCl,mPk','lensing':'yes','P_k_max_1/Mpc':3.0})
    
    cosmo.compute()
    
    #lensed cl until l=l_max
    output = cosmo.lensed_cl(l_max)
    #CHECK THAT THIS IS INDEXED CORRECTLY --- NEED TO CHECK CLASS DOCUMENTATION ON OUTPUT OF LENSED_CL
    ls = output['ell'][2:]
    Cls = output['tt'][2:]
    
    
    #ls = np.arange(l_max+1)
    Dls = ls*(ls+1)*Cls/(2*np.pi)
    
    #clean ups
    cosmo.struct_cleanup()
    cosmo.empty()
    
    return ls, Cls, Dls
    
#------------------------------------------------------------------------------------------------------
#                               MCMC
#------------------------------------------------------------------------------------------------------
    
####TO-DO: add truncation scheme

def initiate(params):
    
    directory = os.getcwd()
    fileName = os.path.join(directory,'planck/planck_tt_spectrum_2018.txt')
    l_data, Dl_data, Dl_data_err_lo, Dl_data_err_hi = np.loadtxt(fileName, unpack = True)
    
    
    l_model, Cl_model, Dl_model = get_power(params)
    
    
    ##TO-DO: remove indexing when you add truncation scheme
    return Dl_model, Dl_data[:len(Dl_model)] #, l_max


### to-do: make l_max set by truncation
def MCMC_run(params, numsteps=10, outFile=None):

    if outFile == None:
        #outFile = 'PHENO_vary_ac_and_frac_fld_naxion='+str(params['n_axion'])+'_log10_ac='+str(params['log10_a_c'])+'.txt'
        outFile = 'PHENO_vary_ac_and_frac_fld_naxion='+'ugh'+'_log10_ac='+str(params['log10_a_c'])+'.txt'

    burn_in_steps = 5

    #Dl_model, Dl_data, l_max = initiate(params)
    Dl_model, Dl_data = initiate(params)

    #starting chain
    JSD_current = JSD(Dl_model, Dl_data)
    p_current = [params['log10_a_c'], params['fraction_fld_ac']] #-3.5 #whatever params you're varying in MCMC
    stdDevs = [0.5, 0.005] #standard deviation for params

    #check if file to write to exists; create it if not
    #if not os.path.isfile(outFile):

    for t in range(numsteps):

        write_params_to_file = False

        #suggest a random value for params from a normal distrib centered on current values
        p_propose = [norm(p_current[0], stdDevs[0]).rvs(), norm(p_current[1], stdDevs[1]).rvs()]
        ##reset params array
        #fullParams = params
        ####TO-DO: write this fxn to use whatever variable param you want. hard-coded for now
        [params['log10_a_c'], params['fraction_fld_ac']] = p_propose

        l_propose, Cl_propose, Dl_propose = get_power(params)

        JSD_propose = JSD(Dl_propose, Dl_data)
        x = JSD_propose/JSD_current

        #Metropolis-Hastings acceptance criterion
        #from https://github.com/AstroHackWeek/AstroHackWeek2015/blob/3e13d786ecb86fd4757c08ab63cfc08135933556/hacks/sklearn-CV-Bayes.py
        if x > np.random.uniform():
            p_current = p_propose
            JSD_current = JSD_propose
            write_params_to_file = True
            
            if t > burn_in_steps:
                with open(outFile, 'a') as fileObject:
                    print('writing to file')
                    line = np.append(p_current, JSD_current)
                    np.savetxt(fileObject,np.transpose(line),delimiter=',',newline = ' ')
                    fileObject.write('\n')
                


    fileObject.close()
    


