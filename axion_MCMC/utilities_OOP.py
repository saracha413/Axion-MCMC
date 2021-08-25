
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from scipy.stats import norm
import math
import os, os.path
from tqdm import tqdm
import time


#universal variables
l_max = 2000 #truncations for l values, set to match first paper
l_min = 90




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



#calculate modal fraction for a power spectrum
def modal(dls):
    modal = dls/np.sum(dls)
    
    return modal
    
#calculate Jensen-Shannon divergence
def JSD(mod_Dl, dat_Dl):
    p, q = modal(mod_Dl), modal(dat_Dl)
    r = 1/2 * (p+q)
    
    #try:
    #    return 1/2 * np.nansum(p*np.log(p/r)) + 1/2 * np.nansum(q*np.log(q/r))
    #except RuntimeWarning:
    #    print('p is ', p)
    #    print('q is ', q)
    #    pass
    #if any(p/r <= 0):
    #    for i in range(len(p)):
    #        if p[i]/r[i] <=0:
    #            print('p is ', p[i], ' for i = ', i)
    #if any(q/r <= 0):
    #    for i in range(len(q)):
    #        if q[i]/r[i] <=0:
    #            print('q is ', q[i], ' for i = ', i)

    Djs = 1/2 * np.nansum(p*np.log(p/r)) + 1/2 * np.nansum(q*np.log(q/r))

    if Djs == 0:
        print('Djs = 0!')
        print('First term is ', np.nansum(p*np.log(p/r)),' and second term is ',  np.nansum(q*np.log(q/r)))

    return Djs
        
    
    
    
#get power spectrum from CLASS for a given parameter set
#currently set to TT only
#modified from https://github.com/lesgourg/class_public/wiki/Python-wrapper
def get_power(params):
        


    #create an instance of CLASS wrapper w/correct params
    cosmo = Class()
    cosmo.set(params)
    #cosmo.set({'output':'tCl,pCl,lCl,mPk','lensing':'yes','P_k_max_1/Mpc':3.0})
    
    cosmo.compute()
    
    #lensed cl until l=l_max
    output = cosmo.lensed_cl(l_max)
    #CHECK THAT THIS IS INDEXED CORRECTLY --- NEED TO CHECK CLASS DOCUMENTATION ON OUTPUT OF LENSED_CL
    ls = output['ell'][l_min:]
    Cls = output['tt'][l_min:]
    
    
    #ls = np.arange(l_max+1)
    Dls = ls*(ls+1)*Cls/(2*np.pi)
    
    #clean ups
    cosmo.struct_cleanup()
    cosmo.empty()
    
    return ls, Cls, Dls
    
####TO-DO: add truncation scheme

def initiate(params):
    
    directory = os.getcwd()
    fileName = os.path.join(directory,'planck/planck_tt_spectrum_2018.txt')
    l_data, Dl_data, Dl_data_err_lo, Dl_data_err_hi = np.loadtxt(fileName, unpack = True)
    
    ####TO-DO: add a_c to params (check axiCLASS)
    #log10_axion_ac = -3.7
    
    l_model, Cl_model, Dl_model = get_power(params)

    
    ##TO-DO: remove indexing when you add truncation scheme
    return Dl_model, Dl_data[l_min-1:l_max] #, l_max

#calculate Gelman-Rubin statistic to test for chain convergence
#take an array of arrays of samples for each chain
def grubin(samples, burn_in_steps, num_chains):
    L = np.zeros(num_chains)
    chain_mean = np.zeros(num_chains) #mean w/in chain
    
    for i in range(num_chains):
        #remove burn-in steps
        samples[i] = samples[burn_in_steps:len(samples[i])-1]
        #get length of each chain
        L[i] = len(samples[i])
        #get mean within chain
        #chain_mean[i] = 1/L[i] * np.sum(samples[i])
        
        
    chain_mean = 1/L * np.sum(samples)
    grand_mean = 1/num_chains * np.sum(chain_mean) #mean of means
    B = L/(num_chains-1) * np.sum(chain_mean - grand_mean) #variance between chains
    s = 1/(L-1) * np.sum(samples-chain_mean) #variance within chains
    
    #for j in range(num_chains):
    #    s[j] = 1/(L-1) * np.sum(samples[i]-chain_mean[i])
    
    W = 1/num_chains * np.sum(s)**2 #weighted variance within chains
    R = ((L-1)/L * W + 1/L * B)/W #Gelman-Rubin statistic
    
    return R

class Chain:

    pars_lcdm = {'omega_b':0.022032,
    'omega_cdm':0.12038,
    'h':0.67556,
    'A_s':2.215e-9,
    'n_s':0.9619,
    'tau_reio':0.0925}
    
    def __init__(self, name, params, num_steps):
        self.name = name
        self.params = params
        self.num_steps = num_steps
        #print('initiated chain with name ', self.name)
        
        
    #def say_hi(self, phrase):
    #    print(phrase)
    
    def MCMC_run(self, numsteps=20, burn_in_steps=0):
    
        params = self.params
        num_steps = self.num_steps
    
        #Dl_model, Dl_data, l_max = initiate(params)
        Dl_model, Dl_data = initiate(params)
        
       #print('Starting chain')
        
        #starting chain
        JSD_current = JSD(Dl_model, Dl_data)
        p_current = [params['log10_axion_ac'], params['log10_fraction_axion_ac'], params['omega_cdm'], params['H0']] #whatever params you're varying in MCMC
        stdDevs = [float(params['log10_axion_ac'])*0.05, float(params['log10_fraction_axion_ac'])*0.05, float(params['omega_cdm'])*0.05, float(params['H0'])*0.05] #standard deviation for params
        stdDevs = [abs(x) for x in stdDevs]
        stdDevs = [str(x) for x in stdDevs] #Python have trouble with type of data, recasting float to str
        
        
        outFile = 'vary_ac_fEDE_wCDM_H0-chain#'+self.name+'.txt'
        
        with open(outFile, 'a') as fileObject:
            line = np.append(p_current, JSD_current)
            line = [float(x) for x in line] #Python have trouble with type of data, recasting str to float
            np.savetxt(fileObject,np.transpose(line),delimiter=',',newline = ' ')
            fileObject.write('\n')
            
        #count steps accepted so can calculate acceptance fraction
        steps_accepted = 0
        
        for t in range(num_steps):

            print('t is , ', t)
            
            #suggest a random value for params from a normal distrib centered on current values
            p_propose = np.random.normal(p_current, stdDevs)
            
            ##reset params array
            #fullParams = params
            ####TO-DO: write this fxn to use whatever variable param you want. hard-coded for now
            params['log10_axion_ac'] = p_propose[0]
            params['log10_fraction_axion_ac'] = p_propose[1]
            params['omega_cdm'] = p_propose[2]
            params['H0'] = p_propose[3]
            
            try:
                l_propose, Cl_propose, Dl_propose = get_power(params)
                    
                JSD_propose = JSD(Dl_propose, Dl_data)
                x = JSD_propose/JSD_current
            except:
                #CLASS didn't converge, keep MCMC going but ensure this step isn't accepted
                x = 0
                pass
                
                
            #if t%50 == 0:
            #    print('Chain ',self.name,' has reached ', str(t/50), ' steps.')
            #Metropolis-Hastings acceptance criterion
            #from https://github.com/AstroHackWeek/AstroHackWeek2015/blob/3e13d786ecb86fd4757c08ab63cfc08135933556/hacks/sklearn-CV-Bayes.py
                
            if x < 1+np.random.uniform():
                p_current = p_propose
                JSD_current = JSD_propose
                steps_accepted = steps_accepted + 1
                
            if t > burn_in_steps:
                with open(outFile, 'a') as fileObject:
                    line = np.append(p_current, JSD_current)
                    line = [float(x) for x in line] #Python have trouble with type of data, recasting str to float
                    np.savetxt(fileObject,np.transpose(line),delimiter=',',newline = ' ')
                    fileObject.write('\n')
        fileObject.close()
        #print('Acceptance fraction is ', steps_accepted/num_steps)
        
