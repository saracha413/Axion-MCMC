
import numpy as np
from classy import Class

#calculate modal fraction for a power spectrum
def modal(dls):
    modal = dls/np.nansum(dls)
    
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

    #print('p is ', p)
    #print('q is ', q)
    #print('r is ', r)
    Djs = 1/2 * np.nansum(p*np.log(p/r)) + 1/2 * np.nansum(q*np.log(q/r))

    #if Djs == 0:
    #    print('Djs = 0!')
    #    print('First term is ', np.nansum(p*np.log(p/r)),' and second term is ',  np.nansum(q*np.log(q/r)))

    return Djs
        
    
    
    
#get power spectrum from CLASS for a given parameter set
#currently set to TT only
#modified from https://github.com/lesgourg/class_public/wiki/Python-wrapper
def get_power(params, l_min, l_max):

    #CLASS gives results in natural units
    #convert to muK^2 to match data
    T_cmb = 2.7255e6 #temp in microkelvins
    #create an instance of CLASS wrapper w/correct params
    cosmo = Class()
    cosmo.set(params)
    #cosmo.set({'output':'tCl,pCl,lCl,mPk','lensing':'yes','P_k_max_1/Mpc':3.0})
    
    cosmo.compute()
    
    #lensed cl until l=l_max
    output = cosmo.raw_cl(l_max)#lensed_cl(l_max)
    ls = output['ell'][l_min:]
    Cls = output['tt'][l_min:]

    Dls = ls*(ls+1)*Cls*T_cmb**2/(2*np.pi)
    
    #clean ups
    cosmo.struct_cleanup()
    cosmo.empty()
    
    return ls, Cls, Dls



#function for ignoring comments in ini files
def ignore_comment(line):
    if '#' in line:
        #save all elements up to the #
        line = line[:line.find('#')]
    if '*' in line:
        line = ''

    return line


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