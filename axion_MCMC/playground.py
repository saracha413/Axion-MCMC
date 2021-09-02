
import numpy as np
from classy import Class

import matplotlib
from matplotlib import pyplot as plt
from utilities_NEW import *



params = read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/Axion-MCMC/axion_MCMC/')
params['n_axion'] = 3


cosmo = Class()
cosmo.set(params)

cosmo.compute()

output = cosmo.lensed_cl(2000)

ls = output['ell'][90:]

Cls = output['tt'][90:]

Dls = ls*(ls+1)*Cls/(2*np.pi)

cosmo.struct_cleanup()
cosmo.empty()

plt.plot(ls, Dls)

