from utilities import * 
from mcmc import *


from classy import Class
import time


def change_starting_place(pars_list):

	new_pars_list = np.zeros(len(pars_list))
	for i in range(len(pars_list)):
		param = float(pars_list[i])
		new_pars_list[i] = np.random.normal(param, abs(param*0.1)) 

	return new_pars_list

if __name__ == '__main__':

	l_min = 90
	l_max = 2000

	model_pars = read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/Axion-MCMC/axion_MCMC/')
	model_pars['n_axion'] = 3

	start = time.time()

	p_propose = [model_pars['log10_axion_ac'], model_pars['log10_fraction_axion_ac'], model_pars['omega_cdm'], model_pars['H0']] #whatever params you're varying in MCMC

	#assign new starting places for each chain
	p_propose = change_starting_place(p_propose)

	model_pars['log10_axion_ac'] = p_propose[0]
	model_pars['log10_fraction_axion_ac'] = p_propose[1]
	model_pars['omega_cdm'] = p_propose[2]
	model_pars['H0'] = p_propose[3]

	#_, _, Dl_propose = get_power(model_pars, l_min, l_max)
	
	end = time.time()
	print(end-start)
	print(np.array([[1,2,3],[1,2,3],[4,5,6]]))