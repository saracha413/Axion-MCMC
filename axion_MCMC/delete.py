from utilities import * 
from mcmc import *

import numpy as np
import multiprocessing as mp
from classy import Class
from classy import CosmoComputationError

import os, os.path, copy

from interruptingcow import timeout


#let 
def change_starting_place(pars_list):

	new_pars_list = np.zeros(len(pars_list))
	for i in range(len(pars_list)):
		param = float(pars_list[i])
		new_pars_list[i] = np.random.normal(param, abs(param*0.1)) 

	return new_pars_list

def func(num_it):

	num_burn_in = 0
	l_min = 90
	l_max = 2000
	n_axion = 3
	log10_axion_ac = -3.531
	log10_fraction_axion_ac = -0.879
	omega_cdm = 0.132
	H0 = 72.81
	print('entered MCMC func')

	l_min, l_max = 90, 2000
	model_pars = read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/Axion-MCMC/axion_MCMC/')
	model_pars['n_axion'] = 3

	#log10_axion_ac = float(model_pars['log10_axion_ac'])
	#log10_fraction_axion_ac = float(model_pars['log10_fraction_axion_ac'])
	#omega_cdm = float(model_pars['omega_cdm'])
	#H0 = float(model_pars['H0'])
	sampling_result = []


	p_current = [model_pars['log10_axion_ac'], model_pars['log10_fraction_axion_ac'], model_pars['omega_cdm'], model_pars['H0']] #whatever params you're varying in MCMC

	#assign new starting places for each chain
	p_current = change_starting_place(p_current)
	stdDevs = [str(abs(float(x)*0.1)) for x in p_current]


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

	for i in range(num_it):

		p_propose = np.random.normal(p_current, stdDevs)

        ####TO-DO: write this fxn to use whatever variable param you want. hard-coded for now
		model_pars['log10_axion_ac'] = p_propose[0]
		model_pars['log10_fraction_axion_ac'] = p_propose[1]
		model_pars['omega_cdm'] = p_propose[2]
		model_pars['H0'] = p_propose[3]



		#log10_axion_ac = np.random.normal(log10_axion_ac, abs(log10_axion_ac*0.1))
		#log10_fraction_axion_ac = np.random.normal(log10_fraction_axion_ac, abs(log10_fraction_axion_ac*0.1))
		#omega_cdm = np.random.normal(omega_cdm, abs(omega_cdm*0.05))
		#H0 = np.random.normal(H0, abs(H0*0.1))
		try:
			with timeout(60, exception=RuntimeError):
				_, _, Dl_propose = get_power(model_pars, l_min, l_max)
		except (RuntimeError, CosmoComputationError):
			pass
			#print('didnt finish within 5 seconds')

		JSD_propose = JSD(Dl_propose, Dl_data)
		x = JSD_propose/JSD_current

		if x < 1+np.random.uniform():
			p_current = p_propose
			JSD_current = JSD_propose


		line = np.append(p_current, JSD_current)
		sampling_result.append(line)

	return sampling_result


if __name__ == '__main__':
	num_chains = 10
	num_steps = 200
	saveFile = True
	fileName = 'Oct-20.txt'

	print('starting code')

	model_pars = read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/Axion-MCMC/axion_MCMC/')
	model_pars['n_axion'] = 3


	pool = mp.Pool(processes=24)
	print('The pool is open! Come swim!')

	n_trials_per_process = [num_steps] * num_chains
	#result = pool.map(func, n_trials_per_process)
	#_, _, Dls = get_power(model_pars, 0, 2000)
	total_sampling_result = pool.map(func, n_trials_per_process)
	#result = pool.map(mcmc, n_trials_per_process)
	pool.close()
	print('The pool is now closed.')

	#if saveFile:
	#	for i in range(num_chains):
	#
	#		chainFileName = fileName.split('.tx')[0] + '#'+str(i)+'.txt'
	#		np.savetxt(total_sampling_result[i],chainFileName)
	#
	#
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
		chain_arr = np.array([log10_axion_ac[i], log10_fraction_axion_ac[i], omega_cdm[i], H0[i], Djs[i]])
		chainFileName = fileName.split('.tx')[0] + '#'+str(i)+'.txt'
		np.savetxt(chainFileName, chain_arr)





#	if saveFile:
#	#create total array to combine all the chain data
#		big_arr = np.zeros((5, len(log10_axion_ac[0])*num_chains))
#
#		big_arr[0] = np.concatenate(([log10_axion_ac[i] for i in range(num_chains)]))
#		big_arr[1] = np.concatenate(([log10_fraction_axion_ac[i] for i in range(num_chains)]))
#		big_arr[2] = np.concatenate(([omega_cdm[i] for i in range(num_chains)]))
#		big_arr[3] = np.concatenate(([H0[i] for i in range(num_chains)]))
#		big_arr[4] = np.concatenate(([Djs[i] for i in range(num_chains)]))
#		#save run to file
#		with open(fileName, 'a') as fileObject:
#			np.savetxt(fileObject, np.transpose(big_arr))
                
           