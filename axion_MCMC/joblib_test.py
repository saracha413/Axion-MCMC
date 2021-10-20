from joblib import Parallel, delayed
import numpy as np
from mcmc import *
import time
from multiprocessing import Pool

if __name__ == '__main__':
	l_min = 90
	l_max = 2000
	num_steps = 10
	num_chains = 3
	num_burn_in = 0
	saveFile = False
	runFromFile = False


	n_axion = 3

	model_pars = read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/Axion-MCMC/axion_MCMC/')
	model_pars['n_axion'] = n_axion

	params = {'num_burn_in': num_burn_in,
              'l_min':  l_min, 'l_max': l_max, 'model_pars':model_pars }
	#get parameters to start chains from

	print('started')
	start = time.time()
	#print(Parallel(n_jobs=1)(delayed(np.sqrt)(i**2) for i in range(10)))
	pool = Pool(processes=num_chains)
	n_trials_per_process = [num_steps] * num_chains
	total_sampling_result = pool.map(mcmc, n_trials_per_process)

	#total_sampling_result = Parallel(n_jobs=num_chains, verbose=10) (delayed(mcmc)(num_steps) for sample_idx in range(num_chains))

	print('finished')
	end = time.time()
	pool.close()
	print('total exec time: ', end-start)