from multiprocessing import Pool
from mcmc_NEW import *
from utilities_NEW import *
import matplotlib
from matplotlib import pyplot as plt
import corner #for triangle plots

if __name__ == '__main__':

	num_walkers = 3
	l_min = 90
	l_max = 2000
	num_steps = 1
	num_burn_in = 0
	name = 'short_test_n=3_'
	n_axion = 3


	pool = Pool()

	pars_array = []



	model_pars = read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/Axion-MCMC/axion_MCMC/')
	model_pars['n_axion'] = n_axion
	#model_pars['log10_axion_ac'] = params['log10_axion_ac']
	#model_pars['log10_fraction_axion_ac'] = params['log10_fraction_axion_ac']
	#model_pars['omega_cdm'] = params['omega_cdm']
	#model_pars['H0'] = params['H0']

	params = {'num_steps': num_steps, 'num_burn_in': num_burn_in,  'name': name, 'l_min':  l_min, 'l_max': l_max, 'model_pars':model_pars }


	#saving param values for less typing
	og_log10_axion_ac = float(params['model_pars']['log10_axion_ac'])
	og_log10_fraction_axion_ac = float(params['model_pars']['log10_fraction_axion_ac'])
	og_omega_cdm = float(params['model_pars']['omega_cdm'])
	og_H0 = float(params['model_pars']['H0'])

	#get the starting points for each chain
	for i in range(num_walkers):

		temp_pars = params
		temp_pars['model_pars']['log10_axion_ac'] = np.random.normal(og_log10_axion_ac, abs(og_log10_axion_ac*0.05))
		temp_pars['model_pars']['log10_fraction_axion_ac'] = np.random.normal(og_log10_fraction_axion_ac, abs(og_log10_fraction_axion_ac*0.05))
		temp_pars['model_pars']['omega_cdm'] = np.random.normal(og_omega_cdm, abs(og_omega_cdm*0.05))
		temp_pars['model_pars']['H0'] = np.random.normal(og_H0, abs(og_H0*0.05))
		temp_pars['name'] = params['name']+'_chain_num_'+str(i)

		pars_array.append(temp_pars)


	print('Starting chains')
	results = pool.map(mcmc, pars_array)
	print('MCMC finished')

	pool.close()

	ndim, nsamples = 4, num_steps

	samples = np.array(results)
	samples = samples.T

	#figure = corner.corner(samples, 
	#                       labels=[r"$f_{EDE}(a_c)$", r"$log_{10}(a_c)$", r"$H_0$", r"$\omega _{CDM}$"], 
	#                      quantiles=[0.16, 0.5, 0.84],
	#                       show_titles=True, title_kwargs={"fontsize": 12})