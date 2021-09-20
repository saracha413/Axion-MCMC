#force MCMC to skip step when it gets stuck 


from mcmc import *
from utilities import *
import matplotlib
import signal
from multiprocessing import Process, TimeoutError, ParamValueError
from classy import Class
from mcmc import *



#instructions for user-defined exceptions from here https://www.programiz.com/python-programming/user-defined-exception
class Error(Exception):
	"""Base class for other exceptions"""
	pass

class ParamValueError(Error):
	"""Raised when strange parameter values prevent AxiCLASS from convering in reasonable time"""
	pass




#borrowed from https://stackoverflow.com/questions/492519/timeout-on-a-function-call
def handler(signum, frame):
	raise ParamValueError

def run_loop(pars):
	model_pars = pars['model_pars']
	l_min = pars['l_min']
	l_max = pars['l_max']

	
	for i in range(10):

		signal.signal(signal.SIGALRM, handler)
		signal.alarm(20)
		try:
			new_pars = model_pars
			new_pars['log10_axion_ac'] = np.random.normal(float(model_pars['log10_axion_ac']), abs(float(model_pars['log10_axion_ac'])*0.2))
			ls, Cls, Dls = get_power(new_pars, l_min, l_max)
		except ParamValueError:
			print('This step took too long! Skipping to next entry.')


if __name__ == '__main__':

	l_min = 90
	l_max = 2000
	num_steps = 2000
	num_burn_in = 0
	name = 'Sept-16_runs'
	n_axion = 3

	num_cores = 10
	processList = []




	model_pars = read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/Axion-MCMC/axion_MCMC/')
	model_pars['n_axion'] = n_axion

	params = {'num_steps': num_steps, 'num_burn_in': num_burn_in,  'name': name, 'l_min':  l_min, 'l_max': l_max, 'model_pars':model_pars }

	#signal.signal(signal.SIGALRM, handler)
	#signal.alarm(20)

	#for i in range(10):
	#	new_pars = model_pars
	#	new_pars['log10_axion_ac'] = np.random.normal(model_pars['log10_axion_ac'], abs(model_pars['log10_axion_ac']*0.2))
	#	ls, Cls, Dls = get_power(new_pars, l_min, l_max)

	pars_list = [params for i in range(num_cores)]

	#saving param values for less typing
	og_log10_axion_ac = float(params['model_pars']['log10_axion_ac'])
	og_log10_fraction_axion_ac = float(params['model_pars']['log10_fraction_axion_ac'])
	og_omega_cdm = float(params['model_pars']['omega_cdm'])
	og_H0 = float(params['model_pars']['H0'])

	for i in range(num_cores):

			pars_list[i]['log10_axion_ac'] = np.random.normal(og_log10_axion_ac, abs(og_log10_axion_ac*0.05))
			pars_list[i]['log10_fraction_axion_ac'] = np.random.normal(og_log10_fraction_axion_ac, abs(og_log10_fraction_axion_ac*0.05))
			pars_list[i]['omega_cdm'] = np.random.normal(og_omega_cdm, abs(og_omega_cdm*0.05))
			pars_list[i]['H0'] = np.random.normal(og_H0, abs(og_H0*0.05))

			pars_list[i]['name'] = name+'#'+str(i)

			#p = Process(target=run_loop, args=(pars_list[i],))
			p = Process(target=mcmc, args=(pars_list[i],))
			processList.append(p)
			p.start()

	for p in processList:
		p.join()



	#try:
	#	run_loop(model_pars)
	#except Exception as exc:
	#	print(exc)


	print('Done!')