###Sara Vannah



import axiclass_mcmc as mcmc
import multiprocessing as mp
from multiprocessing import Pool


if __name__ == '__main__':

    num_cores = mp.cpu_count()

    input_pars = mcmc.read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/axiclass/')
    pars_list = [input_pars for i in range(num_cores)]
    
    print('number of cores is ', num_cores)

    i = 1
    with Pool(processes = num_cores) as pool:
        new_file_name = 'yeetyoot_'+str(i)+'.txt'
        pars_list[i]['log10_axion_ac'] = -2.5+i*0.1
        
        pool.apply_async(mcmc.MCMC_run, args=(pars_list[i], 10, new_file_name))
        print('working?')
        i+=1
        
    pool.close()
    pool.join()
        
    print('All done!')




