###Sara Vannah



import axiclass_mcmc as mcmc
import multiprocessing as mp
from multiprocessing import Process



if __name__ == '__main__':

    processlist = []

    num_cores = 2 #mp.cpu_count()
    
    input_pars = mcmc.read_ini_file('example_axiCLASS.ini', loc='/Users/saravannah/axiclass/')

    pars_list = [input_pars for i in range(num_cores)]

   


    #initiate and start all processes
    for i in range(num_cores):
        #starting at -2.5, increment value of a_c for each walker
        pars_list[i]['log10_axion_ac'] = -2.5+i*0.1
        
        new_file_name = 'weetwoot_'+str(i)+'.txt'
        
        #open the parallel process on processor
        process = Process(target=mcmc.MCMC_run, args=(pars_list[i], 10, new_file_name))
        processlist.append(process)
        process.start()
        
    #join processes
    for process in processlist:
        process.join()

    print('All done!')




