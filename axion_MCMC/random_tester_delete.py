import numpy as np


test_arr = [1,3,5,2]
test_std = [0.1,0.1,0.1,0.1]
test_dict = {'one':1, 'three':3, 'five':5, 'two':2}
test_std_dict = {'one':0.1, 'three':0.1,'five':0.1, 'two':0.1}

param_ranges =  {'log10_axion_ac':[-4.6, -3], 'fEDE': [1e-5,0.15], 'omega_cdm': [0.11, 0.14], 'H0': [65, 80]}

print(np.random.normal(test_arr, test_std))
print(list(param_ranges.values()))