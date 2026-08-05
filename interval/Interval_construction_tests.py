# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 14:07:48 2026

@author: petar.hristov
"""

#%% Imports
import Interval as ival

import numpy as np
#%% Empty constructor 
x = ival.I() #Not allowed; if you want an empty (None) interval specify it explicitly
x = ival.I(None) #but this may not be of any practical use

#%% Two arguments constructor
x = ival.I(1,7)
print(x) #Prints the interval with outward rounding to the specified precision

#%% Change the precision
ival.I.PREC = 3
x = ival.I(1,7)
print(x) #Prints the interval with outward rounding to the specified precision

#%% Iterable constructor - single interval - 1D array like
x = ival.I([1,7])
print('List:\t\t', x) #Prints the interval with outward rounding to the specified precision

x = ival.I((1,7))
print('Tuple:\t\t', x) #Prints the interval with outward rounding to the specified precision

x = ival.I(np.array([1,7]))
print('np.ndarray:\t', x) #Prints the interval with outward rounding to the specified precision

#%% Construct an array-like of intervals
x = ival.interval_array([[1,7], [2,4], [10,12]])
print('Interval array: ', x) #Prints the interval with outward rounding to the specified precision

