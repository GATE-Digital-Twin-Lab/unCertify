# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 15:29:26 2026

@author: petar.hristov
"""

import Interval as ival
from Interval_dependence import *

import numpy as np

#%% Addition tests

#%% Subtraction tests

#%% Multiplication tests 
## Easy to test only for r = {-1,0,1}
intervals = [ival.I(10, 20), ival.I(-5, 1), ival.I(-10,-3)]

for ix in intervals:
    for iy in intervals:
        print(f"Interval x: {ix}, Interval B: {iy}\n")
        for test_r in [1.0, 0.5, 0.0, -0.5, -1.0]:
            res = multiply_dep_2(ix, iy, test_r)
            if test_r == 0.0:
                true = ix*iy
                print(f"r = {test_r:4.1f} -> Computed result = {res}; Exact result = {true}")
            elif test_r == 1.0:
                iz = ival._generate_points_(ix) * ival._generate_points_(iy)
                true = ival.I(np.min(iz), np.max(iz))
                print(f"r = {test_r:4.1f} -> Computed result = {res}; Exact result = {true}")
            elif test_r == -1.0:
                iz = ival._generate_points_(ix) * ival._generate_points_(iy)[-1::-1]
                true = ival.I(np.min(iz), np.max(iz))
                print(f"r = {test_r:4.1f} -> Computed result = {res}; Exact result = {true}")
            else: print(f"r = {test_r:4.1f} -> Computed result = {res}")

#%% Division tests - defined only for non-straddling divisors
## Easy to test only for r = {-1,0,1}
intervals = [ival.I(10, 20), ival.I(-5, 1), ival.I(-10,-3)]

for ix in intervals:
    for iy in intervals:
        print(f"Interval x: {ix}, Interval B: {iy}\n")
        for test_r in [1.0, 0.5, 0.0, -0.5, -1.0]:
            try:
                res = divide_dep_1(ix, iy, test_r)
            except Exception as e:
                print(e)
                continue
            
            if test_r == 0.0:
                true = ix/iy
                print(f"r = {test_r:4.1f} -> Computed result = {res}; Exact result = {true}")
            elif test_r == 1.0:
                iz = ival._generate_points_(ix) / ival._generate_points_(iy)
                true = ival.I(np.min(iz), np.max(iz))
                print(f"r = {test_r:4.1f} -> Computed result = {res}; Exact result = {true}")
            elif test_r == -1.0:
                iz = ival._generate_points_(ix) / ival._generate_points_(iy)[-1::-1]
                true = ival.I(np.min(iz), np.max(iz))
                print(f"r = {test_r:4.1f} -> Computed result = {res}; Exact result = {true}")
            else: print(f"r = {test_r:4.1f} -> Computed result = {res}")
