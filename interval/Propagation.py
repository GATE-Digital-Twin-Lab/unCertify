# -*- coding: utf-8 -*-
"""
I
Incertitude_propagation_methods.py
Incertitude propagation methods test script


Created on Fri May  6 11:29:21 2022

@author: peter
"""

import numpy as np
import Interval as ival

from itertools import product

from typing import Union

CACHE = {'samples_dict':dict(),
         'eval_count':0,
         'saved_signs':None}

# global samples_dict, eval_count, saved_signs
# samples_dict = dict()
# eval_count = 0
# saved_signs = None


def vertex_method(intervals:Union[list, np.ndarray], fun, individual_dims=False,
                  return_samples=False):
    '''Propagate the input incertitude defined via the array-like of interval
    objects in 'intervals', through the numerical model in 'fun' using the
    vertex method.
    If the model in 'fun' accepts a single, list-like argument of size
    len(intervals) for its d-dimensional input, keep individual_dims=False,
    otherwise set it to True. The function in 'fun' must return a single output.
    
    TO DO:
        Vectorised input
        Multidimensional output
    '''
    
    num_comb = 2**len(intervals)
    samples = np.full((num_comb, 1), np.nan)
    
    interval_ends = ival.to_array(intervals)
    
    if individual_dims:
        for i, p in enumerate(product(*interval_ends)):
            samples[i] = fun(*p)
    else:
        for i, p in enumerate(product(*interval_ends)):
            samples[i] = fun(p)
       
    y = ival.I(np.min(samples), np.max(samples))
    
    # dat_ymin = y.index(yint[0])
    # dat_ymax = y.index(yint[1])
    # dat_int = [INPUT[dat_ymin], INPUT[dat_ymax] ]
    if return_samples: return y, samples
    
    return y

def make_extreme_point_propagator(monotonic_function):
    '''This version will have the ability to accept previously computed signs
    to save d+1 model runs'''
    
    # TODO: open an input to specify if signs should be recomputed regardless of cache
    # This is important for SIR and functions whose extremum occurs inside the domain
    
    #TODO: make this a wrapper also for VP and potentially other deterministic UP approaches
    
    def propagate(intervals, out_use=0, cache=True, out_cache=0):
        saved_signs = CACHE['saved_signs']
        intervals = ival.to_array(intervals)
        
        
        if saved_signs is None or len(saved_signs) != len(intervals):
            saved_signs = []
            
            lower_bounds = [low for low, high in intervals]
            baseline = check_cached_samples(lower_bounds, out_use, cache, out_cache)
                        
            for i, (low, high) in enumerate(intervals):
                test_point = lower_bounds.copy()
                test_point[i] = high
                value = check_cached_samples(test_point, out_use, cache, out_cache)
                difference = value - baseline
                saved_signs.append(1 if difference > 0 else (-1 if difference < 0 else 0))
            CACHE['saved_signs'] = saved_signs 
            
        max_point = [
            high if sign > 0 else low
            for (low, high), sign in zip(intervals, saved_signs)
        ]
        min_point = [
            high if sign < 0 else low
            for (low, high), sign in zip(intervals, saved_signs)
        ]
        
        max_value = check_cached_samples(max_point, out_use, cache, out_cache)
        min_value = check_cached_samples(min_point, out_use, cache, out_cache)

        return ival.I(min_value, max_value)

    def check_cached_samples(test_point, out_use, cache, out_cache):
        s_tstpt = str(test_point)
        if not cache or CACHE['samples_dict'].get(s_tstpt) is None:
            value = monotonic_function(test_point)
            if hasattr(value, '__len__'): #Multioutput
                value_use = value[out_use]
                value_cache = value[out_cache]
            else: #Single output
                value_use = value
                value_cache = [value]
                
            CACHE['samples_dict'][s_tstpt] = value_cache #Comment out to turn caching off
            CACHE['eval_count'] += 1
        else:
            try:
                value_use = CACHE['samples_dict'][s_tstpt][out_use]
            except:
                raise(Exception(f'''Incorrect size of cached values.
                                Requested index was {out_use}, but cached
                                values have size {len(CACHE['samples_dict'][s_tstpt])}'''))
        return value_use
    
    return propagate


def samplingMethod(intervals, fun, n, method='montecarlo', endpoints=False):
    from scipy.stats import qmc
    #Propagate the input incertitude defined via 'intervals'
    #(list of lists or ndarray) through the numerical model
    #in 'fun' using a sampling-based method.
    #
    #The number of samples is specified by the integer 'n'.
    #
    #The default 'method' is Monte Carlo, whereby uniform
    #distributions with support on the intervals will be
    #sampled randomly and samples will be passed to 'fun'.
    #
    #Alternatively a Latin hypercube sampling can be chosen using
    #method='lhs'.
    #
    #The 'endpoints' input determines whether the interval edges
    #are included in the sampling plan (endpoints=True). The edges
    #are excluded by default.
    #
    #Assuming 'fun' accepts a single argument (an ndarray)
    #of size len(intervals) and returns a single output.
    
    x = np.array(intervals).reshape(-1,2) #Just in case
    lo = x.T[0]
    hi = x.T[1]
    #print(x)
    if method == 'montecarlo':
        X = np.random.rand(n, x.shape[0])
    elif method == 'lhs':
        sampler = qmc.LatinHypercube(x.shape[0])
        X = sampler.random(n)
        
    if endpoints:
        X = X - np.min(X, axis=0)
        X = X / np.max(X, axis=0)
        
    X = lo + (hi-lo) * X
    
    y = [] #np.zeros(n)
    for i in range(n):
        y.append(fun(X[i]))
        #y[i] = fun(X[i])
        print(f'Completed eval #{i} of {n}...')
   
    exclude ={0}
    y1 =[sublist for sublist in y if not exclude.intersection(sublist)]
    vol_m = [el[1] for el in y1]
    ymin =  [sl for sl in y if sl[1]==min(vol_m)]
    ymax =  [sl for sl in y if sl[1]==max(vol_m)]
    yint = [ymin, ymax]
    
    dat_vol_m_min = next(index for index, value in enumerate(y) if min(vol_m) in value)
    dat_vol_m_max = next(index for index, value in enumerate(y) if max(vol_m) in value)
    dat_int = [X[dat_vol_m_min], X[dat_vol_m_max] ]
    return yint, dat_int

def cauchydeviateMethod(intervals, fun, n):
    from scipy.optimize import brentq
    #Propagate the input incertitude defined via 'intervals'
    #(list of lists or ndarray) through the numerical model
    #in 'fun' using Cauchy deviate method.
    #
    #The number of samples is specified by the integer 'n'.
    #
    #Assuming 'fun' accepts a single argument (an ndarray)
    #of size len(intervals) and returns a single output.
    
    x = np.array(intervals).reshape(-1,2) #Just in case
    lo = x.T[0]
    hi = x.T[1]

    xtilde = (lo + hi) / 2
    Delta = (hi - lo) / 2
    ytilde = np.array(fun(xtilde))
    
    dOut = ytilde.size #Output dimensionality 
    deltaF = np.zeros((n, dOut))
    
    
    X = []
    for k in range(n):
        r = np.random.rand(x.shape[0])
        c = np.tan(np.pi * (r - 0.5))
        K = np.max(c)
        delta = Delta * c / K
        x = xtilde - delta
        X.append(x) #For writing
        deltaF[k] = K * (ytilde - np.array(fun(x)))
        
        #print(f'K:     {K},\nc:     {c},\ndelta: {delta},\nx:     {x},\ndeltaF: {deltaF[k]}')
        print(f'Completed eval #{k} of {n}...')
    
    
    zRoot = []
    for i in range(dOut):
        Z = lambda Del: n/2 - np.sum(1 / (1 + (deltaF.T[i] / Del)**2))
        zRoot.append(brentq(Z, 0.0001, max(deltaF.T[i])/2))
    #Z = lambda Del: n/2 - np.sum(1 / (1 + (deltaF / Del)**2))
    
    return [ytilde - zRoot, ytilde + zRoot], X, deltaF #How do we get to values producing the interval?


def partition_interval(interval, m):
    """Split [a, b] into m equal sub-intervals."""
    a, b = interval
    step = (b - a) / m
    return [[a + i*step, a + (i+1)*step] for i in range(m)]

def cartesian_interval_partitions(arrays, m:list):
    """
    arrays: list of [start, end] intervals
    m: number of partitions per interval
    
    Returns: all m^n combinations of interval partitions
    """
    partitions = [partition_interval(interval, mi) for interval, mi in zip(arrays, m)]
    return list(product(*partitions))

def subinterval_method(intervals:Union[list, np.ndarray], fun, n:Union[int,list],
                       method='vertex',
                       individual_dims=False, return_intervals=False):
    '''Partition the input incertitute in the array-like of interval objects
    ('intervals') into 'n' subinterval and propagate them through the numerical
    model in 'fun' using 'method' (default is vertex) - to become EPP.
    'n' is either an integer or a list with len(n) = len(intervals) to specify
    different number of partitions in each dimension.
    
    About the propagation method
    If the model in 'fun' accepts a single, list-like argument of size
    len(intervals) for its d-dimensional input, keep individual_dims=False,
    otherwise set it to True. The function in 'fun' must return a single output.
    
    To Do:
        Return samples
        Implement an octree mesh
    '''
    
    if type(n) == int: #All inputs have identical division
        n = [n]*len(intervals)
    
    num_int = np.prod(n)
    interval_ends = ival.to_array(intervals)
    partitions = cartesian_interval_partitions(interval_ends, n)
    
    int_out = np.empty(num_int, dtype='object')
    
    match method.lower():
        case 'interval': #Requires a file which understands the Interval class
            for i, partition in enumerate(partitions):
                # int_out[i] = fun(*ival.to_interval(partition))
                int_out[i] = fun(ival.to_interval(partition))
        case 'vertex':
            #Make sure to implement caching
            for i, partition in enumerate(partitions):
                int_out[i] = vertex_method(ival.to_interval(partition), #This is a bit of a work around; improve if possible
                                           fun, individual_dims)
        case 'epp':
            epp_method = make_extreme_point_propagator(fun)
            for i, partition in enumerate(partitions):
                int_out[i] = epp_method(ival.to_interval(partition)) #This is a bit of a work around; improve if possible
            
        case 'sampling':
            pass
        case 'cauchy':
            pass
        case _: raise Exception(f'Propagation method "{method}" is unknown.')
    
    y = ival.env(int_out)
    
    if return_intervals: return y, int_out
    return y