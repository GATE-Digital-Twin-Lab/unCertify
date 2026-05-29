# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 15:22:15 2025

@author: P.Hristov


General note on problems and issues:
    - Passing

"""
from typing import Union

import numpy as np
import scipy.stats as sts
import Interval as ival

import matplotlib.pyplot as plt

STEPS = 100
TRUNC_L = 0.001
TRUNC_R = 0.999
PREC = 4 #Number of decimal places rounded to

round_left = lambda num, decp: np.floor(np.array(num) * 10**decp)/10**decp
round_right = lambda num, decp: np.ceil(np.array(num) * 10**decp)/10**decp
 

class pbox:
    def __init__(self, left:Union[list, tuple, np.ndarray],
                       right:Union[list, np.ndarray]=None,
                       p_left:Union[list, np.ndarray]=None, #Heights - Should be one less than left
                       p_right:Union[list, np.ndarray]=None, #Heights - Should be one less than right
                       n_step:int=STEPS,
                       verbatim_prob:bool=False):
        
        decp = int(np.log10(n_step)) #Number of decimal places to round to based on 
        n_step_l = np.nan
        n_step_r = np.nan
        
        if right is None: #We can still have all four types of numbers
                left = np.array(left)
                r,c = _get_array_dims_(left)
                
                if 0 <= r <= 1 and c == 0: #Scalar
                    right = left.reshape(-1)
                    left = left.reshape(-1)
                elif (r,c) == (2,0): #Interval
                    right = [left[1]] 
                    left = [left[0]]
                elif r > 2 and c == 0: #Distribution
                    right = left #Make sure no simulateneous modifications are happening
                elif c == 2:  #pbox, incl. interval looking ones
                    if isinstance(left[0], ival.Interval):
                        left = _ival_to_array_(left)
                    left = np.array(left) #Enable slicing
                    right = left[:,1]
                    left = left[:,0]
                else: raise(Exception('Unrecognized format for object of type pbox.'))
        else:
            if not hasattr(left, '__iter__') and not hasattr(right, '__iter__'): #Only the scalar case needs handling
                left = [left]
                right = [right]
            
        left = np.sort(left, axis=0)
        right = np.sort(right, axis=0)
            
        if p_right is None:
            if hasattr(p_left, '__iter__'): #Filters None and scalar probs
                p_left = np.array(p_left)    
                if hasattr(p_left[0], '__iter__'): #No intervals and no distributions        
                    p_right = p_left[:,1]
                    p_left = p_left[:,0]
                else:
                    p_right = p_left #Make sure no simulateneous modifications are happening
            else: #This may be from p_left = None -> don't throw an error
                p_right = np.array([round_right(i/len(right), decp) for i in range(1,len(right))]) #np.arrange doesn't contain representation error well
        else:
            if hasattr(p_left, '__iter__') and hasattr(p_left[0], '__iter__'):
                raise(Exception(
                    "You have passed a nested list-like structure for p_left, "
                    "but have also supplied value for p_right. "
                    "This is conflicting behaviour. Specify one or the other."))
            if verbatim_prob:
                n_step_r = pbox._verbatim_prob_calc_(p_right)
            else:
                p_right = round_right(p_right, decp)
                
        if p_left is None: #Not already set above
            p_left = np.array([round_left(i/len(left), 6) for i in range(1,len(left))]) #np.arrange doesn't contain representation error well
        else:
            if verbatim_prob:
                n_step_l = pbox._verbatim_prob_calc_(p_left)
                n_step = max(n_step_l, n_step_r)
            else:
                p_left = round_left(p_left, decp)
        
        
        x_left = np.ones((n_step,1)) * left[0] 
        x_right = np.ones((n_step,1)) * right[0]
        
        inds = np.int32(round_right(p_right*n_step,0)) #Compute quantile bounds with appropriate probability bounds
        for i, ind in enumerate(inds):
            x_left[ind:] = left[i+1]
        
        inds = np.int32(np.round(p_left*n_step,0)) #!!! Round left bound with normal logic to avoid floating-point-related inflation - a stop-gap
        for i, ind in enumerate(inds):
            x_right[ind:] = right[i+1]

        self.n_step = n_step
        self.left = x_left
        self.right = x_right
        self.p = np.arange(0, 1, 1/n_step) #This is assuming 'nice' n_step will be provided
        self.range = pbox.get_range(self) #!!! Use outward rounding here
        self.width = pbox.get_width(self)
        
    def __str__(self):
        #No distinction between an interval, distro, and a p-box for now, and no moment information
        return f'pbox(range=({round_left(self.range[0], 3):0.3f}, {round_right(self.range[1], 3):0.3f}), width={self.width:0.3f})'
        
    def __repr__(self):
        return self.__str__()
    
    #### SOME OPERATORS ####
    def __eq__(self, other):
        if self.left == other.left & self.right == other.right:
            return True
        return False
            
    
    
    #### OTHER METHODS ####
    def plot(self, c='k', ls='-', lw=1, ax=None, label=None):
        cl = cr = c
        if type(c) == list:
            if type(c[0]) == list:
                cl = c[0]
                cr = c[1]
            else: 
                if type(c[0]) == str:
                    cl = c[0]
                    cr = c[1]
                else: cl = cr = c
                
        if not ax:
            ax = plt.subplot()
            
        h = ax.stairs(np.concatenate([self.p, [1]]), np.concatenate(
            [self.left[0], self.left.reshape(self.n_step), self.left[-1]]),
                  baseline=None, color=cl, linestyle=ls, linewidth=lw, label=label)
        ax.stairs(np.concatenate([self.p, [1]]), np.concatenate(
           [self.right[0], self.right.reshape(self.n_step), self.right[-1]]),
                  baseline=None, color=cr, linestyle=ls, linewidth=lw)
        
        ax.plot([self.left[0], self.right[0]],[0,0], color=cr, linestyle=ls, linewidth=lw) #Bottom horizontal line
        ax.plot([self.left[-1], self.right[-1]], [1,1], color=cr, linestyle=ls, linewidth=lw)
        ax.set_ylim((0,1.05))
        
        return h
    
    def to_array(self, condense=False):
        ulr = np.concatenate([self.left, self.right], axis=1)
        if condense:
            ulr, ind = np.unique(ulr, axis=0, return_index=True)
            up = self.p[ind]
            up = np.concatenate([np.delete(up,0), np.array(1).reshape(-1,)])
            return ulr, up    
        return ulr
    
    def to_list(self, condense=False):
        lr = self.to_array(condense=condense)
        if condense:
            p = lr[1]
            lr = lr[0]
            lst = [tuple(e) for e in lr]
            return lst, p
        
        lst = [tuple(e) for e in lr]
        return lst
       
    def to_risk_calc(self, condense=False):
        lr = self.to_array(condense=condense)
        # for lri, pi in zip(lr[0], lr[1]):
        p = lr[1]
        left = lr[0][:,0]
        right = lr[0][:,1]
        
        str_l = f'@({left[0]},0),'
        str_r = f'({right[0]},0),'
        for i in range(len(left)-1):
            str_l += f'({left[i]}, {p[i]}), ({left[i+1]}, {p[i]}),'
            str_r += f'({right[i]}, {p[i]}), ({right[i+1]}, {p[i]}),'
            
        str_l += f'({left[i+1]}, 1),'
        str_r += f'({right[i+1]}, 1)@'
        str_lr = str_l + str_r
        return str_lr
    
    def condense(self, p_condense=None):
        x_o = np.c_[self.left, self.right]
        p_o = self.p
        if not p_condense:
            p_condense = np.arange(0, 1, 1/STEPS)
        xq = condense_parts(x_o, p_o, p_condense)
        
        return pbox(xq, p_left=p_condense)
    
    def ppf(self, p:np.ndarray): #This should return an interval
        if type(p) is float or type(p) is np.float64: p = np.array([p]) #Accommodate single quantiles
        
        p[np.where(p==0)] = 0.01 #To compute the correct indices
        inds = np.int16(p*self.n_step-1)
        
        if p.shape[0] > 1:
            x = []
            for ind in inds:
                # x.append(tuple([*self.left[ind], *self.right[ind]])) # Tuple as a placeholder for an interval
                x.append(ival.I(*self.left[ind], *self.right[ind])) #Interval class
        else: x = ival.I(*self.left[inds][0], *self.right[inds][0])
        return x
    
    def get_cspi(self, alpha=0.05, precision=PREC):
        '''Compute the conservative symmetric probability interval (CSPI) with
        1-alpha coverage probability.'''
        
        return ival.I(self.ppf(alpha/2).leftval, self.ppf(1-alpha/2).rightval,
                      precision=precision)
    
    def mean(self): #Standard bounds-based mean
        return [np.mean(self.left), np.mean(self.right)]
    
    def get_range(self):
        return (*self.left[0], *self.right[-1])
    
    def get_width(self):
        return np.mean(self.right - self.left)
    
    def get_interval(self, left, right=None, precision=PREC):
        if not right:
            if type(left) is float: #Central interval - left is interpreted as content
                right = 0.5 + left/2
                left = 0.5 - left/2
            elif (type(left) is list or type(left) is np.ndarray): #Bounds - take verbatim
                right = left[1]
                left = left[0]
            else:
                raise(Exception('Unsupported type of interval provided.'))    
        else:
            if type(left) is not float or type(right) is not float:
                raise(Exception('Unsupported type of interval provided.'))
        
        l = round_right(self.ppf(left)[0][0], precision)
        r = round_left(self.ppf(right)[0][1], precision)
        
        if r-l < 1e-8: return l
        
        return (l,r)       
    
    def _verbatim_prob_calc_(p):
        from re import split
        p_s = str(p).strip('[]')
        # decp = min(len(max(split('\.|\s', p_s), key=len)), PREC)
        # Test to work with integral probabilities
        spl = split('\.|\s', p_s)
        if len(spl) == len(p): #Integers only
            decp = 0
        else: decp = min(len(max(spl, key=len)), PREC)
        #End test
        n_step = 10**decp
        
        return n_step
 
    
 
def _get_array_dims_(nparray):
    sh = nparray.shape
    if len(sh) == 0: return 0, 0
    if len(sh) == 1: #1-D array; but this may be an array of Interval
        if isinstance(nparray[0], ival.Interval): return sh[0], int(2)
        return sh[0], int(0)
    return sh

def _ival_to_array_(ndarray_of_ival):
    left_right = [[ival.leftval, ival.rightval] for ival in ndarray_of_ival]
 
    return left_right

def prep_focal(*pboxes):
    p = []
    for pbox in pboxes:
        p.append(pbox.to_list(True)[1])
        
    p = np.unique(np.concatenate(p))
    p_ppf = p.copy() #Otherwise p's get modified via pbox.ppf
    return [pbox.ppf(p_ppf) for pbox in pboxes], p

def condense_parts(x_orig:np.ndarray, p_orig:np.ndarray, p_cond:np.ndarray):
    '''Condense a the components of a p-box (quantiles and probabilities)
    to a specified list of probabilities. If the original quantisation 
    is coarser than that specified in p_cond, the resulting components will
    not change.
    
    Use condense_parts to condense the components of a p-box, before
    creating one, for example to reduce memory requirements.
    
    To condense an existing p-box, use the pb.condense method.'''
     
    lq = []
    rq = []
    
    # p_cond = np.array(p_cond)
    p_cond = np.round(p_cond, PREC)
    p_orig = np.round(p_orig, PREC) #Test
    
    for pi in p_cond:
        # lq.append(np.where(np.abs(p_orig[:,0] - pi) <= 1e-6)[0][0])
        # rq.append(np.where(np.abs(p_orig[1:,1] - pi) <= 1e-6)[0][0])
        # lq.append(np.where(p_orig[:,0] <= pi)[0][-1]) # p should never be two-dimensional
        # rq.append(np.where(p_orig[1:,1] >= pi)[0][0])
        lq.append(np.where(p_orig <= pi)[0][-1])
        rq.append(np.where(p_orig[1:] >= pi)[0][0])
    
    xq_0 = np.concatenate([ [x_orig[0,0]], x_orig[lq,0] ])
    xq_1 = np.concatenate([ x_orig[rq,1], [x_orig[-1,1]] ])

    xq = np.vstack([xq_0, xq_1]).T  
    
    return xq


def imp(pboxes):
    left, right, n_step = [], [], []
    
    for pb in pboxes:
        left.append(pb.left)
        right.append(pb.right)
        n_step.append(pb.n_step)
    
    left_imp = np.max(left, axis=0)
    right_imp = np.min(right, axis=0)
    n_step = max(n_step)
    
    if any(left_imp > right_imp):
        res = None #Conisider making this a special p-box
    else:
        res = pbox(left_imp, right_imp, n_step=n_step)
    return res
     
def distro(distro_name, loc, scale, p, trunc=None):
    '''trunc can be None (default = no truncation) or a 2-list of floats containing
    the quantiles for truncation. For distributions with (semi-)infinite support
    appropriate trunction happens automatically, based on the global parameters
    TRUNC_L and TRUNC_R, if not specified otherwise.'''
    
    if not trunc: 
        trunc = [TRUNC_L, TRUNC_R]
    
    trunc_bounded = [0, 1]

    dist_obj = getattr(sts, distro_name)(loc, scale)
    sup = dist_obj.support()
    
    #Only use the truncation for (semi-) infinite distributions, as otherwise a 
    #truncated distribution will have to be created, which does not necessarily
    #exist
    trunc_bounded = np.where(np.isinf(sup), trunc, trunc_bounded) 
        
    # n = len(p)
      
    lb = dist_obj.ppf(trunc_bounded[0])
    rb = dist_obj.ppf(trunc_bounded[1])
    
    left = dist_obj.ppf(p)
    left = np.concatenate([ [lb], left ])
    
    right = dist_obj.ppf(p)
    right = np.concatenate([ right, [rb] ])
    
    # left = dist_obj.ppf(p[:n-1])
    # right = dist_obj.ppf(p[1:])
    
    res = pbox(left, right, p, p, verbatim_prob=True)
    
    return res
    