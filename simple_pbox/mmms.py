# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 09:49:33 2026

@author: P.Hristov
"""
# from pyuncertainnumber import Interval as I #This is too difficult to use
import Interval as ival
import numpy as np
import matplotlib.pyplot as plt

class pbox:
    defaultSigDigs = 6
    defaultNumSteps = 100
    
    
    def __init__(self, left=None, right=None, mean=None, std=None):
        #Allow empty p-boxes to be constructed
        self.left = left
        self.right = right
        self.mean = mean
        self.std = std
        
        if (left is not None) and (right is not None): self.range = self.get_range()
        else: self.range = np.inf
        self.n_step = pbox.defaultNumSteps
        self.p = np.arange(0, 1, 1/self.n_step)
        
    def __str__(self):
        ml = self.mean.leftval
        mr = self.mean.rightval
        sl = self.std.leftval
        sr = self.std.rightval
        return f"P-box (range=[{self.range[0]}, {self.range[1]}], mean=[{ml:0.{self.defaultSigDigs}g}, {mr:0.{self.defaultSigDigs}g}], std=[{sl:0.{self.defaultSigDigs}g}, {sr:0.{self.defaultSigDigs}g}])"
    
    def __repr__(self):
        ml = self.mean.leftval
        mr = self.mean.rightval
        sl = self.std.leftval
        sr = self.std.rightval
        return f"P-box (range=[{self.range[0]}, {self.range[1]}], mean=[{ml:0.{self.defaultSigDigs}g}, {mr:0.{self.defaultSigDigs}g}], std=[{sl:0.{self.defaultSigDigs}g}, {sr:0.{self.defaultSigDigs}g}])"
        
    def get_range(self):
       return (*self.left[0], *self.right[-1])
   
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
    
    
def mmms(minimum, maximum, mean:ival.I, stddev:ival.I):
    #if (nothing(maximum - minimum)) return RandomNbr((minimum + maximum) / 2.0);
    mean = ival.I(mean)
    stddev = ival.I(stddev)
    
    zero = 0.0; one = 1.0
    p = x2 = x3 = x4 = x5 = x6 = rng = maximum - minimum;
    m = constrain(mean, ival.I(minimum, maximum), "(mean)") #Interval
    s = constrain(stddev, ival.envelope(ival.I(0.0), ival.sqrt(ival.abs(rng * rng / 4.0 - (maximum - mean - rng / 2.0)**2))), " (dispersion)");
    ml = (m.leftval - minimum) / rng; sl = s.leftval / rng
    mr = (m.rightval - minimum) / rng; sr = s.rightval / rng

    z = pbox()
    n = z.n_step
    
    u = np.full((n,1), np.nan) #This is upper probably => left quantile bound
    d = np.full((n,1), np.nan) #This is lower probability => right quantile bound
    
    for i in range(n):
        p = i / n
        if p == zero: x2 = zero
        else: x2 = ml - sr * np.sqrt(one / p - one)
        if (ml + p) <= one: x3 = zero
        else:
            x5 = p * p + sl * sl - p;
            if x5 >= zero:
                x4 = one - p + np.sqrt(x5);
                if x4 < ml: x4 = ml;
            else: x4 = ml
            x3 = (p + sl * sl + x4 * x4 - one) / (x4 + p - one);
            
        if (p <= zero) or (p <= (one - ml)): x6 = zero
        else: x6 = (ml - one) / p + one
        u[i] = np.max([x2, x3, x6, zero]) * rng + minimum;

        p = (i + 1) / n #Clever way to update p for the right bound
        if p >= one: x2 = one
        else: x2 = mr + sr * np.sqrt(one / (one / p - one))
        if mr + p >= one: x3 = one
        else:
            x5 = p * p + sl * sl - p
            if x5 >= zero:
                x4 = one - p - np.sqrt(x5)
                if x4 > mr: x4 = mr
            else: x4 = mr
            x3 = (p + sl * sl + x4 * x4 - one) / (x4 + p - one) - one
        
        if ((one - mr) <= p) or (one <= p): x6 = one
        else: x6 = mr / (one - p)
        d[i] = np.min([x2, x3, x6, one]) * rng + minimum;
    
    z = pbox(u,d,m,s) #This is a necessary waste because there are no getter and setter methods
    # z.distrib = RandomNbr::RangeMoments;
    return z

def constrain(a:ival.I, b:ival.I, par):
    c = a - b
    if not c.straddles():
        raise Exception(f"Math Problem: impossible parameter {par}.")
    return ival.imposition(a, b)
    