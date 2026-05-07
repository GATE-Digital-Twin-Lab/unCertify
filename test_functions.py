# -*- coding: utf-8 -*-
"""
Created on Mon May  4 15:08:24 2026

@author: petar.hristov
"""
import numpy as np


def branin(x):
    '''
    x is a 1x2 vector 
    '''
    #% Input scaling
    # x1 = x#x[:,0] #Uniform scales
    # x2 = y#x[:,1]
    
    x1 = 15*x[0]-5
    x2 = 15*x[1]
    # x1 = x[0]
    # x2 = x[1]
    
    #% Mean
    a = 1
    b = 5.1/4/np.pi**2
    c = 5/np.pi
    r = 6
    s = 10
    t = 1/(8*np.pi)
    
    y = (a*(x2-b*x1**2+c*x1-r)**2+s*(1-t)*np.cos(x1)+s)+5*x1

    return y

def ackley(x, a=20, b=0.2, c=2*np.pi):
    d = len(x) #Only accept single array-like x's
    
    sum1 = 0
    sum2 = 0
    for xi in x:
    	sum1 += xi**2
    	sum2 += np.cos(c*xi)
    
    term1 = -a * np.exp(-b*np.sqrt(sum1/d))
    term2 = -np.exp(sum2/d)
    
    return term1 + term2 + a + np.exp(1)

def egg(x):
    '''Domain = [-512, 512]^2'''
    x1 = -512 + 1024*x[0]
    x2 = -512 + 1024*x[1]
    
    term1 = -(x2+47) * np.sin(np.sqrt(np.abs(x2+x1/2+47)))
    term2 = -x1 * np.sin(np.sqrt(np.abs(x1-(x2+47))))
    
    return term1 + term2
    

def piston(x):
    #================================
    #  OUTPUT AND INPUT:
    # 
    #  C = cycle time
    #  x = [M, S, V0, k, P0, Ta, T0]
    #================================
    # M ∈ [30, 60]         	piston weight (kg)
    # S ∈ [0.005, 0.020]	piston surface area (m2)
    # V0 ∈ [0.002, 0.010]	initial gas volume (m3)
    # k ∈ [1000, 5000]	    spring coefficient (N/m)
    # P0 ∈ [90000, 110000] atmospheric pressure (N/m2)
    # Ta ∈ [290, 296]   	ambient temperature (K)
    # T0 ∈ [340, 360]   	filling gas temperature (K)
    #================================
    
    M  = 30 + 30*x[0]
    S  = 0.005 + 0.015*x[1]
    V0 = 0.002 + 0.008*x[2]
    k  = 1000 + 4000*x[3]
    P0 = 90_000 + 20_000*x[4]
    Ta = 290 + 6*x[5]
    T0 = 340 + 20*x[6]
    
    Aterm1 = P0 * S
    Aterm2 = 19.62 * M
    Aterm3 = -k*V0 / S
    A = Aterm1 + Aterm2 + Aterm3
    
    Vfact1 = S / (2*k)
    Vfact2 = np.sqrt(A**2 + 4*k*(P0*V0/T0)*Ta)
    V = Vfact1 * (Vfact2 - A)
    fact1 = M
    fact2 = k + (S**2)*(P0*V0/T0)*(Ta/(V**2))
    
    C = 2 * np.pi * np.sqrt(fact1/fact2)
    
    return C

def wingweight(x):
    #================================
    #  OUTPUT AND INPUT:
    # 
    #  y  = wing weight
    #  x = [Sw, Wfw, A, LamCaps, q, lam, tc, Nz, Wdg, Wp]
    #================================
    # Sw ∈ [150, 200]	    wing area (ft2)
    # Wfw ∈ [220, 300]	    weight of fuel in the wing (lb)
    # A ∈ [6, 10]	        aspect ratio
    # Λ ∈ [-10, 10]        quarter-chord sweep (degrees)
    # q ∈ [16, 45]	        dynamic pressure at cruise (lb/ft2)
    # λ ∈ [0.5, 1]	        taper ratio
    # tc ∈ [0.08, 0.18]	    aerofoil thickness to chord ratio
    # Nz ∈ [2.5, 6] 	    ultimate load factor
    # Wdg ∈ [1700, 2500]   	flight design gross weight (lb)
    # Wp ∈ [0.025, 0.08]	paint weight (lb/ft2)
    #================================
    
    Sw      = 150 + 50*x[0]
    Wfw     = 220 + 80*x[1]
    A       = 6 + 4*x[2]
    LamCaps = -10 + 20*x[3] * np.pi/180
    q       = 16 + 29*x[4]
    lam     = 0.5 + 0.5*x[5]
    tc      = 0.08 + 0.1*x[6]
    Nz      = 2.5 + 3.5*x[7]
    Wdg     = 1700 + 800*x[8]
    Wp      = 0.025 + 0.055*x[9]
    
    fact1 = 0.036 * Sw**0.758 * Wfw**0.0035
    fact2 = (A / np.cos(LamCaps)**2)**0.6
    fact3 = q**0.006 * lam**0.04
    fact4 = (100*tc / np.cos(LamCaps))**(-0.3)
    fact5 = (Nz*Wdg)**0.49
    
    term1 = Sw * Wp
    
    y = fact1*fact2*fact3*fact4*fact5 + term1
    
    return y