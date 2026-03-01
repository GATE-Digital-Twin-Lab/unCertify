# %%
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 14:39:57 2026

@author: P.Hristov and Premjit Saha
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import sympy as sp
import cvxpy as cp

import importlib
from scipy.optimize import root
from scipy.stats import qmc

import AA_zonotope as aa

from Affine_ArithmeticClass import AAZonotopeModel # class fle developed for AA
                                                   # and zonotope plot

importlib.reload(aa)
# %%

#-------------------------------------------------------------------------
#                            ~~~ Some 2-D Examples ~~~
#-------------------------------------------------------------------------

# Problem formulation

# epsilon variable
eps = sp.symbols('eps0', real=True)
vars = [eps]

# interval
V_base = [[-1.0, 2.0]]

# affine variable
x = aa.IAVar_AAVar(eps, V_base[0])


# target 2-D functions
#------------------------------------------

f_target = x**2               # Example-1
# f_target = x**2 - 2*x + 1   # Example-2
# f_target = x*(x+1)          # Example-3

# %%

# %%
#             ~~~~~~~~~~~ Chebyshev minimax approximation ~~~~~~~~~~


model = AAZonotopeModel(vars, f_target, V_base)

c_Ch, delta_Ch = model.chebyshev_fit(n_starts=11, M=100)
center_Ch, G, vertices_Ch, f_rangeCh = model.zonotope(c_Ch, delta_Ch)

print("optimal param (Chebyshev minimax):", c_Ch)
print("optimal error (Chebyshev minimax):", delta_Ch)

print("Zonotope vertices (Chebyshev minimax):")
print(vertices_Ch)

print("f_target range (Chebyshev minimax method):")
print(f_rangeCh)

model.plot_zonotope_2d(
    c_Ch,
    vertices_Ch,
    xlabel="x",
    ylabel="f(x)",
    title="2D Zonotope",
    grid_points=101,
    fill_alpha=0.3
)

# %%
#              ~~~~~~~~~~~ Minimum range approximation ~~~~~~~~~~~

c_lp, delta_lp = model.affine_min_rangeLP(M=50, solver=cp.HIGHS)
center_MR, G, vertices_MR, f_rangeMR = model.zonotope(c_lp, delta_lp)


print("optimal param (Minimum range):", c_lp)
print("optimal error (Minimum range):", delta_lp)

print("Zonotope vertices (Minimum range):")
print(vertices_MR)

print("f_target range (Minimum range method):")
print(f_rangeMR)

model.plot_zonotope_2d(
    c_lp,
    vertices_MR,
    xlabel="x",
    ylabel="f(x)",
    title="2D Zonotope",
    grid_points=101,
    fill_alpha=0.3
)

# %%
#-------------------------------------------------------------------------
#                            ~~~ Some 3-D Examples ~~~
#-------------------------------------------------------------------------


""" x, y = sp.symbols('x y')
f_target = f(x,y)
vars = [x, y] """

# interval parameters

V_base = [[-0.25, 0.5], [-0.25, 0.5]]

N = len(V_base)

# eps vectors
eps = sp.symbols(f'eps0:{N}', real=True)

# vars must be flat list of symbols
vars = list(eps)

# affine variables
x = [
    aa.IAVar_AAVar(eps[i], V_base[i])
    for i in range(N)
]

# Example: 3D functions
#-------------------------------------------------------------

# f_target = x[0]**2 + x[1]**2  # Example-1
# f_target = x[0]**2 + x[1]**2 + sp.Rational(1,2)*x[0] - x[1] # Example-2
f_target = x[0]**2 - x[1]**2 # Example-3

#--------------------------------------------------------------
# For x-label, y-label and z-label, useful in plotting 
# the function

labels =[            
    "f(x,y)",
    "x",
    "y"
]
# %%
#             ~~~~~~~~~~~ Chebyshev minimax approximation ~~~~~~~~~~

model = AAZonotopeModel(vars, f_target, V_base)

c_Ch, delta_Ch = model.chebyshev_fit(n_starts=101, M=1001)
center_Ch, G, vertices_Ch, f_rangeCh = model.zonotope(c_Ch, delta_Ch)

print("optimal param (Chebyshev minimax):", c_Ch)
print("optimal error (Chebyshev minimax):", delta_Ch)

print("Zonotope vertices (Chebyshev minimax):")
print(vertices_Ch)

print("f_target range (Chebyshev minimax method):")
print(f_rangeCh)

model.plot_zonotope_3d(
    vertices=vertices_Ch,
    var_names= labels[1:],
    f_name=labels[0],
    grid_points=101,
    surface_alpha=0.8,
    poly_alpha=0.25
)

# %%
#              ~~~~~~~~~~~ Minimum range approximation ~~~~~~~~~~~

c_lp, delta_lp = model.affine_min_rangeLP(M=1001, solver=cp.HIGHS)
center_MR, G, vertices_MR, f_rangeMR = model.zonotope(c_lp, delta_lp)


print("optimal param (Minimum range):", c_lp)
print("optimal error (Minimum range):", delta_lp)

print("Zonotope vertices (Minimum range):")
print(vertices_MR)

print("f_target range (Minimum range method):")
print(f_rangeMR)

model.plot_zonotope_3d(
    vertices=vertices_MR,
    var_names=labels[1:],
    f_name=labels[0],
    grid_points=101,
    surface_alpha=0.8,
    poly_alpha=0.25
)