# %%
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 14:39:57 2026

@author: P.Hristov and Premjit Saha
"""
# import os

# If you have to change the working directory in Jupyter notebook

# change it to your computer path
#----------------------------------------------------------------
# os.chdir(r"C:\Users\premj\Documents\GATE Institute\Digital Twin Lab\Interval Arithmatic codes")

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import sympy as sp
import cvxpy as cp

import importlib
from scipy.optimize import root
from scipy.stats import qmc

import simple_constraints_cdPS as sc
import AA_zonotope as aa

from Affine_ArithmeticClass import AAZonotopeModel # class fle developed for AA
                                                   # and zonotope plot

importlib.reload(sc)
importlib.reload(aa)

# %% Inputs
# Example, common sense inputs - [a,b] denotes an interval; the codes do not currently work with intervals

# Performance parameters - they are not known precisely during conceptual design
cl_max_to = [1.4, 1.7] #[1.4, 1.7] #Max lift coefficient in t-o configuration
cl_to = [0.8, 1] #[0.8, 1] #Lift coefficient in t-o configuration
cd_to = [0.05, 0.12] #[0.05, 0.12] #Drag coefficient in t-o configuration
cd_min_clean = [0.01, 0.03] #[0.01, 0.03] #Minimum drag coefficient in the clean configuration
mu_r = [0.01, 0.05] #[0.01, 0.05] #Coefficient of rolling resistance

# Design definition parameters
aspect_ratio = [9,12] #[9,12] #This will usually come from planform, span limitations and wing loading

# Design brief specs
groundrun_m = 300
cruisealt_m = 5000
cruisespeed_mps = 100

# Wing loading sweep - uncomment the appropriate
#wing_ld_pa_list = np.arange(700, 1600, 2.5) #Use for visual constraint analysis
wing_ld_pa = 1200 #Use to test particular value

# %%  Pre-processing the interval variables for this specific example
base_params = [
    cl_max_to,
    cl_to,
    cd_to,
    cd_min_clean,
    mu_r,
    aspect_ratio
]

mid_base_params = [0.5 * (p[0] + p[1]) for p in base_params]

x = sp.symbols('x0:6', real=True)

AA_vars = []   # interval parameters to AA parameters

for i, param in enumerate(base_params):
    AA_vars.append(aa.IAVar_AAVar(x[i], param))

# %%  This cell is for plot label purpose
param_names = [
    "C_L,max,TO",
    "C_L,TO",
    "C_D,TO",
    "C_D,min,clean",
    "μ_r",
    "Aspect Ratio"
]

y_names = [
    "Thrust-to-Weight-cruise",
    "Thrust-to-Weight-Take-off"
]

# %% 
#----------------------------------------------------------------------------
#                                ~~~ Example 1 ~~~
#----------------------------------------------------------------------------

#-----------------------------------------------
#  thrust_to_weight_cruise case
#-----------------------------------------------
# Select case
import cases
import importlib
importlib.reload(cases)
vars = x

from cases import case_1

case = case_1 

result = case(sc, wing_ld_pa, cruisealt_m, cruisespeed_mps, vars, AA_vars, base_params)

f_target = result["f_target"]
vars     = result["vars"]
V_base   = result["V_base"]
para_indices = result["param_indices"] 


#-----------------------------------------

labels = [param_names[i] for i in para_indices] # X-axis and Y-axis label name

case_id = para_indices[0]    
labels[0] = y_names[case_id] # f_target  name

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

c_lp, delta_lp = model.affine_min_rangeLP(M=50, solver=cp.HIGHS)
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

# %%
#----------------------------------------------------------------------------
#                                ~~~ Example 2 ~~~
#----------------------------------------------------------------------------

#-----------------------------------------------
#  thrust_to_weight_take_off case
#-----------------------------------------------
importlib.reload(cases)
vars = x

from cases import case_2, case_3, case_4, case_5, case_6, case_7
case = case_3  # change to case_2, case_3 etc.

result = case(sc, wing_ld_pa, vars, AA_vars, base_params, mid_base_params, groundrun_m)

f_target = result["f_target"]
vars     = result["vars"]
V_base   = result["V_base"]
para_indices = result["param_indices"]

#-----------------------------------------

labels = [param_names[i] for i in para_indices] # X-axis and Y-axis label name

case_id = para_indices[0]    
labels[0] = y_names[case_id] # f_target  name

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

# %%

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

c_lp, delta_lp = model.affine_min_rangeLP(M=50, solver=cp.HIGHS)
center_MR, G, vertices_MR, f_rangeMR = model.zonotope(c_lp, delta_lp)


print("optimal param (Minimum range):", c_lp)
print("optimal error (Minimum range)", delta_lp)

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

# %%

#--------------   3-D and 4-D cases  ----------------

importlib.reload(cases)
vars = x

from cases import case_8, case_9, case_10, case_11, case_12
case = case_11  # change to case_8, case_9, ... case_11 etc.  3-D case
                # change to case_12 for 4-D case

result = case(sc, wing_ld_pa, vars, AA_vars, base_params, mid_base_params, groundrun_m)

f_target = result["f_target"]
vars     = result["vars"]
V_base   = result["V_base"]

#-----------------------------------------

labels = [param_names[i] for i in para_indices] # X-axis and Y-axis label name

case_id = para_indices[0]    
labels[0] = y_names[case_id] # f_target  name

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

# %%

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

c_lp, delta_lp = model.affine_min_rangeLP(M=50, solver=cp.HIGHS)
center_MR, G, vertices_MR, f_rangeMR = model.zonotope(c_lp, delta_lp)


print("optimal param (Minimum range):", c_lp)
print("optimal error (Minimum range)", delta_lp)

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





#-----------------------------( Unnecessary parts )-----------------------------------


# # %%
# #--------------------------------------------------------------------
# #                  Chebyshev minimax approximation
# #--------------------------------------------------------------------

# solt = aa.chebyshev_fit_Nd(
#     f_target, # target function
#     vars,     # independent variables
#     V_base,   # Interval varables
#     n_starts=11, # initilization for search of stationary points
#     M=100,    # number of points chosen randomly in the latin hypercube constructed by vars independent variables
#     bounds=(-1.0, 1.0) # bounds for noise parameters, epscilons
# )

# print("optimal param",solt[0])
# print("optimal error",solt[1])
# print("maximum error",solt[-1])

# center_ch, Gmat, vertices_Ch = aa.zonotope_vertices(solt, V_base)

# print("Zonotope vertices (Chebyshev minimax):")
# print(vertices_Ch)
# # %%
# #--------------------------------------------------------------------
# #                     Minimum range approximation
# #--------------------------------------------------------------------

# cvec, delta = aa.affine_min_range_lp(
#     f_target, # target function
#     vars,     # independent variables
#     V_base,   # Interval varables
#     M=50,     # sample number of chosen randomly in the latin hypercube constructed by vars independent variables
#     solver=cp.HIGHS,  # solver chosen for solving LP problem 
#     bounds = (-1.0, 1.0)
# )

# print("cvec =", cvec)
# print("delta =", delta)

# solt = (cvec, delta)

# center_MR, Gmat, vertices_MR = aa.zonotope_vertices(solt, V_base)

# print("Zonotope vertices (Minimum range):")
# print(vertices_MR)

# %%

# x_vals = np.linspace(V_base[0][0], V_base[0][1], 101)
# y_vals = np.linspace(V_base[1][0], V_base[1][1], 101)
# z_vals = np.linspace(-1, 1, 101)

# f_num = sp.lambdify(vars, f_target, modules="numpy")


# X, Y = np.meshgrid(x_vals, y_vals)
# X1, Y1 = np.meshgrid(z_vals, z_vals)
# Z = f_num(X1, Y1)
# # %%
# #-----------------------------------------------------------------------
# #                           Draw the polytope box
# #------------------------------------------------------------------------

# from scipy.spatial import ConvexHull
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# vertices = vertices_Ch   # Chebyshev minimax vertices
# # vertices = vertices_MR   # Minimum range vertices

# #------------------------------

# param_names = [
#     "C_L,max,TO",
#     "C_L,TO",
#     "C_D,TO",
#     "C_D,min,clean",
#     "μ_r",
#     "Aspect Ratio"
# ]

# y_names = [
#     "Thrust-to-Weight-cruise",
#     "Thrust-to-Weight-Take-off"
# ]
# labels = [param_names[i] for i in para_indices]

# case_id = para_indices[0]
# labels[0] = y_names[case_id]

# #------------------------------
# # %%

# hull = ConvexHull(vertices)
# faces = [vertices[simplex] for simplex in hull.simplices]

# fig = plt.figure()
# ax = fig.add_subplot(111, projection="3d")

# poly = Poly3DCollection(
#     faces,
#     alpha=0.25,
#     edgecolor="black"
#     # edgecolor="None"

# )

# ax.add_collection3d(poly)

# # Plot vertices (optional but useful)
# ax.scatter(
#     vertices[:, 0],
#     vertices[:, 1],
#     vertices[:, 2],
#     color="red",
#     s=20
# )

# # Set equal aspect ratio
# mins = vertices.min(axis=0)
# maxs = vertices.max(axis=0)
# ax.set_xlim(mins[0], maxs[0])
# ax.set_ylim(mins[1], maxs[1])
# ax.set_zlim(mins[2], maxs[2])
# # ax.set_box_aspect(maxs - mins)
# # ax.view_init(elev=25, azim=135) 
# ax.set_xlabel(labels[1])
# ax.set_ylabel(labels[2])
# ax.set_zlabel(labels[0])
# ax.plot_surface(X, Y, Z, cmap="viridis", alpha=0.8)
# plt.show()
