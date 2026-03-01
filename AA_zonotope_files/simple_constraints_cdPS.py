# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 16:09:39 2026

@author: petar.hristov and Premjit Saha
"""
import numpy as np
import sympy as sp

def get_conditions_by_alt(altitide_m):
    #These formulas are valid for the troposphere and at no offset
    lapse_rate = -0.0065 #K/m
    T0_K = 288.15
    
    temp_K = T0_K + lapse_rate*altitide_m
    pressure_pa = (temp_K/T0_K)**5.256 * 1.0133e5
    rho_kgpm3 = (temp_K/T0_K)**4.256 * 1.225
    
    temp_c = temp_K - 273.15
    
    return temp_c, rho_kgpm3, pressure_pa

def map_to_sea_level(altitude_m, propulsion='piston'):
        """Altitude corrections, depending on propulsion system (or propulsion system type)"""

        temp_c, rho_kgpm3, pressure_pa = get_conditions_by_alt(altitude_m)
        
        tcorr = 1  # Default correction value (No correction required) for Thrust

        if propulsion == 'piston':
            #Gagg-Ferrar model. Multiply by this to get power at given density
            sigma = rho_kgpm3 / 1.225 # Density ratio
            tcorr = 1.132 * sigma - 0.132 #Why is this assigned to the thrust correction - PH?

        return tcorr

def thrust_to_weight_take_off(wing_load,
                              Cl_max_to, Cl_to, Cd_to,
                              dist_gr,  mu_r, rnwy_elev=0):
    '''This is the averaged T/W required during take-off, as forces in this
    phase are in continuous flux.
    
    The function returns the sea-level static thrust required.'''
    
    rho_to = get_conditions_by_alt(rnwy_elev)[1] #Only standard ISA
    
    ttw = 1.21 * wing_load / rho_to / Cl_max_to / 9.806 / dist_gr +\
        0.5 * (Cd_to/Cl_to + mu_r)
        
    sttw_sl = ttw * map_to_sea_level(rnwy_elev)
    
    return sttw_sl

def induceddragfact(aspectr):
    """Lift induced drag factor k estimate (Cd = Cd0 + K.Cl^2) based on the relationship
        (k = 1 / pi * AR * e_0).

    **Parameters:**
    
    aspectr
        float or ival.I: aspect ratio of the wing
            
    
    **Outputs:**

    induceddragfactor
        float, an estimate for the coefficient of Cl^2 in the drag polar (Cd = Cd0 + K.Cl^2)

    """
    sqrtterm = 4 + aspectr**2
    oswaldeff = 2 / (2 - aspectr + sp.sqrt(sqrtterm))

    return 1.0 / (np.pi * aspectr * oswaldeff)


def thrust_to_weight_cruise(wing_load, cruise_alt, V_cruise, Cd_min, aspectr):
    '''This is the thrust required to cruise at the requested altitude and
    speed, given as sea-level static thrust.'''
    
    rho_cruise = get_conditions_by_alt(cruise_alt)[1]
    q = 0.5 * rho_cruise * V_cruise**2
    
    k = induceddragfact(aspectr)
    
    ttw = q/wing_load*Cd_min + k/q*wing_load

    sttw_sl = ttw * map_to_sea_level(cruise_alt)
    
    return sttw_sl



def thrust_to_power_propeller(tw, airspeed_mps, eta_prop):
    return tw * airspeed_mps / eta_prop 
    # Mapping to sea-level thrust via map_thrust_to_sea_level; no mapping to
    # static thrust or sea-level power for piston engines