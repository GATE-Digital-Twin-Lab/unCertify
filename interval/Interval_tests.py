# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 13:09:51 2026

@author: petar.hristov
"""

#%% Interval power test table
import Interval as ival

#%% Base: negative, not degenerate
x = ival.I(-7.458, -0.1287)

# Real negative exponents
print("\n======= Real negative exponents =======")
# These cases should work
p = -3; print(f"Result of {x}^{p} = {x**p}") # Integral negative power
p = -1/3; print(f"Result of {x}^{p} = {x**p}") # Fractional negative power with odd denominator
p = -0.76; print(f"Result of {x}^{p} = {x**p}") #Decimal negative power - need a test that fraction has odd denominator

#Real positive exponents
print("\n======= Real positive exponents =======")
p = 3; print(f"Result of {x}^{p} = {x**p}") # Integral positive power
p = 1/3; print(f"Result of {x}^{p} = {x**p}") # Fractional posititve power with odd denominator
p = 0.76; print(f"Result of {x}^{p} = {x**p}") #Decimal positive power

# Interval negative exponents
print("\n======= Interval negative exponents =======")
p = ival.I(-5, -3); print(f"Result of {x}^{p} = {x**p}") # Integral negative power
p = ival.I(-12/5, -1/3); print(f"Result of {x}^{p} = {x**p}") # Fractional negative power with odd denominator
p = ival.I(-2.3561, -0.7566); print(f"Result of {x}^{p} = {x**p}") #Decimal negative power

# Interval straddling exponents
print("\n======= Interval straddling exponents =======")
p = ival.I(-5, 3); print(f"Result of {x}^{p} = {x**p}") # Integral negative power
p = ival.I(-12/5, 1/3); print(f"Result of {x}^{p} = {x**p}") # Fractional negative power with odd denominator
p = ival.I(-2.3561, 0.7566); print(f"Result of {x}^{p} = {x**p}") #Decimal negative power

# Interval positive exponents
print("\n======= Interval positive exponents =======")
p = ival.I(3, 5); print(f"Result of {x}^{p} = {x**p}") # Integral negative power
p = ival.I(1/3, 12/5); print(f"Result of {x}^{p} = {x**p}") # Fractional negative power with odd denominator
p = ival.I(0.7566, 2.3561); print(f"Result of {x}^{p} = {x**p}") #Decimal negative power


#%% Cube root of a negative interval
#%% Wing weight also doesn't work




def _clear_cache_():
    global CACHE
    CACHE = {'samples_dict':dict(),
             'eval_count':0,
             'saved_signs':None}
    
    
    

