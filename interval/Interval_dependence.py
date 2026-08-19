# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 09:15:15 2026

@author: petar.hristov
"""

import numpy as np
import Interval as ival

def add_dep_1(x,y,r=0):
    if r >= 0:
        return x + y
    
    if r == -1:
        return ival.I(x.leftval+y.rightval, x.rightval+y.leftval).make_proper()
    
    x1, x2 = x.bounds()
    y1, y2 = y.bounds()
    
    left = ival.I(
        w(x,-r)+y1, x1+w(y,-r)
        ).make_proper()
    
    right = ival.I(
        x2+w(y,1+r), w(x,1+r)+y2
        ).make_proper()
    
    return ival.envelope(left, right)

def subtract_dep_1(x,y,r=0):
    return add_dep_1(x, -y, -r)

def multiply_dep_1(x, y, r=0):
    """
    Calculates the exact math bounds for multiplying any two intervals 
    (positive, negative, or mixed) under a skewed rectangle model for any r.
    
    Tracks both the polygon vertices and the curved peaks along the walls.
    """
    if not (-1.0 <= r <= 1.0):
        raise ValueError("Correlation r must be between -1 and 1.")
        
    x_min, x_max = x.bounds()
    y_min, y_max = y.bounds()
    
    w_x = x.width()
    w_y = y.width()
    
    # Handle zero-width intervals
    if w_x == 0.0 or w_y == 0.0:
        corners = [x_min * y_min, x_min * y_max, x_max * y_min, x_max * y_max]
        return min(corners), max(corners)

    alpha = 1.0 - abs(r)
    
    # Gather specific segments of the boundary walls.
    # Each segment is defined as a line: b = m*a + c, valid for a between [start_x, end_x]
    segments = []
    
    if r >= 0:
        # Constraint in normalized space: -alpha <= u - v <= alpha
        # Line 1 (Upper Bound of Band): v = u + alpha -> b = y_min + ((a - x_min)/w_x + alpha)*w_y
        m1 = w_y / w_x
        c1 = y_min - m1 * x_min + alpha * w_y
        # This line cuts the box. Find the valid 'a' range inside the rectangle
        start_x1 = max(x_min, x_min + (y_min - y_min - alpha * w_y) / m1) # simplified: x_min
        end_x1 = min(x_max, x_min + (y_max - y_min - alpha * w_y) / m1)
        if start_x1 <= end_x1: segments.append((m1, c1, start_x1, end_x1))
        
        # Line 2 (Lower Bound of Band): v = u - alpha -> b = y_min + ((a - x_min)/w_x - alpha)*w_y
        m2 = w_y / w_x
        c2 = y_min - m2 * x_min - alpha * w_y
        start_x2 = max(x_min, x_min + (y_min - y_min + alpha * w_y) / m2)
        end_x2 = min(x_max, x_min + (y_max - y_min + alpha * w_y) / m2)
        if start_x2 <= end_x2: segments.append((m2, c2, start_x2, end_x2))
        
    else:
        # Constraint in normalized space: 1 - alpha <= u + v <= 1 + alpha
        # Line 1 (Upper Bound of Band): v = 1 + alpha - u
        m1 = -w_y / w_x
        c1 = y_min + (1.0 + alpha) * w_y - m1 * x_min
        start_x1 = max(x_min, (y_max - c1) / m1 if m1 != 0 else x_min)
        end_x1 = min(x_max, (y_min - c1) / m1 if m1 != 0 else x_max)
        # Sort because negative slope swaps direction
        start_x1, end_x1 = min(start_x1, end_x1), max(start_x1, end_x1)
        if start_x1 <= end_x1: segments.append((m1, c1, start_x1, end_x1))
        
        # Line 2 (Lower Bound of Band): v = 1 - alpha - u
        m2 = -w_y / w_x
        c2 = y_min + (1.0 - alpha) * w_y - m2 * x_min
        start_x2 = max(x_min, (y_max - c2) / m2 if m2 != 0 else x_min)
        end_x2 = min(x_max, (y_min - c2) / m2 if m2 != 0 else x_max)
        start_x2, end_x2 = min(start_x2, end_x2), max(start_x2, end_x2)
        if start_x2 <= end_x2: segments.append((m2, c2, start_x2, end_x2))
    
    candidates = []
    
    # Evaluate all segment endpoints (the corners)
    for m, c, s_x, e_x in segments:
        candidates.append(s_x * (m * s_x + c))
        candidates.append(e_x * (m * e_x + c))
        
        # Check for parabola peak along the line: f(a) = m*a^2 + c*a
        if m != 0:
            x_peak = -c / (2.0 * m)
            if s_x < x_peak < e_x:
                candidates.append(x_peak * (m * x_peak + c))
                
    # Add outer rectangle box corners that are allowed by the band
    # We do a quick check to see if the corner pairs satisfy the band restriction
    box_corners = [(x_min, y_min), (x_min, y_max), (x_max, y_min), (x_max, y_max)]
    for a, b in box_corners:
        u = (a - x_min) / w_x
        v = (b - y_min) / w_y
        if r >= 0 and -alpha - 1e-9 <= (u - v) <= alpha + 1e-9:
            candidates.append(a * b)
        elif r < 0 and 1.0 - alpha - 1e-9 <= (u + v) <= 1.0 + alpha + 1e-9:
            candidates.append(a * b)

    return ival.I(np.min(candidates), np.max(candidates))

def divide_dep_1(x, y, r=0):
    """
    Calculates the exact math bounds for dividing interval A by interval B (A / B)
    under a skewed rectangle model for any r.
    
    Guaranteed to find true extrema by checking the geometric vertices.
    """
    if not (-1.0 <= r <= 1.0):
        raise ValueError("Correlation r must be between -1 and 1.")
        
    if y.straddles(): #Should there be an exception for dividing excatly the same interval under perfect correlation???
        raise(Exception('Cannot divide with intervals strddling 0.'))
        
    x_min, x_max = x.bounds()
    y_min, y_max = y.bounds()
    
    w_x = x.width()
    w_y = y.width()
    
    alpha = 1.0 - abs(r)
    
    vertices = []
    
    if r >= 0:
        # Wall 1: a = x_min
        for v in [0.0, alpha]:
            y_val = y_min + v * w_y
            if y_min <= y_val <= y_max: vertices.append((x_min, y_val))
        # Wall 2: a = x_max
        for v in [1.0, 1.0 - alpha]:
            y_val = y_min + v * w_y
            if y_min <= y_val <= y_max: vertices.append((x_max, y_val))
        # Wall 3: b = y_min
        for u in [0.0, alpha]:
            x_val = x_min + u * w_x
            if x_min <= x_val <= x_max: vertices.append((x_val, y_min))
        # Wall 4: b = y_max
        for u in [1.0, 1.0 - alpha]:
            x_val = x_min + u * w_x
            if x_min <= x_val <= x_max: vertices.append((x_val, y_max))
    else:
        # Wall 1: a = x_min
        for v in [1.0, 1.0 - alpha]:
            y_val = y_min + v * w_y
            if y_min <= y_val <= y_max: vertices.append((x_min, y_val))
        # Wall 2: a = x_max
        for v in [0.0, alpha]:
            y_val = y_min + v * w_y
            if y_min <= y_val <= y_max: vertices.append((x_max, y_val))
        # Wall 3: b = y_min
        for u in [1.0, 1.0 - alpha]:
            x_val = x_min + u * w_x
            if x_min <= x_val <= x_max: vertices.append((x_val, y_min))
        # Wall 4: b = y_max
        for u in [0.0, alpha]:
            x_val = x_min + u * w_x
            if x_min <= x_val <= x_max: vertices.append((x_val, y_max))

    # Add standard bounding box corners allowed by the correlation band
    box_corners = [(x_min, y_min), (x_min, y_max), (x_max, y_min), (x_max, y_max)]
    for a, b in box_corners:
        u = (a - x_min) / w_x if w_x != 0 else 0.0
        v = (b - y_min) / w_y if w_y != 0 else 0.0
        if r >= 0 and -alpha - 1e-9 <= (u - v) <= alpha + 1e-9:
            vertices.append((a, b))
        elif r < 0 and 1.0 - alpha - 1e-9 <= (u + v) <= 1.0 + alpha + 1e-9:
            vertices.append((a, b))

    # Clean up duplicate vertices
    unique_vertices = list(set((round(a, 9), round(b, 9)) for a, b in vertices))
    
    # Compute DIVISION (a / b) instead of multiplication!
    quotients = [a / b for a, b in unique_vertices]
    
    return ival.I(np.min(quotients), np.max(quotients))
    
def w(x:ival.Interval, p:float):
    x1 = x.leftval
    x2 = x.rightval
    
    if 0 <= p <= 1:
        return p*(x.rightval - x.leftval) + x.leftval
    else: return x1 if p < 0 else x2



#%% 
import jax
import jax.numpy as jnp

# class DepInterval:
    # def __init__(self, interval, name=None, var_map=None):
    #     """
    #     val_range: tuple (min_val, max_val)
    #     name: optional string ID. If None, an auto-incrementing ID is assigned.
    #     var_map: internal dictionary tracking geometric alignment to other variables.
    #              Format: {var_id: directional_sensitivity_coefficient}
    #     """
    #     self.leftval, self.rightval = interval.bounds()
        
    #     self.mid = interval.mid()
    #     self.h = interval.width() / 2.0
        
    #     if name is not None:
    #         self.id = name
    #     else:
    #         DepInterval._id_counter += 1
    #         self.id = f"v{DepInterval._id_counter}"
            
    #     if var_map is not None:
    #         self.var_map = var_map
    #     else:
    #         # Base Variable setup: a variable is perfectly comonotone (1.0) with itself
    #         self.var_map = {self.id: 1.0}
    
class DepInterval(ival.Interval):

    # A global counter to automatically give every new variable a unique tracking ID
    _id_counter = 0
    
    def __init__(self, left, right=None, name=None, var_map=None, **kwargs):
        super().__init__(left, right, **kwargs)
        self.mid = self.mid()
        self.h = self.width() / 2.0
        
        if name is not None:
            self.id = name
        else:
            DepInterval._id_counter += 1
            self.id = f"v{DepInterval._id_counter}"
            
        if var_map is not None:
            self.var_map = var_map
        else:
            # Base Variable setup: a variable is perfectly comonotone (1.0) with itself
            self.var_map = {self.id: 1.0}
                
    

    @classmethod
    def from_bivariate_correlation(cls, x:ival.Interval, y:ival.Interval, r_xy,
                                   name_x="x", name_y="y"):
        """
        Helper method to create two correlated base variables using your 
        geometric correlation factor r_xy.
        """
        # hx = x.width() / 2.0
        # hy = y.width() / 2.0
        
        # Geometrically decompose r_xy into orthogonal tracking coordinates.
        # x points fully along its own tracking axis.
        # y splits its alignment between x's axis and its own unique axis.
        x_map = {name_x: 1.0}
        y_map = {name_x: float(r_xy), name_y: float(np.sqrt(1.0 - r_xy**2))}
        
        x = cls(x, name=name_x, var_map=x_map)
        y = cls(y, name=name_y, var_map=y_map)
        return x, y

    def _combine(self, other, jax_op):
        """
        Core engine combining two nodes using JAX forward-mode JVP AD.
        """
        if not isinstance(other, DepInterval): #An exception must be added here for the Interval class
            # Treat pure numbers as zero-width constant intervals
            other = DepInterval(ival.I(other), name=f"c_{other}", var_map={})
            
        # 1. Identify all unique underlying variables involved in this operation
        all_keys = set(self.var_map.keys()).union(set(other.var_map.keys()))
        
        # 2. Define the pure mathematical operation for JAX to evaluate
        def prim_op(x_val, y_val):
            return jax_op(x_val, y_val)
        
        # 3. Calculate the nominal output midpoint
        z_mid = float(prim_op(self.mid, other.mid))
        
        # 4. Use jax.jvp to find partial derivatives at the center point
        _, df_dx = jax.jvp(prim_op, (self.mid, other.mid), (1.0, 0.0))
        _, df_dy = jax.jvp(prim_op, (self.mid, other.mid), (0.0, 1.0))
        
        df_dx = float(df_dx)
        df_dy = float(df_dy)
        
        # 5. Propagate sensitivities across all underlying variables
        z_map = {}
        h_z_sq = 0.0
        
        for k in all_keys:
            # Read the sensitivity coefficient of each variable component
            s_val = self.var_map.get(k, 0.0) * self.h
            o_val = other.var_map.get(k, 0.0) * other.h
            
            # Combine them using the calculus chain rule
            z_component_sensitivity = df_dx * s_val + df_dy * o_val
            
            # Save the component to track correlations down the line
            if np.abs(z_component_sensitivity) > 1e-12:
                # We normalize back by dividing by the yet-to-be-calculated total h_z
                z_map[k] = z_component_sensitivity
                h_z_sq += z_component_sensitivity ** 2
                
        h_z = float(np.sqrt(h_z_sq))
        
        if h_z > 1e-12:
            # Finish normalizing the tracking coefficients
            for k in z_map:
                z_map[k] /= h_z
        else:
            z_map = {}
            
        # 6. Build the final output node
        return DepInterval(ival.I(z_mid - h_z, z_mid + h_z), var_map=z_map)

    # --- Overloading standard operators for clean expression syntax ---
    def __add__(self, other): return self._combine(other, lambda x, y: x + y)
    def __radd__(self, other): return self._combine(other, lambda x, y: y + x)
    def __sub__(self, other): return self._combine(other, lambda x, y: x - y)
    def __rsub__(self, other): return self._combine(other, lambda x, y: y - x)
    def __mul__(self, other): return self._combine(other, lambda x, y: x * y)
    def __rmul__(self, other): return self._combine(other, lambda x, y: y * x)
    
    # Custom method to support unary operations like sin(x)
    def sin(self):
        return self._combine(0.0, lambda x, y: jnp.sin(x))

    def get_correlation_with(self, other):
        """
        Calculates the exact geometric correlation coefficient between 
        this interval node and any other compound interval node in the library.
        """
        all_keys = set(self.var_map.keys()).union(set(other.var_map.keys()))
        r = 0.0
        for k in all_keys:
            r += self.var_map.get(k, 0.0) * other.var_map.get(k, 0.0)
        return float(np.clip(r, -1.0, 1.0))

    def __repr__(self):
        self = ival.outerBound(self)
        return f"interval([{self.leftval:0.{self.precision}f}, {self.rightval:0.{self.precision}f}])"
        # return f"Interval([{self.leftval:.4f}, {self.max_val:.4f}], mid={self.mid:.2f}, h={self.h:.2f})"
    
    
#%% Woodpile
# def multiply_dep_1(x,y,r=0):
#     if not (-1.0 <= r <= 1.0):
#         raise ValueError("Correlation r must be between -1 and 1.")
    
#     if r==0: return x*y
    
#     both_pos = x.is_pos and y.is_pos  # Positive intervals
#     both_neg = x.is_neg and y.is_neg  # Negative intervals
    
#     if both_pos and r == 1:
#         return ival.I(x.leftval * y.leftval, x.rightval * y.rightval)
    
#     if both_neg and r == 1:
#         return ival.I(x.rightval * y.rightval, x.leftval * y.leftval)
    
#     if both_pos or both_neg and r == -1: #Need to test this
#         xL, xU = x.bounds() 
#         yL, yU = y.bounds() 
#         prod_max = 0.5*(xL + xU*yU/(yU-yL))
        
#         f = lambda xx: xx * yU - (yU-yL)/(xU-xL)*(xx-xL)
        
#         prod_max = f(prod_max) if ival.inside(prod_max, x)\
#                                else np.maximum(f(xL), f(xU))
   
#     #General r and intervals - only complete correlation models


