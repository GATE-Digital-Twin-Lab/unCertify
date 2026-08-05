# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 22:53:01 2026

@author: petar.hristov
"""

#!/usr/bin/env python3

"""
Grid-based enclosure test for interval exponentiation.

Assumptions:
    import ival

    base = ival.I(a, b)
    expo = ival.I(c, d)

    result = base ** expo

You may need to adapt interval_bounds() to match the
endpoint names used by your interval class.
"""

import math
from fractions import Fraction

import numpy as np
import Interval as ival


# ============================================================
# CONFIGURATION
# ============================================================

BASE_GRID_POINTS = 301
EXP_GRID_POINTS = 301
TOL = 1e-12


# ============================================================
# ADAPT THIS IF NECESSARY
# ============================================================

def interval_bounds(iv):
    """
    Return (lower, upper) interval bounds.
    """
    return iv.bounds()


# ============================================================
# REAL POWER EVALUATION
# ============================================================

def real_power(x, y):
    """
    Evaluate x**y in the real domain.

    Raises ValueError when the result is not real.
    """

    # Positive base
    if x > 0:
        return x ** y

    # Zero base
    if x == 0:
        if y > 0:
            return 0.0
        raise ValueError("0 raised to non-positive exponent")

    # Negative base
    if float(y).is_integer():
        return x ** int(y)

    frac = Fraction(y).limit_denominator(1000)

    p = frac.numerator
    q = frac.denominator

    # Odd denominator => real root exists

    if q % 2 == 1:
        val = abs(x) ** (p / q)
        if p % 2:
            return -val
        return val
    raise ValueError("Complex-valued result")


# ============================================================
# ENHANCED GRID
# ============================================================

def enhanced_grid(lo, hi, n, special_points):
    pts = list(np.linspace(lo, hi, n))

    pts.append(lo)
    pts.append(hi)

    for s in special_points:
        if lo <= s <= hi:
            pts.append(float(s))

    return sorted(set(pts))


# ============================================================
# SINGLE TEST
# ============================================================

def test_enclosure(a, b, c, d):
    base = ival.I(a, b)
    expo = ival.I(c, d)

    try:
        result = base ** expo
    except Exception as exc:
        return {
            "status": "OPERATOR_EXCEPTION",
            "exception": repr(exc),
            "base": (a, b),
            "exp": (c, d),
        }

    rlo, rhi = interval_bounds(result)

    xs = enhanced_grid(a, b, BASE_GRID_POINTS, [-10, -1, 0, 1, 10])

    ys = enhanced_grid(c,d,EXP_GRID_POINTS,
        [-10,-5,-4,-3,-2,-1,-2/3,-1/2,-1/3,0,1/3,1/2,2/3,1,2,3,4,5,10,]
    )

    values = []
    for x in xs:
        for y in ys:
            try:
                z = real_power(x, y)
            except Exception:
                continue
            if math.isfinite(z):
                values.append(z)

    if not values:
        return {
            "status": "NO_REAL_VALUES",
            "base": (a, b),
            "exp": (c, d),
            "result": result,
        }

    actual_min = min(values)
    actual_max = max(values)

    enclosure_ok = (
        actual_min >= rlo - TOL
        and actual_max <= rhi + TOL
    )

    width_actual = actual_max - actual_min
    width_interval = rhi - rlo

    return {
        "status": "PASS" if enclosure_ok else "FAIL",
        "base": (a, b),
        "exp": (c, d),
        "result": result,
        "actual_min": actual_min,
        "actual_max": actual_max,
        "computed_min": rlo,
        "computed_max": rhi,
        "actual_width": width_actual,
        "result_width": width_interval,
        "overestimation": width_interval - width_actual,
    }


# ============================================================
# CURATED TEST CASES
# ============================================================

TESTS = [

    # --------------------------------------------------------
    # Degenerate base / exponent
    # --------------------------------------------------------

    (2, 2, 3, 3),
    (2, 2, -3, -3),
    (4, 4, 0.5, 0.5),
    (8, 8, 1 / 3, 1 / 3),
    (-8, -8, 1 / 3, 1 / 3),
    (-2, -2, 2, 2),
    (-2, -2, 3, 3),

    # --------------------------------------------------------
    # Positive intervals
    # --------------------------------------------------------

    (1, 2, 2, 2),
    (1, 2, 3, 3),
    (1, 5, -1, -1),
    (1, 5, -2, 2),
    (1, 10, -3, 3),

    # --------------------------------------------------------
    # Touching zero
    # --------------------------------------------------------

    (0, 2, 2, 2),
    (0, 2, 3, 3),
    (0, 2, -1, -1),
    (0, 2, -2, 2),
    (0, 4, 0.5, 0.5),

    # --------------------------------------------------------
    # Spanning zero
    # --------------------------------------------------------

    (-1, 1, 2, 2),
    (-1, 1, 3, 3),
    (-2, 3, 2, 2),
    (-2, 3, 3, 3),
    (-5, 5, 2, 2),
    (-5, 5, 3, 3),
    (-5, 5, 4, 4),

    # --------------------------------------------------------
    # Negative intervals
    # --------------------------------------------------------

    (-4, -1, 2, 2),
    (-4, -1, 3, 3),
    (-4, -1, -1, -1),
    (-8, -1, 1 / 3, 1 / 3),

    # --------------------------------------------------------
    # Exponent intervals
    # --------------------------------------------------------

    (1, 10, 0, 5),
    (1, 10, -5, 0),
    (1, 10, -5, 5),

    # --------------------------------------------------------
    # Negative bases with exponent intervals
    # --------------------------------------------------------

    (-5, -1, 2, 4),
    (-5, -1, 1, 5),

    # --------------------------------------------------------
    # Fractional exponent intervals
    # --------------------------------------------------------

    (1, 100, 0.1, 0.9),
    (1, 16, 0.25, 0.75),

    # --------------------------------------------------------
    # Extreme mixed cases
    # --------------------------------------------------------

    (-10, 10, -3, 3),
    (-100, 100, 2, 3),
    (-100, 100, -2, 2),

    # --------------------------------------------------------
    # Very small values
    # --------------------------------------------------------

    (1e-12, 1e-10, 2, 2),
    (1e-12, 1e-10, -1, -1),

    # --------------------------------------------------------
    # Very large values
    # --------------------------------------------------------

    (1e10, 1e12, 2, 2),
    (1e10, 1e12, -1, -1),

]

#%%
# ============================================================
# MAIN
# ============================================================

def main():

    failures = 0

    print("=" * 80)
    print("INTERVAL EXPONENTIATION GRID ENCLOSURE TEST")
    print("=" * 80)

    for idx, (a, b, c, d) in enumerate(TESTS, start=1):

        report = test_enclosure(a, b, c, d)

        status = report["status"]

        if status == "PASS":

            print(
                f"[PASS] #{idx:02d} "
                f"base=[{a},{b}] "
                f"exp=[{c},{d}]"
            )

        else:

            failures += 1

            print()
            print("=" * 80)
            print(f"FAILURE #{failures}")
            print("=" * 80)

            print(f"base      = [{a}, {b}]")
            print(f"exponent  = [{c}, {d}]")
            print(f"status    = {status}")

            for k, v in report.items():
                if k != "status":
                    print(f"{k:16s}: {v}")

            print()

    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total tests : {len(TESTS)}")
    print(f"Failures    : {failures}")
    print("=" * 80)


if __name__ == "__main__":
    main()
    
    
#%% Debug
a,b,c,d = TESTS[36]
report = test_enclosure(a, b, c, d)

# ival.I(-5,-1) ** ival.I(2,4)