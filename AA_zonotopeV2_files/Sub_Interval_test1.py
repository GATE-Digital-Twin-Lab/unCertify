# ===========================================================================
# Example 1 — Branin function
#
#   f(x₁, x₂) = a(x₂ − b·x₁² + c·x₁ − r)² + s(1−t)cos(x₁) + s
#   Input domain: x ∈ [0, 1]²  (scaled internally to the standard Branin domain)
#   Known global minima ≈ 0.397 at three locations
# ===========================================================================

# %%
import numpy as np
from itertools import product

from Affine_ArithmeticClassV2 import AffineZonotope

def branin(x):
    """Branin function — accepts 1-D ndarray x of shape (2,)."""
    x1 = 15 * x[0] - 5
    x2 = 15 * x[1]

    a = 1
    b = 5.1 / (4 * np.pi**2)
    c = 5  / np.pi
    r = 6
    s = 10
    t = 1 / (8 * np.pi)

    return a * (x2 - b * x1**2 + c * x1 - r)**2 + s * (1 - t) * np.cos(x1) + s



def split_interval(lower, upper, splits_per_dim):
    """
    lower, upper      : bound lists
    splits_per_dim    : list of ints, one split count per dimension
                        e.g. [3, 2] → 3 splits in dim-0, 2 splits in dim-1
    """
    n_dims = len(lower)  # infer dimension from l_branin length
    assert len(splits_per_dim) == n_dims, (
        f"Expected {n_dims} split values, got {len(splits_per_dim)}"
    )

    # Build subinterval edges per dimension
    dim_intervals = []
    for d in range(n_dims):
        grid = np.linspace(lower[d], upper[d], splits_per_dim[d] + 1)
        intervals = [(grid[i], grid[i+1]) for i in range(splits_per_dim[d])]
        dim_intervals.append(intervals)

    # Cartesian product → all subinterval combinations
    subintervals = []
    for combo in product(*dim_intervals):
        sub_lower = [c[0] for c in combo]
        sub_upper = [c[1] for c in combo]
        subintervals.append((sub_lower, sub_upper))

    return subintervals


# --- Original bounds stay untouched ---
l_branin = [0.0, 0.0]
u_branin = [1.0, 1.0]

# results = {}
# m=1; n=0
# splits_per_dim = [m+1, n+1]
# subintervals = split_interval(l_branin, u_branin, splits_per_dim)

# # Reset accumulators for each (m, n) configuration
# lo_Ch, hi_Ch = np.inf, -np.inf
# lo_MR, hi_MR = np.inf, -np.inf

# for ind, (lb, ub) in enumerate(subintervals):

#     # --- Chebyshev fit ---
#     az_branin = AffineZonotope(branin, lb, ub)
#     az_branin.fit_chebyshev(solver='trust-constr', verbose=False)
#     az_branin.build_zonotope()
#     dum1 = az_branin.f_range

#     lo_Ch = min(lo_Ch, dum1[0])
#     hi_Ch = max(hi_Ch, dum1[1])

#     # --- Min-range fit ---
#     az_branin.fit_min_range(solver='CLARABEL')
#     az_branin.build_zonotope()
#     dum2 = az_branin.f_range

#     lo_MR = min(lo_MR, dum2[0])
#     hi_MR = max(hi_MR, dum2[1])

#     # Store result for this (m, n) configuration
#     results[(m+1, n+1)] = {
#         'f_rangeCh' : [lo_Ch, hi_Ch],
#         'f_rangeMR' : [lo_MR, hi_MR],
#     }

#     print(f"splits={splits_per_dim}  |  "
#         f"Chebyshev: [{lo_Ch:.6f}, {hi_Ch:.6f}]  |  "
#         f"Min-range: [{lo_MR:.6f}, {hi_MR:.6f}]")

dim1 = 15
dim2 = 15

# Store results per split configuration
results = {}

for m in range(dim1):
    for n in range(dim2):

        splits_per_dim = [m+1, n+1]
        subintervals = split_interval(l_branin, u_branin, splits_per_dim)

        # Reset accumulators for each (m, n) configuration
        lo_Ch, hi_Ch = np.inf, -np.inf
        lo_MR, hi_MR = np.inf, -np.inf

        for ind, (lb, ub) in enumerate(subintervals):

            # --- Chebyshev fit ---
            az_branin = AffineZonotope(branin, lb, ub)
            # az_branin.fit_chebyshev(solver='trust-constr', verbose=False)
            az_branin.fit_chebyshev(tol= 1e-8, solver='SLSQP', verbose=False)

            az_branin.build_zonotope()
            dum1 = az_branin.f_range

            lo_Ch = min(lo_Ch, dum1[0])
            hi_Ch = max(hi_Ch, dum1[1])

            # --- Min-range fit ---
            # az_branin.fit_min_range(solver='CLARABEL')
            az_branin.fit_min_range(solver='OSQP')
            az_branin.build_zonotope()
            dum2 = az_branin.f_range

            lo_MR = min(lo_MR, dum2[0])
            hi_MR = max(hi_MR, dum2[1])

        # Store result for this (m, n) configuration
        results[(m+1, n+1)] = {
            'f_rangeCh' : [lo_Ch, hi_Ch],
            'f_rangeMR' : [lo_MR, hi_MR],
        }

        print(f"splits={splits_per_dim}  |  "
              f"Chebyshev: [{lo_Ch:.6f}, {hi_Ch:.6f}]  |  "
              f"Min-range: [{lo_MR:.6f}, {hi_MR:.6f}]")
        


# %%
import pickle

# Save after loop completes
with open('branin_results.pkl', 'wb') as f:
    pickle.dump(results, f)
print(f"\nSaved {len(results)} results to branin_results.pkl")
# %%
