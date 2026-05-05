
# %%

import numpy as np
from itertools import product

from Affine_ArithmeticClassV2 import AffineZonotope


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
 
    return a * (x2 - b * x1**2 + c * x1 - r)**2 + s * (1 - t) * np.cos(x1) + s + 5*x1


# --- Original bounds stay untouched ---
l_branin = [0.0, 0.0]
u_branin = [1.0, 1.0]

# %%

# Store results per split configuration
results_2n = {}
dim1 = 8

for n in range(dim1+1):

    splits_per_dim = [2**n, 2**n]
    subintervals = split_interval(l_branin, u_branin, splits_per_dim)

    # Reset accumulators for each (m, n) configuration
    lo_Ch, hi_Ch = np.inf, -np.inf
    lo_MR, hi_MR = np.inf, -np.inf

    total_f_evals_Ch = 0
    total_f_evals_MR = 0

    for ind, (lb, ub) in enumerate(subintervals):

        # --- Chebyshev fit ---
        az_Ch = AffineZonotope(branin, lb, ub)
        az_Ch.fit_chebyshev(tol=1e-8, solver='SLSQP', verbose=False)
        az_Ch.build_zonotope()
        dum1 = az_Ch.f_range

        lo_Ch = min(lo_Ch, dum1[0])
        hi_Ch = max(hi_Ch, dum1[1])
        total_f_evals_Ch += az_Ch.sip_stats['total_f_evals']

        # --- Min-range fit (fresh object to avoid state contamination) ---
        az_MR = AffineZonotope(branin, lb, ub)
        az_MR.fit_min_range(solver='OSQP', verbose=False)
        az_MR.build_zonotope()
        dum2 = az_MR.f_range

        lo_MR = min(lo_MR, dum2[0])
        hi_MR = max(hi_MR, dum2[1])
        total_f_evals_MR += az_MR.mr_stats['total_f_evals']

    # FIX: key by actual split count (2**n), not n+1
    results_2n[n+1] = {
        'splits_per_dim': splits_per_dim,
        'f_rangeCh'     : [lo_Ch, hi_Ch],
        'f_rangeMR'     : [lo_MR, hi_MR],
    }

    print(f"splits={splits_per_dim}  |  "
        f"Chebyshev: [{lo_Ch:.6f}, {hi_Ch:.6f}]"
        f"total function count={total_f_evals_Ch}   |"
        f"Min-range: [{lo_MR:.6f}, {hi_MR:.6f}]"
        f"total function count={total_f_evals_MR}")
    
# %%
import pickle

# Save after loop completes
with open('branin_results_2n.pkl', 'wb') as f:
    pickle.dump(results_2n, f)
print(f"\nSaved {len(results_2n)} results to branin_results_2n.pkl")

# %%

N = 2**dim1  

# Initialize 2D grids with NaN for empty off-diagonal cells
lo_Ch_grid = np.full((N, N), np.nan)
lo_MR_grid = np.full((N, N), np.nan)
hi_Ch_grid = np.full((N, N), np.nan)
hi_MR_grid = np.full((N, N), np.nan)

# Fill only the diagonal from results_2n
for n in range(dim1+1):
    k = n + 1  # your dict key

    lo_Ch_grid[2**n-1, 2**n-1] = results_2n[k]['f_rangeCh'][0]
    hi_Ch_grid[2**n-1, 2**n-1] = results_2n[k]['f_rangeCh'][1]
    lo_MR_grid[2**n-1, 2**n-1] = results_2n[k]['f_rangeMR'][0]
    hi_MR_grid[2**n-1, 2**n-1] = results_2n[k]['f_rangeMR'][1]

tick_labels = [str(i+1) for i in range(N)]

# %%
import matplotlib.pyplot as plt

# --- Plot ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Chebyshev lower bound
im1 = axes[0].imshow(lo_Ch_grid, origin='lower', aspect='auto', cmap='viridis')
# axes[0].set_xscale('log')
# axes[0].set_yscale('log')
axes[0].set_title('Lower Bound — Chebyshev ($lo_{Ch}$)', fontsize=13)
axes[0].set_xlabel('Splits dim 2 (n)')
axes[0].set_ylabel('Splits dim 1 (m)')
axes[0].set_xticks(range(N)); axes[0].set_xticklabels(tick_labels)
axes[0].set_yticks(range(N)); axes[0].set_yticklabels(tick_labels)
plt.colorbar(im1, ax=axes[0])

# Min-range lower bound
im2 = axes[1].imshow(lo_MR_grid, origin='lower', aspect='auto', cmap='plasma')
# axes[1].set_xscale('log')
# axes[1].set_yscale('log')
axes[1].set_title('Lower Bound — Min-Range ($lo_{MR}$)', fontsize=13)
axes[1].set_xlabel('Splits dim 2 (n)')
axes[1].set_ylabel('Splits dim 1 (m)')
axes[1].set_xticks(range(N)); axes[1].set_xticklabels(tick_labels)
axes[1].set_yticks(range(N)); axes[1].set_yticklabels(tick_labels)
plt.colorbar(im2, ax=axes[1])

plt.suptitle('Branin Function — Lower Bound vs Split Configuration', fontsize=14)
plt.tight_layout()
plt.savefig('branin_lower_bounds.png', dpi=150, bbox_inches='tight')
plt.show()

# %%

# dim1 = 256
# dim2 = 256

# # Store results per split configuration
# results = {}

# for m in range(dim1):
#     for n in range(dim2):

#         splits_per_dim = [m+1, n+1]
#         subintervals = split_interval(l_branin, u_branin, splits_per_dim)

#         # Reset accumulators for each (m, n) configuration
#         lo_Ch, hi_Ch = np.inf, -np.inf
#         lo_MR, hi_MR = np.inf, -np.inf

#         for ind, (lb, ub) in enumerate(subintervals):

#             # --- Chebyshev fit ---
#             az_Ch = AffineZonotope(branin, lb, ub)
#             az_Ch.fit_chebyshev(tol=1e-8, solver='SLSQP', verbose=False)
#             az_Ch.build_zonotope()
#             dum1 = az_Ch.f_range

#             lo_Ch = min(lo_Ch, dum1[0])
#             hi_Ch = max(hi_Ch, dum1[1])
#             total_f_evals_Ch += az_Ch.sip_stats['total_f_evals']

#             # --- Min-range fit (fresh object to avoid state contamination) ---
#             az_MR = AffineZonotope(branin, lb, ub)
#             az_MR.fit_min_range(solver='OSQP', verbose=False)
#             az_MR.build_zonotope()
#             dum2 = az_MR.f_range

#             lo_MR = min(lo_MR, dum2[0])
#             hi_MR = max(hi_MR, dum2[1])
#             total_f_evals_MR += az_MR.mr_stats['total_f_evals']

#         # Store result for this (m, n) configuration
#         results[(m+1, n+1)] = {
#             'splits_per_dim': splits_per_dim,
#             'f_rangeCh' : [lo_Ch, hi_Ch],
#             'f_rangeMR' : [lo_MR, hi_MR],
#         }

    

# # %%
# import pickle

# # Save after loop completes
# with open('branin_results_nm.pkl', 'wb') as f:
#     pickle.dump(results, f)
# print(f"\nSaved {len(results)} results to branin_results_nm.pkl")

# # %%

# # --- Load results ---
# with open('branin_results_nm.pkl', 'rb') as f:
#     results = pickle.load(f)

# # --- Extract grid dimensions ---
# keys   = list(results.keys())
# dim1   = max(k[0] for k in keys)   
# dim2   = max(k[1] for k in keys)   

# # --- Build 2D arrays for lo_Ch and lo_MR ---
# lo_Ch_grid = np.zeros((dim1, dim2))
# lo_MR_grid = np.zeros((dim1, dim2))
# hi_Ch_grid = np.zeros((dim1, dim2))
# hi_MR_grid = np.zeros((dim1, dim2))

# for (m, n), val in results.items():
#     lo_Ch_grid[m-1, n-1] = val['f_rangeCh'][0]   # lo
#     hi_Ch_grid[m-1, n-1] = val['f_rangeCh'][1]   # hi
#     lo_MR_grid[m-1, n-1] = val['f_rangeMR'][0]   # lo
#     hi_MR_grid[m-1, n-1] = val['f_rangeMR'][1]   # hi

# tick_labels = [str(i+1) for i in range(dim1)]

# # %%
# import matplotlib.pyplot as plt

# # --- Plot ---
# fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# # Chebyshev lower bound
# im1 = axes[0].imshow(lo_Ch_grid, origin='lower', aspect='auto', cmap='viridis')
# axes[0].set_title('Lower Bound — Chebyshev ($lo_{Ch}$)', fontsize=13)
# axes[0].set_xlabel('Splits dim 2 (n)')
# axes[0].set_ylabel('Splits dim 1 (m)')
# axes[0].set_xticks(range(dim2)); axes[0].set_xticklabels(tick_labels)
# axes[0].set_yticks(range(dim1)); axes[0].set_yticklabels(tick_labels)
# plt.colorbar(im1, ax=axes[0])

# # Min-range lower bound
# im2 = axes[1].imshow(lo_MR_grid, origin='lower', aspect='auto', cmap='plasma')
# axes[1].set_title('Lower Bound — Min-Range ($lo_{MR}$)', fontsize=13)
# axes[1].set_xlabel('Splits dim 2 (n)')
# axes[1].set_ylabel('Splits dim 1 (m)')
# axes[1].set_xticks(range(dim2)); axes[1].set_xticklabels(tick_labels)
# axes[1].set_yticks(range(dim1)); axes[1].set_yticklabels(tick_labels)
# plt.colorbar(im2, ax=axes[1])

# plt.suptitle('Branin Function — Lower Bound vs Split Configuration', fontsize=14)
# plt.tight_layout()
# plt.savefig('branin_lower_bounds.png', dpi=150, bbox_inches='tight')
# plt.show()

# # %%
