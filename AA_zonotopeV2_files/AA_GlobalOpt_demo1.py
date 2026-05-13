# %%
"""
demo_AA_zonotope.py
===================
Demonstrates the AffineZonotope class on two benchmark functions:

  1. Branin function   — a standard nonlinear test function on [0,1]²
  2. f(x) = x₁² − x₂² — a saddle function on [−1,1]²

For each function both approximation methods are tested:
  - Chebyshev minimax fit  (fit_chebyshev)
  - Minimum-range LP fit   (fit_min_range)

Workflow for each case:
  AffineZonotope(f, l, u)  →  fit_*()  →  build_zonotope()  →  plot / inspect
"""

import numpy as np

from Affine_ArithmeticClassV2 import AffineZonotope

# ===========================================================================
# Helper: print zonotope results
# ===========================================================================

def print_results(label, az):
    """Print the key outputs of a fitted and built AffineZonotope."""
    print(f"\n{'─' * 55}")
    print(f"  {label}")
    print(f"{'─' * 55}")
    print(f"  Output range : {az.f_range}")
    print(f"  Center       : {az.center}")
    print(f"  G matrix     :\n{az.G}")
    print(f"  delta         : {az.t_opt}")
    print(f"  alpha        : {az.a_opt}")




# ===========================================================================
# Example 1 — Branin function
#
#   f(x₁, x₂) = a(x₂ − b·x₁² + c·x₁ − r)² + s(1−t)cos(x₁) + s
#   Input domain: x ∈ [0, 1]²  (scaled internally to the standard Branin domain)
#   Known global minima ≈ 0.397 at three locations
# ===========================================================================

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


l_branin = [0.0, 0.0]
u_branin = [1.0, 1.0]

# ---------------------------------------------------------------------------
# 1a. Branin — Chebyshev minimax fit
#     Finds the affine h(x) = aᵀx + b minimising the worst-case error
#     max_{x ∈ [l,u]} |f(x) − h(x)|  via a SIP cutting-plane loop.
#     Outer solver: SLSQP (well-suited for linear objective + nonlinear constraints)
# ---------------------------------------------------------------------------
az_branin = AffineZonotope(branin, l_branin, u_branin)

az_branin.fit_chebyshev(
    solver='SLSQP',
    verbose=True,        # set True to watch SIP iterations
)
az_branin.build_zonotope()

print_results("Branin — Chebyshev minimax fit", az_branin)

az_branin.plot_zonotope_3d(
    var_names=["x_1", "x_2"],
    f_name="Branin function",
)

# ---------------------------------------------------------------------------
# 1b. Branin — Minimum-range LP fit
#     Finds alpha* minimising Σ|x1_i · alpha_i| via LP (CVXPY / OSQP),
#     then computes the tightest gamma ± delta enclosure of f(x) − alpha*ᵀx.
#     OSQP is fast and accurate for small LP dimensions.
# ---------------------------------------------------------------------------
az_branin.fit_min_range(
    solver='OSQP',
    verbose=True,
)
az_branin.build_zonotope()

print_results("Branin — Minimum-range LP fit", az_branin)

az_branin.plot_zonotope_3d(
    var_names=["x_1", "x_2"],
    f_name="Branin function",
)


# ===========================================================================
# Example 2 — Saddle function  f(x) = x₁² − x₂²
#
#   Symmetric saddle centred at the origin.
#   Gradient is zero at (0,0) — a challenging case for LP gradient bounds.
#   Input domain: x ∈ [−1, 1]²
# ===========================================================================

def f_saddle(x):
    """Saddle function — accepts 1-D ndarray x of shape (2,)."""
    return x[0]**2 - x[1]**2


l_saddle = [-1.0, -1.0]
u_saddle = [ 1.0,  1.0]

# ---------------------------------------------------------------------------
# 2a. Saddle — Chebyshev minimax fit
# ---------------------------------------------------------------------------
az_saddle = AffineZonotope(f_saddle, l_saddle, u_saddle)

az_saddle.fit_chebyshev(
    #solver='SLSQP',
    verbose=True,
)
az_saddle.build_zonotope()

# Only change what you need — rest stay at _DE_DEFAULTS_CHEB
az_saddle.fit_chebyshev(
    tol=1e-12,
    max_iter=500,
    verbose=True,
    de_kwargs={'seed': 7, 'maxiter': 2000}
)

print_results("Saddle f = x₁²−x₂² — Chebyshev minimax fit", az_saddle)

az_saddle.plot_zonotope_3d(
    var_names=["x_1", "x_2"],
    f_name="x₁² − x₂²",
)

# ---------------------------------------------------------------------------
# 2b. Saddle — Minimum-range QP/LP fit
# ---------------------------------------------------------------------------
az_saddle.fit_min_range(
    solver = 'CLARABEL',
    verbose=True,
    de_kwargs={'seed': 7, 'maxiter': 2000}
)
az_saddle.build_zonotope()

print_results("Saddle f = x₁²−x₂² — Minimum-range using QP/LP fit", az_saddle)

az_saddle.plot_zonotope_3d(
    var_names=["x_1", "x_2"],
    f_name="x₁² − x₂²",
)


# Only change seed and popsize — rest stay at _DE_DEFAULTS_MR
az2 = AffineZonotope(f_saddle, l_saddle, u_saddle)
az2.fit_min_range(
    grad_eps=1e-6,
    de_kwargs={'seed': 19, 'popsize': 30}
)

# %%
repr(az2)
print(az2)

az2.build_zonotope()

print_results("Saddle f = x₁²−x₂² — Minimum-range using QP/LP fit", az2)
# %%
