"""
AA_zonotope.py
==============
Affine Arithmetic zonotope builder.
 
Workflow
--------
1.  Instantiate with a target function and its input bounds.
2.  Choose an approximation method:
      - fit_chebyshev()  →  Chebyshev / SIP minimax affine fit
      - fit_min_range()  →  Minimum-range affine fit
3.  Build the zonotope:
      - build_zonotope()
 
All intermediate and final results are stored as instance attributes
(see attribute table in __init__).
"""
 
import itertools
import math
 
import cvxpy as cp
import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
from scipy.optimize import approx_fprime, differential_evolution, minimize
 
sig_dig = 12  # global significant digits for outward rounding

# ---------------------------------------------------------------------------
# Default differential-evolution settings shared by both fit methods
# ---------------------------------------------------------------------------
_DE_DEFAULTS_CHEB = dict(
    seed=42,
    tol=1e-8,
    atol=1e-8,
    maxiter=1000,
    popsize=15,
    mutation=(0.5, 1.0),
    recombination=0.7,
)

_DE_DEFAULTS_MR = dict(
    seed=42,
    tol=1e-6,
    atol=1e-6,
    maxiter=1000,
    popsize=15,
    mutation=(0.5, 1.0),
    recombination=0.7,
)


class AffineZonotope:
    """
    Constructs an affine-arithmetic zonotope for a nonlinear scalar function
    over a box domain, using either a Chebyshev minimax or a minimum-range
    LP affine approximation.

    Parameters
    ----------
    f : callable
        Target function.  Must accept a 1-D ndarray x of shape (n,) and
        return a scalar float.
    l : array-like, shape (n,)
        Lower bounds of the input box.
    u : array-like, shape (n,)
        Upper bounds of the input box.

    Attributes set after __init__
    -----------------------------
    f, l, u, n   — inputs as stored
    x0            — box midpoint,   shape (n,)
    x1            — box half-widths, shape (n,)
    box           — list of (lo, hi) tuples for scipy optimisers

    Attributes set after fit_chebyshev / fit_min_range
    ---------------------------------------------------
    a_opt         — slope coefficients,  shape (n,)
    b_opt         — intercept (scalar)
    t_opt         — minimax error / delta (scalar)

    Attributes set after build_zonotope
    ------------------------------------
    center        — zonotope center,     shape (n+1,)
    G             — generator matrix,    shape (n+1, n+1)
    vertices      — all 2^(n+1) vertices, shape (2^(n+1), n+1)
    f_range       — [f_min, f_max] output enclosure
    """

    # -----------------------------------------------------------------------
    # Construction
    # -----------------------------------------------------------------------

    def __init__(self, f_target: callable, l, u):
        self.f  = f_target
        self.l  = np.asarray(l, dtype=float)
        self.u  = np.asarray(u, dtype=float)
        self.n  = len(self.l)

        # Derived geometry
        self.x0  = 0.5 * (self.l + self.u)      # midpoint
        self.x1  = 0.5 * (self.u - self.l)      # half-widths
        self.box = list(zip(self.l, self.u))

        # Placeholders — populated by fit_* and build_zonotope
        self.a_opt    = None
        self.b_opt    = None
        self.t_opt    = None
        self.center   = None
        self.G        = None
        self.vertices = None
        self.f_range  = None

    # -----------------------------------------------------------------------
    # Private helpers
    # -----------------------------------------------------------------------

    def _f_vec(self, x: np.ndarray) -> float:
        """Scalar wrapper: accepts a vector, returns float."""
        return float(self.f(x))

    def _merge_de_kwargs(self, defaults: dict, overrides: dict | None) -> dict:
        """Merge user-supplied de_kwargs on top of method-specific defaults."""
        merged = defaults.copy()
        if overrides:
            merged.update(overrides)
        return merged
    

    def _oB(self, x, left):
        """Round x outward: left bound rounds down, right bound rounds up."""
        str_rep   = f'{x:0.{sig_dig}g}'
        x_print   = float(str_rep)
 
        prec      = math.floor(math.log10(abs(x) + 1e-100) + 1.0) - sig_dig
        least_sig = 10.0 ** prec * 0.5
 
        if left:
            if x_print > x and not (abs(x_print - x) < 1e-100):
                x -= least_sig
        else:
            if x_print < x and not (abs(x_print - x) < 1e-100):
                x += least_sig
 
        return round(x, sig_dig)
 
    def _outer_bound(self, lo, hi):
        """Apply outward rounding to an [lo, hi] pair."""
        return self._oB(lo, left=True), self._oB(hi, left=False)
    
    # -----------------------------------------------------------------------
    # Fit method 1 — Chebyshev / SIP minimax
    # -----------------------------------------------------------------------

    def fit_chebyshev(
        self,
        tol: float       = 1e-6,
        max_iter: int    = 1000,
        de_kwargs: dict  = None,
        solver: str      = 'trust-constr',
        verbose: bool    = True,
    ):
        """
        Solve  min_{a,b,t} t  s.t.  |f(x) − aᵀx − b| ≤ t  ∀ x ∈ [l, u]
        via a Sequential / Cutting-plane (SIP) outer loop with a
        Differential Evolution inner oracle.

        Parameters
        ----------
        tol      : convergence tolerance (inner worst-case ≤ t + tol)
        max_iter : maximum SIP iterations
        de_kwargs: dict — overrides for differential_evolution (inner oracle)
        verbose  : print iteration log
        solver   : algorithm implemented in outer loop of SIP
                   default: 'trust-constr', or user can put 'SLSQP'

        Returns
        -------
        self  (for method chaining)
        """
        de_kw = self._merge_de_kwargs(_DE_DEFAULTS_CHEB, de_kwargs)
        n     = self.n

        # -- Inner oracle: maximise |f(x) − aᵀx − b| over the box ----------
        def _inner_solve(a, b):
            def neg_abs_res(x):
                return -np.abs(self._f_vec(x) - np.dot(a, x) - b)

            res = differential_evolution(neg_abs_res, self.box, **de_kw)
            return res.x, -res.fun   # (x_worst, worst_residual)

        # -- SIP outer loop --------------------------------------------------
        active_pts = [self.x0.copy()]

        for it in range(max_iter):

            # Outer LP/NLP: minimise t subject to current active constraints
            def objective(z):
                return z[-1]

            def _make_con(xk):
                def con(z):
                    a, b, t = z[:n], z[n], z[n + 1]
                    return t - np.abs(self._f_vec(xk) - np.dot(a, xk) - b)
                return {'type': 'ineq', 'fun': con}

            constraints = [_make_con(xk) for xk in active_pts]
            z0          = np.zeros(n + 2)

            res_outer = minimize(
                objective, z0,
                method=solver,
                #method='SLSQP',
                constraints=constraints,
            )

            a_opt = res_outer.x[:n]
            b_opt = res_outer.x[n]
            t_opt = res_outer.x[n + 1]

            # Inner: find global worst-case point
            x_worst, val_worst = _inner_solve(a_opt, b_opt)

            if verbose:
                print(
                    f"[Chebyshev] iter {it + 1:3d}: "
                    f"t = {t_opt:.8f}, "
                    f"worst_res = {val_worst:.8f}, "
                    f"x_worst = {x_worst}"
                )

            # Convergence
            if val_worst <= t_opt + tol:
                if verbose:
                    print("[Chebyshev] Converged — global optimum certified.\n")
                break

            active_pts.append(x_worst)

        _, t_cheb = self._outer_bound(0.0, float(t_opt))

        # Store results
        self.a_opt = a_opt
        self.b_opt = b_opt
        self.t_opt = t_cheb

        return self

    # -----------------------------------------------------------------------
    # Fit method 2 — Minimum-range 
    # -----------------------------------------------------------------------

    def fit_min_range(
        self,
        solver: str     = None,
        grad_eps: float = 1e-8,
        de_kwargs: dict = None,
    ):
        """
        Minimum-range affine approximation via LP.

        Finds alpha* minimising  Σ |x1_i · alpha_i|  subject to gradient
        bounds, then computes the tightest enclosure of  g(x) = f(x) − alpha*ᵀx
        globally over the box.

        Parameters
        ----------
        solver   : CVXPY solver string (None → CVXPY default)
        grad_eps : finite-difference step for numerical gradient
        de_kwargs: dict — overrides for differential_evolution (global search)

        Returns
        -------
        self  (for method chaining)

        Extra attributes set
        --------------------
        c_lp     : ndarray shape (n+1,)  — [gamma, alpha*]
        delta_lp : float                 — half-width of g-residual enclosure
        """
        de_kw = self._merge_de_kwargs(_DE_DEFAULTS_MR, de_kwargs)
        x0, x1, box = self.x0, self.x1, self.box
        dim          = self.n

        # -- Step 1: gradient bounds at corners ------------------------------
        g_low  = approx_fprime(self.l, self._f_vec, grad_eps)
        g_high = approx_fprime(self.u, self._f_vec, grad_eps)

        d_min = np.minimum(g_low, g_high)
        d_max = np.maximum(g_low, g_high)

        # -- Step 2: LP — minimise Σ |x1_i · alpha_i| -----------------------
        alpha = cp.Variable(dim)
        z     = cp.Variable(dim)

        lp_constraints = [
            z >=  cp.multiply(x1, alpha),
            z >= -cp.multiply(x1, alpha),
            z >= 0,
            alpha >= d_min,
            alpha <= d_max,
        ]

        prob = cp.Problem(cp.Minimize(cp.sum(z)), lp_constraints)
        prob.solve(solver=solver)

        alpha_star = alpha.value   # shape (dim,)

        # -- Step 3: global min/max of g(x) = f(x) − alpha*ᵀx --------------
        def g_scalar(x):
            return self._f_vec(x) - float(alpha_star @ x)

        res_min = differential_evolution(g_scalar,           box, **de_kw)
        res_max = differential_evolution(lambda x: -g_scalar(x), box, **de_kw)

        g_global_min = res_min.fun
        g_global_max = -res_max.fun

        # Cross-check at all 2^n corners of the box
        corners_eps = np.array(list(itertools.product([-1, 1], repeat=dim)))
        X_corners   = np.array([x0 + np.diag(x1) @ e for e in corners_eps])
        g_corners   = np.array([g_scalar(xc) for xc in X_corners])

        g_lo = min(g_global_min, g_corners.min())
        g_hi = max(g_global_max, g_corners.max())

        # -- Step 4: gamma and delta -----------------------------------------
        gamma    = 0.5 * (g_hi + g_lo)
        _, delta = self._outer_bound(0.0, 0.5 * (g_hi - g_lo))  # round delta UP
 
        # Store MR-specific results
        self.c_mr     = np.hstack((gamma, alpha_star))
        self.delta_mr = float(delta)
 
        # Map to common interface used by build_zonotope
        self.a_opt = alpha_star
        self.b_opt = gamma
        self.t_opt = float(delta)
 
        return self

    # -----------------------------------------------------------------------
    # Zonotope construction
    # -----------------------------------------------------------------------

    def build_zonotope(self):
        """
        Construct the zonotope from the stored affine approximation.
        Outward rounding is always applied via self._outer_bound.
 
        Must be called after fit_chebyshev() or fit_min_range().
 
        Returns
        -------
        self  (for method chaining)
 
        Attributes set
        --------------
        center   : ndarray, shape (n+1,)
        G        : ndarray, shape (n+1, n+1)
        vertices : ndarray, shape (2^(n+1), n+1)
        f_range  : [f_min, f_max]
        """

        if self.a_opt is None:
            raise RuntimeError(
                "No affine fit found. Call fit_chebyshev() or fit_min_range() first."
            )
        
        x0, x1 = self.x0, self.x1
        dim     = self.n
        c0      = self.b_opt
        c_lin   = self.a_opt   # shape (n,)
 
        # Zonotope center: [x0, f_affine(x0)]
        y0          = c0 + c_lin @ x0
        self.center = np.hstack([x0, y0])   # shape (n+1,)
 
        # Generator matrix G, shape (n+1, n+1)
        G = np.zeros((dim + 1, dim + 1))
 
        for i in range(dim):
            G[i, i] = x1[i]               # input generator
 
            val = c_lin[i] * x1[i]        # output component
            if val >= 0:
                _, G[-1, i] = self._outer_bound(0.0, val)   # positive: round UP
            else:
                G[-1, i], _ = self._outer_bound(val, 0.0)   # negative: round DOWN
 
        G[-1, -1] = self.t_opt            # error generator
 
        self.G = G

        # All 2^(n+1) vertices
        eps           = np.array(list(itertools.product([-1, 1], repeat=dim + 1)))
        self.vertices = np.array([self.center + G @ e for e in eps])

        # Output range with outward rounding
        f_lo, f_hi   = self._outer_bound(
            float(self.vertices[:, -1].min()),
            float(self.vertices[:, -1].max()),
        )
        self.f_range = [f_lo, f_hi]
 
        return self

    # -----------------------------------------------------------------------
    # Plotting
    # -----------------------------------------------------------------------
 
    def plot_zonotope_3d(
        self,
        var_names=None,
        f_name=None,
        grid_points=40,
        surface_alpha=0.75,
        poly_alpha=0.15,
    ):
        """
        Plot a 3D zonotope/polytope with the true function surface.
        Uses self.vertices — call build_zonotope() first.
 
        Parameters
        ----------
        var_names     : list[str]  — input axis labels (default ["x_1", "x_2"])
        f_name        : str        — output axis label (default "f(x_1, x_2)")
        grid_points   : int        — surface resolution (default 40)
        surface_alpha : float      — surface transparency (default 0.4)
        poly_alpha    : float      — polytope transparency (default 0.15)
        """
        if self.vertices is None:
            raise RuntimeError("Call build_zonotope() before plotting.")
        if self.vertices.shape[1] != 3:
            raise ValueError(
                f"3D plot requires n=2. Got vertices shape {self.vertices.shape}."
            )
 
        if var_names is None:
            var_names = ["x_1", "x_2"]
        if f_name is None:
            f_name = "f(x_1, x_2)"
 
        vertices = self.vertices
 
        # -- Surface grid over the real x/y domain ---------------------------
        x_vals = np.linspace(vertices[:, 0].min(), vertices[:, 0].max(), grid_points)
        y_vals = np.linspace(vertices[:, 1].min(), vertices[:, 1].max(), grid_points)
        X, Y   = np.meshgrid(x_vals, y_vals)
 
        # Evaluate f point-by-point using the stored _f_vec wrapper
        f_vect = np.vectorize(lambda xi, yi: self._f_vec(np.array([xi, yi])))
        Z      = f_vect(X, Y)
 
        # -- Convex hull of zonotope vertices ---------------------------------
        hull  = ConvexHull(vertices)
        faces = [vertices[simplex] for simplex in hull.simplices]
 
        # -- Plot -------------------------------------------------------------
        fig = plt.figure()
        ax  = fig.add_subplot(111, projection="3d")
 
        ax.add_collection3d(
            Poly3DCollection(faces, alpha=poly_alpha, edgecolor="black")
        )
 
        ax.scatter(
            vertices[:, 0], vertices[:, 1], vertices[:, 2],
            color="red", s=20, label="Vertices"
        )
 
        ax.plot_surface(X, Y, Z, cmap="viridis", alpha=surface_alpha)
 
        mins, maxs = vertices.min(axis=0), vertices.max(axis=0)
        ax.set_xlim(mins[0], maxs[0])
        ax.set_ylim(mins[1], maxs[1])
        ax.set_zlim(mins[2], maxs[2])
 
        ax.set_xlabel(var_names[0])
        ax.set_ylabel(var_names[1])
        ax.set_zlabel(f_name)
 
        plt.tight_layout()
        plt.show()
 
    def plot_zonotope_2d(
        self,
        xlabel="x",
        ylabel="f(x)",
        title="2D Zonotope",
        grid_points=200,
        fill_alpha=0.2,
    ):
        """
        Plot a 2D zonotope/polytope with the true function and affine approximation.
        Uses self.vertices, self.b_opt, self.a_opt — call build_zonotope() first.
 
        Parameters
        ----------
        xlabel     : str   — x-axis label (default "x")
        ylabel     : str   — y-axis label (default "f(x)")
        title      : str   — plot title
        grid_points: int   — curve resolution (default 200)
        fill_alpha : float — zonotope fill transparency (default 0.2)
        """
        if self.vertices is None:
            raise RuntimeError("Call build_zonotope() before plotting.")
        if self.vertices.shape[1] != 2:
            raise ValueError(
                f"2D plot requires n=1. Got vertices shape {self.vertices.shape}."
            )
 
        vertices = np.asarray(self.vertices)
        cvec     = np.array([self.b_opt, self.a_opt[0]])  # [intercept, slope]
 
        # -- Order vertices counter-clockwise for a clean polygon ------------
        centroid = vertices.mean(axis=0)
        angles   = np.arctan2(
            vertices[:, 1] - centroid[1],
            vertices[:, 0] - centroid[0],
        )
        order          = np.argsort(angles)
        vertices_ord   = vertices[order]
        vertices_closed = np.vstack([vertices_ord, vertices_ord[0]])   # close polygon
 
        # -- Evaluate true function over x domain ----------------------------
        x_vals = np.linspace(vertices[:, 0].min(), vertices[:, 0].max(), grid_points)
        f_vect = np.vectorize(lambda xi: self._f_vec(np.array([xi])))
        y_true = f_vect(x_vals)
 
        # -- Affine approximation --------------------------------------------
        y_aff = cvec[0] + cvec[1] * x_vals
 
        # -- Plot ------------------------------------------------------------
        plt.figure(figsize=(6, 6))
 
        plt.plot(x_vals, y_true, "r-",  label="True function")
        plt.plot(x_vals, y_aff,  "b--", label="Affine approx")
 
        plt.plot(
            vertices_closed[:, 0], vertices_closed[:, 1],
            color="black", linewidth=2, label="Zonotope boundary"
        )
 
        plt.scatter(vertices_ord[:, 0], vertices_ord[:, 1], color="gray", zorder=3)
 
        plt.fill(vertices_ord[:, 0], vertices_ord[:, 1], alpha=fill_alpha)
 
        plt.scatter(
            centroid[0], centroid[1],
            color="black", marker="x", s=100, label="Center"
        )
 
        plt.axis("equal")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
 
    # -----------------------------------------------------------------------
    # Convenience
    # -----------------------------------------------------------------------
 
    def __repr__(self):
        fitted = self.a_opt is not None
        zono   = self.center is not None
        return (
            f"AffineZonotope(n={self.n}, fitted={fitted}, zonotope_built={zono})"
        )