import numpy as np
import sympy as sp
import cvxpy as cp
import itertools
from scipy.optimize import root
from scipy.stats import qmc

def IAVar_AAVar(var1,var2):
    x0 = (var2[0] + var2[1])/2
    x1 = (var2[1]-var2[0])/2
    AAvar = x0 + x1*var1
    return AAvar

def chebyshev_fit_Nd(
    f_target,      # target function
    vars,          # independent variables
    V_base,        # Interval varables
    n_starts=None, # number of initial guess points in N-D latin hypercube
    M=None,        # number of points chosen to calculate f-(c_0 + c_i*x_i)
    bounds=None    # bounds for noise parameters, epscilons
    ):
    """
    Chebyshev minimax fit for a N-D symbolic function.


    Parameters
    ----------
    f_target : sympy.Expr
        Scalar function f(x)
    vars : list[sympy.Symbol]
        Variables [x1, ..., xN]
    V_base : ndarray, shape (N, 2)
        Lower/upper bounds per dimension
    n_starts : int
        number of initial guess points in N-D latin hypercube
        to calculate stationary points
    M : int
        Number of LHS samples (without addition of corners)


    Returns
    -------
    c_val : ndarray, shape (N+1,)
        [alpha*, gamma]
    E_val : float
        Residual half-width
    max_err : float
        Residual half-width = h_infinity norm
    """

    # --------------------------------------------------
    # 1. Symbolic gradient
    # --------------------------------------------------
    # grad = [
    #     sp.diff(f_target, x[3]),
    #     sp.diff(f_target, x[5])
    # ]

    grad = [sp.diff(f_target, v) for v in vars]

    grad_f = sp.lambdify(vars, grad, modules="numpy")
    
    def grad_f_vec(z):
        return np.asarray(grad_f(*z), dtype=float).ravel()
    

    # --------------------------------------------------
    # 2. Multi-start root finding
    # --------------------------------------------------
    xmin, xmax = bounds
    dim = len(vars)

    initial_guesses = np.random.uniform(
    low=xmin,
    high=xmax,
    size=(n_starts, dim)
    )
    
    solutions = []
    
    for x0 in initial_guesses:
        sol = root(grad_f_vec, x0, method="hybr")

        if sol.success:
            solutions.append(sol.x)

    solutions = np.array(solutions)

    if solutions.size == 0:
        solutions = solutions.reshape(0, dim)

    # --------------------------------------------------
    # 3. Boundary vertices
    # --------------------------------------------------
    eps = np.array(list(itertools.product([-1, 1], repeat=dim)))

    all_points = np.vstack((eps, solutions))
    all_points = np.unique(np.round(all_points, 18), axis=0)

    # --------------------------------------------------
    # 4. Latin Hypercube samples
    # --------------------------------------------------
    
    sampler = qmc.LatinHypercube(d=dim)
    samples = sampler.random(n=M)

    lower = np.ones(dim) * xmin
    upper = np.ones(dim) * xmax

    samples_scaled = qmc.scale(samples, lower, upper)

    all_sample = np.vstack((all_points, samples_scaled))

    # V_base = np.asarray(V_base, dtype=float)

    # x0 = 0.5 * (V_base[:, 0] + V_base[:, 1])
    # x1 = 0.5 * (V_base[:, 1] - V_base[:, 0])

    # samples = np.array([
    #     x0 + np.diag(x1) @ e for e in sample_set
    # ])

    # all_sample = np.unique(np.round(sample_set, 18), axis=0)
    # --------------------------------------------------
    # 5. Target function evaluation
    # --------------------------------------------------
    f_target_f = sp.lambdify(vars, f_target, modules="numpy")
    y = np.asarray(f_target_f(*all_sample.T), dtype=float).ravel()

    # --------------------------------------------------
    # 6. Build V matrix
    # --------------------------------------------------

    V = []

    for ind in all_sample:
        row = [1.0]  # constant term
        for j, var in enumerate(V_base):
            xi = IAVar_AAVar(ind[j], var)
            row.append(xi)
        
        V.append(row)

    V = np.array(V, dtype=float)

    # --------------------------------------------------
    # 7. Chebyshev minimax optimization
    # --------------------------------------------------
    m, p = V.shape

    c = cp.Variable(p)
    E = cp.Variable(nonneg=True)

    constraints = [
        V @ c - y <= E,
        V @ c - y >= -E
    ]

    problem = cp.Problem(cp.Minimize(E), constraints)
    problem.solve()

    c_val = c.value
    E_val = E.value
    max_err = np.max(np.abs(V @ c_val - y))

    return c_val, E_val, max_err

def affine_min_range_lp(
    f_target,      # target function
    vars,          # independent variables
    V_base,        # Interval varables
    M=None,        # number of points chosen to calculate f-(c_0 + c_i*x_i)
    solver=None,   # solver chosen for solving LP problem 
    bounds = None  # bounds for noise parameters, epscilons
    ):
    """
    Compute minimum-range affine approximation via LP (N-D, fully general).

    Parameters
    ----------
    f_target : sympy.Expr
        Scalar function f(x)
    vars : list[sympy.Symbol]
        Variables [x1, ..., xN]
    V_base : ndarray, shape (N, 2)
        Lower/upper bounds per dimension
    M : int
        Number of LHS samples (without addition of corners)
    solver : cvxpy solver
        LP solver (HIGHS recommended)

    Returns
    -------
    cvec : ndarray, shape (N+1,)
        [alpha*, gamma]
    delta : float
        Residual half-width
    """

    # -----------------------
    # Basic geometry
    # -----------------------
    xmin, xmax = bounds

    V_base = np.asarray(V_base, dtype=float)
    dim = V_base.shape[0]

    x0 = 0.5 * (V_base[:, 0] + V_base[:, 1])
    x1 = 0.5 * (V_base[:, 1] - V_base[:, 0])

    # -----------------------
    # Gradient bounds
    # -----------------------
    grad = [sp.diff(f_target, v) for v in vars]
    grad_f = sp.lambdify(vars, grad, modules="numpy")

    g_low  = np.asarray(grad_f(*(xmin * np.ones(dim))), dtype=float)
    g_high = np.asarray(grad_f(*(xmax * np.ones(dim))), dtype=float)

    d_min = np.minimum(g_low, g_high)
    d_max = np.maximum(g_low, g_high)

    # -----------------------
    # LP: minimize generator contribution
    # -----------------------
    alpha = cp.Variable(dim)
    z = cp.Variable(dim)

    objective = cp.Minimize(cp.sum(z))

    constraints = [
        z >=  cp.multiply(x1, alpha),
        z >= -cp.multiply(x1, alpha),
        z >= 0,
        alpha >= d_min,
        alpha <= d_max
    ]

    prob = cp.Problem(objective, constraints)
    prob.solve(solver=solver)

    alpha_star = alpha.value

    # -----------------------
    # Sample domain (corners + LHS)
    # -----------------------

    sampler = qmc.LatinHypercube(d=dim)
    samples = sampler.random(n=M)

    lower = np.ones(dim) * xmin
    upper = np.ones(dim) * xmax

    samples_scaled = qmc.scale(samples, lower, upper)

    eps_corners = np.array(list(itertools.product([-1, 1], repeat=dim)))

    eps_all = np.vstack((samples_scaled, eps_corners))

    samples = np.array([
        x0 + np.diag(x1) @ e for e in eps_all
    ])

    # -----------------------
    # Residual computation
    # -----------------------
    f_fun = sp.lambdify(vars, f_target, modules="numpy")

    f_vals = f_fun(*eps_all.T)


    g_vals = f_vals - samples @ alpha_star

    gamma = 0.5 * (g_vals.max() + g_vals.min())
    delta = 0.5 * (g_vals.max() - g_vals.min())

    # -----------------------
    # Outputs
    # -----------------------
    cvec = np.hstack((gamma, alpha_star))

    return cvec, delta


def zonotope_vertices(solt, V_base):
    """
    Construct an (x, y) zonotope from Chebyshev coefficients.

    Parameters
    ----------
    solt : tuple or list
        solt[0] = coefficient vector c (length dim+1)
        solt[1] = Chebyshev error bound E
    V_base : array-like, shape (dim, 2)
        Lower/upper bounds for each input variable

    Returns
    -------
    cvec : ndarray, shape (dim+1,)
        Zonotope center
    Gmat : ndarray, shape (dim+1, dim+1)
        Generator matrix
    vertices : ndarray, shape (2^(dim+1), dim+1)
        Zonotope vertices
    """

    # ---- extract solution ----
    c_vals = np.asarray(solt[0], dtype=float)
    E0 = float(solt[1])

    dim = len(c_vals) - 1

    # ---- domain center and domain radius ----
    V_base = np.asarray(V_base, dtype=float)

    x0 = (V_base[:, 0] + V_base[:, 1]) / 2
    x1 = (V_base[:, 1] - V_base[:, 0]) / 2

    # ---- coefficients ----
    c0 = c_vals[0]
    c_lin = c_vals[1:]

    # ---- center vector ----
    y0 = c0 + c_lin @ x0
    cvec = np.hstack([x0, y0])

    # ---- generator matrix ----
    Gmat = np.zeros((dim + 1, dim + 1))

    for i in range(dim):
        Gmat[i, i] = x1[i]
        Gmat[-1, i] = c_lin[i] * x1[i]

    Gmat[-1, -1] = E0

    # ---- zonotope vertices ----
    eps = np.array(list(itertools.product([-1, 1], repeat=dim + 1)))
    vertices = np.array([cvec + Gmat @ e for e in eps])

    return cvec, Gmat, vertices