import numpy as np
from itertools import product
from scipy.fft import dct
from numpy.polynomial.chebyshev import chebval
from numpy.polynomial.chebyshev import chebvander


# --------------------------------------------------
# Chebyshev-Lobatto nodes
# --------------------------------------------------

def chebyshev_lobatto_nodes(n, a, b):

    k = np.arange(n + 1)

    xi = np.cos(np.pi * k / n)

    x = 0.5 * (a + b) + 0.5 * (b - a) * xi

    return x

# --------------------------------------------------
# N-dimensional tensor-product Chebyshev grid
# --------------------------------------------------


def chebyshev_grid_nd(nvec, combo):

    d = len(nvec)

    if len(combo) != d:
        raise ValueError(
            f"Dimension mismatch: nvec has {d} entries "
            f"but combo has {len(combo)} intervals."
        )

    # Nodes for each dimension
    nodes = []

    for k in range(d):

        a, b = combo[k]

        nk = nvec[k]

        nodes_k = chebyshev_lobatto_nodes(
            nk,
            a,
            b
        )

        nodes.append(nodes_k)

    # Tensor-product grid
    points = np.array(
        list(product(*nodes))
    )

    return nodes, points

# --------------------------------------------------
# Evaluate function at all grid points
# --------------------------------------------------

def evaluate_function_on_grid(func, points):

    return np.array([
        func(point)
        for point in points
    ])

def reshape_grid_values(F, nvec):

    shape = tuple(
        n + 1
        for n in nvec
    )

    return F.reshape(shape)

# --------------------------------------------------
# Calculate N-D Chebyshev coefficients
# --------------------------------------------------

def chebyshev_coefficients(F):

    d = F.ndim

    C = F.copy()

    # Degree in each dimension
    degrees = [
        F.shape[axis] - 1
        for axis in range(d)
    ]

    # DCT-I in every dimension
    for axis in range(d):
        C = dct(C, type=1, axis=axis)

    # Normalization and endpoint weights
    for axis in range(d):

        n = degrees[axis]

        # DCT-I normalization
        C /= n

        # Endpoint weights
        w = np.ones(n + 1)
        w[0] = 0.5
        w[-1] = 0.5

        # Reshape so weights act only along current axis
        shape = [1] * d
        shape[axis] = n + 1

        C *= w.reshape(shape)

    return C


def evaluate_chebyshev_nd(C, x, combo):

    x = np.asarray(x, dtype=float)

    d = len(combo)

    if len(x) != d:
        raise ValueError(
            "Dimension of x and combo do not match."
        )

    # Map each physical coordinate [a_k,b_k]
    # to Chebyshev coordinate [-1,1]
    xi = np.empty(d)

    for k, (a, b) in enumerate(combo):

        xi[k] = (
            2.0 * x[k] - (a + b)
        ) / (b - a)

    result = C

    # Evaluate one Chebyshev dimension at a time
    for k in range(d):

        result = chebval(
            xi[k],
            result,
            tensor=False
        )

    return result    

def chebyshev_approximation(func, nvec, combo):

    # 1. Generate Chebyshev grid
    nodes, points = chebyshev_grid_nd(
        nvec,
        combo
    )

    # 2. Evaluate function
    F = evaluate_function_on_grid(
        func,
        points
    )

    # 3. Reshape
    F_grid = reshape_grid_values(
        F,
        nvec
    )

    # 4. Calculate coefficients
    C = chebyshev_coefficients(
        F_grid
    )

    return nodes, points, F_grid, C