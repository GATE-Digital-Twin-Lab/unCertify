import numpy as np
from itertools import product
from scipy.special import comb
from scipy.special import roots_legendre


# ============================================================
# NODE GENERATION
# ============================================================

def chebyshev_lobatto_nodes(n, a, b):
    """
    Chebyshev-Lobatto nodes on the physical interval [a,b].

    Parameters
    ----------
    n : int
        Polynomial degree. Returns n+1 nodes.
    a, b : float
        Physical interval.

    Returns
    -------
    x : ndarray
        Chebyshev-Lobatto nodes on [a,b].
    """

    if n < 0:
        raise ValueError("n must be non-negative.")

    if b <= a:
        raise ValueError("Require b > a.")

    if n == 0:
        return np.array([(a + b) / 2.0])

    k = np.arange(n + 1)

    # Chebyshev-Lobatto nodes on [-1,1]
    xi = np.cos(np.pi * k / n)

    # Map [-1,1] -> [a,b]
    x = (
        0.5 * (a + b)
        + 0.5 * (b - a) * xi
    )

    return x


def legendre_lobatto_nodes(n, a, b):
    """
    Legendre-Lobatto nodes on the physical interval [a,b].

    The nodes consist of the two endpoints -1 and 1 together
    with the roots of the derivative of the Legendre polynomial
    P_n(x).

    Parameters
    ----------
    n : int
        Polynomial degree. Returns n+1 nodes.
    a, b : float
        Physical interval.

    Returns
    -------
    x : ndarray
        Legendre-Lobatto nodes on [a,b].
    """

    if n < 0:
        raise ValueError("n must be non-negative.")

    if b <= a:
        raise ValueError("Require b > a.")

    if n == 0:
        return np.array([(a + b) / 2.0])

    if n == 1:
        xi = np.array([-1.0, 1.0])

    else:
        # Roots of P'_n(x) are the interior
        # Legendre-Lobatto nodes.
        #
        # roots_legendre(n-1) gives roots of P_(n-1),
        # which are NOT the desired roots, so we obtain
        # the roots of P'_n using the companion relation.
        #
        # For numerical robustness, use the eigenvalue
        # construction for Legendre-Lobatto nodes.

        from numpy.polynomial.legendre import Legendre

        Pn = Legendre.basis(n)

        dPn = Pn.deriv()

        interior = dPn.roots()

        xi = np.concatenate(
            ([-1.0], interior, [1.0])
        )

    # Map [-1,1] -> [a,b]
    x = (
        0.5 * (a + b)
        + 0.5 * (b - a) * xi
    )

    return x


def bernstein_nodes(n, a, b):
    """
    Bernstein nodes on the physical interval [a,b].

        x_k = a + (b-a) k/n

    Parameters
    ----------
    n : int
        Polynomial degree.
    a, b : float
        Physical interval.

    Returns
    -------
    x : ndarray
        Bernstein nodes on [a,b].
    """

    if n < 0:
        raise ValueError("n must be non-negative.")

    if b <= a:
        raise ValueError("Require b > a.")

    if n == 0:
        return np.array([(a + b) / 2.0])

    k = np.arange(n + 1)

    x = (
        a
        + (b - a) * k / n
    )

    return x


# ============================================================
# NODE DISPATCHER
# ============================================================

def get_nodes(n, a, b, node_type="bernstein"):
    """
    Generate interpolation nodes according to node_type.

    Parameters
    ----------
    n : int
        Polynomial degree.
    a, b : float
        Physical interval.
    node_type : str
        Node family:

            "bernstein"
            "chebyshev"
            "legendre"

    Returns
    -------
    ndarray
        n+1 interpolation nodes on [a,b].
    """

    node_type = node_type.lower()

    if node_type == "bernstein":

        return bernstein_nodes(
            n, a, b
        )

    elif node_type == "chebyshev":

        return chebyshev_lobatto_nodes(
            n, a, b
        )

    elif node_type == "legendre":

        return legendre_lobatto_nodes(
            n, a, b
        )

    else:

        raise ValueError(
            "Unknown node_type. Choose "
            "'bernstein', 'chebyshev', or 'legendre'."
        )


# ============================================================
# BERNSTEIN GRID
# ============================================================

def bernstein_grid_nd(
    nvec,
    combo,
    node_type="bernstein"
):
    """
    Construct a tensor-product interpolation grid.

    Parameters
    ----------
    nvec : list or tuple
        Polynomial degree in each dimension.

    combo : list of tuples
        Physical subdomain:

            [(a1,b1), (a2,b2), ..., (ad,bd)]

    node_type : str
        Interpolation node family.

    Returns
    -------
    nodes : list
        One node array for each dimension.

    points : ndarray
        Tensor-product grid points.
    """

    d = len(nvec)

    if len(combo) != d:
        raise ValueError(
            "Length of combo must equal len(nvec)."
        )

    nodes = []

    for k in range(d):

        n = nvec[k]

        a, b = combo[k]

        nodes_k = get_nodes(
            n,
            a,
            b,
            node_type=node_type
        )

        nodes.append(nodes_k)

    points = np.array(
        list(product(*nodes))
    )

    return nodes, points


# ============================================================
# FUNCTION VALUES ON GRID
# ============================================================

def evaluate_function_on_grid(func, points):
    """
    Evaluate a scalar function at all grid points.
    """

    return np.array([
        func(point)
        for point in points
    ])


def reshape_grid_values(F, nvec):
    """
    Reshape flattened grid values into tensor-product form.
    """

    shape = tuple(
        n + 1
        for n in nvec
    )

    return F.reshape(shape)


# ============================================================
# BERNSTEIN BASIS
# ============================================================

def bernstein_basis(k, n, x):
    """
    Evaluate the k-th Bernstein basis polynomial:

        B_{k,n}(x)
        = C(n,k) x^k (1-x)^(n-k)

    Parameters
    ----------
    k : int
        Basis index.
    n : int
        Polynomial degree.
    x : float or ndarray
        Point(s) in [0,1].

    Returns
    -------
    float or ndarray
        Bernstein basis value.
    """

    if k < 0 or k > n:
        raise ValueError(
            "Require 0 <= k <= n."
        )

    return (
        comb(n, k)
        * x**k
        * (1.0 - x)**(n - k)
    )


# def bernstein_basis(n, x):
#     """
#     Construct all Bernstein basis polynomials of degree n
#     using the Bernstein recurrence.

#     Returns
#     -------
#     basis : list
#         basis[k] = B_{k,n}(x)
#     """

#     # Degree zero
#     basis = [np.ones_like(x, dtype=float)]

#     for degree in range(1, n + 1):

#         new_basis = []

#         # B_{0,degree}
#         new_basis.append(
#             (1.0 - x) * basis[0]
#         )

#         # B_{k,degree}, 0 < k < degree
#         for k in range(1, degree):

#             value = (
#                 (1.0 - x) * basis[k]
#                 + x * basis[k - 1]
#             )

#             new_basis.append(value)

#         # B_{degree,degree}
#         new_basis.append(
#             x * basis[-1]
#         )

#         basis = new_basis

#     return basis


# ============================================================
# BERNSTEIN MATRIX
# ============================================================

def bernstein_matrix(n, nodes):
    """
    Construct the Bernstein interpolation matrix.

        A[i,k] = B_{k,n}(nodes[i])

    Parameters
    ----------
    n : int
        Bernstein polynomial degree.
    nodes : ndarray
        Nodes in [0,1].

    Returns
    -------
    A : ndarray
        Bernstein basis matrix.
    """

    nodes = np.asarray(
        nodes,
        dtype=float
    )

    A = np.empty(
        (len(nodes), n + 1)
    )

    for i, x in enumerate(nodes):

        for k in range(n + 1):

            A[i, k] = bernstein_basis(
                k,
                n,
                x
            )

    return A



# ============================================================
# BERNSTEIN MATRIX
# ============================================================

# def bernstein_matrix(n, nodes):
#     """
#     Construct the Bernstein interpolation matrix.

#         A[i,k] = B_{k,n}(nodes[i])

#     Parameters
#     ----------
#     n : int
#         Bernstein polynomial degree.

#     nodes : ndarray
#         Nodes in [0,1].

#     Returns
#     -------
#     A : ndarray
#         Bernstein basis matrix.
#     """

#     nodes = np.asarray(
#         nodes,
#         dtype=float
#     )

#     A = np.empty(
#         (len(nodes), n + 1)
#     )

#     for i, x in enumerate(nodes):

#         basis = bernstein_basis(
#             n,
#             x
#         )

#         for k in range(n + 1):

#             A[i, k] = basis[k]

#     return A


# ============================================================
# BERNSTEIN COEFFICIENTS
# ============================================================

def bernstein_coefficients(
    F,
    nvec,
    combo,
    node_type="bernstein"
):
    """
    Compute tensor-product Bernstein coefficients.

    The function values F are assumed to correspond to the
    tensor-product grid generated by get_nodes().

    Parameters
    ----------
    F : ndarray
        Function values on the grid.

        Shape:
            (n1+1, n2+1, ..., nd+1)

    nvec : list or tuple
        Polynomial degree in each dimension.

    combo : list of tuples
        Physical subdomain.

    node_type : str
        Interpolation node family:

            "bernstein"
            "chebyshev"
            "legendre"

    Returns
    -------
    C : ndarray
        Tensor of Bernstein coefficients.
    """

    F = np.asarray(
        F,
        dtype=float
    )

    d = len(nvec)

    if len(combo) != d:
        raise ValueError(
            "Length of combo must equal len(nvec)."
        )

    expected_shape = tuple(
        n + 1
        for n in nvec
    )

    if F.shape != expected_shape:
        raise ValueError(
            f"Expected F with shape {expected_shape}, "
            f"but received {F.shape}."
        )

    # --------------------------------------------------------
    # Generate interpolation nodes
    # --------------------------------------------------------

    nodes = []

    for k in range(d):

        n = nvec[k]

        a, b = combo[k]

        nodes_k = get_nodes(
            n,
            a,
            b,
            node_type=node_type
        )

        nodes.append(nodes_k)

    # --------------------------------------------------------
    # Map physical nodes -> [0,1]
    #
    # Bernstein basis is always defined on [0,1].
    # --------------------------------------------------------

    t_nodes = []

    for k in range(d):

        a, b = combo[k]

        t = (
            nodes[k] - a
        ) / (b - a)

        t_nodes.append(t)

    # --------------------------------------------------------
    # Solve interpolation system dimension by dimension
    # --------------------------------------------------------

    C = F.copy()

    for axis in range(d):

        n = nvec[axis]

        A = bernstein_matrix(
            n,
            t_nodes[axis]
        )

        # Bring the current dimension to axis 0.
        C_axis = np.moveaxis(
            C,
            axis,
            0
        )

        original_shape = C_axis.shape

        # Flatten remaining dimensions.
        C_flat = C_axis.reshape(
            n + 1,
            -1
        )

        # Solve:
        #
        # A C = F
        #
        C_flat = np.linalg.solve(
            A,
            C_flat
        )

        C_axis = C_flat.reshape(
            original_shape
        )

        # Move axis back.
        C = np.moveaxis(
            C_axis,
            0,
            axis
        )

    return C


# ============================================================
# 1-D BERNSTEIN EVALUATION
# ============================================================

def evaluate_bernstein_1d(c, x):
    """
    Evaluate a 1-D Bernstein polynomial:

        p(x) = sum_k c_k B_{k,n}(x)

    Parameters
    ----------
    c : ndarray
        Bernstein coefficients.
    x : float
        Evaluation point in [0,1].

    Returns
    -------
    float
        Polynomial value.
    """

    c = np.asarray(
        c,
        dtype=float
    )

    n = len(c) - 1

    result = 0.0

    for k in range(n + 1):

        result += (
            c[k]
            * bernstein_basis(k, n, x)
        )

    return result


# ============================================================
# MAP PHYSICAL DOMAIN -> [0,1]^d
# ============================================================

def map_to_bernstein(x, combo):
    """
    Map a physical point x from combo to [0,1]^d.

        t_k = (x_k-a_k)/(b_k-a_k)

    Parameters
    ----------
    x : array-like
        Physical evaluation point.

    combo : list of tuples
        Physical subdomain.

    Returns
    -------
    t : ndarray
        Point mapped to [0,1]^d.
    """

    x = np.asarray(
        x,
        dtype=float
    )

    d = len(combo)

    if len(x) != d:
        raise ValueError(
            "Dimension of x does not match combo."
        )

    t = np.empty(d)

    for k, (a, b) in enumerate(combo):

        if b <= a:
            raise ValueError(
                f"Invalid interval [{a},{b}] "
                f"in dimension {k}."
            )

        t[k] = (
            x[k] - a
        ) / (b - a)

    return t


# ============================================================
# ND BERNSTEIN EVALUATION
# ============================================================

def evaluate_bernstein_nd(
    C,
    x,
    combo
):
    """
    Evaluate a tensor-product Bernstein polynomial.

    Parameters
    ----------
    C : ndarray
        Tensor of Bernstein coefficients.

    x : array-like
        Physical evaluation point.

    combo : list of tuples
        Physical subdomain.

    Returns
    -------
    float
        Polynomial value.
    """

    x = np.asarray(
        x,
        dtype=float
    )

    d = len(combo)

    if len(x) != d:
        raise ValueError(
            "Dimension of x does not match combo."
        )

    # Map physical point -> [0,1]^d.
    t = map_to_bernstein(
        x,
        combo
    )

    result = C

    # Evaluate one dimension at a time.
    for axis in range(d):

        result = np.apply_along_axis(
            evaluate_bernstein_1d,
            axis,
            result,
            t[axis]
        )

    return result


# ============================================================
# COMPLETE BERNSTEIN APPROXIMATION
# ============================================================

def bernstein_approximation(
    func,
    nvec,
    combo,
    node_type="bernstein"
):
    """
    Construct a tensor-product Bernstein polynomial
    approximation of func over combo.

    Parameters
    ----------
    func : callable
        Function to approximate.

    nvec : list or tuple
        Polynomial degree in each dimension.

    combo : list of tuples
        Physical subdomain.

    node_type : str
        Interpolation node family:

            "bernstein"
            "chebyshev"
            "legendre"

    Returns
    -------
    nodes : list
        Interpolation nodes in each dimension.

    points : ndarray
        Tensor-product grid points.

    F_grid : ndarray
        Function values on the grid.

    C : ndarray
        Bernstein coefficients.
    """

    # --------------------------------------------------------
    # Construct grid
    # --------------------------------------------------------

    nodes, points = bernstein_grid_nd(
        nvec,
        combo,
        node_type=node_type
    )

    # --------------------------------------------------------
    # Evaluate original function
    # --------------------------------------------------------

    F = evaluate_function_on_grid(
        func,
        points
    )

    # --------------------------------------------------------
    # Reshape values into tensor form
    # --------------------------------------------------------

    F_grid = reshape_grid_values(
        F,
        nvec
    )

    # --------------------------------------------------------
    # Compute Bernstein coefficients
    # --------------------------------------------------------

    C = bernstein_coefficients(
        F_grid,
        nvec,
        combo,
        node_type=node_type
    )

    return (
        nodes,
        points,
        F_grid,
        C
    )

