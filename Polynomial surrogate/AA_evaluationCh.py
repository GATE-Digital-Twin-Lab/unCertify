import numpy as np
from math import comb

from Affine_ArithmeticClassV3 import AffineArray

def map_to_chebyshev(X, combo):
    """
    Map an AffineArray from a physical subdomain to [-1, 1]^d.

    Parameters
    ----------
    X : AffineArray
        AA representation of the variables.

    combo : list of tuple
        Subdomain intervals:
            [(a1, b1), (a2, b2), ..., (ad, bd)]

    Returns
    -------
    xi : AffineArray
        Variables mapped from the physical subdomain to [-1, 1]^d.
    """

    xi = []

    for k, (a, b) in enumerate(combo):

        if b <= a:
            raise ValueError(
                f"Invalid interval in dimension {k}: {(a, b)}"
            )

        # Physical interval [a,b] -> Chebyshev interval [-1,1]
        xi_k = (2.0 * X[k] - (a + b)) / (b - a)

        xi.append(xi_k)

    return xi


def polylist_Ch(x, n, cheb=False):
    """
    Construct the Chebyshev basis

        [T_0(x), T_1(x), ..., T_n(x)]

    for an AA variable x.

    Parameters
    ----------
    x : AffineScalar
        AA variable.

    n : int
        Polynomial degree.

    cheb : bool, optional
        Passed to AffineScalar operations.

    Returns
    -------
    basis : list
        Chebyshev basis evaluated using affine arithmetic.
    """

    if n < 0:
        raise ValueError("Polynomial degree must be non-negative.")

    basis = [None] * (n + 1)

    # T_0(x) = 1
    # basis[0] = AffineArray.from_intervals([(1.0, 1.0)])
    basis[0] = 1.0


    if n == 0:
        return basis

    # T_1(x) = x
    basis[1] = x

    # Recurrence:
    #
    # T_k(x) = 2*x*T_{k-1}(x) - T_{k-2}(x)

    for k in range(2, n + 1):

        basis[k] = (
            2.0 * x * basis[k - 1]
            - basis[k - 2]
        )

    return basis


def chebyshev_basis_AA(xi, nvec, cheb=False):
    """
    Construct the tensor-product Chebyshev basis using AA.

    For d dimensions and degrees nvec = [n1,...,nd],

        basis[i1,...,id]
        =
        T_i1(xi[0]) * ... * T_id(xi[d-1])

    Parameters
    ----------
    xi : list
        Chebyshev-mapped AA variables.

    nvec : list
        Polynomial degree in each dimension.

    cheb : bool, optional
        Passed to AffineArray operations.

    Returns
    -------
    base_mat : ndarray
        Object array containing the AA basis functions.
    """

    d = len(nvec)

    if len(xi) != d:
        raise ValueError(
            "Length of xi must match length of nvec."
        )

    # Build one-dimensional basis in every dimension
    bases = []

    for k in range(d):

        basis_k = polylist_Ch(
            xi[k],
            nvec[k],
            cheb=cheb
        )

        bases.append(basis_k)

    # Create tensor-product basis
    shape = tuple(n + 1 for n in nvec)

    base_mat = np.empty(shape, dtype=object)

    for index in np.ndindex(shape):

        term = bases[0][index[0]]

        for k in range(1, d):

            term = term * bases[k][index[k]]

        base_mat[index] = term

    return base_mat


def evaluate_chebyshev_AA(C, X, combo, nvec, cheb=False):
    """
    Evaluate a Chebyshev polynomial with AA arithmetic.

    Parameters
    ----------
    C : ndarray
        Chebyshev coefficient tensor.

    X : AffineArray
        AA variables in the physical domain.

    combo : list of tuple
        Physical subdomain.

    nvec : list
        Polynomial degree in each dimension.

    cheb : bool, optional
        Passed to AffineArray operations.

    Returns
    -------
    p : AffineScalar
        AA enclosure of the Chebyshev polynomial.
    """

    xi = map_to_chebyshev(X, combo)

    base_mat = chebyshev_basis_AA(
        xi,
        nvec,
        cheb=cheb
    )

    p = None

    for index in np.ndindex(base_mat.shape):

        term = C[index] * base_mat[index]

        if p is None:
            p = term
        else:
            p = p + term

    return p


def residual_eval_AA(
    X,
    C,
    combo,
    nvec,
    func_AA,
    cheb=False
):
    """
    Evaluate f(X), p(X), and the residual

        r(X) = f(X) - p(X)

    using affine arithmetic.

    Parameters
    ----------
    X : AffineArray
        AA representation of the variables.

    C : ndarray
        Chebyshev coefficient tensor.

    combo : list of tuple
        Physical subdomain.

    nvec : list
        Polynomial degree in each dimension.

    func_AA : callable
        AA implementation of the function being approximated.

    cheb : bool, optional
        Passed to AA operations.

    Returns
    -------
    f : AffineArray / AffineScalar
        AA evaluation of the original function.

    p : AffineArray / AffineScalar
        AA evaluation of the polynomial approximation.

    residual : AffineArray / AffineScalar
        AA evaluation of f - p.
    """

    # Original function
    f = func_AA(
        X,
        cheb=cheb
    )

    # Chebyshev polynomial
    p = evaluate_chebyshev_AA(
        C,
        X,
        combo,
        nvec,
        cheb=cheb
    )

    # Residual
    residual = f - p

    return f, p, residual


