import numpy as np
from scipy.special import comb

from Affine_ArithmeticClassV3 import AffineArray


# ============================================================
# Bernstein basis as AffineArray
# ============================================================

# def polylist_Bernstein(x, n, cheb=False):

#     """
#     Construct the Bernstein basis polynomials

#         B_{k,n}(x) = C(n,k) x^k (1-x)^(n-k)

#     using AffineArray arithmetic.

#     Parameters
#     ----------
#     x : AffineArray
#         Affine variable defined over [0,1].

#     n : int
#         Bernstein polynomial degree.

#     cheb : bool
#         Passed to AffineArray operations.

#     Returns
#     -------
#     basis : list
#         List containing

#             [B_0,n(x), ..., B_n,n(x)]
#     """

#     basis = []

#     if n == 0:
#         return [1.0]

#     # ------------------------------------------
#     # k = 0
#     # B_0,n(x) = (1-x)^n
#     # ------------------------------------------

#     B0 = (
#         (1.0 - x).pow(
#             n,
#             cheb=cheb
#         )
#     )

#     basis.append(B0)

#     # ------------------------------------------
#     # k = 1,...,n-1
#     # ------------------------------------------

#     for k in range(1, n):

#         Bk = (
#             comb(n, k)
#             * x.pow(
#                 k,
#                 cheb=cheb
#             )
#             * (1.0 - x).pow(
#                 n - k,
#                 cheb=cheb
#             )
#         )

#         basis.append(Bk)

#     # ------------------------------------------
#     # k = n
#     # B_n,n(x) = x^n
#     # ------------------------------------------

#     Bn = x.pow(
#         n,
#         cheb=cheb
#     )

#     basis.append(Bn)

#     return basis

def polylist_Bernstein(x, n, cheb=False):

    """
    Construct all Bernstein basis polynomials of degree n
    using the Bernstein recurrence and AffineArray arithmetic.

    Returns
    -------
    basis : list
        [B_0,n(x), ..., B_n,n(x)]
    """

    # --------------------------------------------------------
    # Degree zero
    # --------------------------------------------------------

    if n == 0:
        return [
            AffineArray.from_intervals(
                [(1.0, 1.0)]
            )[0]
        ]

    # --------------------------------------------------------
    # Start with degree zero:
    #
    # B_0,0(x) = 1
    # --------------------------------------------------------

    basis = [
        AffineArray.from_intervals(
            [(1.0, 1.0)]
        )[0]
    ]

    one_minus_x = 1.0 - x

    # --------------------------------------------------------
    # Recursively construct degree 1,...,n
    # --------------------------------------------------------

    for degree in range(1, n + 1):

        new_basis = []

        # ----------------------------------------------------
        # B_0,degree
        # ----------------------------------------------------

        B0 = (
            one_minus_x
            * basis[0]
        )

        new_basis.append(B0)

        # ----------------------------------------------------
        # B_k,degree
        # ----------------------------------------------------

        for k in range(1, degree):

            Bk = (
                one_minus_x * basis[k]
                + x * basis[k - 1]
            )

            new_basis.append(Bk)

        # ----------------------------------------------------
        # B_degree,degree
        # ----------------------------------------------------

        Bn = (
            x * basis[-1]
        )

        new_basis.append(Bn)

        basis = new_basis

    return basis


# ============================================================
# Map physical AA variable -> [0,1]
# ============================================================

def map_to_bernstein(X_sub, combo):

    """ Map an AffineArray X defined on the physical subdomain 
    combo = [(a1,b1), ..., (ad,bd)] 
    to Bernstein coordinates xi in [0,1]^d: xi_k = (X_k - a_k) / (b_k - a_k) 
    X is a single AffineArray containing all input variables. 
    """

    # Number of variables represented by X

    d = X_sub.x0.size

    if len(combo) != d:
        raise ValueError(
            "Dimension of X_sub and combo do not match."
        )

    xi = []

    for k in range(d):

        a, b = combo[k]

        if b == a:
            raise ValueError(
                f"Zero-width interval in dimension {k}."
            )

        # Extract the kth affine variable
        xi_k = (
            X_sub[k] - a
        ) / (b - a)

        xi.append(xi_k)

    return xi


# ============================================================
# Tensor-product Bernstein basis as AffineArray
# ============================================================

def bernstein_basis_AA(
    t,
    nvec,
    cheb=False
):
    """
    Construct the tensor-product Bernstein basis using
    AffineArray arithmetic.

    Parameters
    ----------
    t : list
        AffineArray variables in [0,1].

    nvec : list or tuple
        Polynomial degree in each dimension.

    cheb : bool
        Passed to AffineArray operations.

    Returns
    -------
    base_mat : ndarray(dtype=object)
        Tensor containing the tensor-product Bernstein basis.

    Example
    -------
    For d=2 and nvec=[2,2]:

        base_mat[i,j]
        =
        B_i,2(t0) * B_j,2(t1)
    """

    d = len(nvec)

    if len(t) != d:
        raise ValueError(
            "Length of t must equal len(nvec)."
        )

    # --------------------------------------------------------
    # Construct 1-D bases
    # --------------------------------------------------------

    bases = []

    for k in range(d):

        bases.append(
            polylist_Bernstein(
                t[k],
                nvec[k],
                cheb=cheb
            )
        )

    # --------------------------------------------------------
    # Construct tensor-product basis
    # --------------------------------------------------------

    shape = tuple(
        n + 1
        for n in nvec
    )

    base_mat = np.empty(
        shape,
        dtype=object
    )

    for index in np.ndindex(shape):

        # First dimension
        term = bases[0][index[0]]

        # Remaining dimensions
        for k in range(1, d):

            term = (
                term
                * bases[k][index[k]]
            )

        base_mat[index] = term

    return base_mat


# ============================================================
# Bernstein polynomial evaluation using AA
# ============================================================

def evaluate_bernstein_AA(
    C,
    X,
    combo,
    nvec,
    cheb=False
):
    """
    Evaluate a tensor-product Bernstein polynomial using
    Affine Arithmetic.

    Parameters
    ----------
    C : ndarray
        Bernstein coefficient tensor.

    X : AffineArray
        Physical-domain affine variables.

    combo : list of tuples
        Physical subdomain.

    nvec : list or tuple
        Polynomial degree in each dimension.

    cheb : bool
        Passed to AffineArray operations.

    Returns
    -------
    p_AA : AffineArray
        Affine enclosure of the Bernstein polynomial.
    """

    # --------------------------------------------------------
    # Map physical variables -> [0,1]
    # --------------------------------------------------------

    t = map_to_bernstein(
        X,
        combo
    )

    # --------------------------------------------------------
    # Construct tensor-product basis
    # --------------------------------------------------------

    base_mat = bernstein_basis_AA(
        t,
        nvec,
        cheb=cheb
    )

    # --------------------------------------------------------
    # Polynomial evaluation
    #
    # p(t) = sum C_i * B_i(t)
    # --------------------------------------------------------

    p = None

    for index in np.ndindex(
        base_mat.shape
    ):

        term = (
            C[index]
            * base_mat[index]
        )

        if p is None:
            p = term
        else:
            p = p + term

    return p


# ============================================================
# Residual evaluation
# ============================================================

def residual_eval_Bernstein_AA(
    X,
    C,
    combo,
    nvec,
    func_AA,
    cheb=False
):
    """
    Evaluate the original function, Bernstein polynomial,
    and their residual using Affine Arithmetic.

    The residual is

        r(x) = f(x) - p_B(x)

    Two different residual bounds are returned:

    1. residual_AA

       Direct AA evaluation of

           f_AA - p_AA

    2. residual_Bernstein

       Residual bound obtained from the convex-hull property
       of the Bernstein polynomial.

       Since

           p(x) = sum_i C_i B_i(x)

       and

           B_i(x) >= 0,
           sum_i B_i(x) = 1,

       we have

           min(C_i) <= p(x) <= max(C_i).

    Parameters
    ----------
    X : AffineArray
        Physical-domain affine variables.

    C : ndarray
        Bernstein coefficient tensor.

    combo : list of tuples
        Physical subdomain.

    nvec : list or tuple
        Polynomial degree in each dimension.

    func_AA : callable
        Function implemented using AffineArray operations.

    cheb : bool
        Passed to AffineArray operations.

    Returns
    -------
    dict
        Dictionary containing:

            f_AA
            p_AA
            p_Bernstein_coeff
            residual_AA
            residual_Bernstein
    """

    # --------------------------------------------------------
    # Original function
    # --------------------------------------------------------

    f_AA = func_AA(
        X,
        cheb=cheb
    )

    # --------------------------------------------------------
    # Bernstein polynomial evaluated using AA
    # --------------------------------------------------------

    p_AA = evaluate_bernstein_AA(
        C,
        X,
        combo,
        nvec,
        cheb=cheb
    )

    # --------------------------------------------------------
    # Direct AA residual
    # --------------------------------------------------------

    residual_AA = (
        f_AA
        - p_AA
    )

    # --------------------------------------------------------
    # Bernstein coefficient convex-hull bound
    # --------------------------------------------------------

    C_values = np.asarray(
        C,
        dtype=float
    )

    p_lo = np.min(
        C_values
    )

    p_hi = np.max(
        C_values
    )

    # --------------------------------------------------------
    # Function interval
    # --------------------------------------------------------

    f_lo, f_hi = f_AA.interval

    # --------------------------------------------------------
    # Residual:
    #
    # r = f - p
    #
    # Lower bound:
    #
    #     f_lo - p_hi
    #
    # Upper bound:
    #
    #     f_hi - p_lo
    # --------------------------------------------------------

    f_AA_Int = AffineArray.from_intervals([(f_lo, f_hi)])
    f_p = AffineArray.from_intervals([(p_lo, p_hi)])

    r_bernstein = f_AA_Int[0] - f_p[0]

    residual_lo, residual_hi = r_bernstein.interval


    # residual_lo = (
    #     f_lo - p_hi
    # )

    # residual_hi = (
    #     f_hi - p_lo
    # )

    # Construct an AffineArray containing this interval.
    # residual_Bernstein = (
    #     AffineArray.from_intervals(
    #         [
    #             (
    #                 residual_lo,
    #                 residual_hi
    #             )
    #         ]
    #     )
    # )

    residual_Bernstein = (
        residual_lo,
        residual_hi
    )

    # Bernstein polynomial coefficient bound
    # p_Bernstein_coeff = (
    #     AffineArray.from_intervals(
    #         [
    #             (
    #                 p_lo,
    #                 p_hi
    #             )
    #         ]
    #     )
    # )

    p_Bernstein_coeff = (
        p_lo,
        p_hi
    )

    return {
        "f_AA": f_AA,
        "p_AA": p_AA,
        "p_Bernstein_coeff": p_Bernstein_coeff,
        "residual_AA": residual_AA,
        "residual_Bernstein": residual_Bernstein
    }

