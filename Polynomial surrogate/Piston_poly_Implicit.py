import numpy as np
from itertools import product
from scipy.special import comb

from Affine_ArithmeticClassV3 import AffineArray

from BernsteinFiles import(
    bernstein_approximation,
    bernstein_grid_nd,
    reshape_grid_values,
    bernstein_coefficients
)

from AA_evaluationBr import (
    residual_eval_Bernstein_AA,
    evaluate_bernstein_AA
)

from subintervalization import (
    split_interval,
    hull_intervals
)

def scale(x, lo, hi):
    return lo + (hi - lo)*x

def piston_AB_numpy(x):
    """
    Numerical evaluation of A(x) and B(x)
    such that C(x) = A(x) / B(x).
    """

    M  = scale(x[0], 30, 60)
    S  = scale(x[1], 0.005, 0.02)
    V0 = scale(x[2], 0.002, 0.01)
    k  = scale(x[3], 1000, 5000)
    P0 = scale(x[4], 90000, 110000)
    Ta = scale(x[5], 290, 296)
    T0 = scale(x[6], 340, 360)

    # --------------------------------------------------------
    # Piston volume equation
    # --------------------------------------------------------

    A0 = (
        P0 * S
        + 19.62 * M
        - k * V0 / S
    )

    V = (
        S / (2 * k)
        * (
            np.sqrt(
                A0**2
                + 4 * k * (P0 * V0) / T0 * Ta
            )
            - A0
        )
    )

    # --------------------------------------------------------
    # Implicit C = A / B
    # --------------------------------------------------------

    A = 2 * np.pi * np.sqrt(M)

    B = np.sqrt(
        k
        + S**2
        * (P0 * V0) / T0
        * Ta
        / V**2
    )

    return A, B


# def bernstein_design_matrix(points, nvec, combo):
#     """
#     Construct the tensor-product Bernstein design matrix.

#     Parameters
#     ----------
#     points : array_like, shape (m, d)
#         Sample points in the physical subdomain.

#     nvec : tuple/list
#         Polynomial degree in each dimension.

#     combo : list of (lo, hi)
#         Physical bounds of the current subdomain.

#     Returns
#     -------
#     Q : ndarray, shape (m, ncoeff)
#         Bernstein design matrix.
#     """

#     points = np.asarray(points, dtype=float)

#     d = len(nvec)

#     # Multi-index of Bernstein coefficients
#     indices = list(product(*[
#         range(n + 1) for n in nvec
#     ]))

#     Q = np.empty(
#         (len(points), len(indices)),
#         dtype=float
#     )

#     for i, x in enumerate(points):

#         # ----------------------------------------------------
#         # Map physical coordinates x -> local Bernstein
#         # coordinates t in [0, 1]^d
#         # ----------------------------------------------------
#         t = np.empty(d, dtype=float)

#         for dim in range(d):

#             lo, hi = combo[dim]

#             t[dim] = (
#                 x[dim] - lo
#             ) / (
#                 hi - lo
#             )

#         # ----------------------------------------------------
#         # Tensor-product Bernstein basis
#         # ----------------------------------------------------
#         for j, idx in enumerate(indices):

#             value = 1.0

#             for dim in range(d):

#                 k = idx[dim]
#                 n = nvec[dim]

#                 value *= (
#                     comb(n, k)
#                     * t[dim]**k
#                     * (1.0 - t[dim])**(n - k)
#                 )

#             Q[i, j] = value

#     return Q

# def implicit_bernstein_coefficients(
#     A_values,
#     B_values,
#     Q,
# ):
#     """
#     Estimate Bernstein coefficients from an implicit relation

#         A(x) - B(x) p(x) = 0

#     by solving

#         min_c || A - B * (Q @ c) ||_2.

#     Parameters
#     ----------
#     A_values : array_like, shape (m,)
#         Values of A at the sample points.

#     B_values : array_like, shape (m,)
#         Values of B at the sample points.

#     Q : ndarray, shape (m, ncoeff)
#         Bernstein design matrix.

#     Returns
#     -------
#     c : ndarray, shape (ncoeff,)
#         Estimated Bernstein coefficients.
#     residuals : ndarray
#         Least-squares residual information.
#     rank : int
#         Numerical rank of the weighted design matrix.
#     singular_values : ndarray
#         Singular values of the weighted design matrix.
#     """

#     A_values = np.asarray(A_values, dtype=float)
#     B_values = np.asarray(B_values, dtype=float)
#     Q = np.asarray(Q, dtype=float)

#     M = B_values[:, None] * Q

#     c, residuals, rank, singular_values = np.linalg.lstsq(
#         M,
#         A_values,
#         rcond=None
#     )

#     return c, residuals, rank, singular_values

# def piston_bernstein_coefficients(
#     nvec,
#     combo,
#     node_type="bernstein"
# ):

#     # --------------------------------------------------------
#     # Generate interpolation grid
#     # --------------------------------------------------------

#     nodes, points = bernstein_grid_nd(
#         nvec,
#         combo,
#         node_type=node_type
#     )

#     # --------------------------------------------------------
#     # Evaluate A and B at every node
#     # --------------------------------------------------------

#     A_values = []
#     B_values = []

#     for point in points:

#         A, B = piston_AB_numpy(point)

#         A_values.append(A)
#         B_values.append(B)

#     A_values = np.asarray(A_values)
#     B_values = np.asarray(B_values)

#     # --------------------------------------------------------
#     # Explicit C values ONLY for coefficient construction
#     # --------------------------------------------------------

#     C_values = A_values / B_values

#     # --------------------------------------------------------
#     # Reshape tensor grid
#     # --------------------------------------------------------

#     A_grid = reshape_grid_values(
#         A_values,
#         nvec
#     )

#     B_grid = reshape_grid_values(
#         B_values,
#         nvec
#     )

#     C_grid = reshape_grid_values(
#         C_values,
#         nvec
#     )

#     # --------------------------------------------------------
#     # Bernstein coefficients
#     # --------------------------------------------------------

#     Q_matrix = bernstein_design_matrix(
#     points,
#     nvec,
#     combo
#     )

#     C_coeff, residuals, rank, singular_values = (
#         implicit_bernstein_coefficients(
#             A_values,
#             B_values,
#             Q_matrix
#         )
#     )

#     return (
#         nodes,
#         points,
#         A_grid,
#         B_grid,
#         C_grid,
#         C_coeff
#     )


# def piston_bernstein_coefficients(
#     nvec,
#     combo,
#     node_type="bernstein"
# ):
#     """
#     Compute Bernstein coefficients for the implicit piston model

#         A(x) - B(x) p(x) = 0.

#     At the interpolation nodes,

#         p(x) = A(x) / B(x).

#     The resulting values are then converted to Bernstein
#     coefficients dimension-by-dimension.
#     """

#     # --------------------------------------------------------
#     # Generate interpolation grid
#     # --------------------------------------------------------

#     nodes, points = bernstein_grid_nd(
#         nvec,
#         combo,
#         node_type=node_type
#     )

#     # --------------------------------------------------------
#     # Evaluate A and B
#     # --------------------------------------------------------

#     A_values = []
#     B_values = []

#     for point in points:

#         A, B = piston_AB_numpy(point)

#         A_values.append(A)
#         B_values.append(B)

#     A_values = np.asarray(A_values)
#     B_values = np.asarray(B_values)

#     # --------------------------------------------------------
#     # Explicit C values at interpolation nodes
#     # --------------------------------------------------------

#     C_values = A_values / B_values

#     # --------------------------------------------------------
#     # Reshape tensor grid
#     # --------------------------------------------------------

#     # A_grid = reshape_grid_values(
#     #     A_values,
#     #     nvec
#     # )

#     # B_grid = reshape_grid_values(
#     #     B_values,
#     #     nvec
#     # )

#     C_grid = reshape_grid_values(
#         C_values,
#         nvec
#     )

#     # --------------------------------------------------------
#     # Compute Bernstein coefficients
#     # --------------------------------------------------------

#     C_coeff = bernstein_coefficients(
#         C_grid,
#         nvec,
#         combo,
#         node_type=node_type
#     )

#     return (
#         nodes,
#         points,
#         C_grid,
#         C_coeff
#     )


def piston_bernstein_coefficients(
    nvec,
    combo,
    node_type="bernstein"
):

    # --------------------------------------------------------
    # Generate interpolation grid
    # --------------------------------------------------------

    nodes, points = bernstein_grid_nd(
        nvec,
        combo,
        node_type=node_type
    )

    # --------------------------------------------------------
    # Evaluate A and B
    # --------------------------------------------------------

    A_values = []
    B_values = []

    for point in points:

        A, B = piston_AB_numpy(point)

        A_values.append(A)
        B_values.append(B)

    A_values = np.asarray(A_values)
    B_values = np.asarray(B_values)

    # --------------------------------------------------------
    # Explicit C values
    # --------------------------------------------------------

    C_values = A_values / B_values

    # --------------------------------------------------------
    # Reshape tensor grids
    # --------------------------------------------------------

    A_grid = reshape_grid_values(
        A_values,
        nvec
    )

    B_grid = reshape_grid_values(
        B_values,
        nvec
    )

    C_grid = reshape_grid_values(
        C_values,
        nvec
    )

    # --------------------------------------------------------
    # Bernstein coefficients
    # --------------------------------------------------------

    C_coeff = bernstein_coefficients(
        C_grid,
        nvec,
        combo,
        node_type=node_type
    )

    B_coeff = bernstein_coefficients(
        B_grid,
        nvec,
        combo,
        node_type=node_type
    )

    return (
        nodes,
        points,
        A_grid,
        B_grid,
        C_grid,
        B_coeff,
        C_coeff
    )


# def piston_AB_AA(x, cheb=False):

#     M  = scale(x[0], 30, 60)
#     S  = scale(x[1], 0.005, 0.02)
#     V0 = scale(x[2], 0.002, 0.01)
#     k  = scale(x[3], 1000, 5000)
#     P0 = scale(x[4], 90000, 110000)
#     Ta = scale(x[5], 290, 296)
#     T0 = scale(x[6], 340, 360)

#     # --------------------------------------------------------
#     # Piston volume equation
#     # --------------------------------------------------------

#     A0 = (
#         P0 * S
#         + 19.62 * M
#         - k * V0 / S
#     )

#     V = (
#         S / (2 * k)
#         * (
#             (
#                 A0.pow(2, cheb=cheb)
#                 + 4 * k * (P0 * V0) / T0 * Ta
#             ).sqrt(cheb=cheb)
#             - A0
#         )
#     )

#     # --------------------------------------------------------
#     # Implicit A and B
#     # --------------------------------------------------------

#     A = (
#         2 * np.pi
#         * M.sqrt(cheb=cheb)
#     )

#     B = (
#         k
#         + S.pow(2, cheb=cheb)
#         * (P0 * V0) / T0
#         * Ta
#         / V.pow(2, cheb=cheb)
#     ).sqrt(cheb=cheb)

#     return A, B


def piston_A_AA(x, cheb=False):

    M = scale(x[0], 30, 60)

    A = (
        2 * np.pi
        * M.sqrt(cheb=cheb)
    )

    return A


# def residual_piston_Bernstein_AA(
#     X,
#     C,
#     B,
#     combo,
#     nvec,
#     cheb=False
# ):

#     # --------------------------------------------------------
#     # Evaluate A using Affine Arithmetic
#     # --------------------------------------------------------

#     A_AA = piston_A_AA(
#         X,
#         cheb=cheb
#     )

#     # --------------------------------------------------------
#     # Evaluate p1(x) ≈ C(x)
#     # --------------------------------------------------------

#     p1_AA = evaluate_bernstein_AA(
#         C,
#         X,
#         combo,
#         nvec,
#         cheb=cheb
#     )

#     # --------------------------------------------------------
#     # Evaluate p2(x) ≈ B(x)
#     # --------------------------------------------------------

#     p2_AA = evaluate_bernstein_AA(
#         B,
#         X,
#         combo,
#         nvec,
#         cheb=cheb
#     )

#     # --------------------------------------------------------
#     # Implicit residual
#     #
#     # δ(x) = A(x) - p1(x) p2(x)
#     # --------------------------------------------------------

#     delta_AA = (
#         A_AA
#         - p1_AA * p2_AA
#     )

#     return {
#         "A_AA": A_AA,
#         "p1_AA": p1_AA,
#         "p2_AA": p2_AA,
#         "delta_AA": delta_AA
#     }

def residual_piston_Bernstein_AA(
    X,
    C,
    B,
    combo,
    nvec,
    cheb=False
):

    # --------------------------------------------------------
    # Evaluate A using Affine Arithmetic
    # --------------------------------------------------------

    A_AA = piston_A_AA(
        X,
        cheb=cheb
    )

    # --------------------------------------------------------
    # Construct independent affine variables for p1 and p2
    # --------------------------------------------------------

    C_interval = (
        np.min(C),
        np.max(C)
    )

    B_interval = (
        np.min(B),
        np.max(B)
    )

    P = AffineArray.from_intervals(
        list(combo) + [C_interval, B_interval]
    )

    n = len(combo)

    p1_AA = P[n]
    p2_AA = P[n+1]

    # --------------------------------------------------------
    # Implicit residual
    # --------------------------------------------------------

    delta_AA = (
        A_AA
        - p1_AA * p2_AA
    )

    return {
        "A_AA": A_AA,
        "p1_AA": p1_AA,
        "p2_AA": p2_AA,
        "delta_AA": delta_AA
    }


def piston_implicit_subinterval(
    nvec,
    splits_per_dim,
    cheb=False
):

    subdomains = split_interval(
        splits_per_dim
    )

    results = []

    for combo in subdomains:

        # ====================================================
        # 1. Numerical Bernstein approximation
        # ====================================================

        (
            nodes,
            points,
            A_grid,
            B_grid,
            C_grid,
            B_coeff,
            C_coeff
        ) = piston_bernstein_coefficients(
            nvec,
            combo,
            node_type="bernstein"
        )


        # ====================================================
        # 2. Construct AA variables
        # ====================================================

        X = AffineArray.from_intervals(
            combo
        )

        # ====================================================
        # 3. Evaluate implicit residual
        # ====================================================

        result = residual_piston_Bernstein_AA(
            X,
            C_coeff,
            B_coeff,
            combo,
            nvec,
            cheb=cheb
        )

        p1_AA = result["p1_AA"]
        p2_AA = result["p2_AA"]
        delta_AA = result["delta_AA"]

        # ====================================================
        # 4. AA polynomial range
        # ====================================================

        p2_AA_interval = p2_AA.interval

        # ====================================================
        # 5. AA implicit residual range
        # ====================================================

        delta_interval = delta_AA.interval

        results.append({

            "combo": combo,

            "B_coeff": B_coeff,
            "p1_coeff": C_coeff,

            "p1_AA": p1_AA,
            "p2_AA": p2_AA,

            "p1_AA_interval": p1_AA.interval,
            "p2_AA_interval": p2_AA_interval,

            "delta_AA": delta_AA,
            "delta_interval": delta_interval
        })

    # ============================================================
    # 6. Global hulls over all subdomains
    # ============================================================

    p1_coeff_intervals = [
        (np.min(r["p1_coeff"]), np.max(r["p1_coeff"]))
        for r in results
    ]

    p1_AA_intervals = [
        r["p1_AA_interval"]
        for r in results
    ]

    p2_AA_intervals = [
        r["p2_AA_interval"]
        for r in results
    ]

    delta_intervals = [
        r["delta_interval"]
        for r in results
    ]

    global_p1_coeff = hull_intervals(
        p1_coeff_intervals
    )

    global_p1_AA = hull_intervals(
        p1_AA_intervals
    )

    global_p2_AA = hull_intervals(
        p2_AA_intervals
    )

    global_delta = hull_intervals(
        delta_intervals
    )

    return (
        results,
        global_p1_coeff,
        global_p1_AA,
        global_p2_AA,
        global_delta
    )


#--------------------------------------------------------------------------
#   User defined nvec and split values
#---------------------------------------------------------------------------

nvec = [2] * 7

split_values = [1,2,4]

#===========================================================

for splits in split_values:

    splits_per_dim = [splits] * 7

    (
        results,
        global_p1_coeff,
        global_p1_AA,
        global_p2_AA,
        global_delta
    ) = piston_implicit_subinterval(
        nvec=nvec,
        splits_per_dim=splits_per_dim
    )

    print(f"\n{'=' * 60}")
    print(f"Splits = {splits} x 7")
    print(f"Number of subdomains = {len(results)}")
    print(f"{'=' * 60}")

    print("Global p1 coefficient hull :", global_p1_coeff)
    print("Global p1_AA hull         :", global_p1_AA)
    print("Global p2_AA hull         :", global_p2_AA)
    print("Global delta hull         :", global_delta)



