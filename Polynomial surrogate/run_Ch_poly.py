import numpy as np

from Affine_ArithmeticClassV3 import AffineArray

from ChebyshevFiles import (
    chebyshev_approximation
)

from subintervalization import (
    split_interval,
    interval_width,
    hull_intervals
)

from AA_evaluationCh import (
    residual_eval_AA
)

from test_functions import (
    Branin_numpy,
    test_BraninAA,
    Ackley_numpy,
    test_AckleyAA,
    Eggholder_numpy,
    test_EggholderAA
)


# ============================================================
# Configuration
# ============================================================

nvec = [3, 3]

cheb = False

# Selected benchmark problem

# fun_numpy = Eggholder_numpy
# fun_AA = test_EggholderAA

fun_numpy = Ackley_numpy
fun_AA = test_AckleyAA
# ============================================================
# Subintervalization test
# ============================================================

def test_subintervalization(
    func_numpy,
    func_AA,
    splits_per_dim,
    nvec,
    cheb=False
):
    """
    Perform local Chebyshev approximation and AA residual
    evaluation over all subdomains.

    Parameters
    ----------
    func_numpy : callable
        NumPy implementation of the function.

    func_AA : callable
        Affine arithmetic implementation of the function.

    splits_per_dim : list
        Number of subintervals in each dimension.

    nvec : list
        Chebyshev polynomial degree in each dimension.

    cheb : bool
        Whether to use Chebyshev-aware AA operations.

    Returns
    -------
    results : list of dict
        Results for every subdomain.
    """

    # --------------------------------------------------------
    # Generate subdomains
    # --------------------------------------------------------

    subdomains = split_interval(splits_per_dim)

    results = []

    # --------------------------------------------------------
    # Process every subdomain
    # --------------------------------------------------------

    for combo in subdomains:

        # Create AA representation of this subdomain
        X_sub = AffineArray.from_intervals(combo)

        # ----------------------------------------------------
        # Local Chebyshev approximation
        # ----------------------------------------------------

        nodes, points, F_grid, C = chebyshev_approximation(
            func_numpy,
            nvec,
            combo
        )

        # ----------------------------------------------------
        # AA evaluation of f, p and residual
        # ----------------------------------------------------

        f_AA, p_AA, residual_AA = residual_eval_AA(
            X_sub,
            C,
            combo,
            nvec,
            func_AA,
            cheb=cheb
        )

        # ----------------------------------------------------
        # Store local results
        # ----------------------------------------------------

        results.append(
            {
                "combo": combo,

                "f_interval": f_AA.interval,

                "p_interval": p_AA.interval,

                "residual_interval": residual_AA.interval
            }
        )

    return results


# ============================================================
# Run subintervalization experiment
# ============================================================

for splits in [1, 2, 4, 8, 16, 32, 64, 128, 256]:

    results = test_subintervalization(
        fun_numpy,
        fun_AA,
        [splits, splits],
        nvec,
        cheb=cheb
    )

    # --------------------------------------------------------
    # Collect local intervals
    # --------------------------------------------------------

    f_intervals = [
        r["f_interval"]
        for r in results
    ]

    p_intervals = [
        r["p_interval"]
        for r in results
    ]

    residual_intervals = [
        r["residual_interval"]
        for r in results
    ]

    # --------------------------------------------------------
    # Global hull
    # --------------------------------------------------------

    f_hull = hull_intervals(
        f_intervals
    )

    p_hull = hull_intervals(
        p_intervals
    )

    residual_hull = hull_intervals(
        residual_intervals
    )

    # --------------------------------------------------------
    # Print
    # --------------------------------------------------------

    print(
        f"\nSplits = {splits} x {splits}"
    )

    print(
        "Number of subdomains =",
        len(results)
    )

    print(
        "Global f hull        =",
        f_hull
    )

    print(
        "Global p hull        =",
        p_hull
    )

    print(
        "Global residual hull =",
        residual_hull
    )

    print(
        "Width(f)             =",
        f_hull[1] - f_hull[0]
    )

    print(
        "Width(p)             =",
        p_hull[1] - p_hull[0]
    )

    print(
        "Width(residual)      =",
        residual_hull[1] - residual_hull[0]
    )

