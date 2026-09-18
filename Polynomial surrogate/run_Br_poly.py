import numpy as np

from Affine_ArithmeticClassV3 import AffineArray

from BernsteinFiles import bernstein_approximation
from AA_evaluationBr import residual_eval_Bernstein_AA

from subintervalization import (
    split_interval,
    hull_intervals
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
# USER SETTINGS
# ============================================================

# Selected benchmark problem

# fun_numpy = Eggholder_numpy
# fun_AA = test_EggholderAA

fun_numpy = Ackley_numpy
fun_AA = test_AckleyAA

# Polynomial degree in each dimension
nvec = [3, 3]

# Bernstein interpolation nodes:
#
# "bernstein"  -> Bernstein nodes
# "chebyshev"  -> Chebyshev-Lobatto nodes
# "legendre"   -> Legendre-Lobatto nodes
#
node_type = "chebyshev"

# Whether AffineArray special operations should use
# the Chebyshev-enhanced evaluation.
cheb = False

# Domain subdivision levels
split_values = [1,2,4,8,16,32,64,128,256]
# split_values = [1,2,4,8,16,32]


# ============================================================
# TEST ONE SUBINTERVALIZATION
# ============================================================

def test_subintervalization(
    func_numpy,
    func_AA,
    splits_per_dim,
    nvec,
    node_type="bernstein",
    cheb=False
):
    """
    Construct a Bernstein polynomial approximation on every
    subdomain and evaluate the function, polynomial, and
    residual using Affine Arithmetic.

    Parameters
    ----------
    func_numpy : callable
        Numerical version of the function.

    func_AA : callable
        AffineArray version of the function.

    splits_per_dim : list or tuple
        Number of subdivisions in each dimension.

    nvec : list or tuple
        Bernstein polynomial degree in each dimension.

    node_type : str
        Interpolation node family.

    cheb : bool
        Passed to AffineArray operations.

    Returns
    -------
    results : list of dict
        Results for every subdomain.
    """

    # --------------------------------------------------------
    # Generate all subdomains
    # --------------------------------------------------------

    subdomains = split_interval(
        splits_per_dim
    )

    results = []

    # --------------------------------------------------------
    # Process every subdomain
    # --------------------------------------------------------

    for combo in subdomains:

        # ----------------------------------------------------
        # Affine variables for the current subdomain
        # ----------------------------------------------------

        X_sub = AffineArray.from_intervals(
            combo
        )

        # ----------------------------------------------------
        # Construct Bernstein approximation
        # ----------------------------------------------------

        (
            nodes,
            points,
            F_grid,
            C
        ) = bernstein_approximation(
            func_numpy,
            nvec,
            combo,
            node_type=node_type
        )

        # ----------------------------------------------------
        # Evaluate f, p, and residual using AA
        # ----------------------------------------------------

        AA_results = residual_eval_Bernstein_AA(
            X_sub,
            C,
            combo,
            nvec,
            func_AA,
            cheb=cheb
        )

        # ----------------------------------------------------
        # Store results
        # ----------------------------------------------------

        results.append(
            {
                "combo": combo,

                "nodes": nodes,

                "points": points,

                "F_grid": F_grid,

                "coefficients": C,

                "f_interval":
                    AA_results["f_AA"].interval,

                "p_Bernstein_coeff":
                    AA_results[
                        "p_Bernstein_coeff"
                    ],

                "p_AA":
                    AA_results[
                        "p_AA"
                    ].interval,

                "residual_AA":
                    AA_results[
                        "residual_AA"
                    ].interval,

                "residual_Bernstein":
                    AA_results[
                        "residual_Bernstein"
                    ]
            }
        )

    return results


# ============================================================
# RUN BERNSTEIN EXPERIMENT
# ============================================================

def run_experiment(
    func_numpy,
    func_AA,
    nvec,
    split_values,
    node_type="bernstein",
    cheb=False
):
    """
    Run the Bernstein approximation experiment for all
    requested subintervalization levels.

    Returns
    -------
    all_results : dict
        Results indexed by number of splits.
    """

    all_results = {}

    for splits in split_values:

        results = test_subintervalization(
            func_numpy,
            func_AA,
            [splits, splits],
            nvec,
            node_type=node_type,
            cheb=cheb
        )

        all_results[splits] = results

    return all_results


# ============================================================
# PRINT RESULTS
# ============================================================

def print_results(
    all_results,
    split_values
):
    """
    Print global interval hulls and widths.
    """

    for splits in split_values:

        results = all_results[splits]

        # ----------------------------------------------------
        # Collect intervals
        # ----------------------------------------------------

        f_intervals = [
            r["f_interval"]
            for r in results
        ]

        p_coeff_intervals = [
            r["p_Bernstein_coeff"]
            for r in results
        ]

        p_AA_intervals = [
            r["p_AA"]
            for r in results
        ]

        residual_AA_intervals = [
            r["residual_AA"]
            for r in results
        ]

        residual_Bernstein_intervals = [
            r["residual_Bernstein"]
            for r in results
        ]

        # ----------------------------------------------------
        # Global hulls
        # ----------------------------------------------------

        f_hull = hull_intervals(
            f_intervals
        )

        p_coeff_hull = hull_intervals(
            p_coeff_intervals
        )

        p_AA_hull = hull_intervals(
            p_AA_intervals
        )

        residual_AA_hull = hull_intervals(
            residual_AA_intervals
        )

        residual_Bernstein_hull = hull_intervals(
            residual_Bernstein_intervals
        )

        # ----------------------------------------------------
        # Widths
        # ----------------------------------------------------

        f_width = (
            f_hull[1]
            - f_hull[0]
        )

        p_coeff_width = (
            p_coeff_hull[1]
            - p_coeff_hull[0]
        )

        p_AA_width = (
            p_AA_hull[1]
            - p_AA_hull[0]
        )

        residual_AA_width = (
            residual_AA_hull[1]
            - residual_AA_hull[0]
        )

        residual_Bernstein_width = (
            residual_Bernstein_hull[1]
            - residual_Bernstein_hull[0]
        )

        # ----------------------------------------------------
        # Print
        # ----------------------------------------------------

        print(
            f"\nSplits = {splits} x {splits}"
        )

        print(
            "Number of subdomains =",
            len(results)
        )

        print(
            "Global f hull              =",
            f_hull
        )

        print(
            "Global Bernstein coeff hull =",
            p_coeff_hull
        )

        print(
            "Global p_AA hull            =",
            p_AA_hull
        )

        print(
            "Global residual AA hull     =",
            residual_AA_hull
        )

        print(
            "Global residual Bernstein hull =",
            residual_Bernstein_hull
        )

        print(
            "Width(f)                    =",
            f_width
        )

        print(
            "Width(Bernstein coeff)      =",
            p_coeff_width
        )

        print(
            "Width(p_AA)                 =",
            p_AA_width
        )

        print(
            "Width(residual AA)          =",
            residual_AA_width
        )

        print(
            "Width(residual Bernstein)   =",
            residual_Bernstein_width
        )


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    print(
        "\nBernstein Polynomial Approximation"
    )


    print(
        "Polynomial degree =",
        nvec
    )

    print(
        "Node type         =",
        node_type
    )

    print(
        "AA Chebyshev mode =",
        cheb
    )

    # --------------------------------------------------------
    # Run experiment
    # --------------------------------------------------------

    all_results = run_experiment(
        fun_numpy,
        fun_AA,
        nvec,
        split_values,
        node_type=node_type,
        cheb=cheb
    )

    # --------------------------------------------------------
    # Print results
    # --------------------------------------------------------

    print_results(
        all_results,
        split_values
    )
