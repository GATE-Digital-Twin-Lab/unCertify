
import numpy as np
from itertools import product

def split_interval(splits_per_dim):

    d = len(splits_per_dim)

    # Original domain: [0, 1]^d
    dim_intervals = []

    for i in range(d):

        grid = np.linspace(
            0.0,
            1.0,
            splits_per_dim[i] + 1
        )

        intervals = [
            (grid[j], grid[j + 1])
            for j in range(splits_per_dim[i])
        ]

        dim_intervals.append(intervals)

    # Cartesian product of intervals
    return list(product(*dim_intervals))


def interval_width(x):
    lo, hi = x.interval
    return hi - lo

def hull_intervals(intervals):

    lower = min(
        lo for lo, hi in intervals
    )

    upper = max(
        hi for lo, hi in intervals
    )

    return (lower, upper)