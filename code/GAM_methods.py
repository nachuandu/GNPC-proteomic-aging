"""Natural cubic spline roughness penalty used by the fixed-EDF fit."""
import numpy as np


def penalty_matrix(basis, knots):
    # Two-point Gauss-Legendre is exact for squared piecewise-linear f''.
    midpoint = (knots[1:] + knots[:-1]) / 2
    halfwidth = np.diff(knots) / 2
    nodes = np.column_stack([midpoint - halfwidth / np.sqrt(3), midpoint + halfwidth / np.sqrt(3)]).ravel()
    weights = np.repeat(halfwidth, 2)
    d2 = basis(nodes, 2)
    return d2.T @ (weights[:, None] * d2)
