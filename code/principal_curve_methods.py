"""Custom single principal curve, not the full Slingshot algorithm.

Numerical functions extracted from the working trajectory implementation.
"""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.stats import spearmanr


def unit_range(x):
    span = np.ptp(x)
    if not np.isfinite(span) or span <= 0:
        raise ValueError('Degenerate ordering.')
    return (x - np.min(x)) / span


def project_polyline(points, curve):
    """Exact nearest segment projections on a piecewise-linear curve.

    Distances use the input geometry; both coordinates have the same weight.
    Arc length is continuous within segments (not a point-index rank).
    """
    delta = np.diff(curve, axis=0)
    lengths = np.linalg.norm(delta, axis=1)
    keep = np.r_[True, lengths > 1e-14]
    curve = curve[keep]
    if len(curve) < 2:
        raise ValueError('Degenerate principal curve.')
    starts = curve[:-1]
    delta = np.diff(curve, axis=0)
    squared_lengths = np.einsum('ij,ij->i', delta, delta)
    lengths = np.sqrt(squared_lengths)
    offsets = np.r_[0., np.cumsum(lengths)]
    arc = np.empty(len(points))
    projections = np.empty_like(points)
    d2 = np.empty(len(points))
    for left in range(0, len(points), 256):
        p = points[left:left+256]
        displacement = p[:, None, :] - starts[None, :, :]
        fractions = np.einsum('nsi,si->ns', displacement, delta) / squared_lengths
        fractions = np.clip(fractions, 0, 1)
        residual = displacement - fractions[:, :, None] * delta[None, :, :]
        distances = np.einsum('nsi,nsi->ns', residual, residual)
        nearest = np.argmin(distances, axis=1)
        row = np.arange(len(p))
        frac = fractions[row, nearest]
        arc[left:left+len(p)] = offsets[nearest] + frac * lengths[nearest]
        projections[left:left+len(p)] = starts[nearest] + frac[:, None] * delta[nearest]
        d2[left:left+len(p)] = distances[row, nearest]
    return arc, projections, d2, curve


def natural_spline_fit(x, Y, n_knots):
    """Fixed-complexity natural cubic regression spline; no age-based tuning."""
    knots = np.unique(np.quantile(x, np.linspace(0, 1, n_knots)))
    if len(knots) < 4:
        raise ValueError('Too few unique knot locations.')
    basis = CubicSpline(knots, np.eye(len(knots)), bc_type='natural', extrapolate=False)
    design = basis(x)
    coefficients, _, rank, _ = np.linalg.lstsq(design, Y, rcond=None)
    if rank != len(knots):
        raise ValueError('Rank-deficient spline design.')
    return basis, coefficients, knots


def fit_principal_curve(xy, n_knots=6, n_grid=1001, max_iter=100, tolerance=1e-6):
    """Alternating spline smoothing and orthogonal projection, PCA initialized.

    A common centering/scaling constant is used for numerical conditioning.
    Independent axis scaling is forbidden: native distances are preserved.
    Age and contributor are not arguments and cannot influence the curve.
    """
    xy = np.asarray(xy, dtype=np.float64)
    if xy.ndim != 2 or xy.shape[1] != 2 or not np.isfinite(xy).all():
        raise ValueError('Expected a finite n x 2 coordinate matrix.')
    center = xy.mean(axis=0)
    scale = float(np.max(np.ptp(xy, axis=0)))
    if scale <= 0:
        raise ValueError('No coordinate variation.')
    Z = (xy-center)/scale
    _, _, vh = np.linalg.svd(Z, full_matrices=False)
    direction = vh[0]
    dominant = int(np.argmax(np.abs(direction)))
    if direction[dominant] < 0:
        direction = -direction
    score = Z @ direction
    u = unit_range(score)
    previous_projection = score[:, None] * direction
    history = []
    converged = False
    consecutive = 0
    for iteration in range(1, max_iter+1):
        basis, coefficients, knots = natural_spline_fit(u, Z, n_knots)
        grid = np.linspace(u.min(), u.max(), n_grid)
        curve = basis(grid) @ coefficients
        arc, projection, d2, curve = project_polyline(Z, curve)
        new_u = unit_range(arc)
        max_change = float(np.max(np.abs(new_u-u)))
        projection_change = float(np.max(np.linalg.norm(projection-previous_projection, axis=1)))
        history.append(dict(iteration=iteration, mse_scaled=float(d2.mean()),
                            max_ptime_change=max_change,
                            max_projection_change_scaled=projection_change))
        print(f'Principal curve iteration {iteration}: delta_pt={max_change:.3g}, delta_projection={projection_change:.3g}', flush=True)
        if max_change < tolerance and projection_change < tolerance:
            consecutive += 1
        else:
            consecutive = 0
        u = new_u
        previous_projection = projection
        if consecutive >= 3:
            converged = True
            break
    native_curve = curve*scale+center
    native_projection = projection*scale+center
    # Independent re-projection verifies that saved ptime really is native arc length.
    raw_arc, saved_projection, saved_d2, native_curve = project_polyline(xy, native_curve)
    np.testing.assert_allclose(saved_projection, native_projection, atol=scale*1e-9)
    raw_pt = unit_range(raw_arc)
    curve_arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(native_curve, axis=0), axis=1))]
    curve_pt = (curve_arc-raw_arc.min())/np.ptp(raw_arc)
    return raw_pt, native_curve, saved_projection, np.sqrt(saved_d2), curve_pt, {
        'converged': converged, 'iterations': len(history), 'max_iterations': max_iter,
        'tolerance': tolerance, 'history': history, 'spline_knots': knots.tolist(),
        'number_of_knots': n_knots, 'evaluation_grid_points': n_grid,
        'common_coordinate_scale': scale, 'coordinate_center': center.tolist(),
        'native_curve_total_length': float(curve_arc[-1]),
        'native_observed_arc_span': float(np.ptp(raw_arc))}


def orient(raw_pt, ages):
    rho = float(spearmanr(ages, raw_pt).statistic)
    if not np.isfinite(rho):
        raise ValueError('Cannot orient a constant/uninformative ptime.')
    flip = rho < 0
    return 1-raw_pt if flip else raw_pt.copy(), flip, rho


def test_algorithm():
    # Curved synthetic data: ordering, row-permutation, scale and rotation checks.
    rng = np.random.default_rng(91)
    s = np.linspace(-1, 1, 220)
    xy = np.column_stack([s, .35*s*s]) + rng.normal(0, .008, size=(len(s), 2))
    p, curve, projected, distance, cp, report = fit_principal_curve(xy)
    assert report['converged'], report
    assert abs(spearmanr(p, s).statistic) > .99
    assert p.min() == 0 and p.max() == 1 and np.isfinite(distance).all()
    order = rng.permutation(len(s))
    p2 = fit_principal_curve(xy[order])[0]
    np.testing.assert_allclose(p[order], p2, atol=1e-6)
    theta = .71
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    p3 = fit_principal_curve((xy @ rotation)*3.7+10)[0]
    assert min(np.max(np.abs(p-p3)), np.max(np.abs(p-(1-p3)))) < 1e-5
    arc, pr, dist, _ = project_polyline(np.array([[.2, 2.], [1.8, -1.]]), np.array([[0., 0.], [1., 0.], [2., 0.]]))
    np.testing.assert_allclose(arc, [.2, 1.8])
    np.testing.assert_allclose(dist, [4., 1.])
    near_line = np.column_stack([s, 1e-5*np.sin(s*3)])
    pn, cn, _, _, _, rn = fit_principal_curve(near_line)
    assert rn['converged'] and abs(spearmanr(pn, s).statistic) > .999
    # A unit test of orientation does not feed age back into curve fitting.
    oriented, flip, _ = orient(np.linspace(0, 1, 20), np.arange(20)[::-1])
    assert flip and np.all(np.diff(oriented) < 0)
    print('SELF-TEST PASSED: continuous segment projection; row-order, rigid-rotation and common-scale invariance; near-line behavior.', flush=True)
