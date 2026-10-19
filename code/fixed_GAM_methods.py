"""Unconstrained natural cubic smoothing spline with fixed total EDF."""
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline, make_smoothing_spline
from scipy.linalg import eigh
from scipy.optimize import brentq


def fixed_gam(age,pt,methods,total_edf=6):
    lo,span=float(age.min()),float(np.ptp(age))
    u=(age-lo)/span
    knots,inverse,counts=np.unique(u,return_inverse=True,return_counts=True)
    means=np.bincount(inverse,weights=pt)/counts
    basis=CubicSpline(knots,np.eye(len(knots)),bc_type="natural",extrapolate=False)
    penalty=methods.penalty_matrix(basis,knots)
    invroot=1/np.sqrt(counts)
    eigen,vectors=eigh(invroot[:,None]*penalty*invroot[None,:])
    tolerance=1e-10*max(1,eigen[-1])
    assert eigen[0]>=-tolerance and np.max(np.abs(eigen[:2]))<tolerance and eigen[2]>0
    eigen[:2]=0
    z=vectors.T@(np.sqrt(counts)*means)
    lam=10**brentq(lambda loglam:np.sum(1/(1+10**loglam*eigen))-total_edf,-16,12,xtol=1e-12)
    shrink=1/(1+lam*eigen)
    edf=float(shrink.sum())
    np.testing.assert_allclose(edf,total_edf,atol=1e-8)
    coef=invroot*(vectors@(shrink*z))
    grid=np.linspace(lo,lo+span,1001)
    predicted=basis((grid-lo)/span)@coef
    fitted=coef[inverse]
    independent=make_smoothing_spline(knots,means,w=counts.astype(float),lam=lam)
    independent_pred=independent((grid-lo)/span)
    np.testing.assert_allclose(predicted,independent_pred,atol=3e-7,rtol=2e-6)
    np.testing.assert_allclose(fitted,independent(u),atol=3e-7,rtol=2e-6)
    report={"method":"Gaussian identity-link GAM / natural cubic smoothing spline",
        "effective_degrees_of_freedom":edf,"lambda_normalized_age":float(lam),
        "complexity_choice":"Same user-requested fixed total EDF=6 as earlier figures; exploratory, not GCV-selected",
        "knots_normalized":knots.tolist(),"coefficients":coef.tolist(),
        "age_center":lo,"age_scale":span,"RSS":float(np.sum((pt-fitted)**2)),
        "independent_scipy_max_abs_difference":float(np.max(np.abs(predicted-independent_pred))),
        "no_monotonicity_convexity_or_concavity_constraint":True,"no_prediction_clipping":True,
        "no_CI_or_confirmatory_nonlinearity_test":True,"no_aging_rate_inference":True}
    return pd.DataFrame({"age":grid,"fitted_pseudotime":predicted}),report
