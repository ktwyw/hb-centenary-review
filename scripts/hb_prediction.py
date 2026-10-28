"""Exact, no-slip, laminar circular-pipe forward model in SI units."""
import numpy as np
from scipy.optimize import brentq


def pipe_flux(wall_stress, radius, yield_stress, consistency, index):
    """Positive-direction volume flow [m^3/s] at nonnegative wall stress [Pa]."""
    if (not np.all(np.isfinite([radius,yield_stress,consistency,index]))
            or min(radius, consistency, index) <= 0 or yield_stress < 0):
        raise ValueError('Invalid HB parameters or radius')
    w = np.asarray(wall_stress, dtype=float)
    if np.any(~np.isfinite(w)) or np.any(w < 0):
        raise ValueError('Finite nonnegative wall stress is required')
    flowing = w > yield_stress
    d = np.where(flowing, w-yield_stress, 0.)
    den = np.where(flowing, w, 1.)
    # Dimensionless factors avoid powers of large dimensional stresses.
    phi = np.where(flowing, yield_stress/den, 0.)
    z = d/den
    q = np.pi*radius**3*(d/consistency)**(1/index)*z*(
        z*z/(3+1/index)+2*phi*z/(2+1/index)+phi*phi/(1+1/index))
    return np.where(flowing, q, 0.)


def pressure_gradient(flow, radius, yield_stress, consistency, index):
    """Required positive pressure gradient [Pa/m] for flow [m^3/s]."""
    if not np.isfinite(flow) or flow <= 0:
        raise ValueError('Positive imposed flow is required')
    high = max(1.,yield_stress*2)
    while pipe_flux(high,radius,yield_stress,consistency,index) < flow:
        high *= 2
    w=brentq(lambda w:pipe_flux(w,radius,yield_stress,consistency,index)-flow,
             yield_stress,high,xtol=1e-11,rtol=1e-12)
    return 2*w/radius
