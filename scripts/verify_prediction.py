"""Independent checks of the SI-unit pipe-flow helper (hb_prediction.py) used by figS3_pipe_uncertainty.py.

Uses numerical integration of the local constitutive law and classical limits;
this is mathematical verification, not experimental validation.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from hb_prediction import pipe_flux, pressure_gradient


def main():
    comparisons = 0
    worst = 0.
    for n in (.35, .7, 1., 1.4):
        for radius in (.00775, .01):
            ty, consistency = 6., .8
            for ratio in (1.001, 2., 10.):
                wall = ratio*ty
                plug = radius*ty/wall
                # Independent radial moment of the local shear rate.
                numerical = np.pi*quad(
                    lambda r: ((wall*r/radius-ty)/consistency)**(1/n)*r**2,
                    plug, radius, epsabs=1e-28, epsrel=1e-10)[0]
                exact = float(pipe_flux(wall, radius, ty, consistency, n))
                error = abs(exact/numerical-1.)
                worst = max(worst, error)
                np.testing.assert_allclose(exact, numerical, rtol=2e-8, atol=0.)
                recovered = pressure_gradient(numerical, radius, ty, consistency, n)
                np.testing.assert_allclose(recovered, 2*wall/radius, rtol=2e-9)
                comparisons += 1
            with np.errstate(all='raise'):
                stopped = pipe_flux(np.array([0., ty/2, ty]), radius, ty, consistency, n)
            np.testing.assert_array_equal(stopped, np.zeros(3))

    for radius in (.001, .01, .1):
        for gradient in (10., 1000.):
            viscosity = .8
            reference = np.pi*radius**4*gradient/(8*viscosity)
            np.testing.assert_allclose(pipe_flux(gradient*radius/2, radius, 0., viscosity, 1.),
                                       reference, rtol=1e-12, atol=0.)
            np.testing.assert_allclose(pressure_gradient(reference, radius, 0., viscosity, 1.),
                                       gradient, rtol=1e-9)

    result = {'status': 'PASS', 'direct_radial_quadratures': comparisons,
              'independent_flux_to_gradient_checks': comparisons,
              'Newtonian_forward_and_inverse_cases': 6,
              'unyielded_vectors_checked': 8,
              'maximum_quadrature_relative_error': worst,
              'scope': 'Positive-direction, laminar, no-slip circular pipe in SI units.'}
    root = Path(__file__).resolve().parent.parent
    (root/'results').mkdir(exist_ok=True)
    (root/'results'/'prediction_verification.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
