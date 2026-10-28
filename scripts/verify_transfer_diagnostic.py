"""Check the three-rate (triplet) estimates of n, K and sigma_y (Supplementary Note S8) and the
variance of the three-rate exponent against exact laws and finite differences.

Synthetic algebra checks only; no claim of experimental validation.
"""
from pathlib import Path
import json
import numpy as np


def exponent(stress, ratio):
    d = np.diff(stress)
    return np.log(d[1] / d[0]) / np.log(ratio)


def main():
    cases = 0
    for n in (0.2, 0.5, 1.0, 1.5):
        for ratio in (2.0, 4.0, 10.0):
            for g in (0.1, 1.0, 10.0):
                rates = g * ratio ** np.arange(3)
                stress = 10 + 5 * rates**n
                nhat = exponent(stress, ratio)
                d1 = stress[1] - stress[0]
                khat = d1 / (g**nhat * (ratio**nhat - 1))
                tyhat = stress[0] - d1 / (ratio**nhat - 1)
                np.testing.assert_allclose([nhat, khat, tyhat], [n, 5, 10], rtol=1e-10)
                np.testing.assert_allclose(exponent(stress + 35, ratio), n, rtol=1e-10)
                cases += 1

    # Independently differentiate n(s0,s1,s2) and propagate stress covariance.
    s = np.array([15., 20., 30.])
    h = 1e-4
    grad = np.array([(exponent(s + h * e, 4) - exponent(s - h * e, 4)) / (2*h)
                     for e in np.eye(3)])
    sigma = 0.1
    numerical_variance = sigma**2 * grad @ grad
    d1, d2 = np.diff(s)
    formula_variance = 2*sigma**2 / np.log(4)**2 * (1/d1**2 + 1/d2**2 + 1/(d1*d2))
    np.testing.assert_allclose(numerical_variance, formula_variance, rtol=1e-8)

    # A changing rate dependence must not pass merely because it is monotonic.
    rates = 4.**np.arange(6)
    crossover = 10 + 5*rates**0.5 + 0.05*rates
    estimates = [exponent(crossover[i:i+3], 4) for i in range(4)]
    assert max(estimates) - min(estimates) > 0.1
    result = {'status': 'PASS', 'exact_reconstruction_and_offset_cases': cases,
              'variance_finite_difference_check': 'PASS',
              'synthetic_crossover_detected': True,
              'illustrative_exponent_standard_uncertainty': float(np.sqrt(formula_variance)),
              'scope': 'Exact algebra and first-order propagation; not an experimental validation.'}
    root = Path(__file__).resolve().parent.parent
    (root / 'results').mkdir(exist_ok=True)
    (root / 'results' / 'transfer_verification.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
