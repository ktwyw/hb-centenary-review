"""Fig. S3 of the supplementary material (figS3_pipe.pdf): Supplementary pipe calculations.

(a) Pressure-gradient consequences of the three rounded HB triplets of
    Magnon and Cayeux, Fluids 6, 157 (2021), doi:10.3390/fluids6040157,
    in their 15.5-mm tube (exact no-slip laminar pipe solution).
(b) Synthetic repeated experiments: joint HB fits from three calibration
    designs propagated to pipe pressure. The broad and high-rate designs
    use seed 20260920 (draws in a fixed order); the third
    design truncates the usable rates at 1 s^-1, as for a material whose
    critical shear rate is about 1 s^-1 so that steady data below it are
    unavailable.
Outputs: figures/figS3_pipe.pdf, data/pipe_supplement_results.json and the csv files of
fitted triplets and pressure ratios (data/*_window_fits.csv, data/*_window_pressure.csv).
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares
from hbstyle import C, plt, save, label_panels
from hb_prediction import pressure_gradient, pipe_flux

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'

# ---------------- (a) published triplets
sets = [('Rotational', 1.198, .2717, .6389),
        ('Pipe, nonlinear fit', 1.030, .2560, .6490),
        ('Pipe, corrected method', .909, .2720, .6410)]
R_TUBE = .0155 / 2
flow_lmin = np.logspace(-2, 1, 120)
q_tube = flow_lmin * 1e-3 / 60
grads = [np.array([pressure_gradient(v, R_TUBE, ty, K, n) for v in q_tube]) for _, ty, K, n in sets]
rel = [100 * (gr / grads[0] - 1) for gr in grads]
wall_ref = grads[0] * R_TUBE / 2


def apparent_wall_rate(q, r):
    return 4 * q / (np.pi * r**3)


pub = {'source_doi': '10.3390/fluids6040157',
       'sets': [dict(zip(['method', 'yield_stress_Pa', 'K_Pa_s_n', 'n'], r)) for r in sets],
       'yield_differences_percent': [100 * (r[1] / sets[0][1] - 1) for r in sets],
       'pressure_difference_percent_at_Q_L_min': {}}
for qv in (0.01, 0.1, 1.0, 10.0):
    i = int(np.argmin(abs(np.log(flow_lmin / qv))))
    pub['pressure_difference_percent_at_Q_L_min'][f'{qv:g}'] = {
        'Q_L_min': float(flow_lmin[i]),
        'apparent_wall_rate_s_inv': float(apparent_wall_rate(q_tube[i], R_TUBE)),
        'nonlinear_fit': float(rel[1][i]), 'corrected_method': float(rel[2][i])}
pub['max_abs_pressure_difference_percent'] = {'nonlinear_fit': float(np.max(abs(rel[1]))),
                                              'corrected_method': float(np.max(abs(rel[2])))}

# ---------------- (b) synthetic repeated experiments
rng = np.random.default_rng(20260920)
truth = np.array([10., 5., .5])   # sigma_y, A at 1 s^-1, n
noise_sd = 2.
repeats = 240
designs = {'Broad window': np.logspace(-1, 3, 24),
           'High-rate window': np.logspace(2, 3, 24),
           'Cutoff at 1 s$^{-1}$': np.logspace(0, 3, 24)}
radius = .01
wall_rates = np.logspace(-2, 3, 55)
true_w = truth[0] + truth[1] * wall_rates**truth[2]
flow = pipe_flux(true_w, radius, *truth)
true_grad = 2 * true_w / radius
summary = {}
bands = {}
for name, rates in designs.items():
    fits = []
    for k in range(repeats):
        y = truth[0] + truth[1] * rates**truth[2] + rng.normal(0, noise_sd, len(rates))
        residual = lambda t: (t[0] + t[1] * rates**t[2] - y) / noise_sd
        cands = [least_squares(residual, [10., 5., n0], bounds=([0., 1e-8, .1], [200., 200., 1.4]),
                               xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=1500) for n0 in (.25, .5, .85)]
        conv = [r for r in cands if r.success]
        if not conv:
            raise RuntimeError(f'No converged fit: {name} {k}')
        fits.append(min(conv, key=lambda r: r.cost).x)
    fits = np.array(fits)
    G = np.array([[pressure_gradient(q, radius, *t) for q in flow] for t in fits])
    r = G / true_grad
    lo, med, hi = np.quantile(r, [.025, .5, .975], axis=0)
    bands[name] = (lo, med, hi, fits)
    tag = {'Broad window': 'broad_window', 'High-rate window': 'high_rate_window',
           'Cutoff at 1 s$^{-1}$': 'cutoff_1s_window'}[name]
    np.savetxt(DATA / f'{tag}_fits.csv', fits, delimiter=',', header='sigma_y_Pa,A_Pa_at_1_per_s,n', comments='')
    np.savetxt(DATA / f'{tag}_pressure.csv', np.column_stack([flow * 6e4, wall_rates, lo, med, hi]), delimiter=',',
               header='Q_L_min,true_wall_rate_s_inv,q025_G_over_true,q50_G_over_true,q975_G_over_true', comments='')
    at = {}
    for wr in (0.01, 0.1, 1.0, 10.0, 1000.0):
        i = int(np.argmin(abs(np.log(wall_rates / wr))))
        at[f'{wr:g}'] = [float(lo[i]), float(hi[i])]
    summary[tag] = {'rate_min_s_inv': float(rates.min()), 'rate_max_s_inv': float(rates.max()),
                    'parameter_q025_q50_q975': np.quantile(fits, [.025, .5, .975], axis=0).tolist(),
                    'zero_yield_boundary_fits': int(np.sum(fits[:, 0] < 1e-5)),
                    'pressure_ratio_95_range_at_true_wall_rate': at}

# ---------------- figure
fig, ax = plt.subplots(1, 2, figsize=(6.8, 2.8))
for rr, (lab, *_), c, ls in zip(rel, sets, [C['black'], C['blue'], C['red']], ['-', '--', '-.']):
    ax[0].semilogx(flow_lmin, rr, ls, color=c, label=lab)
ax[0].set(xlabel=r'$Q$ (L min$^{-1}$)', ylabel=r'$\Delta G/G_{\rm rot}$ (%)')
ax[0].legend(loc='lower right', fontsize=8)
for (name, (lo, med, hi, _)), c, ls in zip(bands.items(), [C['blue'], C['red'], C['green']], ['-', '--', ':']):
    ax[1].loglog(wall_rates, lo, ls, color=c, lw=1.2, label=name, gid=f'{name}: lower bound of 95% range')
    ax[1].loglog(wall_rates, hi, ls, color=c, lw=1.2, gid=f'{name}: upper bound of 95% range')
ax[1].axhline(1, color='0.5', lw=0.6, gid='G/G_true = 1')
ax[1].set(xlabel=r'True wall shear rate (s$^{-1}$)', ylabel=r'$G/G_{\rm true}$, 95% range', ylim=(0.02, 4))
ax[1].legend(loc='lower right', fontsize=8)
label_panels(ax, x=-0.2)
fig.tight_layout()
save(fig, "figS3_pipe")
(DATA / 'pipe_supplement_results.json').write_text(json.dumps({
    'published_triplets': pub,
    'synthetic': {'seed': 20260920, 'repeats': repeats, 'observations_per_dataset': 24,
                  'truth_sigma_y_A_n': truth.tolist(), 'reference_rate_s_inv': 1.0,
                  'noise': 'independent Gaussian stress, SD 2 Pa', 'radius_m': radius,
                  'fit_bounds': [[0, 1e-8, .1], [200, 200, 1.4]], 'results': summary}}, indent=2))
print(json.dumps(pub, indent=2))
print(json.dumps(summary, indent=2))
