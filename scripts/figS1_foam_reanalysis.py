"""Fig. S1 of the supplementary material (figS1_foam.pdf): Supplementary reanalysis of the SLES/betaine foam local flow curve.

Data: 152 marker centers digitized from Fig. 3b of Ovarlez, Krishan and
Cohen-Addad, EPL 91, 68005 (2010) (data/foam_sles_fig3b_digitized.csv).

Error model (errors in variables). An observed local rate g_i and stress s_i
arise from a latent true rate x_i on the model branch sigma = f(x):
    g_i = x_i + e_g,  Var(e_g) = (delta0*x)^2 + eps_g^2   (relative error plus
                                                         an absolute rate floor
                                                         set by MRI resolution)
    s_i = f(x_i) + e_s, Var(e_s) = vs  (vs estimated by maximum likelihood)
For the untruncated models the latent rate is eliminated with the standard
effective-variance approximation, V_i = vs + f'(g_i)^2 Var(e_g)(g_i).
For the critical-shear-rate model the steady branch exists only for x >= gc,
so the latent rate is profiled numerically on a grid restricted to x >= gc.

Models: Herschel-Bulkley (HB), Bingham, Casson and HB-gc (HB branch with
steady flow impossible below a critical shear rate gc; sigma_y is then the
extrapolated intercept and sigma_c = sigma_y + K gc^n the stopping stress).
Outputs: data/foam_reanalysis_results.json, figures/figS_foam.pdf.
"""
from pathlib import Path
import csv
import json
import numpy as np
from scipy.optimize import minimize
from hbstyle import plt, save, C, label_panels

ROOT = Path(__file__).resolve().parents[1]
DELTA0, EPS_G = 0.05, 0.01          # baseline rate-error model
CHI95 = 3.841458820694124            # chi^2_1 0.95 quantile
SOURCE_FIT = (32.0, 12.0, 0.46)      # Ovarlez et al. (2010), Fig. 3 caption
SOURCE_GC_BOUND = 0.04               # upper bound from local velocity profiles
DENKOV_GC = 1.5                      # Denkov et al. prediction quoted by source (Table 1)


def load():
    with (ROOT / 'data/foam_sles_fig3b_digitized.csv').open() as f:
        rows = list(csv.DictReader(f))
    g = np.array([float(r['shear_rate_s_inv']) for r in rows])
    s = np.array([float(r['stress_Pa']) for r in rows])
    grp = np.array([r['group'] for r in rows])
    return g, s, grp


def f_model(name, x, p):
    if name == 'HB':
        return p[0] + p[1] * x**p[2]
    if name == 'Bingham':
        return p[0] + p[1] * x
    if name == 'Casson':
        return (np.sqrt(p[0]) + np.sqrt(p[1] * x))**2
    raise ValueError(name)


def df_model(name, x, p):
    if name == 'HB':
        return p[1] * p[2] * x**(p[2] - 1)
    if name == 'Bingham':
        return p[1] * np.ones_like(x)
    if name == 'Casson':
        return (np.sqrt(p[0]) + np.sqrt(p[1] * x)) * np.sqrt(p[1] / x)
    raise ValueError(name)


def var_rate(x, delta0=DELTA0, eps=EPS_G):
    return (delta0 * x)**2 + eps**2


def nll_ev(name, theta, g, s, delta0=DELTA0, eps=EPS_G):
    """Negative log-likelihood, effective-variance approximation.
    theta = model parameters followed by log(stress SD)."""
    p, ls = theta[:-1], theta[-1]
    if name == 'HB' and not (p[0] >= 0 and p[1] > 0 and 0.02 < p[2] < 1.5):
        return 1e12
    if name in ('Bingham', 'Casson') and not (p[0] >= 0 and p[1] > 0):
        return 1e12
    V = np.exp(2 * ls) + df_model(name, g, p)**2 * var_rate(g, delta0, eps)
    r = s - f_model(name, g, p)
    return 0.5 * np.sum(r * r / V + np.log(2 * np.pi * V))


START = {'HB': [[30, 16, 0.36, 0.5], [32, 12, 0.46, 0.5], [25, 20, 0.3, 0.5]],
         'Bingham': [[40, 2, 0.5], [30, 5, 1.0]],
         'Casson': [[30, 5, 0.5], [20, 20, 1.0]]}


def fit_ev(name, g, s, starts=None, fixed=None, **kw):
    """Fit by Nelder-Mead from several starts. `fixed` = (index, value) holds one
    model parameter fixed for profile likelihood."""
    starts = starts or START[name]
    best = None
    for st in starts:
        st = np.array(st, float)
        if fixed is not None:
            i, v = fixed
            free = [k for k in range(len(st)) if k != i]
            def obj(th, free=free, i=i, v=v):
                full = np.empty(len(free) + 1); full[free] = th; full[i] = v
                return nll_ev(name, full, g, s, **kw)
            x0 = st[free]
        else:
            obj = lambda th: nll_ev(name, th, g, s, **kw)
            x0 = st
        res = minimize(obj, x0, method='Nelder-Mead',
                       options={'xatol': 1e-9, 'fatol': 1e-10, 'maxiter': 40000, 'maxfev': 40000})
        res = minimize(obj, res.x, method='Nelder-Mead',
                       options={'xatol': 1e-10, 'fatol': 1e-11, 'maxiter': 40000, 'maxfev': 40000})
        if best is None or res.fun < best.fun:
            best = res
            if fixed is not None:
                full = np.empty(len(st)); full[free] = res.x; full[i] = v
                best_full = full
            else:
                best_full = res.x
    return best_full, best.fun


# ---------------- critical-shear-rate model (grid-profiled latent rate) -------------
XGRID = np.logspace(-4, 2, 900)


def nll_gc(theta, g, s, gc, delta0=DELTA0, eps=EPS_G):
    sy, K, n, ls = theta
    if not (sy >= 0 and K > 0 and 0.02 < n < 1.5):
        return 1e12
    x = XGRID[XGRID >= gc]
    vg = var_rate(x, delta0, eps)
    vs = np.exp(2 * ls)
    fx = sy + K * x**n
    ll = (-0.5 * (g[:, None] - x[None, :])**2 / vg[None, :] - 0.5 * np.log(2 * np.pi * vg)[None, :]
          - 0.5 * (s[:, None] - fx[None, :])**2 / vs - 0.5 * np.log(2 * np.pi * vs))
    return -np.sum(ll.max(axis=1))


def fit_gc(g, s, gc, start, **kw):
    obj = lambda th: nll_gc(th, g, s, gc, **kw)
    res = minimize(obj, start, method='Nelder-Mead',
                   options={'xatol': 1e-7, 'fatol': 1e-8, 'maxiter': 6000, 'maxfev': 6000})
    res = minimize(obj, res.x, method='Nelder-Mead',
                   options={'xatol': 1e-8, 'fatol': 1e-9, 'maxiter': 6000, 'maxfev': 6000})
    return res.x, res.fun


def gc_profile(g, s, gcs, start, **kw):
    out = []
    st = np.array(start, float)
    for gc in gcs:
        th, f = fit_gc(g, s, gc, st, **kw)
        out.append((gc, f, th))
        st = th
    return out


def upper_bound(gcs, dev):
    """Smallest gc at which the profile deviance first exceeds CHI95 (log interpolation)."""
    for k in range(1, len(gcs)):
        if dev[k] > CHI95 >= dev[k - 1]:
            a, b = np.log(gcs[k - 1]), np.log(gcs[k])
            t = (CHI95 - dev[k - 1]) / (dev[k] - dev[k - 1])
            return float(np.exp(a + t * (b - a)))
    return None


def profile_interval(name, idx, grid, g, s, best_nll, start):
    devs = []
    for v in grid:
        _, f = fit_ev(name, g, s, starts=[start], fixed=(idx, v))
        devs.append(2 * (f - best_nll))
    devs = np.array(devs)
    inside = grid[devs <= CHI95]
    return devs, (float(inside.min()), float(inside.max()))


def run():
    g, s, grp = load()
    N = len(g)
    results = {'data': {'n': int(N), 'series': int(len(np.unique(grp))),
                        'rate_min_s_inv': float(g.min()), 'rate_max_s_inv': float(g.max()),
                        'points_below_1.5_s_inv': int(np.sum(g < DENKOV_GC)),
                        'points_below_0.04_s_inv': int(np.sum(g < SOURCE_GC_BOUND))},
               'error_model': {'rate_relative_sd': DELTA0, 'rate_floor_sd_s_inv': EPS_G,
                               'stress_sd': 'maximum likelihood, common to all points'},
               'source': {'fit_sigma_y_K_n': SOURCE_FIT, 'gc_upper_bound_s_inv': SOURCE_GC_BOUND,
                          'denkov_predicted_gc_s_inv': DENKOV_GC}}
    # ---- untruncated models
    fits = {}
    for name in ['HB', 'Bingham', 'Casson']:
        th, f = fit_ev(name, g, s)
        k = len(th)
        rmse = float(np.sqrt(np.mean((f_model(name, g, th[:-1]) - s)**2)))
        fits[name] = th
        results[name] = {'parameters': th[:-1].tolist(), 'stress_sd_Pa': float(np.exp(th[-1])),
                         'nll': float(f), 'k': k, 'AIC': float(2 * f + 2 * k),
                         'BIC': float(2 * f + k * np.log(N)), 'rmse_Pa': rmse}
    hb = fits['HB']
    nll_hb = results['HB']['nll']
    # source fit evaluated under the same error model (stress SD re-optimized)
    def nll_src(ls):
        return nll_ev('HB', np.array([*SOURCE_FIT, ls[0]]), g, s)
    r = minimize(nll_src, [0.5], method='Nelder-Mead')
    results['source_fit_under_error_model'] = {
        'deviance_vs_HB_ML': float(2 * (r.fun - nll_hb)),
        'rmse_Pa': float(np.sqrt(np.mean((f_model('HB', g, SOURCE_FIT) - s)**2)))}
    # ---- profile intervals for sigma_y and n (independent-point likelihood)
    sy_grid = np.linspace(24, 37, 53)
    n_grid = np.linspace(0.29, 0.44, 61)
    dev_sy, ci_sy = profile_interval('HB', 0, sy_grid, g, s, nll_hb, hb)
    dev_n, ci_n = profile_interval('HB', 2, n_grid, g, s, nll_hb, hb)
    results['HB']['profile_95_sigma_y_Pa'] = ci_sy
    results['HB']['profile_95_n'] = ci_n
    # ---- joint (sigma_y, n) profile region for the figure
    SYg, NNg = np.meshgrid(np.linspace(24, 38, 57), np.linspace(0.28, 0.50, 56))
    joint = np.zeros_like(SYg)
    for i in range(SYg.shape[0]):
        for j in range(SYg.shape[1]):
            def obj(t, sy=SYg[i, j], nn=NNg[i, j]):
                return nll_ev('HB', np.array([sy, t[0], nn, t[1]]), g, s)
            # K start from linearization at the data median
            K0 = max(np.median(s - SYg[i, j]) / np.median(g)**NNg[i, j], 1e-3)
            rr = minimize(obj, [K0, hb[-1]], method='Nelder-Mead',
                          options={'xatol': 1e-8, 'fatol': 1e-9, 'maxiter': 4000})
            joint[i, j] = 2 * (rr.fun - nll_hb)
    # ---- series block bootstrap (seven acquisition series)
    rng = np.random.default_rng(20261025)
    labels = np.unique(grp)
    boot = []
    for _ in range(400):
        pick = rng.choice(labels, size=len(labels), replace=True)
        idx = np.concatenate([np.where(grp == q)[0] for q in pick])
        th, _ = fit_ev('HB', g[idx], s[idx], starts=[hb])
        boot.append(th[:-1])
    boot = np.array(boot)
    results['HB']['series_bootstrap_95_sigma_y_Pa'] = np.quantile(boot[:, 0], [.025, .975]).tolist()
    results['HB']['series_bootstrap_95_n'] = np.quantile(boot[:, 2], [.025, .975]).tolist()
    results['HB']['series_bootstrap_replicates'] = 400
    # ---- withheld-series prediction for HB (one-sentence note only)
    held = np.empty(N)
    for q in labels:
        tr = grp != q
        th, _ = fit_ev('HB', g[tr], s[tr], starts=[hb])
        held[~tr] = f_model('HB', g[~tr], th[:-1])
    results['HB']['withheld_series_rmse_Pa'] = float(np.sqrt(np.mean((held - s)**2)))
    # ---- critical shear rate profile
    start = hb.copy()
    th0, f0 = fit_gc(g, s, 1e-4, start)              # gc -> 0 reference on the same grid
    gcs = np.logspace(-3, np.log10(3), 40)
    prof = gc_profile(g, s, gcs, th0)
    fmin = min(f0, min(p[1] for p in prof))
    dev_gc = np.array([2 * (p[1] - fmin) for p in prof])
    ub = upper_bound(gcs, dev_gc)
    results['HB_gc'] = {'reference_nll_gc_to_0': float(f0),
                        'gc_grid_s_inv': gcs.tolist(), 'profile_deviance': dev_gc.tolist(),
                        'upper_bound_95_s_inv': ub}
    th_d, f_d = fit_gc(g, s, DENKOV_GC, th0)
    results['HB_gc']['deviance_at_denkov_1.5'] = float(2 * (f_d - fmin))
    results['HB_gc']['parameters_at_denkov_1.5'] = th_d[:3].tolist()
    # sensitivity of the bound to the rate-error model
    sens = []
    for d0, eps in [(0.05, 0.003), (0.05, 0.03), (0.02, 0.01), (0.10, 0.01)]:
        t0, ff0 = fit_gc(g, s, 1e-4, th0, delta0=d0, eps=eps)
        pr = gc_profile(g, s, gcs, t0, delta0=d0, eps=eps)
        fm = min(ff0, min(p[1] for p in pr))
        dv = np.array([2 * (p[1] - fm) for p in pr])
        _, fd = fit_gc(g, s, DENKOV_GC, t0, delta0=d0, eps=eps)
        sens.append({'rate_relative_sd': d0, 'rate_floor_sd_s_inv': eps,
                     'upper_bound_95_s_inv': upper_bound(gcs, dv),
                     'deviance_at_denkov_1.5': float(2 * (fd - fm)), 'deviance': dv.tolist()})
    results['HB_gc']['sensitivity'] = sens
    # rate floor needed before gc = 1.5 s^-1 becomes compatible (deviance < CHI95)
    floor_needed = None
    for eps in np.logspace(-2, 0.5, 26):
        t0, ff0 = fit_gc(g, s, 1e-4, th0, eps=eps)
        _, fd = fit_gc(g, s, DENKOV_GC, t0, eps=eps)
        if 2 * (fd - min(ff0, fd)) < CHI95:
            floor_needed = float(eps); break
    results['HB_gc']['rate_floor_sd_needed_for_1.5_s_inv'] = floor_needed
    results['HB_gc']['nll_min'] = float(fmin)
    # ---- figure
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.6))
    marks = ['o', 'D', 's', 'v', '<', '>', '^']
    for q, m in zip(labels, marks):
        t = grp == q
        ax[0].semilogx(g[t], s[t], m, ms=3, mfc='white', mec='0.35', mew=0.5,
                       gid=f'digitized points, marker group {q}')
    xx = np.logspace(-2.2, 1.6, 400)
    for name, c, ls in [('HB', C['black'], '-'), ('Bingham', C['blue'], '--'), ('Casson', C['red'], ':')]:
        ax[0].semilogx(xx, f_model(name, xx, fits[name][:-1]), color=c, ls=ls, label=name)
    ax[0].axvline(DENKOV_GC, color='0.45', lw=0.7, ls='-.', gid='predicted critical shear rate')
    ax[0].text(DENKOV_GC * 0.8, 88, r'predicted $\dot\gamma_c$', fontsize=8, ha='right', va='top', rotation=90, color='0.3')
    ax[0].set(xlim=(6e-3, 40), ylim=(20, 90), xlabel=r'$\dot\gamma$ (s$^{-1}$)', ylabel=r'$\sigma$ (Pa)')
    ax[0].legend(loc='upper left', fontsize=8)
    ax[1].semilogx(gcs, dev_gc, 'k-', label='0.01', gid='profile deviance, rate floor 0.01 1/s')
    for sv, ls in zip(sens[:2], ['--', ':']):
        ax[1].semilogx(gcs, sv['deviance'], color='k', ls=ls, label=f"{sv['rate_floor_sd_s_inv']:g}",
                       gid=f"profile deviance, rate floor {sv['rate_floor_sd_s_inv']:g} 1/s")
    ax[1].axhline(CHI95, color='0.5', lw=0.7, gid='95% threshold (chi-squared, 1 dof)')
    ax[1].axvline(SOURCE_GC_BOUND, color='0.5', lw=0.7, ls='-.', gid='bound on critical shear rate reported by the source')
    ax[1].set(ylim=(0, 40), xlim=(1e-3, 3), xlabel=r'Critical shear rate $\dot\gamma_c$ (s$^{-1}$)',
              ylabel='Profile deviance')
    ax[1].legend(loc='upper left', fontsize=8, title=r'$\epsilon_g$ (s$^{-1}$)', title_fontsize=8, handlelength=1.6); ax[1].set_yticks([10, 20, 30, 40])
    cs = ax[2].contour(SYg, NNg, joint, levels=[2.2789, 5.9915], colors='k', linestyles=['-', '--'], linewidths=0.8)
    ax[2].plot(hb[0], hb[2], 'k+', ms=6, gid='maximum-likelihood estimate (sigma_y, n)')
    ax[2].plot(boot[:, 0], boot[:, 2], '.', ms=1.5, color='0.6', zorder=0, gid='block-bootstrap fits (sigma_y, n)')
    ax[2].plot(SOURCE_FIT[0], SOURCE_FIT[2], 's', mfc='white', mec='k', ms=4, gid='fit reported by the source')
    ax[2].annotate('source fit', (SOURCE_FIT[0], SOURCE_FIT[2]), xytext=(-34, 6), textcoords='offset points', fontsize=8)
    ax[2].set(xlim=(24, 38), ylim=(0.28, 0.50), xlabel=r'$\sigma_y$ (Pa)', ylabel=r'$n$')
    fig.tight_layout(w_pad=2.2)
    label_panels(ax, x=-0.3)
    # panel (c): the joint profile-deviance surface whose contours are drawn
    fig.hb_extra = {"joint_deviance_grid": [{"sigma_y_Pa": SYg[i, j], "n": NNg[i, j], "profile_deviance": joint[i, j]}
                                            for i in range(SYg.shape[0]) for j in range(SYg.shape[1])]}
    save(fig, "figS1_foam")
    results['figure_notes'] = {'panel_c_contours': 'joint profile deviance 2.28 (68%) and 5.99 (95%), two parameters',
                               'panel_c_gray_points': 'series block-bootstrap HB fits'}
    (ROOT / 'data/foam_reanalysis_results.json').write_text(json.dumps(results, indent=2))
    short = {k: {kk: vv for kk, vv in v.items() if kk not in ('gc_grid_s_inv', 'profile_deviance', 'sensitivity')}
             if isinstance(v, dict) else v for k, v in results.items()}
    print(json.dumps(short, indent=2))
    print('sensitivity', [(x['rate_relative_sd'], x['rate_floor_sd_s_inv'], x['upper_bound_95_s_inv'],
                           round(x['deviance_at_denkov_1.5'], 1)) for x in sens])


if __name__ == '__main__':
    run()
