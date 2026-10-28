"""Fig. 7 of the article (fig07_fitting.pdf): Ill-conditioning of the three-parameter HB fit (synthetic data)."""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from hbstyle import *

rng = np.random.default_rng(1)
ty, K, n = 10.0, 2.0, 0.5
gd = np.logspace(-2, 3, 26)
tau_meas = hb_stress(gd, ty, K, n) * (1 + 0.05 * rng.standard_normal(gd.size))
model = lambda g, a, b, c: a + b * g**c
bounds = ([0, 1e-6, 0.05], [50, 100, 2]); p0 = [5, 1, 0.7]
p_log, _ = curve_fit(lambda g, a, b, c: np.log(model(g, a, b, c)), gd, np.log(tau_meas), p0=p0, bounds=bounds)
mask = gd >= 10.0
p_high, _ = curve_fit(lambda g, a, b, c: np.log(model(g, a, b, c)), gd[mask], np.log(tau_meas[mask]), p0=p0, bounds=bounds)
p_lin, _ = curve_fit(model, gd, tau_meas, p0=p0, bounds=bounds)
fits = [(p_log, C["black"], "-", "(i) log, all rates"),
        (p_high, C["red"], "--", "(ii) log, high rates"),
        (p_lin, C["blue"], ":", "(iii) linear, all rates")]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=DOUBLE)
gfine = np.logspace(-2.2, 3.2, 300)
ax1.axvspan(1e-2 / 1.6, 10.0, color=C["lightgrey"], alpha=0.35, lw=0, gid="rates excluded in fit (ii)")
ax1.loglog(gd, tau_meas, "^", ms=4, mfc="none", mec="k", mew=0.7, label="synthetic data", zorder=5)
for p, col, ls, lab in fits:
    ax1.loglog(gfine, model(gfine, *p), color=col, ls=ls,
               label=lab)
ax1.text(0.011, 22, "excluded in (ii)", fontsize=9.5, color=C["grey"], va="bottom")
ax1.set_xlabel("$\\dot\\gamma$ (s$^{-1}$)"); ax1.set_ylabel("$\\sigma$ (Pa)")
ax1.set_ylim(6, 600); ax1.legend(loc="upper left", fontsize=9.5)

TY = np.linspace(0, 20, 161); NN = np.linspace(0.2, 1.0, 161)
RMS = np.zeros((NN.size, TY.size)); w = 1.0 / tau_meas
for i, nn in enumerate(NN):
    X = gd**nn
    for j, tt in enumerate(TY):
        yv = tau_meas - tt
        Kopt = max(np.sum(w**2 * X * yv) / np.sum(w**2 * X * X), 1e-9)
        RMS[i, j] = np.sqrt(np.mean(((yv - Kopt * X) * w)**2)) * 100
cs = ax2.contourf(TY, NN, RMS, levels=np.linspace(4, 32, 15), cmap="viridis_r", extend="max")
ax2.contour(TY, NN, RMS, levels=[6, 8, 12], colors="white", linewidths=0.5)
cb = fig.colorbar(cs, ax=ax2, pad=0.02); cb.set_label("RMS residual (%)"); cb.ax.minorticks_off()
ax2.plot(ty, n, "*", color=C["orange"], ms=12, mec="k", mew=0.4, label="true", zorder=3)
for (p, col, ls, lab), mk, mcol in zip(fits, ["o", "s", "D"], ["white", C["red"], C["cyan"]]):
    ax2.plot(p[0], p[2], mk, color=mcol, ms=4.5, mec="k", mew=0.5, zorder=4, label=lab.split(")")[0] + ")")
ax2.set_xlabel("$\\sigma_y$ (Pa)"); ax2.set_ylabel("$n$"); ax2.minorticks_off()
ax2.legend(loc="lower right", fontsize=9.5, labelcolor="white")
label_panels([ax1, ax2])
fig.tight_layout(w_pad=1.6)
# panel (b) is a filled contour map: export the underlying grid
fig.hb_extra = {"rms_residual_grid": [{"sigma_y_Pa": tt, "n": nn, "rms_residual_percent": RMS[i, j]}
                                      for i, nn in enumerate(NN) for j, tt in enumerate(TY)]}
save(fig, "fig07_fitting")

from pathlib import Path
import json
(Path(__file__).resolve().parents[1]/'data/fitting_illustration.json').write_text(json.dumps({'seed':1,'true_parameters':[ty,K,n],'log_full':p_log.tolist(),'log_high':p_high.tolist(),'linear_full':p_lin.tolist()},indent=2))
