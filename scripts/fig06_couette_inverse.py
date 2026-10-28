"""Fig. 6 of the article (fig06_couette_inverse.pdf): Partially yielded Couette gap and Newtonian-conversion bias."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import curve_fit
from hbstyle import *

tau_y, K, n, Ri = 10.0, 2.0, 0.5, 1.0
def omega(tau_b, kappa):
    Ro = kappa * Ri
    if tau_b <= tau_y: return 0.0
    r_max = min(Ri * np.sqrt(tau_b / tau_y), Ro)
    return quad(lambda rr: ((tau_b * (Ri / rr)**2 - tau_y) / K)**(1 / n) / rr, Ri, r_max, limit=200)[0]
model = lambda g, a, b, c: a + b * g**c

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=DOUBLE)
tb = np.linspace(tau_y * 1.0001, 80, 300)
kappas = [1.068, 1.5, 2.0]; cols = [C["black"], C["red"], C["blue"]]; lss = ["-", "--", ":"]; mks = ["o", "x", "s"]
for kappa, col, ls in zip(kappas, cols, lss):
    frac = np.clip((Ri * np.sqrt(tb / tau_y) - Ri) / (kappa * Ri - Ri), 0, 1)
    ax1.plot(tb / tau_y, frac, color=col, ls=ls, label=f"$R_o/R_i={kappa}$")
ax1.set_xlabel("$\\sigma_b/\\sigma_y$"); ax1.set_ylabel("$(r_y-R_i)/(R_o-R_i)$")
ax1.set_xlim(1, 6.2); ax1.set_ylim(0, 1.05); ax1.legend(loc="lower right", bbox_to_anchor=(1.0, 0.09))
ax1.text(4.15, 0.45, "fully yielded for\n$\\sigma_b/\\sigma_y\\geq(R_o/R_i)^2$", fontsize=9.5)

gfine = np.logspace(-1, 3, 200)
ax2.loglog(gfine, hb_stress(gfine, tau_y, K, n), color="k", lw=1.6, label="true flow curve")
fit_records = []
for kappa, col, mk in zip(kappas, cols, mks):
    Ro = kappa * Ri
    tbs = np.logspace(np.log10(tau_y * 1.02), np.log10(tau_y * 40), 40)
    Om = np.array([omega(t, kappa) for t in tbs]); g_app = 2 * Om * Ro**2 / (Ro**2 - Ri**2); ok = g_app > 0
    p, _ = curve_fit(lambda g, a, b, c: np.log(model(g, a, b, c)), g_app[ok], np.log(tbs[ok]), p0=[8, 2, 0.5], bounds=([0, 1e-6, 0.05], [50, 100, 2]))
    fit_records.append({'radius_ratio': kappa, 'tau_y_Pa': float(p[0]),
                        'K_Pa_s_to_n': float(p[1]), 'n': float(p[2]),
                        'apparent_rate_range_per_s': [float(g_app[ok].min()), float(g_app[ok].max())]})
    ax2.loglog(g_app[ok], tbs[ok], mk, ms=3.5, color=col, mfc="none" if mk != "x" else col, mew=0.7,
               label=f"$R_o/R_i={kappa}$")
ax2.set_xlabel("$\\dot\\gamma_{\\rm app}$ (s$^{-1}$)"); ax2.set_ylabel("$\\sigma_b$ (Pa)")
ax2.set_xlim(1e-1, 1e3); ax2.set_ylim(8, 150); ax2.legend(loc="upper left", fontsize=9.5)
label_panels([ax1, ax2])
fig.tight_layout(w_pad=2.2)
(Path(__file__).resolve().parents[1] / 'data/couette_illustration.json').write_text(
    json.dumps({'evidence_type': 'Calculated Couette conversion and fit illustration',
                'objective': 'Unweighted squared log-stress residuals',
                'true_triplet': [tau_y, K, n], 'fits': fit_records}, indent=2) + '\n')
save(fig, "fig06_couette_inverse")
