"""Fig. S7 of the supplementary material (figS7_friction_factor.pdf): Fanning friction factor vs Metzner-Reed Reynolds number.

Laminar line f_F = 16/Re_MR (Metzner & Reed, AIChE J. 1, 434, 1955); turbulent curves from the
Dodge-Metzner correlation with its published coefficients (Dodge & Metzner, AIChE J. 5, 189, 1959).
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from hbstyle import *

def dodge_metzner(Re, n):
    return brentq(lambda ff: 1 / np.sqrt(ff) - (4 / n**0.75) * np.log10(Re * ff**(1 - n / 2)) + 0.4 / n**1.2, 1e-5, 1.0)

fig, ax = plt.subplots(figsize=(3.4, 2.6))
Re_lam = np.logspace(2, np.log10(2100), 50)
ax.loglog(Re_lam, 16 / Re_lam, color="k", lw=1.4, label="laminar, $f_F=16/\\mathrm{Re_{MR}}$")
Re_t = np.logspace(np.log10(3000), 5.5, 80)
for n, col, ls in zip([1.0, 0.8, 0.6, 0.4], [C["grey"], C["blue"], C["red"], C["purple"]], ["-", "--", "-.", ":"]):
    ax.loglog(Re_t, [dodge_metzner(R, n) for R in Re_t], color=col, ls=ls, label=f"turbulent, $n$′ = {n}")
ax.axvspan(2100, 3000, color=C["lightgrey"], alpha=0.5, lw=0, gid="illustrative transition band")
ax.text(np.sqrt(2100 * 3000), 1.25e-3, "transition", rotation=90,
        ha="center", va="bottom", fontsize=9.5, color=C["grey"])
slope_triangle(ax, 2e2, 4e-2, 0.7, -1.0, "$1$", side="below", run_label="$-1$")
ax.set_xlabel("$\\mathrm{Re_{MR}}$"); ax.set_ylabel("$f_F$")
ax.set_xlim(1e2, 3e5); ax.set_ylim(1e-3, 3e-1); ax.set_xticks([1e3, 1e4, 1e5]); ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.0), fontsize=8.25, handlelength=1.6, borderaxespad=0.3)
save(fig, "figS7_friction_factor")
