"""Fig. 4 of the article (fig04_flow_curves.pdf): Flow curves of the Herschel-Bulkley model and its limiting cases.
(a) Linear axes: shear stress vs shear rate. (b) Log-log apparent viscosity with
slope triangles for the low-rate (-1) and high-rate (n-1) asymptotes."""
import numpy as np
import matplotlib.pyplot as plt
from hbstyle import *

gdot = np.linspace(0, 20, 400)
gd_log = np.logspace(-3, 3, 400)
cases = [
    ("Newtonian",                       0.0, 1.0, 1.0, C["grey"],   "-"),
    ("power law, $n=0.5$",              0.0, 3.0, 0.5, C["green"],  ":"),
    ("Bingham",                         8.0, 1.0, 1.0, C["blue"],   "--"),
    ("Herschel–Bulkley, $n=0.5$",       8.0, 3.0, 0.5, C["black"],  "-"),
    ("Herschel–Bulkley, $n=1.5$",       8.0, 0.3, 1.5, C["red"],    "-."),
]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=DOUBLE)
for lab, ty, K, n, col, ls in cases:
    ax1.plot(gdot, hb_stress(gdot, ty, K, n), ls, color=col, label=lab)
    ax2.loglog(gd_log, hb_viscosity(gd_log, ty, K, n), ls, color=col, label=lab)
ax1.axhline(8.0, color="k", lw=0.5, ls=(0, (1, 2)), gid="yield stress sigma_y")
ax1.text(19.5, 8.6, "$\\sigma_y$", ha="right", va="bottom", fontsize=9.5)
ax1.set_xlabel("$\\dot\\gamma$ (s$^{-1}$)")
ax1.set_ylabel("$\\sigma$ (Pa)")
ax1.set_xlim(0, 20); ax1.set_ylim(0, 52); ax1.set_xticks(range(0, 21, 5))
ax1.legend(loc="upper left")
slope_triangle(ax2, 1e-2, 2.0e3, 1.0, -1.0, "$1$", side="above", run_label="$-1$")
slope_triangle(ax2, 1e2, 0.15, 0.8, -0.5, "$1-n$", side="below", run_label="$-1$")
ax2.set_xlabel("$\\dot\\gamma$ (s$^{-1}$)")
ax2.set_ylabel("$\\eta=\\sigma/\\dot\\gamma$ (Pa s)")
ax2.set_ylim(1e-2, 1e5)
label_panels([ax1, ax2])
fig.tight_layout(w_pad=2.2)
save(fig, "fig04_flow_curves")
