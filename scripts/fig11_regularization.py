"""Fig. 11 of the article (fig11_regularization.pdf): Regularization of the Herschel-Bulkley model. (a) Flow curves;
(b) apparent viscosity, exact vs Papanastasiou (three m) and bi-viscosity."""
import numpy as np
import matplotlib.pyplot as plt
from hbstyle import *

tau_y, K, n = 10.0, 2.0, 0.5
gd = np.logspace(-5, 2, 600)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=DOUBLE)
ax1.loglog(gd, hb_stress(gd, tau_y, K, n), color="k", lw=1.6, label="exact HB")
ax2.loglog(gd, hb_viscosity(gd, tau_y, K, n), color="k", lw=1.6, label="exact HB")
for m, col, ls in zip([10, 100, 1000], [C["green"], C["blue"], C["purple"]], [":", "--", "-."]):
    eta = papanastasiou_viscosity(gd, tau_y, K, n, m)
    ax1.loglog(gd, eta * gd, color=col, ls=ls, label=f"Papanastasiou, $m={m}$ s")
    ax2.loglog(gd, eta, color=col, ls=ls, label=f"Papanastasiou, $m={m}$ s")
eta_bv, gc = biviscosity(gd, tau_y, K, n, eta_r=3e3)
ax1.loglog(gd, eta_bv * gd, color=C["red"], ls=(0, (5, 1, 1, 1)), label="bi-viscosity, $\\eta_r=3\\times10^3$ Pa s")
ax2.loglog(gd, eta_bv, color=C["red"], ls=(0, (5, 1, 1, 1)), label="bi-viscosity, $\\eta_r=3\\times10^3$ Pa s")
ax2.axvline(gc, color=C["red"], lw=0.5, ls=(0, (1, 2)), gid="bi-viscosity switching rate gamma_dot_r")
ax2.text(gc * 1.4, 2.5e-1, "$\\dot\\gamma_r$", color=C["red"], fontsize=9.5)
ax1.axhline(tau_y, color="k", lw=0.5, ls=(0, (1, 2)), gid="yield stress sigma_y")
ax1.text(2e-5, tau_y * 1.3, "$\\sigma_y$", fontsize=9.5)
ax1.set_xlabel("$\\dot\\gamma$ (s$^{-1}$)"); ax1.set_ylabel("$\\sigma$ (Pa)")
ax1.set_ylim(1e-2, 1e2)
ax2.set_xlabel("$\\dot\\gamma$ (s$^{-1}$)"); ax2.set_ylabel("$\\eta$ (Pa s)")
ax2.set_ylim(1e-1, 1e7)
ax2.annotate("residual\ndivergence of\n$K\\dot\\gamma^{\\,n-1}$ ($n<1$)", xy=(3e-5, 4.0e2), xytext=(1.0e-5, 2.0e-1),
             fontsize=9.5, arrowprops=dict(arrowstyle="->", lw=0.5, shrinkA=2), ha="left", va="bottom")
label_panels([ax1, ax2])
fig.tight_layout(w_pad=2.2)
h, l = ax2.get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.0), handlelength=2.4)
save(fig, "fig11_regularization")
