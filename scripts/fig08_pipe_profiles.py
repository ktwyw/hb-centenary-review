"""Fig. 8 of the article (fig08_pipe_profiles.pdf): Laminar pipe flow of a Herschel-Bulkley fluid.
(a) Schematic of the stress distribution and unyielded plug. (b,c) Velocity
profiles normalized by the mean velocity for several phi = tau_y/tau_w."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from hbstyle import *


def hb_pipe_profile(r, phi, n, R=1.0, tau_w=1.0, K=1.0):
    tau_y = phi * tau_w
    excess = np.clip(tau_w * r / R - tau_y, 0.0, None)
    pref = (n / (n + 1.0)) * (R / tau_w) * K**(-1.0 / n)
    return pref * ((tau_w - tau_y)**((n + 1) / n) - excess**((n + 1) / n))


def mean_velocity(phi, n, R=1.0, tau_w=1.0, K=1.0):
    Q = (np.pi * R**3 * n * (tau_w / K)**(1 / n) * (1 - phi)**((n + 1) / n)
         * ((1 - phi)**2 / (3 * n + 1) + 2 * phi * (1 - phi) / (2 * n + 1) + phi**2 / (n + 1)))
    return Q / (np.pi * R**2)


r = np.linspace(-1, 1, 601)
phis = [0.0, 0.2, 0.5, 0.8]
fig = plt.figure(figsize=(6.6, 2.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.0], wspace=0.5)
ax0 = fig.add_subplot(gs[0]); ax1 = fig.add_subplot(gs[1]); ax2 = fig.add_subplot(gs[2])

ax0.set_xlim(0, 1); ax0.set_ylim(-1.3, 1.25); ax0.axis("off")
ax0.add_patch(Rectangle((0.05, -1), 0.9, 2, fc="#dcebf5", ec="none"))
ax0.add_patch(Rectangle((0.05, -0.45), 0.9, 0.9, fc="#f7e1d3", ec=C["red"], lw=0.7))
ax0.plot([0.05, 0.95], [1, 1], color="k", lw=1.2); ax0.plot([0.05, 0.95], [-1, -1], color="k", lw=1.2)
ax0.plot([0.05, 0.95], [0, 0], color="k", lw=0.5, ls="-.")
for y in np.linspace(-0.9, 0.9, 7):
    L = 0.35 * abs(y)
    if L > 0.02:
        ax0.add_patch(FancyArrowPatch((0.15, y), (0.15 + L, y), arrowstyle="->", mutation_scale=6, lw=0.6, color=C["blue"]))
ax0.text(0.5, 0.14, "plug, $|\\sigma|<\\sigma_y$", ha="center", va="center", fontsize=9.5)
ax0.text(0.5, 0.7, "$|\\sigma|>\\sigma_y$", ha="center", fontsize=9.5)
ax0.text(0.5, -0.78, "$|\\sigma|>\\sigma_y$", ha="center", fontsize=9.5)
ax0.text(0.95, 1.03, "$r=R$", fontsize=9.5, ha="right", va="bottom")
ax0.text(0.94, 0.47, "$r=r_p$", fontsize=9.5, ha="right", va="bottom", color=C["red"])
ax0.text(0.5, -1.05, "arrows: $|\\sigma|=\\sigma_w|r|/R$", fontsize=9.5, color=C["blue"], ha="center", va="top")

for ax, n in [(ax1, 0.5), (ax2, 1.0)]:
    for phi, col, ls in zip(phis, [C["black"], C["red"], C["blue"], C["green"]], ["-", "--", "-.", ":"]):
        u = hb_pipe_profile(np.abs(r), phi, n); um = mean_velocity(phi, n)
        ax.plot(u / um, r, color=col, ls=ls, label=f"$\\phi={phi}$")
    ax.set_xlabel("$u/\\bar u$"); ax.set_ylabel("$r/R$")
    ax.set_xlim(0, 2.2); ax.set_ylim(-1, 1)
    ax.text(1.35, -0.92, f"$n={n}$", fontsize=9.5)
h, l = ax1.get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=4, bbox_to_anchor=(0.64, 0.0), handlelength=2.4)
label_panels([ax0], x=-0.15); label_panels([ax1, ax2], labels=["(b)", "(c)"], x=-0.32)
save(fig, "fig08_pipe_profiles")
