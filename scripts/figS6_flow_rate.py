"""Fig. S6 of the supplementary material (figS6_flow_rate.pdf): Closed-form pipe-flow relations for an HB fluid."""
import numpy as np
import matplotlib.pyplot as plt
from hbstyle import *

def Qstar(phi, n):
    return n * (1 - phi)**((n + 1) / n) * ((1 - phi)**2 / (3 * n + 1) + 2 * phi * (1 - phi) / (2 * n + 1) + phi**2 / (n + 1))
def plug_over_mean(phi, n):
    return (n / (n + 1)) * (1 - phi)**((n + 1) / n) / Qstar(phi, n)

phi = np.linspace(0, 0.999, 400)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=DOUBLE)
for n, col, ls in zip([0.3, 0.5, 0.7, 1.0, 1.5], PALETTE, LS):
    mk = dict(marker="o", ms=2.6, mfc="white", markevery=0.12) if n == 1.5 else {}
    ax1.plot(phi, Qstar(phi, n), color=col, ls=ls, label=f"$n={n}$", **mk)
    ax2.plot(phi, plug_over_mean(phi, n), color=col, ls=ls, label=f"$n={n}$", **mk)
ax1.set_xlabel("$\\phi=\\sigma_y/\\sigma_w$"); ax1.set_ylabel("$Q/[\\pi R^3(\\sigma_w/K)^{1/n}]$")
ax1.set_xlim(0, 1); ax1.set_ylim(0, 0.3); ax1.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); ax1.legend(loc="upper right")
zero_as_integer(ax1, xdec=1, ydec=2)
ax2.plot(phi, phi, color=C["grey"], lw=0.8, ls="-", label="$r_p/R=\\phi$")
ax2.set_xlabel("$\\phi=\\sigma_y/\\sigma_w$"); ax2.set_ylabel("$u_p/\\bar u,\\; r_p/R$")
ax2.set_xlim(0, 1); ax2.set_ylim(0, 2.8); ax2.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); ax2.legend(loc="upper right", ncol=2, fontsize=9.5)
zero_as_integer(ax2, xdec=1, ydec=1)
ax2.text(0.76, 0.44, "$r_p/R$", fontsize=9.5, rotation=22, color=C["grey"])
ax2.text(0.04, 0.55, "$u_p/\\bar u$", fontsize=9.5)
label_panels([ax1, ax2])
fig.tight_layout(w_pad=2.2)
save(fig, "figS6_flow_rate")
