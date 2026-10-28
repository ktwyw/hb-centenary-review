"""Fig. 12 of the article (fig12_hbp_pipe.pdf): Convergence of the HB-Papanastasiou model in pipe flow.

Scaled units: tau_y = K = R = 1, so the rate scale gdot_y = (tau_y/K)^(1/n) = 1 and the
dimensionless regularization parameter M = m*gdot_y is numerically equal to m here."""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid
from hbstyle import *

tau_y, K, n, R, phi = 1.0, 1.0, 0.5, 1.0, 0.5
tau_w = tau_y / phi
r = np.linspace(0, R, 2001); tau_r = tau_w * r / R

def gdot_hbp(t, m):
    return 0.0 if t <= 0 else brentq(lambda g: tau_y * (1 - np.exp(-m * g)) + K * g**n - t, 0.0, 1e6)
def profile(g):
    return cumulative_trapezoid(g[::-1], r[::-1], initial=0)[::-1] * -1
g_ex = np.where(tau_r > tau_y, (np.clip(tau_r - tau_y, 0, None) / K)**(1 / n), 0.0)
u_ex = profile(g_ex)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=DOUBLE)
ax1.plot(u_ex / u_ex.max(), r, color="k", lw=1.6, label="exact")
for m, col, ls in zip([10, 100, 1000, 10000], [C["green"], C["blue"], C["purple"], C["red"]], [":", "--", "-.", (0, (5, 1, 1, 1))]):
    g = np.array([gdot_hbp(t, m) for t in tau_r])
    ax1.plot(profile(g) / u_ex.max(), r, color=col, ls=ls, label=f"$M={m:g}$")
ax1.axhline(phi, color=C["grey"], lw=0.5, ls=(0, (1, 2)), gid="ideal plug radius r_p/R = phi")
ax1.text(0.03, phi + 0.02, "$r_p/R=\\phi$", fontsize=9.5, color=C["grey"])
ax1.set_xlabel("$u/u_{\\max}$"); ax1.set_ylabel("$r/R$")
ax1.set_xlim(0, 1.35); ax1.set_ylim(0, 1); ax1.set_xticks([0, 0.25, 0.5, 0.75, 1.0, 1.25]); ax1.legend(loc="lower left")
zero_as_integer(ax1, xdec=2, ydec=1)

m_grid = np.logspace(0.5, 5, 40)
for eps, col, ls in [(1e-2, C["blue"], ":"), (1e-3, C["red"], "--"), (1e-4, C["purple"], "-.")]:
    rp = np.minimum(1.0, (tau_y * (-np.expm1(-m_grid * eps)) + K * eps**n) / tau_w)
    ax2.semilogx(m_grid, rp, color=col, ls=ls, label=f"$\\dot\\gamma/\\dot\\gamma_y<{eps:g}$")
ax2.axhline(phi, color="k", lw=1.0, label="exact, $\\sigma<\\sigma_y$")
ax2.set_xlabel("$M=m\\dot\\gamma_y$"); ax2.set_ylabel("$r_p^{\\rm app}/R$")
ax2.set_ylim(0, 0.98); ax2.legend(loc="upper left")
zero_as_integer(ax2, ydec=1)
label_panels([ax1, ax2])
fig.tight_layout(w_pad=2.2)
save(fig, "fig12_hbp_pipe")
