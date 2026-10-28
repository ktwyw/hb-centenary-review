"""Fig. S4 of the supplementary material (figS4_mullineux.pdf): Mullineux's one-dimensional reduction of the HB least-squares fit.

A complete implementation of the method for fitting measured flow curves (Python package with tests,
confidence intervals and model comparison, plus a browser app) is available at
https://github.com/ktwyw/herschel-bulkley-determinant-fit
"""
import numpy as np
import matplotlib.pyplot as plt
from hbstyle import *

rng = np.random.default_rng(1)
ty, K, n = 10.0, 2.0, 0.5
gd = np.logspace(-2, 3, 26)
tau = hb_stress(gd, ty, K, n) * (1 + 0.05 * rng.standard_normal(gd.size))

def reduced(nn, w):
    """Profile objective F(n): least squares in (sigma_y, K) at fixed n, with sigma_y >= 0.

    If the unconstrained optimum has sigma_y < 0, the constrained optimum lies on the boundary
    sigma_y = 0 and K is refitted there (Supplementary Note S6 B)."""
    A = np.column_stack([np.ones_like(gd), gd**nn]) * w[:, None]; b = tau * w
    coef, *_ = np.linalg.lstsq(A, b, rcond=None)
    if coef[0] < 0:
        k = np.dot(A[:, 1], b) / np.dot(A[:, 1], A[:, 1])
        coef = np.array([0.0, k])
    return np.sum((A @ coef - b)**2), coef[0], coef[1]

NN = np.linspace(0.2, 1.0, 400)
out = {lab: np.array([reduced(nn, w) for nn in NN]) for lab, w in
       [("absolute weighting", np.ones_like(tau)), ("relative weighting", 1 / tau)]}
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=DOUBLE)
for (lab, res), col, ls in zip(out.items(), [C["blue"], C["black"]], ["--", "-"]):
    F = res[:, 0] / res[:, 0].min(); i = np.argmin(res[:, 0])
    ax1.semilogy(NN, F, color=col, ls=ls, label=lab)
    ax1.plot(NN[i], 1.0, "o", color=col, ms=3.5, gid=f"{lab}: minimum n*")
    if col == C["black"]:
        ax1.text(NN[i] + 0.05, 0.85, f"$n^*={NN[i]:.2f}$", fontsize=9.5, color=col, ha="left", va="top")
    else:
        ax1.text(NN[i] - 0.03, 0.9, f"$n^*={NN[i]:.2f}$", fontsize=9.5, color=col, ha="right", va="top")
ax1.axvline(n, color=C["grey"], lw=0.5, ls=(0, (1, 2)), gid="true n")
ax1.text(n - 0.015, 40, "true $n$", fontsize=9.5, color=C["grey"], va="center", ha="right")
ax1.set_xlabel("$n$"); ax1.set_ylabel("$F(n)/\\min F$")
ax1.set_xlim(0.2, 1.0); ax1.set_ylim(0.4, 300); ax1.legend(loc="upper right")

res = out["relative weighting"]
ax2.plot(NN, res[:, 1], color="k", label="$\\sigma_y(n)$")
ax2b = ax2.twinx()
ax2b.plot(NN, res[:, 2], color=C["red"], ls="--", label="$K(n)$")
ax2b.set_yscale("log"); ax2b.set_ylabel("$K(n)$ (Pa s$^{n}$)", color=C["red"])
ax2b.tick_params(axis="y", colors=C["red"], which="both")
ax2.set_ylabel("$\\sigma_y(n)$ (Pa)"); ax2.set_xlabel("$n$")
ax2.set_xlim(0.2, 1.0); ax2.set_ylim(0, 14); ax2.set_xticks([0.4, 0.6, 0.8, 1.0])
ax2.axvline(n, color=C["grey"], lw=0.5, ls=(0, (1, 2)), gid="true n")
ax2.plot(n, ty, "*", color=C["red"], ms=10, mec="k", mew=0.4, zorder=5, gid="true (n, sigma_y)")
ax2.text(0.30, 10.2, "$\\sigma_y$", fontsize=9.5); ax2.text(0.23, 3.0, "$K$", fontsize=9.5, color=C["red"])
label_panels([ax1, ax2])
fig.tight_layout(w_pad=2.6)
save(fig, "figS4_mullineux")
