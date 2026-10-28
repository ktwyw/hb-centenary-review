"""Fig. 5 of the article (fig05_microstructure.pdf): Microstructural origins of Herschel-Bulkley behavior (schematic).

(a)-(c) Three classes of yield-stress material: a soft-particle glass (jammed
    deformable particles above random close packing), an attractive colloidal
    gel (percolated network of attractive particles), and a suspension of
    non-colloidal particles in a yield-stress matrix.
(d) Theory schematic for dimensionless excess stress: the n=1/2 elastohydrodynamic
    scaling of Seth et al. (2011) for soft-particle glasses. No experimental or synthetic observations are plotted.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from hbstyle import *

rng = np.random.default_rng(3)
fig = plt.figure(figsize=(7.0, 4.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.25], hspace=0.55, wspace=0.25)
axs = [fig.add_subplot(gs[0, i]) for i in range(3)]
axd = fig.add_subplot(gs[1, :])


def packing(ax, n_side, jitter, deform, color, r=0.5):
    """Draw a jittered hexagonal packing; deform>0 draws slightly polygonal (squashed) particles."""
    pts = []
    for i in range(n_side):
        for j in range(n_side):
            x = i + 0.5 * (j % 2) + jitter * rng.standard_normal()
            y = j * 0.87 + jitter * rng.standard_normal()
            pts.append((x, y))
    for (x, y) in pts:
        if deform:
            ang = np.linspace(0, 2 * np.pi, 7)[:-1] + rng.uniform(0, 1)
            rad = r * (1 + 0.08 * rng.standard_normal(6))
            ax.add_patch(Polygon(np.c_[x + rad * np.cos(ang), y + rad * np.sin(ang)],
                                 closed=True, fc=color, ec="k", lw=0.5, alpha=0.85))
        else:
            ax.add_patch(Circle((x, y), r, fc=color, ec="k", lw=0.5, alpha=0.85))
    return pts


# (a) soft-particle glass: dense, deformed, touching
ax = axs[0]
packing(ax, 6, 0.02, True, C["cyan"], r=0.55)
ax.text(0.5, -0.04, "Athermal jammed particles\n(microgels, emulsions, foams)", transform=ax.transAxes, ha="center", va="top", fontsize=9.5)
# (b) attractive gel: sparse network
ax = axs[1]
pts = []
for _ in range(40):
    pts.append((rng.uniform(0, 5.5), rng.uniform(0, 4.5)))
pts = np.array(pts)
for p in pts:
    d = np.linalg.norm(pts - p, axis=1)
    for q in pts[(d < 1.3) & (d > 0)]:
        ax.plot([p[0], q[0]], [p[1], q[1]], color=C["orange"], lw=1.0, alpha=0.6, zorder=1)
for p in pts:
    ax.add_patch(Circle(p, 0.22, fc=C["orange"], ec="k", lw=0.5, zorder=2))
ax.text(0.5, -0.04, "Attractive colloidal gel\n(clays, cement, flocculated sludge)", transform=ax.transAxes, ha="center", va="top", fontsize=9.5)
# (c) suspension in yield-stress matrix
ax = axs[2]
ax.add_patch(plt.Rectangle((-0.5, -0.5), 6.5, 5.5, fc=C["cyan"], alpha=0.25, ec="none"))
placed = []
while len(placed) < 14:
    x, y, r = rng.uniform(0.3, 5.2), rng.uniform(0.3, 3.85), rng.uniform(0.3, 0.55)
    if all((x - a) ** 2 + (y - b) ** 2 > (r + s + 0.08) ** 2 for a, b, s in placed):
        placed.append((x, y, r))
for x, y, r in placed:
    ax.add_patch(Circle((x, y), r, fc=C["grey"], ec="k", lw=0.5))
ax.text(2.75, 4.65, "yield-stress matrix", ha="center", fontsize=9.5, color=C["blue"])
ax.text(0.5, -0.04, "Non-colloidal particles in a\nyield-stress matrix (concrete)", transform=ax.transAxes, ha="center", va="top", fontsize=9.5)
for ax in axs:
    ax.set_xlim(-0.5, 6); ax.set_ylim(-0.5, 5); ax.set_aspect("equal"); ax.axis("off")

# (d) theory schematic
x = np.logspace(-6, -1, 200)
axd.loglog(x, np.sqrt(x), color="k", lw=1.4, label="soft-particle glass, elastohydrodynamic scaling, $n=1/2$ (Seth et al. 2011)")
slope_triangle(axd, 1e-4, 5e-3, 1.0, 0.5, "$1/2$", side="below")
axd.set_xlabel("$\\eta_s\\dot\\gamma/G_0$")
axd.set_ylabel("$(\\sigma-\\sigma_y)/\\sigma_y$")
axd.set_xlim(1e-6, 1e-1); axd.set_ylim(5e-4, 1)
axd.legend(loc="upper left")

for ax, lab in zip(axs, ["(a)", "(b)", "(c)"]):
    ax.text(-0.02, 1.0, lab, transform=ax.transAxes, fontsize=9.5, style="italic", va="bottom")
axd.text(-0.105, 1.03, "(d)", transform=axd.transAxes, fontsize=9.5, style="italic", va="bottom")
save(fig, "fig05_microstructure")
