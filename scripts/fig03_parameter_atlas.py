"""Fig. 3 of the article (fig03_parameter_atlas.pdf): one law, many materials.

Herschel-Bulkley yield stress and flow index reported for materials from several
fields, taken from studies cited in the review (data/hb_parameter_atlas.csv; the
source table, section or data file is given for every row). Values are those
reported by the sources (or by this review's reanalyses); they are not refitted.
The 1926 rubber solutions have sigma_y = 0 and are shown in a separate strip;
fresh concrete (de Larrard et al. 1998) is shown by its reported exponent range
only, because its yield stresses are not compiled here.
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from hbstyle import C, save, label_panels

ROOT = Path(__file__).resolve().parents[1]
rows = list(csv.DictReader(open(ROOT / "data" / "hb_parameter_atlas.csv")))

style = {  # class: (color, marker, legend label)
    "Drilling and well fluids": (C["orange"], "D", "Drilling and well fluids"),
    "Foams and aphrons": (C["cyan"], "o", "Foams and aphrons"),
    "Emulsions": (C["blue"], "s", "Concentrated emulsions"),
    "Carbopol and polymer gels": (C["purple"], "^", "Carbopol and polymer gels"),
    "Clay suspensions": (C["red"], "v", "Clay suspension (kaolin)"),
}

fig = plt.figure(figsize=(6.6, 3.0), layout="constrained")
gs = fig.add_gridspec(1, 3, width_ratios=[0.16, 1.0, 0.22], wspace=0.02)
ax0 = fig.add_subplot(gs[0]); ax = fig.add_subplot(gs[1], sharey=ax0); ax2 = fig.add_subplot(gs[2], sharey=ax0)
YL = (0.15, 1.95)

for cls, (col, mk, lab) in style.items():
    first = True
    for r in [r for r in rows if r["class"] == cls]:
        s, n = float(r["sigma_y_Pa"]), float(r["n"])
        if r["sigma_y_max_Pa"]:
            s2 = float(r["sigma_y_max_Pa"])
            lw, z = (3.2, 4) if cls == "Emulsions" else (1.4, 3)
            ax.plot([s, s2], [n, n], color=col, lw=lw, solid_capstyle="butt", zorder=z)
            ax.plot([s, s2], [n, n], ls="none", marker="|", ms=7, color=col, mew=1.2, zorder=z)
            continue
        ax.plot(s, n, ls="none", marker=mk, ms=5.5, mfc=col, mec="k", mew=0.4, alpha=0.95, zorder=5,
                label=lab if first else None)
        first = False

# the same drilling fluid fitted over two shear-rate windows
d = [r for r in rows if r["source_key"] == "SaasenYtrehus2020"]
ax.annotate("", xy=(0.2, float(d[1]["n"])), xytext=(0.2, float(d[0]["n"])),
            arrowprops=dict(arrowstyle="-|>", color=C["orange"], lw=0.8, shrinkA=4, shrinkB=4))
ax.text(0.26, 0.885, "same fluid,\nnarrower fit window", fontsize=8.5, color="0.25", va="center")
ax.axhline(1.0, color="0.55", lw=0.6, ls="--")
ax.text(700, 1.03, "Bingham ($n=1$)", fontsize=8.5, color="0.35", ha="right", va="bottom")
ax.set_xscale("log"); ax.set_xlim(0.1, 1000); ax.set_ylim(*YL)
ax.set_xlabel(r"Yield stress $\sigma_y$ (Pa)")
ax.legend(loc="upper left", fontsize=8.5, handletextpad=0.3, borderaxespad=0.3)
plt.setp(ax.get_yticklabels(), visible=False)

# left strip: sigma_y = 0 (1926 rubber solutions, reanalyzed)
for r in [r for r in rows if r["class"] == "1926 rubber solutions"]:
    ax0.plot(0, float(r["n"]), ls="none", marker="o", ms=5.5, mfc="white", mec="k", mew=0.8)
ax0.set_xlim(-1, 1); ax0.set_xticks([0]); ax0.set_xticklabels(["0"])
ax0.set_ylabel("Flow index $n$")
ax0.text(0, 0.95, "1926\nrubber\nsolutions", ha="center", va="bottom", fontsize=8.5)

# right strip: fresh concrete, exponent range only
ax2.plot([0, 0], [1.10, 1.88], color=C["grey"], lw=5, solid_capstyle="butt")
ax2.set_xlim(-1, 1); ax2.set_xticks([])
ax2.text(0, 1.06, "fresh\nconcrete", ha="center", va="top", fontsize=8.5)
ax2.set_xlabel(r"$\sigma_y$ not shown", fontsize=8.5)
plt.setp(ax2.get_yticklabels(), visible=False)
for a in (ax0, ax2):
    a.tick_params(axis="x", length=0)
# the points are rows of data/hb_parameter_atlas.csv: export them with their sources
fig.hb_table = [{k: r[k] for k in r} for r in rows]
# the right strip is a range of reported exponents, not a row of the atlas
fig.hb_table.append({**{k: "" for k in rows[0]}, "material": "Fresh concrete (range of reported exponents)",
                     "class": "fresh concrete (right strip)", "n": "1.10-1.88",
                     "fit": "range of n reported for HB fits; yield stresses not compiled",
                     "source_key": "deLarrard1998", "source_location": "Sec. III B of the review"})
save(fig, "fig03_parameter_atlas")
