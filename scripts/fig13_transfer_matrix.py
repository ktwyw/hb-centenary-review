"""Fig. 13 of the article (fig13_transfer_matrix.pdf): Overview of the experimental tests of Table XII: which change relative to
the calibration each study examined, and its main finding. Encodes the findings summarized in Sec. IX (Table XII);
no data are computed here.

Finding codes: 'T' predicted with the shear-calibrated parameters;
'A' predicted after an extension, a fitted coefficient or a parameter change;
'F' physics not captured by the baseline model.
"""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from hbstyle import C, save

COLS = ["Local\nflow curve", "Flow\ntype", "Geometry,\nscale", "Walls,\ngap", "History,\ntime"]
DX = 1.45
# (group, label, {col_index: outcome})
CASES = [
 ("A", "Emulsions, wide-gap MRI (Ovarlez 2008)",        {0: "T"}),
 ("A", "Carbopol gel, MRI and cone–plate (Coussot 2009)", {0: "T"}),
 ("A", "Foam, local MRI (Ovarlez, Krishan 2010)",                 {0: "T", 4: "T"}),
 ("B", "Emulsion in microchannels (Goyon 2008)",        {3: "A"}),
 ("B", "Adhesive emulsion, smooth narrow gap (Bécu 2006)", {3: "F", 4: "F"}),
 ("C", "Clay-loaded emulsion (Fall 2010)",              {4: "F"}),
 ("C", "Bentonite, rest-time series (Møller 2006)",     {4: "F"}),
 ("D", "Rotation plus squeeze (Ovarlez, Barral 2010)",           {1: "T"}),
 ("D", "Normal stresses (de Cagny 2019)",               {1: "F"}),
 ("D", "Capillary thinning, onset (Martinie 2013)",     {1: "T"}),
 ("D", "Slipping plates (Zhang 2018)",                  {1: "F"}),
 ("D", "Gravity filaments, onset (Geffrault 2021)",     {1: "T"}),
 ("D", "Gravity filaments, flow curve (Geffrault 2021)", {1: "A"}),
 ("D", "Cross-slot, EVP model (Varchanis 2020)",        {1: "A"}),
 ("E", "Sphere drag and arrest (Tabuteau 2007)",        {2: "T"}),
 ("E", "Flow field around sphere (Putz 2008)",          {2: "F"}),
 ("F", "Pipe rheometer (Magnon 2021)",                  {2: "A", 3: "A"}),
 ("F", "Free-surface channel, kaolin (Chambon 2014)",   {2: "T"}),
 ("F", "Free-surface channel, Carbopol (Chambon 2014)", {2: "A"}),
 ("G", "Bead packings, pressure drop (Chevalier 2013)", {2: "A"}),
 ("G", "Bead packing, NMR field (Chevalier 2014)",      {2: "F"}),
 ("G", "Pillar array (Waisbord 2019)",                  {2: "A", 3: "A", 4: "A"}),
]
# Colour repeats the fill code (filled / half-filled / open), so the figure also reads in grayscale.
# Okabe-Ito hues: blue = predicted, orange = predicted after extension, vermilion = not captured.
DARK_ORANGE = "#B87A00"   # darker edge for the Okabe-Ito orange, for contrast on white
STYLE = {"T": dict(marker="o", ms=7.5, mfc=C["blue"], mec=C["blue"], mew=1.0),
         "A": dict(marker="o", ms=7.5, mfc=C["orange"], mec=DARK_ORANGE, mew=1.0, fillstyle="left",
                   markerfacecoloralt="white"),
         "F": dict(marker="o", ms=7.5, mfc="white", mec=C["red"], mew=1.4)}

fig, ax = plt.subplots(figsize=(6.6, 6.1))
n = len(CASES)
prev = None
for i, (g, lab, outs) in enumerate(CASES):
    y = n - 1 - i
    if g != prev:
        ax.axhline(y + 0.5, color=C["lightgrey"], lw=0.5)
        ax.text(-0.8, y, g, ha="center", va="center", fontsize=9, fontweight="bold")
        prev = g
    for j in range(len(COLS)):
        ax.plot(j*DX, y, marker=".", color=C["lightgrey"], ms=2)
    for j, o in outs.items():
        ax.plot(j*DX, y, linestyle="none", **STYLE[o])
    ax.text(-1.1, y, lab, ha="right", va="center", fontsize=8.5)
ax.set_xlim(-1.0, (len(COLS) - 1) * DX + 0.6)
ax.set_ylim(-0.6, n - 0.4)
ax.set_xticks([j*DX for j in range(len(COLS))])
ax.set_xticklabels(COLS, fontsize=8.5)
ax.xaxis.tick_top()
ax.tick_params(axis="x", length=0, pad=2)
ax.set_yticks([])
ax.minorticks_off()
for s in ax.spines.values():
    s.set_visible(False)
h = [Line2D([], [], linestyle="none", **STYLE["T"]), Line2D([], [], linestyle="none", **STYLE["A"]),
     Line2D([], [], linestyle="none", **STYLE["F"])]
ax.legend(h, ["predicted with shear parameters", "predicted after extension or refit", "not captured by the baseline model"],
          loc="upper center", bbox_to_anchor=(0.4, -0.005), ncol=3, fontsize=8.5, handletextpad=0.3, columnspacing=1.2)
fig.subplots_adjust(left=0.44, right=0.99, top=0.91, bottom=0.05)
# the figure is a matrix of findings: export it as a table
_codes = {"T": "predicted with shear parameters", "A": "predicted after extension or refit",
          "F": "not captured by the baseline model"}
fig.hb_table = [{"group": g, "study": lab, "what changes": COLS[j].replace("\n", " "), "finding": _codes[o]}
                for g, lab, outs in CASES for j, o in sorted(outs.items())]
save(fig, "fig13_transfer_matrix")
