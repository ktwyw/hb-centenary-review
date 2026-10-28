"""Fig. 2 of the article (fig02_usage.pdf): a century of use of the Herschel-Bulkley model (OpenAlex counts).

(a) Works per year whose title or abstract (filled) or full text (open) contains
    "Herschel-Bulkley", 1950-2025 (2026 is a partial year and omitted).
(b) Works with "Herschel-Bulkley" in title or abstract by application area,
    grouping OpenAlex's automatic primary topics (data/hb_bibliometrics.json).
(c) Share of the works naming the HB, Bingham and Casson models (title or abstract)
    by OpenAlex research field, grouped into six classes (data/hb_bibliometrics.json,
    "model_comparison").
"""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from hbstyle import C, save, label_panels

ROOT = Path(__file__).resolve().parents[1]
d = json.loads((ROOT / "data" / "hb_bibliometrics.json").read_text())
ta = np.array([r for r in d["title_abstract_by_year"] if 1950 <= r[0] <= 2025])
ft = np.array([r for r in d["fulltext_by_year"] if 1950 <= r[0] <= 2025])
topics = dict(d["top_topics_title_abstract"])
groups = {
    "Fluid mechanics and computation": ["Rheology and Fluid Dynamics Studies", "Fluid Dynamics and Mixing",
        "Fluid Dynamics Simulations and Interactions", "Fluid Dynamics and Heat Transfer", "Heat Transfer and Optimization",
        "Hydraulic flow and structures"],
    "Drilling, oil and gas": ["Drilling and Well Engineering", "Petroleum Processing and Analysis",
        "Enhanced Oil Recovery Techniques", "Hydraulic Fracturing and Reservoir Analysis"],
    "Foods and biopolymers": ["Polysaccharides Composition and Applications", "Food composition and properties",
        "Proteins in Food Systems"],
    "Construction materials, grouts": ["Concrete and Cement Materials Research",
        "Innovations in Concrete and Construction Materials", "Grouting, Rheology, and Soil Mechanics"],
    "Mudflows, tailings, slurries": ["Landslides and related hazards", "Tailings Management and Properties",
        "Coal Combustion and Slurry Processing"],
    "MR and ER fluids": ["Vibration Control and Rheological Fluids"],
    "Blood and biomedicine": ["Blood properties and coagulation"],
    "Lubricants and greases": ["Tribology and Lubrication Engineering", "Lubricants and Their Additives"],
    "3D printing": ["Additive Manufacturing and 3D Printing Technologies"],
}
vals = {k: sum(topics.get(t, 0) for t in v) for k, v in groups.items()}
order = sorted(vals, key=vals.get)

fig = plt.figure(figsize=(6.6, 4.9), layout="constrained")
top, bottom = fig.subfigures(2, 1, height_ratios=[3.1, 1.6], hspace=0.02)
ax1, ax2 = top.subplots(1, 2, gridspec_kw={"width_ratios": [1.0, 1.0]})
ax3 = bottom.subplots(1, 1)
ax1.semilogy(ft[:, 0], ft[:, 1], "o", ms=3, mfc="white", mec=C["black"], mew=0.7, label="full text")
ax1.semilogy(ta[:, 0], ta[:, 1], "s", ms=3, color=C["blue"], label="title or abstract")
ax1.set_xlim(1948, 2028); ax1.set_ylim(0.7, 2000)
ax1.set_xlabel("Year"); ax1.set_ylabel("Publications per year")
ax1.legend(loc="upper left")
y = np.arange(len(order))
ax2.barh(y, [vals[k] for k in order], color=C["lightgrey"], edgecolor=C["black"], linewidth=0.5, height=0.65,
         gid="works with Herschel-Bulkley in title or abstract")
for yi, k in zip(y, order):
    ax2.text(vals[k] + 25, yi, f"{vals[k]}", va="center", fontsize=8.5)
ax2.set_yticks(y); ax2.set_yticklabels(order, fontsize=8.5)
ax2.set_xlabel("Publications, 1930–2026")
ax2.set_xlim(0, 2050)
ax2.minorticks_off(); ax2.tick_params(axis="y", length=0)
# (c) field composition of the works naming each model
mc = d["model_comparison"]["models"]
classes = [
    ("Engineering", ["Engineering"], "#4d4d4d"),
    ("Chemical engineering", ["Chemical Engineering"], "#a6a6a6"),
    ("Agricultural and biological sciences", ["Agricultural and Biological Sciences"], C["green"]),
    ("Medicine and health", ["Medicine", "Nursing", "Health Professions", "Pharmacology, Toxicology and Pharmaceutics",
                             "Dentistry", "Neuroscience", "Immunology and Microbiology"], C["red"]),
    ("Environmental and Earth sciences", ["Environmental Science", "Earth and Planetary Sciences"], C["cyan"]),
    ("Other fields", None, "#d9c7e8"),
]
rows = [("Herschel–Bulkley", "HB"), ("Bingham", "Bingham"), ("Casson", "Casson")]
for yi, (lab, key) in enumerate(rows):
    byf = dict(mc[key]["by_field"]); tot = sum(byf.values()); used = 0; left = 0.0
    for name, members, col in classes:
        cnt = sum(byf.get(m, 0) for m in members) if members else tot - used
        used += cnt if members else 0
        w = 100.0 * cnt / tot
        bars = ax3.barh(yi, w, left=left, color=col, edgecolor="white", linewidth=0.6, height=0.62,
                        label=name if yi == 0 else None)
        for p in bars:
            p.set_gid(f"{lab}: {name} (% of works)")
        if w >= 6:
            ax3.text(left + w / 2, yi, f"{w:.0f}", ha="center", va="center", fontsize=8.5,
                     color="white" if col in ("#4d4d4d", C["red"], C["green"]) else "black")
        left += w
ax3.set_yticks(range(len(rows))); ax3.set_yticklabels([r[0] for r in rows], fontsize=8.5)
ax3.invert_yaxis(); ax3.set_xlim(0, 100); ax3.set_xlabel("Share of works by research field (%)")
ax3.tick_params(axis="y", length=0); ax3.minorticks_off()
ax3.legend(loc="upper center", bbox_to_anchor=(0.5, -0.42), ncol=3, fontsize=8.5, handlelength=1.2, columnspacing=1.2)
label_panels([ax1, ax2, ax3], x=-0.22)

save(fig, "fig02_usage")
