"""Fig. 1 of the article (fig01_hb1926.pdf): the 1926 capillary data of Herschel and Bulkley in modern variables.

Data: Table III of Herschel & Bulkley, Kolloid-Z. 39, 291 (1926), Bingham
viscometers Nr. 3 (d = 0.1018 cm, l = 10.0 cm) and Nr. 2 (d = 0.0689 cm,
l = 14.6 cm), 30 C, transcribed from the printed page (data/hb1926_table3.csv).
P is the corrected driving pressure in g(force)/cm^2, q the flow rate in cm^3/s.

(a) The 1926 representation: log q against log P with the printed fits
    q = (d^4/(L I)) (P - K)^n of Table IV.
(b) Modern representation: wall shear stress sigma_w = P d / (4 L) against the
    true wall shear rate from the Weissenberg-Rabinowitsch-Mooney correction,
    for the two instruments; no slip assumed. Curves are Herschel-Bulkley fits
    (sigma_y >= 0) to the pooled data of the two instruments at matched
    concentration.
"""
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from hbstyle import *
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
G = 980.665          # cm s^-2 ; 1 gf/cm^2 = 98.0665 Pa
DIM = {"B3": (0.1018, 10.0), "B2": (0.0689, 14.6)}
raw = np.genfromtxt(ROOT / "data" / "hb1926_table3.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
t4 = np.genfromtxt(ROOT / "data" / "hb1926_table4.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")

def series():
    out = {}
    for inst, c, P, q in raw:
        out.setdefault((inst, float(c)), []).append((float(P), float(q)))
    return {k: np.array(sorted(v)) for k, v in out.items()}

S = series()
res = {"note": "Wall stress and Rabinowitsch-Mooney wall rate from Table III of Herschel & Bulkley (1926); no slip assumed.", "series": [], "pooled_fits": []}
mod = {}
for (inst, c), arr in S.items():
    d, L = DIM[inst]
    P, q = arr[:, 0], arr[:, 1]
    sw = P * 98.0665 * d / (4 * L)                      # Pa
    gapp = 32 * q / (np.pi * d**3)                       # s^-1
    sl = np.polyfit(np.log(gapp), np.log(sw), 1)[0]      # n' (local, constant fit)
    # local slope from quadratic fit in log space
    cfit = np.polyfit(np.log(gapp), np.log(sw), 2)
    nprime = 2 * cfit[0] * np.log(gapp) + cfit[1]
    gw = gapp * (3 * nprime + 1) / (4 * nprime)
    mod[(inst, c)] = (sw, gw)
    res["series"].append({"instrument": inst, "conc_pct": c, "n_points": int(len(P)),
                          "sigma_w_range_Pa": [round(float(sw.min()), 2), round(float(sw.max()), 2)],
                          "gdot_w_range_per_s": [round(float(gw.min()), 0), round(float(gw.max()), 0)],
                          "nprime_mean": round(float(sl), 3)})

def hbfit(sw, gw):
    def r(p):
        ty, lk, n = p
        return np.log(ty + np.exp(lk) * gw**n) - np.log(sw)
    best = None
    for ty0 in [0.0, 0.5, 2.0]:
        sol = least_squares(r, [ty0, np.log(0.02), 0.8], bounds=([0, -20, 0.2], [np.min(sw), 5, 1.5]))
        if best is None or sol.cost < best.cost:
            best = sol
    ty, lk, n = best.x
    rms = float(np.sqrt(np.mean(best.fun**2)))
    # power-law fit (sigma_y = 0) for comparison
    pl = np.polyfit(np.log(gw), np.log(sw), 1)
    return dict(sigma_y=float(ty), K=float(np.exp(lk)), n=float(n), rms_log=rms,
                powerlaw_K=float(np.exp(pl[1])), powerlaw_n=float(pl[0]))

pairs = [(("B3", 0.693), ("B2", 0.692)), (("B3", 0.993), ("B2", 1.004)), (("B3", 1.21), ("B2", 1.19))]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=DOUBLE)
cols = [C["black"], C["red"], C["blue"]]
for i, (a, b) in enumerate(pairs):
    for j, key in enumerate((a, b)):
        if key not in S:
            continue
        arr = S[key]
        mk = "o" if key[0] == "B3" else "s"
        mfc = cols[i] if key[0] == "B3" else "white"
        ax1.loglog(arr[:, 0], arr[:, 1], mk, color=cols[i], mfc=mfc, ms=3.5, mew=0.7,
                   gid=f"{key[1]:.3f}%, {'Nr. 3' if key[0]=='B3' else 'Nr. 2'}: data of 1926 Table III")
        row = [r for r in t4 if r["instrument"] == ("Bingham 3" if key[0] == "B3" else "Bingham 2") and abs(r["conc_pct"] - key[1]) < 1e-3]
        if row:
            n26, I26, K26 = float(row[0]["exponent_n"]), float(row[0]["I"]), float(row[0]["K_g_per_cm2"])
            d, L = DIM[key[0]]
            Pg = np.logspace(np.log10(max(arr[:, 0].min(), K26 + 1)), np.log10(arr[:, 0].max()), 50)
            ax1.loglog(Pg, d**4 / (L * I26) * (Pg - K26)**n26, ["-", "--"][j], color=cols[i], lw=0.8,
                       gid=f"{key[1]:.3f}%, {'Nr. 3' if key[0]=='B3' else 'Nr. 2'}: 1926 fit, Table IV constants")
        sw, gw = mod[key]
        ax2.loglog(gw, sw, mk, color=cols[i], mfc=mfc, ms=3.5, mew=0.7,
                   label=f"{key[1]:.3f}%, {'Nr. 3' if key[0]=='B3' else 'Nr. 2'}")
    sw = np.concatenate([mod[a][0], mod[b][0]]); gw = np.concatenate([mod[a][1], mod[b][1]])
    f = hbfit(sw, gw)
    f.update({"pair": [list(a), list(b)]})
    res["pooled_fits"].append(f)
    gg = np.logspace(np.log10(gw.min()), np.log10(gw.max()), 60)
    ax2.loglog(gg, f["sigma_y"] + f["K"] * gg**f["n"], "-", color=cols[i], lw=0.8,
               gid=f"HB fit to the pooled data, {a[1]:.3f}% and {b[1]:.3f}%")
ax1.set_xlabel("$P$ (g cm$^{-2}$)")
ax1.set_ylabel("$q$ (cm$^{3}$ s$^{-1}$)")
ax2.set_xlabel("$\\dot\\gamma_w$ (s$^{-1}$)")
ax2.set_ylabel("$\\sigma_w$ (Pa)")
ax2.legend(loc="lower right", fontsize=8.5, ncol=1)
from matplotlib.ticker import NullFormatter
for ax in (ax1, ax2):
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
ax2.set_yticks([3, 10, 30, 70]); ax2.set_yticklabels(["3", "10", "30", "70"])
ax1.set_xticks([20, 50, 100, 200, 500]); ax1.set_xticklabels(["20", "50", "100", "200", "500"])
label_panels([ax1, ax2])
fig.tight_layout(w_pad=2.2)
save(fig, "fig01_hb1926")
(ROOT / "data" / "hb1926_reanalysis.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res["pooled_fits"], indent=1))
for s in res["series"]:
    print(s)
