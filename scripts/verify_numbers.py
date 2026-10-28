"""Check that numbers quoted in the article match the stored analysis results.

Run from the repository root after make_all_figures.py:  python scripts/verify_numbers.py
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
D = lambda f: json.loads((ROOT / "data" / f).read_text())
checks = []

def check(label, value, target, tol):
    ok = abs(value - target) <= tol
    checks.append((label, value, target, ok))

# Sec. II: 1926 data (Fig. 1)
h = D("hb1926_reanalysis.json")
for f, n in zip(h["pooled_fits"], [0.81, 0.74, 0.68]):
    check(f"1926 pooled fit sigma_y {f['pair'][0][1]}%", f["sigma_y"], 0.0, 1e-6)
    check(f"1926 pooled fit n {f['pair'][0][1]}%", f["n"], n, 0.005)
    assert f["rms_log"] < 0.035
r = [100 * f["rms_log"] for f in h["pooled_fits"]]
check("1926 pooled-fit rms deviation min %", min(r), 1.7, 0.05)
check("1926 pooled-fit rms deviation max %", max(r), 3.2, 0.05)

# Note S2: foam reanalysis
fo = D("foam_reanalysis_results.json")
check("foam sigma_y", fo["HB"]["parameters"][0], 30.3, 0.05)
check("foam n", fo["HB"]["parameters"][2], 0.36, 0.005)
check("foam sigma_y 95% low", fo["HB"]["profile_95_sigma_y_Pa"][0], 28.0, 0.05)
check("foam sigma_y 95% high", fo["HB"]["profile_95_sigma_y_Pa"][1], 32.3, 0.06)
check("foam n 95% low", fo["HB"]["profile_95_n"][0], 0.325, 0.001)
check("foam n 95% high", fo["HB"]["profile_95_n"][1], 0.395, 0.001)
check("dAIC Casson", fo["Casson"]["AIC"] - fo["HB"]["AIC"], 90, 1)
check("dAIC Bingham", fo["Bingham"]["AIC"] - fo["HB"]["AIC"], 270, 1)
check("gc upper bound", fo["HB_gc"]["upper_bound_95_s_inv"], 0.047, 0.0005)

# Note S3: normal stresses
ns = D("normal_stress_reconstruction.json")["materials"]
inc = [m["equivalent_stress_increment_at_intercept_percent"] for m in ns]
check("invariant increase min %", min(inc), 0.39, 0.005)
check("invariant increase max %", max(inc), 0.56, 0.005)
r = [m["psi2_over_psi1"] for m in ns]
check("psi2/psi1 max", max(r), -0.49, 0.005)
check("psi2/psi1 min", min(r), -0.86, 0.005)
n1 = [m["N1_over_sigma_at_intercept"] for m in ns]
check("N1/sigma_y min", min(n1), 0.17, 0.005)
check("N1/sigma_y max", max(n1), 0.20, 0.005)
# tensor identity used for the invariant
rng = np.random.default_rng(1)
for s, a, b in rng.normal(0, 100, (50, 3)):
    t = np.diag([(2 * a + b) / 3, (-a + b) / 3, (-a - 2 * b) / 3]); t[0, 1] = t[1, 0] = s
    assert abs(np.sqrt(np.sum(t * t) / 2) - np.sqrt(s * s + (a * a + a * b + b * b) / 3)) < 1e-9

# Note S5: pipe
p = D("pipe_supplement_results.json")
y = p["published_triplets"]["yield_differences_percent"]
check("pipe yield diff 1 %", y[1], -14.0, 0.05)
check("pipe yield diff 2 %", y[2], -24.1, 0.05)
d = p["published_triplets"]["pressure_difference_percent_at_Q_L_min"]
check("pressure diff at 10 L/min (max abs)", max(abs(d["10"]["nonlinear_fit"]), abs(d["10"]["corrected_method"])), 1.3, 0.05)
check("pressure diff at ~1 L/min low", abs(d["1"]["nonlinear_fit"]), 5.9, 0.05)
check("pressure diff at ~1 L/min high", abs(d["1"]["corrected_method"]), 6.8, 0.05)
check("pressure diff at 0.01 L/min low", abs(d["0.01"]["nonlinear_fit"]), 13, 0.3)
check("pressure diff at 0.01 L/min high", abs(d["0.01"]["corrected_method"]), 20, 0.5)
s = p["synthetic"]["results"]
for key, lo, hi in [("broad_window", 0.84, 1.14), ("high_rate_window", 0.03, 2.22), ("cutoff_1s_window", 0.74, 1.22)]:
    rr = s[key]["pressure_ratio_95_range_at_true_wall_rate"]["0.01"]
    check(f"{key} 95% low", rr[0], lo, 0.006); check(f"{key} 95% high", rr[1], hi, 0.006)
check("zero-yield boundary fits", s["high_rate_window"]["zero_yield_boundary_fits"], 30, 0)


# Sec. III: bibliometrics (OpenAlex, Oct 2026)
b = D("hb_bibliometrics.json")
check("works, title/abstract", b["title_abstract_total"], 4883, 0)
check("works, full text", b["fulltext_total"], 12300, 0)
check("1926 paper citations", b["hb1926_citations_total"], 1513, 0)
ftb = np.array(b["fulltext_by_year"]); m = (ftb[:, 0] >= 1995) & (ftb[:, 0] <= 2025)
check("full-text doubling time 1995-2025 (yr)", np.log(2) / np.polyfit(ftb[m, 0], np.log(ftb[m, 1]), 1)[0], 5.5, 0.05)
check("share of full-text mentions since 2000 (%)", 100 * ftb[ftb[:, 0] >= 2000, 1].sum() / ftb[:, 1].sum(), 97, 0.5)
check("min full-text mentions/yr 2023-2025", min(v for y, v in b["fulltext_by_year"] if 2023 <= y <= 2025), 966, 0)
c = [v for y, v in b["hb1926_citations_by_year"] if 2020 <= y <= 2026]
check("min citations/yr since 2020", min(c), 101, 0)


# Sec. III B: breadth of use compared with the Bingham and Casson models (OpenAlex)
mc = b["model_comparison"]["models"]
for key, works, share, nf, nsub in [("HB", 4883, 28, 22, 84), ("Bingham", 5179, 15, 18, 75), ("Casson", 3732, 16, 19, 54)]:
    check(f"{key} works", mc[key]["works"], works, 0)
    check(f"{key} share outside engineering (%)", mc[key]["share_outside_engineering_pct"], share, 0.5)
    check(f"{key} fields", mc[key]["fields"], nf, 0)
    check(f"{key} subfields", mc[key]["subfields"], nsub, 0)
ag = {k: dict(mc[k]["by_field"]).get("Agricultural and Biological Sciences", 0) for k in mc}
check("HB agri/bio works", ag["HB"], 458, 0); check("Bingham agri/bio works", ag["Bingham"], 65, 0); check("Casson agri/bio works", ag["Casson"], 101, 0)


# Sec. III B / Fig. 3: parameter atlas (values as reported by the cited sources)
import csv
at = [r for r in csv.DictReader(open(ROOT / "data" / "hb_parameter_atlas.csv")) if r["class"] != "1926 rubber solutions"]
sy = [float(r["sigma_y_Pa"]) for r in at] + [float(r["sigma_y_max_Pa"]) for r in at if r["sigma_y_max_Pa"]]
check("atlas min yield stress (Pa)", min(sy), 0.2, 1e-9); check("atlas max yield stress (Pa)", max(sy), 309, 1e-9)
nn = [float(r["n"]) for r in at if r["source_key"] != "SaasenYtrehus2020"]
check("atlas min n (excl. drilling fluid of Saasen)", min(nn), 0.24, 1e-9); check("atlas max n (same)", max(nn), 0.71, 1e-9)
ng = [float(r["n"]) for r in at if r["class"] in ("Emulsions", "Carbopol and polymer gels") or r["material"].startswith("SLES")]
check("foams/emulsions/gels min n", min(ng), 0.31, 0.011); check("foams/emulsions/gels max n", max(ng), 0.64, 0.011)

w = max(len(c[0]) for c in checks)
for label, v, t, ok in checks:
    print(f"{label:<{w}}  {v:10.4g}  quoted {t:<8g} {'OK' if ok else 'MISMATCH'}")
bad = [c for c in checks if not c[3]]
print(f"\n{len(checks) - len(bad)}/{len(checks)} quoted numbers match")
raise SystemExit(1 if bad else 0)
