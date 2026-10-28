"""Load every CSV in data/tables/ and check row counts and a few values against the article.

Run from anywhere: python scripts/check_tables.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "data" / "tables"

# number of rows of every table, as in the article and its supplement
EXPECTED_ROWS = {
    "tableI_materials.csv": 12,
    "tableII_field_reviews.csv": 14,
    "tableIII_milestones.csv": 12,
    "tableIV_dimensionless_groups.csv": 6,
    "tableV_notation.csv": 11,
    "tableVI_predicted_exponents.csv": 8,
    "tableVII_yield_stress_protocols.csv": 7,
    "tableVIII_verification_problems.csv": 9,
    "tableIX_error_budget.csv": 5,
    "tableX_models_compared.csv": 10,
    "tableXI_applicability_coordinates.csv": 10,
    "tableXII_transfer_matrix.csv": 21,
    "tableXIII_predicted_vs_measured.csv": 7,
    "tableSI_foam_model_scores.csv": 3,
    "tableSII_normal_stress_coefficients.csv": 3,
    "tableSIII_pipe_triplets.csv": 3,
    "tableSIV_onset_exponents.csv": 5,
    "tableSV_regression_checks.csv": 16,
    "tableSVI_unit_conversions.csv": 7,
    "tableSVII_digital_record.csv": 7,
}


def load(name):
    return pd.read_csv(TABLES / name, comment="#", keep_default_na=False)


def main():
    files = sorted(p.name for p in TABLES.glob("*.csv"))
    assert set(files) == set(EXPECTED_ROWS), set(files) ^ set(EXPECTED_ROWS)
    tables = {}
    for name in EXPECTED_ROWS:
        df = load(name)
        first = (TABLES / name).read_text(encoding="utf-8").splitlines()[0]
        assert first.startswith("# Table "), f"{name}: missing caption comment"
        print(f"{name:45s} {len(df):3d} rows x {df.shape[1]} columns")
        assert len(df) == EXPECTED_ROWS[name], f"{name}: {len(df)} rows, expected {EXPECTED_ROWS[name]}"
        tables[name] = df

    # Table I: 12 material classes, all with at least one citation key
    t1 = tables["tableI_materials.csv"]
    assert len(t1) == 12 and (t1["bibkeys"] != "").all()
    assert t1["material"].iloc[0] == "Concentrated emulsions"

    # Table X: HB has three parameters and the Bingham and power-law models are its special cases
    t10 = tables["tableX_models_compared.csv"].set_index("model")
    assert int(t10.loc["Herschel–Bulkley", "p"]) == 3 and int(t10.loc["Herschel–Bulkley", "year"]) == 1926
    assert t10.loc["Bingham", "relation_to_hb"] == "Special (n=1)"
    assert t10.loc["Ostwald–de Waele", "relation_to_hb"] == "Special (σ_y=0)"

    # Table XII: seven groups A-G, 21 rows
    t12 = tables["tableXII_transfer_matrix.csv"]
    assert list(t12["group"].unique()) == list("ABCDEFG")

    # Table SI: HB is the reference model (ΔAIC = ΔBIC = 0), Bingham the worst
    s1 = tables["tableSI_foam_model_scores.csv"].set_index("model")
    assert s1.loc["HB", "delta_aic"] == 0 and s1.loc["HB", "delta_bic"] == 0
    assert s1.loc["Bingham", "delta_aic"] == 271 and s1.loc["HB", "k"] == 4
    assert abs(s1.loc["HB", "rmse_pa"] - 1.77) < 1e-12

    # Table SII: the ratio column equals psi2/psi1 to the printed two decimals
    s2 = tables["tableSII_normal_stress_coefficients.csv"]
    ratio = s2["psi2_x1e3_per_pa"] / s2["psi1_x1e3_per_pa"]
    assert (abs(ratio.round(2) - s2["psi2_over_psi1"]) < 1e-9).all()

    # Table SIII: rotational yield stress 1.198 Pa; the printed differences follow from the triplets
    s3 = tables["tableSIII_pipe_triplets.csv"].set_index("identification_route")
    rot = s3.loc["Rotational rheometry (corrected)", "sigma_y_pa"]
    assert rot == 1.198
    for route, printed in (("Pipe, nonlinear optimization", -14.0), ("Pipe, corrected method", -24.1)):
        d = 100 * (s3.loc[route, "sigma_y_pa"] - rot) / rot
        assert abs(d - printed) < 0.05, (route, d, printed)
        assert float(s3.loc[route, "delta_sigma_y_over_sigma_y_rot_percent"]) == printed

    # Table SV: 16 analytical limits, as stated in Sec. VII of the article
    assert len(tables["tableSV_regression_checks.csv"]) == 16

    # Table SVI: lb_f/100 ft^2 to Pa
    s6 = tables["tableSVI_unit_conversions.csv"]
    assert "0.4788 Pa" in s6["si_equivalent"].values

    print("all table checks passed")


if __name__ == "__main__":
    main()
