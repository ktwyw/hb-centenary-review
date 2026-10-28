"""Regenerate every figure of the article (Figs. 1-13) and of its supplementary material (Figs. S1-S7).

Run from the repository root:  python scripts/make_all_figures.py
Each script is run in its own process; the run stops at the first error.
"""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
MAIN = ["fig01_hb1926", "fig02_usage", "fig03_parameter_atlas", "fig04_flow_curves", "fig05_microstructure", "fig06_couette_inverse", "fig07_fitting", "fig08_pipe_profiles", "fig09_transport_exponents", "fig10_slip_network", "fig11_regularization", "fig12_hbp_pipe", "fig13_transfer_matrix"]
NOTES = ["figS1_foam_reanalysis", "figS2_normal_stress", "figS3_pipe_uncertainty", "figS4_mullineux", "figS5_uncertainty", "figS6_flow_rate", "figS7_friction_factor"]
for name in MAIN + NOTES:
    print(f"Generating {name}", flush=True)
    subprocess.run([sys.executable, str(ROOT / "scripts" / f"{name}.py")], cwd=ROOT, check=True)
