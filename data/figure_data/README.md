# Plotted data of every figure

One CSV file per figure of the article (Figs. 1–13) and of its supplementary material (Figs. S1–S7),
written automatically by `scripts/figdata.py` each time a figure script runs, so the files always
match the published figures. Schematic drawings (Fig. 5a–c and Fig. 8a) contain no data and are not
exported.

## Columns (long format, one row per plotted point)

| Column | Meaning |
|---|---|
| `panel` | panel letter as printed, e.g. `(a)`; `(b) right axis` for a second y axis; `(single panel)` |
| `x_label`, `y_label` | axis labels as typeset (Matplotlib mathtext); `[log]` marks a logarithmic axis |
| `series` | legend label, or a description of the curve (e.g. `0.693%, Nr. 3: data of 1926 Table III`) |
| `kind` | `line`, `markers`, `bar`, `band`, `errorbar`, `segment`, `hline` (horizontal reference line, value in `y`), `vline` (value in `x`) |
| `index` | point number within the series |
| `x`, `y` | data coordinates, in the units of the axis labels |
| `x_end`, `y_end` | second corner of a bar or band, or end of a segment |
| `category` | category label of a bar (e.g. the application area in Fig. 2b) |

Read a file with, for example,

```python
import pandas as pd
df = pd.read_csv("data/figure_data/fig08_pipe_profiles.csv", comment="#")
df[(df.panel == "(b)") & (df.series == r"$\phi=0.5$")].plot(x="x", y="y")
```

Two figures are better described by a table, and two contain maps drawn as contours:

- `fig03_parameter_atlas.csv` — the rows of `data/hb_parameter_atlas.csv` plotted in Fig. 3, with the
  source article and the table or section the values come from; the last row is the range of flow
  indices of fresh concrete drawn as the right-hand strip.
- `fig13_transfer_matrix.csv` — the matrix of Fig. 13: group, study, what changes between
  calibration and test, and the finding (Table XII of the article).
- `fig07_fitting_rms_residual_grid.csv` — the relative RMS residual over the (σ_y, n) plane drawn as
  filled contours in Fig. 7b.
- `figS1_foam_joint_deviance_grid.csv` — the joint profile deviance over the (σ_y, n) plane whose
  68% and 95% contours are drawn in Fig. S1c.

## Origin of the plotted values

| File | Figure | Origin of the values |
|---|---|---|
| `fig01_hb1926.csv` | Fig. 1 | (a) data transcribed from Table III of Herschel & Bulkley (1926) and their Table IV fits; (b) the same data in wall variables, and HB fits made for this review |
| `fig02_usage.csv` | Fig. 2 | OpenAlex counts (queried 25 October 2026; queries in `data/hb_bibliometrics.json`) |
| `fig03_parameter_atlas.csv` | Fig. 3 | values as reported by the cited sources, or by this review's reanalyses (one row per point, bar or strip) |
| `fig04_flow_curves.csv` | Fig. 4 | computed from the model (illustrative parameters given in the caption) |
| `fig05_microstructure.csv` | Fig. 5d | the n = 1/2 scaling of Seth et al. (2011) |
| `fig06_couette_inverse.csv` | Fig. 6 | computed (illustrative calculation) |
| `fig07_fitting.csv`, `fig07_fitting_rms_residual_grid.csv` | Fig. 7 | synthetic data and fits made for this review |
| `fig08_pipe_profiles.csv` | Fig. 8 | exact laminar pipe solution |
| `fig09_transport_exponents.csv` | Fig. 9 | analytical calculations |
| `fig10_slip_network.csv` | Fig. 10 | analytical calculations |
| `fig11_regularization.csv` | Fig. 11 | regularized and exact constitutive laws |
| `fig12_hbp_pipe.csv` | Fig. 12 | numerical pipe solutions of the HB–Papanastasiou law |
| `fig13_transfer_matrix.csv` | Fig. 13 | findings of the cited experimental studies, as summarized in Table XII |
| `figS1_foam.csv`, `figS1_foam_joint_deviance_grid.csv` | Fig. S1 | points **digitized** from Fig. 3b of Ovarlez, Krishan & Cohen-Addad (2010), and fits made for this review |
| `figS2_normal_stress.csv` | Fig. S2 | **computed from published coefficients** (fits of de Cagny et al., 2019) |
| `figS3_pipe.csv` | Fig. S3 | (a) **computed from the published parameter triplets** of Magnon & Cayeux (2021); (b) synthetic calibration ensembles |
| `figS4_mullineux.csv` | Fig. S4 | synthetic data |
| `figS5_uncertainty.csv` | Fig. S5 | synthetic data |
| `figS6_flow_rate.csv` | Fig. S6 | closed-form pipe-flow relations |
| `figS7_friction_factor.csv` | Fig. S7 | **computed from published correlations** (Metzner & Reed, 1955; Dodge & Metzner, 1959) |

Values digitized from, or computed from coefficients printed in, published articles remain the
results of those articles: cite the original source when you reuse them.
