# Data

All files are inputs to, or outputs of, the scripts in `scripts/`. None is a new experimental measurement.

The folder [`figure_data/`](figure_data/README.md) holds the plotted data of every figure of the article and of its supplementary material, one CSV file per figure, written when the figure is made. The folder [`tables/`](tables/README.md) holds every numbered table of the article (I–XIII) and of the supplementary material (SI–SVII) as CSV files.

| File(s) | Content | Generator |
|---|---|---|
| `hb1926_table3.csv` | Flow rate vs corrected pressure for the two Bingham viscometers (Nr. 3, d = 0.1018 cm, L = 10.0 cm; Nr. 2, d = 0.0689 cm, L = 14.6 cm), 30 °C, transcribed from Table III of Herschel & Bulkley, Kolloid-Z. 39, 291 (1926). | transcription |
| `hb1926_table4.csv` | Fitted constants (n, I, K) of Table IV of the same paper. | transcription |
| `hb1926_reanalysis.json` | Wall stress and wall shear rate (Weissenberg–Rabinowitsch–Mooney, no slip) and HB fits to the pooled capillary data at matched concentration (article Fig. 1). | `fig01_hb1926.py` |
| `hb_bibliometrics.json` | OpenAlex counts (queried 25 October 2026) of works mentioning "Herschel-Bulkley" by year (title/abstract and full text), by primary topic and by source, and citations of the 1926 Kolloid-Z. paper; query strings included (article Fig. 2, Sec. III). | OpenAlex API (queries stored in the file); read by `fig02_usage.py` |
| `foam_sles_fig3b_digitized.csv` | 152 local (shear rate, stress) points digitized from Fig. 3b of Ovarlez, Krishan & Cohen-Addad, EPL 91, 68005 (2010) (arXiv:1007.4133v2). Plot-derived values. | `digitize_foam.py` (needs a separately obtained copy of the source PDF) |
| `foam_digitization_provenance.json`, `foam_selection_audit.json` | Axis calibration, selected marker paths and source-file digest for the digitization. `digitize_foam.py` also writes a check image, `foam_extraction_overlay.png`, cropped from the source figure; it is excluded from the repository (`.gitignore`) because it reproduces part of a published figure. | `digitize_foam.py` |
| `foam_reanalysis_results.json` | HB, Bingham, Casson and critical-shear-rate fits with errors in stress and rate; profile and bootstrap intervals (Supplementary Note S2). | `figS1_foam_reanalysis.py` |
| `normal_stress_reconstruction.json`, `normal_stress_curves.csv` | Consequences of the published fits N_i = ψ_i σ² of de Cagny et al., J. Rheol. 63, 285 (2019) (Supplementary Note S3). | `figS2_normal_stress.py` |
| `pipe_supplement_results.json` | Pipe-flow predictions from the triplets of Magnon & Cayeux, Fluids 6, 157 (2021), and the synthetic calibration-window study (Supplementary Note S5). | `figS3_pipe_uncertainty.py` |
| `broad_window_*.csv`, `high_rate_window_*.csv`, `cutoff_1s_window_*.csv` | Synthetic fitted triplets and pressure-ratio ensembles (seed 20260920). | `figS3_pipe_uncertainty.py` |
| `hb_parameter_atlas.csv` | Herschel–Bulkley yield stress, flow index and (where reported) consistency for 26 materials, taken from cited studies or from this review's reanalyses; each row names the source and the table, section or data file the values come from (article Fig. 3). Values are as reported, not refitted. | transcription |
| `couette_illustration.json`, `fitting_illustration.json` | Coefficients for article Figs. 6 and 7 (synthetic). | `fig06_couette_inverse.py`, `fig07_fitting.py` |

Values digitized from, or computed from coefficients printed in, published articles are facts reported by those articles; cite the original sources when reusing them.
