# A century of the Herschel–Bulkley model (1926–2026)

Code and data for

> Y. Wang, "A century of the Herschel–Bulkley model (1926–2026): Widespread use, community achievements, predictive capability, and engineering implications," submitted to *Journal of Rheology* (2026).

Yanwei Wang · Department of Chemical and Materials Engineering, School of Engineering, and Institute of New Materials and Energy Technologies, Nazarbayev University, Astana, Kazakhstan · yanwei.wang@nu.edu.kz · ORCID [0000-0002-8488-9833](https://orcid.org/0000-0002-8488-9833)

This repository contains everything named in the article's Data Availability statement:

- **Python source for every figure** of the article (Figs. 1–13) and of its supplementary material (Figs. S1–S7), one script per figure, named after the figure it draws.
- **Data** read and written by the scripts, with the origin of each file in [`data/README.md`](data/README.md). Values digitized from published figures and curves computed from published coefficients are marked as such there and in the figure captions.
- **The plotted data of every figure** of the article and of its supplementary material, one CSV file per figure, in [`data/figure_data/`](data/figure_data/README.md), and **every numbered table** of both documents as CSV files in [`data/tables/`](data/tables/README.md).
- **Eighteen tutorial notebooks** that walk through the calculations of the review and of its Supplementary Notes, written for students and newcomers: [`notebooks/`](notebooks/README.md).
- **Bibliographic counts** of Sec. III of the article (OpenAlex), with the query strings: [`data/hb_bibliometrics.json`](data/hb_bibliometrics.json).
- **Checks**: `verify_numbers.py` recomputes 65 numbers quoted in the text from the stored data; `verify_equations.py` and two further scripts test the analytical results (their reports are written to `results/`); `check_tables.py` loads the table files and checks spot values.

The Supplementary Notes (S1–S19) are published with the article by the journal as supplementary material and are not reproduced here; the scripts and data for their figures and worked examples are included. No new experimental measurements are reported. Source articles are cited, not redistributed.

## Figures

| Figure | Script | Content |
|---|---|---|
| Fig. 1 | [`fig01_hb1926.py`](scripts/fig01_hb1926.py) | The 1926 capillary data of Herschel and Bulkley in modern variables |
| Fig. 2 | [`fig02_usage.py`](scripts/fig02_usage.py) | Publications per year, by application area, and by research field compared with the Bingham and Casson models (OpenAlex) |
| Fig. 3 | [`fig03_parameter_atlas.py`](scripts/fig03_parameter_atlas.py) | Yield stress and flow index for materials from several fields (`data/hb_parameter_atlas.csv`) |
| Fig. 4 | [`fig04_flow_curves.py`](scripts/fig04_flow_curves.py) | HB flow curves and limiting cases |
| Fig. 5 | [`fig05_microstructure.py`](scripts/fig05_microstructure.py) | Microstructural origins of HB behavior (schematic) |
| Fig. 6 | [`fig06_couette_inverse.py`](scripts/fig06_couette_inverse.py) | Partially yielded Couette gap and inversion bias |
| Fig. 7 | [`fig07_fitting.py`](scripts/fig07_fitting.py) | Ill-conditioning of the three-parameter fit |
| Fig. 8 | [`fig08_pipe_profiles.py`](scripts/fig08_pipe_profiles.py) | Laminar pipe flow: stress and velocity profiles |
| Fig. 9 | [`fig09_transport_exponents.py`](scripts/fig09_transport_exponents.py) | Near-onset transport exponents |
| Fig. 10 | [`fig10_slip_network.py`](scripts/fig10_slip_network.py) | Wall slip and parallel channels |
| Fig. 11 | [`fig11_regularization.py`](scripts/fig11_regularization.py) | Regularized constitutive laws |
| Fig. 12 | [`fig12_hbp_pipe.py`](scripts/fig12_hbp_pipe.py) | Convergence of the HB–Papanastasiou model in pipe flow |
| Fig. 13 | [`fig13_transfer_matrix.py`](scripts/fig13_transfer_matrix.py) | Overview of experimental tests of the model beyond the calibration flow |
| Fig. S1 | [`figS1_foam_reanalysis.py`](scripts/figS1_foam_reanalysis.py) | Worked example: a published local foam flow curve |
| Fig. S2 | [`figS2_normal_stress.py`](scripts/figS2_normal_stress.py) | Normal-stress consequences of published fits |
| Fig. S3 | [`figS3_pipe_uncertainty.py`](scripts/figS3_pipe_uncertainty.py) | Published pipe triplets and the synthetic calibration-window study |
| Fig. S4 | [`figS4_mullineux.py`](scripts/figS4_mullineux.py) | One-dimensional reduction of the least-squares fit |
| Fig. S5 | [`figS5_uncertainty.py`](scripts/figS5_uncertainty.py) | Constrained yield-stress profiles |
| Fig. S6 | [`figS6_flow_rate.py`](scripts/figS6_flow_rate.py) | Closed-form pipe-flow relations |
| Fig. S7 | [`figS7_friction_factor.py`](scripts/figS7_friction_factor.py) | Fanning friction factor |

Each script writes `figures/<name>.pdf` (vector, as used in the article), a PNG preview, and `data/figure_data/<name>.csv` with every plotted series; `<name>` is the script name, except that `figS1_foam_reanalysis.py` writes `figS1_foam` and `figS3_pipe_uncertainty.py` writes `figS3_pipe`. Shared style and helpers are in `hbstyle.py` and `hb_prediction.py`. All figures use the colour-blind-safe Okabe–Ito palette (Okabe & Ito 2008) with distinct line styles, so they also read in grayscale.

## Reproduce

```sh
python -m pip install -r requirements.txt
python scripts/make_all_figures.py        # all figures and data/*.json (about 5 min, mostly the foam reanalysis)
python scripts/verify_numbers.py          # numbers quoted in the article vs. data/*.json
python scripts/verify_equations.py        # analytical identities and limits
python scripts/check_tables.py            # the table CSVs in data/tables/ load and match spot values
python scripts/verify_prediction.py
python scripts/verify_transfer_diagnostic.py
```

Tested with Python 3.13, NumPy 2.5, SciPy 1.18 and Matplotlib 3.11 (versions in `requirements.txt`). Random seeds are fixed in the scripts; regenerated numbers agree to optimizer tolerance, not bit for bit.

## Related repository

[**herschel-bulkley-determinant-fit**](https://github.com/ktwyw/herschel-bulkley-determinant-fit) is a tested Python package (`hb_mullineux`) and a browser app ([live](https://ktwyw.github.io/herschel-bulkley-determinant-fit/web/herschel_bulkley_fit.html)) for fitting the Herschel–Bulkley model to your own flow-curve data by the determinant method of Mullineux (2008). The flow index is found first as the root of a one-variable function, and the yield stress and consistency then follow from a linear solve. The package reports confidence intervals and goodness-of-fit statistics, and tests the fit against the Newtonian, power-law and Bingham models. It is the practical companion to the discussion of parameter identifiability in Sec. VI D of the article and to Fig. S4 ([`figS4_mullineux.py`](scripts/figS4_mullineux.py)), which use the same reduction of the least-squares problem.

## Sources of published data

- Herschel & Bulkley, Kolloid-Z. **39**, 291 (1926): Tables III and IV, transcribed.
- Ovarlez, Krishan & Cohen-Addad, EPL **91**, 68005 (2010): Fig. 3b, digitized (`digitize_foam.py` requires a separately obtained copy of arXiv:1007.4133v2; the PDF is not redistributed).
- de Cagny et al., J. Rheol. **63**, 285 (2019): fitted normal-stress coefficients.
- Magnon & Cayeux, Fluids **6**, 157 (2021): published parameter triplets.
- OpenAlex (https://openalex.org), queried 25 October 2026: publication counts; query strings in `data/hb_bibliometrics.json`.
- Metzner & Reed, AIChE J. **1**, 434 (1955), and Dodge & Metzner, AIChE J. **5**, 189 (1959): laminar and turbulent friction-factor correlations with their published coefficients (Fig. S7).
- Seth et al., Nat. Mater. **10**, 838 (2011): the n = 1/2 elastohydrodynamic scaling drawn in Fig. 5(d).
- Fig. 13 encodes the findings of the experimental studies listed in Table XII of the article; no data are computed.
- Herschel–Bulkley parameters compiled in `data/hb_parameter_atlas.csv` (article Fig. 3): each row names its source article and the table, section or file the values come from.

Every figure that redraws or reanalyzes published data cites the original article (and figure or table) in its caption.

## Licences

* Code, including the code cells of the notebooks: MIT licence (`LICENSE`).
* Data generated in this work: the CSV and JSON data files, the figures and the explanatory text of the notebooks are released under CC BY 4.0 (`DATA_LICENSE.md`).
* Third-party data: values transcribed or digitized from published articles, and curves computed from published coefficients, remain the results of those articles; the conditions of the original sources apply (see `DATA_LICENSE.md` and `data/README.md`).

## Citation

If you use this code or data, please cite the paper above and this repository (`CITATION.cff`):

> Y. Wang, Code and data for "A century of the Herschel–Bulkley model (1926–2026)", GitHub repository, https://github.com/ktwyw/hb-centenary-review (2026).

## Contact

Yanwei Wang, Nazarbayev University: [yanwei.wang@nu.edu.kz](mailto:yanwei.wang@nu.edu.kz), ORCID [0000-0002-8488-9833](https://orcid.org/0000-0002-8488-9833).
