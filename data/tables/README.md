# Tables of the article and of its supplementary material

One CSV file per numbered table of the article (Tables I–XIII) and of its supplementary material
(Tables SI–SVII). Each file transcribes one table row per CSV row, with columns named after the
table headers (ASCII snake_case). The first line is a comment, `# Table X of <document>: <first
sentence of the caption>`; `scripts/check_tables.py` loads every file (`pandas.read_csv(path,
comment="#")`) and checks row counts and spot values.

Conventions:

- Values are as printed in the article. Nothing has been refitted, rounded or corrected; numbers
  keep the printed digits (e.g. `1.030`, `0.2560`) with an ASCII minus sign, and equation and
  section numbers are those of the article and supplement as printed (e.g. `Eq. (19)`, `Sec. IX E`).
- LaTeX has been converted to plain text or Unicode: `σ_y` (yield stress), `γ̇` (shear rate),
  `γ̇_c` (critical shear rate), `Pa s^n`, `s^-1`, `≤`, `≥`, `∝`, `≈`, `→`, en dash `–` for ranges.
  Subscripts and superscripts use `_` and `^`; a multi-token script is parenthesized, as in
  `(σ_w-σ_y)^(1+1/n)`. Multi-item cells are `; `-separated.
- Citations are written as printed in the article, `Author (year)` or `(Author, year)`, with the
  same `et al.` and year-letter conventions (e.g. `Coussot et al., 2002b`); a `bibkeys` column lists
  the corresponding keys of `references.bib` in the same order. The Table III cell printed as
  "Mahmoudabadbozchelou et al. (2021, 2022)" is expanded to one citation per key.
- Cited values remain the results of the cited works: Tables SII and SIII reproduce coefficients
  published by de Cagny et al. (2019) and Magnon & Cayeux (2021); Table SI, the reanalysis of this
  review (`scripts/figS1_foam_reanalysis.py`); the other tables are qualitative. Cite the original
  sources when reusing cited values.

| File | Table | Rows | Columns |
|---|---|---|---|
| `tableI_materials.csv` | Table I, complex fluids described with the HB model | 12 | material, microstructure, fields_of_use, examples, bibkeys |
| `tableII_field_reviews.csv` | Table II, recent field reviews | 14 | review, field, emphasis_relevant_to_hb, bibkeys |
| `tableIII_milestones.csv` | Table III, milestones | 12 | period, achievement, representative_contributions, sec, bibkeys |
| `tableIV_dimensionless_groups.csv` | Table IV, dimensionless groups | 6 | name, definition, context |
| `tableV_notation.csv` | Table V, principal notation | 11 | symbol, meaning, si_unit, remarks |
| `tableVI_predicted_exponents.csv` | Table VI, predicted exponents and discriminating tests | 8 | theory_or_mechanism, predicted_exponent, discriminating_test, bibkeys |
| `tableVII_yield_stress_protocols.csv` | Table VII, yield-stress protocols | 7 | protocol, observation, threshold_convention, main_ambiguity, bibkeys |
| `tableVIII_verification_problems.csv` | Table VIII, verification problems for HB computations | 9 | problem, reference, quantities_to_compare, principal_test, bibkeys |
| `tableIX_error_budget.csv` | Table IX, error budget | 5 | error_source, diagnostic_comparison, required_report |
| `tableX_models_compared.csv` | Table X, yield-stress models compared with HB | 10 | model, equation_1d_sigma_gt_sigma_y, year, p, relation_to_hb, typical_fields, note (the footnote on Sisko's law) |
| `tableXI_applicability_coordinates.csv` | Table XI, coordinates of HB applicability | 10 | comparison, diagnostic_quantity, measurement_and_consequence |
| `tableXII_transfer_matrix.csv` | Table XII, experimental tests beyond the calibration flow | 21 | group, study, material, calibration_to_test, main_finding, note (the footnote on the Carbopol hair gel), bibkeys |
| `tableXIII_predicted_vs_measured.csv` | Table XIII, what a shear calibration predicts and what else must be measured | 7 | target, what_a_shear_calibration_supports, additional_measurement |
| `tableSI_foam_model_scores.csv` | Table SI, model comparison for the SLES foam local flow curve | 3 | model, k, parameters, varsigma_s_pa, rmse_pa, delta_aic, delta_bic |
| `tableSII_normal_stress_coefficients.csv` | Table SII, HB and normal-stress coefficients of de Cagny et al. (2019) | 3 | material, sigma_y_pa, k_pa_s_n, n, psi1_x1e3_per_pa, psi2_x1e3_per_pa, psi2_over_psi1 |
| `tableSIII_pipe_triplets.csv` | Table SIII, HB triplets of Magnon & Cayeux (2021) | 3 | identification_route, sigma_y_pa, k_pa_s_n, n, delta_sigma_y_over_sigma_y_rot_percent (empty for the rotational reference, printed "—") |
| `tableSIV_onset_exponents.csv` | Table SIV, near-onset exponents | 5 | geometry, leading_response, reason |
| `tableSV_regression_checks.csv` | Table SV, analytical limits used as regression checks | 16 | test, expected_result, error_exposed |
| `tableSVI_unit_conversions.csv` | Table SVI, unit conversions | 7 | quantity, unit, si_equivalent |
| `tableSVII_digital_record.csv` | Table SVII, suggested digital record | 7 | record_component, contents |

Table XII and Fig. 13 describe the same studies: `figure_data/fig13_transfer_matrix.csv` holds the
cells of the matrix drawn in Fig. 13 (group, study, what changes between calibration and test, and
the symbol class), one row per filled cell, whereas `tableXII_transfer_matrix.csv` holds the text of
the table (calibration → test, main finding). In Table XII the study of Geffrault et al. (2021) is
one row with the finding "Mixed"; Fig. 13 draws it as two rows, one for the yield onset and one for
the flow curve, as its caption states.

Table V lists the principal notation only; the complete list of symbols, acronyms and
abbreviations is Supplementary Note S19 of the supplementary material.
