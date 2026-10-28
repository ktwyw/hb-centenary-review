# Tutorial notebooks

Eighteen Jupyter notebooks that work through the calculations behind the review step by step. Each one explains the physics before computing it, uses the notation
of the article, reproduces the corresponding figures, checks its numbers against the plotted data in
[`../data/figure_data/`](../data/figure_data/README.md), and ends with exercises (with worked
solutions), key points and references. The notebooks are stored with their outputs, so they can be
read on GitHub without running anything.

| Notebook | Topic | In the article |
|---|---|---|
| [`01_herschel_bulkley_1926`](01_herschel_bulkley_1926.ipynb) | The 1926 capillary data in modern variables: wall stress, the Weissenberg–Rabinowitsch–Mooney correction, refitting, and why the 1926 relation was a correlation | Sec. II, Fig. 1 |
| [`02_flow_curves_and_dimensionless_groups`](02_flow_curves_and_dimensionless_groups.ipynb) | The HB law and its limits, apparent viscosity, the tensorial form and yield criterion, and the dimensionless groups with a worked pipe example | Sec. IV, Fig. 4, Table IV |
| [`03_fitting_and_identifiability`](03_fitting_and_identifiability.ipynb) | Fitting σ_y, K and n: the residual valley, Mullineux's reduction, weighting, constraints and profile intervals | Sec. VI D, Figs. 7, S4, S5; Note S6 |
| [`04_couette_rheometry`](04_couette_rheometry.ipynb) | Concentric-cylinder rheometry: partly sheared gaps and the bias of the Newtonian conversion | Sec. VI B, Fig. 6 |
| [`05_pipe_and_channel_flow`](05_pipe_and_channel_flow.ipynb) | Laminar pipe and channel flow: profile, plug, flow rate, onset exponents, the convex-duality reading, and a pumping example | Sec. VII A–C, Figs. 8, 9, S6; Note S9 |
| [`06_wall_slip_and_parallel_channels`](06_wall_slip_and_parallel_channels.ipynb) | Wall slip and apparent thresholds; a bundle of capillaries | Sec. VII D, Fig. 10 |
| [`07_regularization_and_the_plug`](07_regularization_and_the_plug.ipynb) | Papanastasiou and bi-viscosity regularizations, the parameter M, and locating the plug | Sec. VII F–G, Figs. 11, 12 |
| [`08_turbulent_pipe_flow`](08_turbulent_pipe_flow.ipynb) | Metzner–Reed Reynolds number, the Fanning friction factor and a pressure-drop calculation | Fig. S7; Notes S9, S12, S15 |
| [`09_foam_local_flow_curve`](09_foam_local_flow_curve.ipynb) | A worked example with published (digitized) data: the local flow curve of a foam and its critical shear rate | Note S2, Fig. S1 |
| [`10_normal_stresses_and_elastoviscoplasticity`](10_normal_stresses_and_elastoviscoplasticity.ipynb) | Normal stresses, the yield invariant and the Saramito model | Sec. VIII C, Notes S3, S4, Fig. S2 |
| [`11_pipe_rheometry_and_calibration_design`](11_pipe_rheometry_and_calibration_design.ipynb) | Pipe rheometry, disagreeing parameter sets that agree in pressure loss, and calibration design | Sec. VI E–F, Note S5, Fig. S3 |
| [`12_literature_parameters_and_tests`](12_literature_parameters_and_tests.ipynb) | Working with the review's literature data: growth of use, field composition, the parameter atlas, and the tests of the model | Secs. III, IX, Figs. 2, 3, 13 |
| [`13_thin_layers_and_free_surface_arrest`](13_thin_layers_and_free_surface_arrest.ipynb) | Thin-layer flux and the pseudo-plug, arrested shapes and the slump-flow yield-stress estimate with its sensitivity, plastocapillary and Bond numbers, and why finite-time stopping does not prove a yield stress | Secs. VII E, IX E–F; Notes S11 B, S15, S16 |
| [`14_thixotropy_and_structural_kinetics`](14_thixotropy_and_structural_kinetics.ipynb) | Simple versus thixotropic yield-stress fluids: the worked structural-kinetics model, rest recovery, rate steps, sweeps and creep bifurcation, why a descending flow curve cannot identify the kinetics, and the foam as the simple case | Secs. V B, VIII B, IX B; Notes S2, S13 |
| [`15_extensional_flows_and_the_yield_criterion`](15_extensional_flows_and_the_yield_criterion.ipynb) | The tensorial HB law in NumPy, uniaxial and planar extension, von Mises versus Tresca over the Lode angle, the measured extensional yield ratios, and normal stresses in the yield invariant | Secs. IV A–B, IV E, IX C; Note S3, Fig. S2 |
| [`16_curve_shape_diagnostics_and_apparent_exponents`](16_curve_shape_diagnostics_and_apparent_exponents.ipynb) | The three-rate diagnostic with its propagated variance, applied to exact, solvent-modified, noisy and 1926 data; solvent, geometry and fitting-window sources of apparent exponents | Secs. III B, V C, VII C, Fig. 9; Note S8, Table SIV |
| [`17_verifying_hb_calculations`](17_verifying_hb_calculations.ipynb) | A verification checklist for HB code: the sixteen analytical identities of Table SV, each with the error it exposes and a deliberately wrong implementation it rejects; the cancellation problem near onset; the error budget and the numerical verification sequence | Secs. VII H–I, Tables VIII, IX, Fig. 12; Notes S9 E, S11, S18 |
| [`18_units_and_worked_engineering_examples`](18_units_and_worked_engineering_examples.ipynb) | Units and conversions (oilfield, cgs, dial readings, the factor 60^n in K) with a unit-safe parameter helper; the slurry-pipeline, slump-flow and settling-cutting examples, the pipeline redone in oilfield units; drilling practice, the fitting window and Table XIII | Sec. X, Table XIII; Notes S15, S19 |

A suggested path for a first reading is 02 → 01 → 05 → 03 → 04 → 07; the others can be read in any order. Notebooks 13–18 go further into topics that the Supplementary Notes develop (thin layers, thixotropy, extensional flows, curve-shape diagnostics, code verification, units and engineering examples).

## Running the notebooks

```sh
python -m pip install -r requirements.txt     # from the repository root
jupyter lab notebooks/
```

or, to re-run all of them from the command line (about 2–3 min in total):

```sh
cd notebooks
jupyter nbconvert --to notebook --execute --inplace *.ipynb
```

The notebooks read data from `../data/` and helper functions from `../scripts/`; they never overwrite the
published figures or data. Calculations that go beyond the article (cross-checks, extra examples) are
labelled as such in the text.
