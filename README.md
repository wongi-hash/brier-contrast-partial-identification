# Sharp Partial Identification for Survival Model Comparison Without Target Outcomes

Reproducibility materials for **Sharp Partial Identification for Survival Model
Comparison Without Target Outcomes** (submitted to *Biometrics*; earlier working
title: *Comparing Survival Prediction Models without Target Outcomes: Sharp
Partial Identification of Direct Brier-Risk Contrasts under Conditional Outcome
Shift*). Versions 1.2.0 and 1.2.1 add the analyses introduced for the Biometrics
manuscript; `BIOMETRICS_NUMBERING.md` maps every Biometrics figure and table to
its committed file and generator.

The repository follows the terminology and analysis structure used in the
manuscript. It contains statistical code, public-data pipelines, synthetic
fixtures, aggregate results, manuscript tables, and Figures 1–4 and S1–S2. It
does not contain institutional patient-level data, patient identifiers, dates,
individual predictions, or medical images.

## Repository map

| Manuscript component | Public materials |
|---|---|
| Controlled-geometry illustration | `code/simulations/generate_controlled_geometry.py` |
| NSCLC application-scale comparator and mirror study | `code/simulations/application_scale_comparator/` |
| Targeted primary-implementation stress study | `code/simulations/primary_implementation_coverage/` |
| Source-size by censoring coverage factorial (Web Table 8) | `code/simulations/coverage_factorial/` |
| Decision-regime study | `code/simulations/decision_regime/` |
| NSCLC outcome-free evaluation | `code/nsclc/` |
| Rotterdam–GBSG supporting survival analysis | `code/rotterdam_gbsg/` (robustness values: `05_robustness_values.py`) |
| Figures and tables | `build_figures.py`, `build_tables.py` |
| One-step endpoint estimator (Web Appendix A.7; Web Tables 16–17) | `code/onestep/` |
| Scalar versus profile-varying diagnostic κ (Web Table 1) | `code/scalar_diagnostics/` |
| Biometrics figures and Web Tables 1, 16, 17 | `build_biometrics_figures.py`, `build_biometrics_tables.py`, `figures/biometrics/`, `tables/biometrics/` |

## NSCLC data boundary

The manuscript analysis compares a clinical reference predictor with a
source-trained PET radiomic candidate. The statistical analysis code is public
and can be exercised with the included synthetic fixture. Real institutional
inputs are restricted and are not redistributed; therefore the NSCLC numerical
results are auditable from committed aggregate outputs but cannot be fully
rerun from this public repository alone. The required outcome-free input schema
and explicit restricted-input error are documented in
`code/nsclc/input_schema.md`.

## Quick validation

```bash
python -m pip install -r requirements.txt
python code/simulations/generate_controlled_geometry.py
python code/nsclc/run_outcome_free_evaluation.py --synthetic --bootstrap 9 --out _check/nsclc
python code/simulations/application_scale_comparator/run_all.py --smoke --output _check/comparator
python build_tables.py
python build_figures.py
python code/onestep/run_onestep_simulation.py smoke
python code/onestep/run_onestep_applications.py --nsclc-synthetic --bootstrap 9 --out _check/onestep/applications.csv
python build_biometrics_tables.py
python build_biometrics_figures.py
python validate_release.py
```

`build_figures.py` writes regenerated artwork to `build/figures/` by default.
It does not overwrite the committed reference figures. To intentionally stage
release artwork in another directory, set `SIM_REPRO_FIGDIR` explicitly. The
release validator verifies every byte listed in `SHA256SUMS.txt`; it does not
accept file presence alone as an integrity check.

Full application-scale comparator regeneration:

```bash
python code/simulations/application_scale_comparator/run_all.py --full --output _full/comparator
```

Public Rotterdam–GBSG supporting analysis:

```bash
cd code/rotterdam_gbsg
python run_public_pipeline.py --bootstrap 999
```

**R requirements:** R 4.x and the `survival` package are required. Version
`3.8-3` of `survival` was used for the archived analysis. The pipeline invokes
the export stage as:

```bash
Rscript 00_export_survival_data.R
```

The remaining stages are orchestrated by `run_public_pipeline.py`. A future
`survival` version should not be substituted for release verification without
recording the version and checking the locked outputs.

Add `--open-outcomes` only for the explicitly retrospective outcome comparison.
The supporting analysis is not presented as prospectively blinded clinical
validation.

## Numerical checks

- primary-implementation simultaneous endpoint coverage: `0.859 / 0.847 /
  0.943` for target sample sizes 60, 80, and 686;
- corresponding outer-envelope coverage: `0.861 / 0.847 / 0.946`;
- application-scale direct/separate identified-set width ratios at target size
  80 and sensitivity parameter 0.2: `3.44 / 3.44 / 3.62 / 3.43`;
- controlled-geometry width ratios: `1.00 / 1.22 / 1.61 / 2.46 / 5.73`;
- NSCLC: DEFER across all prespecified primary and SRDO sensitivity
  configurations;
- one-step versus plug-in (Web Table 17): identical decisions in all 60
  application cells; Rotterdam–GBSG seed 62, γ = 0.2 one-step envelope
  `[-0.0486, -0.0010]`;
- one-step simulation (Web Table 16): outer coverage `0.863 / 0.844` (law A,
  NSCLC scale, plug-in / one-step), `0.957 / 0.957` (law B, Rotterdam–GBSG
  scale).

## Exact reproduction of archived bootstrap draws (v1.2.0)

In v1.1.0 the public runners derived bootstrap seeds from the renamed public
labels, whereas the committed NSCLC and Rotterdam–GBSG results had been
computed with the seeds of the archived analysis. Rerunning v1.1.0 therefore
reproduced the identified sets and the decisions (verified for all 20 Rotterdam–GBSG cells) but not every envelope digit
(for example, Rotterdam–GBSG seed 62, γ = 0.2 gave `[-0.0493, -0.0010]` instead
of the committed `[-0.0486, -0.0017]`). Version 1.2.0 stores the archived
integer seeds in `code/inference/archived_bootstrap_seeds.json`; with B = 999
the public runners now reproduce `results/nsclc/outcome_free_evaluation.csv`
and `tables/Table_S4_Rotterdam_GBSG_decisions.csv` to floating-point precision.
For other values of B the label-derived seeds are used as before.

## Restricted-input supplementary analyses

Supplementary items S7, S8, S11, and S14 use the restricted NSCLC inputs.
Their committed, non-individual-level summaries can be audited, but these
items are not fully rerunnable from the public repository alone. This is the
same data boundary described above; it is not a claim that the institutional
inputs are publicly available. Where a summary is not present in `results/` or
`tables/`, the corresponding item must not be described as publicly
regenerable.

## Data and software terms

The public TCIA NSCLC-Radiogenomics and R `survival` datasets remain governed
by their original terms. The MIT license applies to code in this repository,
not to external or restricted data.

See `RUN_ORDER.md` for code-to-result traceability and
`VALIDATION_REPORT.md` for validation details.

## Archive

- Repository: https://github.com/wongi-hash/brier-contrast-partial-identification
- Zenodo concept DOI: https://doi.org/10.5281/zenodo.22807124
- Version 1.0.0 DOI: https://doi.org/10.5281/zenodo.22807125
- Version 1.1.0 DOI: https://doi.org/10.5281/zenodo.23028694
- Version 1.2.0 DOI: https://doi.org/10.5281/zenodo.23750710
- Version 1.2.1 DOI: https://doi.org/10.5281/zenodo.23789365
