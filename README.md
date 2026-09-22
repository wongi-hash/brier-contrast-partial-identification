# Direct Brier-risk contrast partial identification

Reproducibility materials for **Comparing Survival Prediction Models without
Target Outcomes: Sharp Partial Identification of Direct Brier-Risk Contrasts
under Conditional Outcome Shift**.

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
| Decision-regime study | `code/simulations/decision_regime/` |
| NSCLC outcome-free evaluation | `code/nsclc/` |
| Rotterdam–GBSG supporting survival analysis | `code/rotterdam_gbsg/` |
| Figures and tables | `build_figures.py`, `build_tables.py` |

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
python validate_release.py
```

Full application-scale comparator regeneration:

```bash
python code/simulations/application_scale_comparator/run_all.py --full --output _full/comparator
```

Public Rotterdam–GBSG supporting analysis:

```bash
cd code/rotterdam_gbsg
python run_public_pipeline.py --bootstrap 999
```

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
  configurations.

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

The next version will be published only after final validation and explicit
approval.
