# Validation report

## Version 1.2.0 additions (checked 2026-10-28)

- **Exact reproduction of archived bootstrap draws.** With the archived seeds in
  `code/inference/archived_bootstrap_seeds.json`, the unchanged public runners
  reproduced all 80 rows of `results/nsclc/outcome_free_evaluation.csv` (restricted
  inputs) and all 20 rows of `tables/Table_S4_Rotterdam_GBSG_decisions.csv`.
  The maximum absolute difference was below 3e-16, and every decision agreed.
  The Rotterdam–GBSG check used analysis inputs identical to those of the archived
  analysis. R was not available on the checking host, so the R export stage was not
  rerun.
- **One-step simulation.** Re-executing committed replications reproduced the stored
  records exactly. `results/onestep/simulation_summary.csv` is regenerated
  bit-for-bit from `results/onestep/simulation_replicates/`.
- **One-step applications.** `results/onestep/applications.csv` reproduces the
  committed plug-in envelopes for all 60 primary cells, with a maximum absolute
  difference below 2e-16. One-step and plug-in decisions agree in all 60 cells.
- **Web Table 1.** `results/scalar_diagnostics/web_table_1.csv` matches the
  manuscript Web Table 1, and |width ratio − κ| ≤ 0.002 in all 12 settings.
- **Biometrics artwork.** `build_biometrics_figures.py` regenerates Figure 4 and
  Web Figure 2 pixel-identically to `figures/biometrics/`.
- **Smoke checks.** The following code paths completed:
  - `run_onestep_simulation.py smoke`;
  - `run_onestep_applications.py --nsclc-synthetic`;
  - `run_scalar_kappa_robustness.py --nsclc-synthetic`;
  - `run_outcome_free_evaluation.py --synthetic`;
  - `build_biometrics_tables.py`;
  - `build_biometrics_figures.py`.

## Version 1.1.0

The release candidate is checked at four levels.

## Numerical invariants

- Primary-implementation simultaneous endpoint coverage is 0.859, 0.847, and
  0.943 at target sample sizes 60, 80, and 686.
- Corresponding outer-envelope coverage is 0.861, 0.847, and 0.946.
- Application-scale direct/separate identified-set width ratios reproduce the
  manuscript values.
- The controlled-geometry illustration reproduces ratios 1.00, 1.22, 1.61,
  2.46, and 5.73.
- NSCLC outcome-free results contain 40 primary and 40 SRDO sensitivity
  configurations, all with decision DEFER.
- In the Rotterdam–GBSG supporting analysis, seed 62 at sensitivity parameter
  0.2 gives ADOPT CANDIDATE for direct identification and DEFER for
  separate-risk subtraction.

## Executable paths

- The NSCLC runner completes with the public synthetic fixture and refuses a
  target file containing outcome columns.
- The application-scale comparator smoke run completes through generation and
  aggregation.
- The public Rotterdam–GBSG pipeline reconstructs the analysis from the R
  `survival` data.
- Figure and table builders regenerate all cited files from committed
  aggregates.
- Figure regeneration writes to `build/figures/` by default and does not
  overwrite the committed reference artwork.

## Release integrity

- `validate_release.py` verifies every entry in `SHA256SUMS.txt` before it
  evaluates numerical or privacy checks.
- A negative-control edit to `README.md` produced the expected checksum
  mismatch and stopped validation; restoring the file returned the release
  audit to PASS.
- `build/` is excluded from the release ledger, so environment-dependent font
  and graphics metadata in newly rendered figures cannot silently change or
  validate the committed reference figures.

## External runtime

The Rotterdam–GBSG pipeline requires R 4.x and the `survival` package. Version
3.8-3 was used for the archived analysis and is recorded in the top-level
README. R was not installed on the Windows host used for the 2026-09-29
checksum/path audit; that audit therefore relied on the previously validated
public-data pipeline and committed aggregates rather than rerunning R.

## Public-data boundary

No institutional patient-level file, patient identifier, date, individual
prediction, raw image, segmentation, credential, or secret is included. Real
NSCLC inputs remain restricted. Public scripts use repository-relative paths;
no local user directory or institutional filesystem path is retained.

## Terminology audit

Public paths, documentation, code labels, aggregate outputs, and table/figure
builders use the terminology in the manuscript. Development task labels,
dated working-folder names, and superseded cohort-comparison suffixes are
excluded by the fail-closed release validator.

The exact automated checks are implemented in `validate_release.py`.
