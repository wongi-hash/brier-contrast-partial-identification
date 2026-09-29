# Validation report

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
