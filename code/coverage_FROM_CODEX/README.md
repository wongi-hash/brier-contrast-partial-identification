# Conditional-censoring coverage and supporting analyses

This directory contains the executed simulation package used for the manuscript's
conditional Cox--Breslow coverage assessment and related supplementary analyses.
No patient-level data are included.

## Locked coverage results

The committed summaries were produced with 1,000 Monte Carlo replications and
199 joint pairs-bootstrap draws per replication on the prespecified grid
`gamma = {0, 0.1, 0.2, 0.3, 0.5}`.

| Setting | Two-sided simultaneous endpoint coverage | Outer-envelope coverage |
|---|---:|---:|
| `n_source=114`, `n_target=60` | 0.859 | 0.861 |
| `n_source=114`, `n_target=80` | 0.847 | 0.847 |
| `n_source=516`, `n_target=686` | 0.943 | 0.946 |

These values are stored in `FINAL_RESULTS_3_ROWS.csv`. The first two rows also
appear in `task_e_primary_coverage_summary.csv`; the third appears in
`task_e2_rotterdam_scale_coverage_summary.csv`.

## Fully synthetic reruns

The coverage studies use simulated source and target outcomes only and can be
rerun without restricted data:

```bash
python run_task_e_parallel.py
python run_task_e2_parallel.py
```

The parallel launchers use the active Python interpreter. Direct single-process
runs are also possible with `RUN_TASK_E_PRIMARY_COVERAGE.py` and
`RUN_TASK_E2_ROTTERDAM_SCALE_COVERAGE.py`, but are slower.

The decision-regime simulation is defined by
`RUN_SHARP_CONTRAST_EXTENSION_SIM.py` and `CONFIG_LOCKED.json`.

## Scripts requiring restricted inputs

`RUN_TASK_D_HORIZON_SENSITIVITY.py` and `RUN_TASK_F_DIAGNOSTICS.py` are supplied
for transparency, but their exact application reruns require restricted or
externally obtained inputs that are not redistributed. Their committed derived
summary is `horizon_sensitivity.csv`.

Relevant environment-variable overrides are documented in the scripts. Keep
restricted inputs outside the public repository (the default placeholder is
`restricted_inputs/`, which is ignored by Git).

## Provenance

- `CONFIG_EXECUTED.json`: executed coverage specification.
- `ENVIRONMENT_PROVENANCE.json`: software environment used for the locked run.
- `MANIFEST_TASK_E.json` and `MANIFEST_TASK_E2.json`: study-level provenance.
- `SHA256SUMS.txt`: hashes for the public files in this directory.

The public copies normalize machine-specific absolute paths. This normalization
does not alter algorithms, configurations, or numerical results.
