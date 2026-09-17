# Public-package validation report

Validation date: 2026-09-17

## Structural checks

- 11 Python source files parsed successfully with `ast.parse`.
- 15 JSON files parsed successfully.
- 7 CSV files loaded successfully.
- No committed CSV contains patient identifiers or individual time/event rows;
  all CSVs are simulation aggregates or configuration-level derived results.
- No credentials, API tokens, or machine-specific absolute paths were detected.
- Generated `__pycache__` and `*.pyc` files are excluded by `.gitignore` and by
  the public manifest.

## Synthetic execution smoke tests

- `RUN_TASK_E_PRIMARY_COVERAGE.py` imported with the locked defaults
  (`reps=1000`, `bootstrap=199`, `gamma={0,0.1,0.2,0.3,0.5}`).
- A reduced-bootstrap synthetic replication completed and returned both
  two-sided simultaneous endpoint coverage and outer-envelope coverage fields.
- The generic NSCLC downstream inference function completed on synthetic frozen
  predictions after reducing the bootstrap count for the smoke test.

## Locked numerical checks

The committed conditional-censoring coverage summaries reproduce:

| `n_source` | `n_target` | Two-sided endpoint coverage | Outer-envelope coverage |
|---:|---:|---:|---:|
| 114 | 60 | 0.859 | 0.861 |
| 114 | 80 | 0.847 | 0.847 |
| 516 | 686 | 0.943 | 0.946 |

For the application-scale comparator at `n_target=80` and `gamma=0.2`, the
committed aggregate reproduces the following direct-versus-separate ratios and
true-contrast containment values:

| Scenario | Identified-set width ratio | Confidence-envelope width ratio | Direct true-contrast containment |
|---|---:|---:|---:|
| No shift | 3.44 | 1.59 | 0.980 |
| Z-varying in-class shift | 3.44 | 1.58 | 0.982 |
| Low overlap | 3.62 | 1.47 | 0.891 |
| Misspecified source outcome regression | 3.43 | 1.58 | 0.968 |

## Application reproducibility boundary

The repository provides the generic NSCLC downstream sharp-contrast runner and
configuration-level derived results. It does not redistribute institutional
patient data, frozen patient-level prediction files, or upstream radiomic
feature matrices. Consequently, the public package verifies the inference code
and committed summaries but does not claim a fully self-contained reconstruction
of the institutional feature-building pipeline.
