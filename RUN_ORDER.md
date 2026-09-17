# Run order — script → manuscript table / result

The conditional-censoring coverage package is in `code/coverage_FROM_CODEX/`.
Its public paths have been normalized; locked configurations and results are
unchanged.

## Tables / numerical results
| Manuscript object | Result file | Producing code |
|---|---|---|
| Table 1 (comparator widths + true-contrast containment) | `final_sim/COMPARATOR_6080_AGGREGATED.csv` | `final_sim/RUN_COMPARATOR_APPMATCHED_SIM_V2.py` → `final_sim/AGGREGATE_COMPARATOR_6080.py` |
| Table S6 (m_S-misspecification containment) | `final_sim/COMPARATOR_6080_AGGREGATED.csv` | comparator study (same as above) |
| Table S4 (Rotterdam→GBSG decisions, seed × γ) | `direct_separate_harmonized_20260902/harmonized_direct_separate.csv` | `code/RUN_DIRECT_SEPARATE_HARMONIZATION_C4_20260902.py` |
| §4.6 / Table 3 (NSCLC direct vs separate) | `ctdate_60_80_results/direct_separate_59_80_vs_60_80.csv` | derived application result; generic downstream runner: `code/RUN_NSCLC_HONEST_SPLIT_SHARP_CONTRAST.py` |
| Table S5 (decision-regime frequencies) | `pending_BDS7_20260907_results/B/AGG_*.json` | `code/coverage_FROM_CODEX/RUN_SHARP_CONTRAST_EXTENSION_SIM.py` |
| Table S9 (γ-grid two-sided **and** outer-envelope coverage: 0.859/0.847/0.943 and 0.861/0.847/0.946) | `code/coverage_FROM_CODEX/FINAL_RESULTS_3_ROWS.csv` | `RUN_TASK_E_PRIMARY_COVERAGE.py` + `RUN_TASK_E2_ROTTERDAM_SCALE_COVERAGE.py`; parallel drivers in the same directory |
| Table S10 (horizon sensitivity) | `code/coverage_FROM_CODEX/horizon_sensitivity.csv` | `code/coverage_FROM_CODEX/RUN_TASK_D_HORIZON_SENSITIVITY.py` |

## Coverage tables (once Codex files are in place)
```bash
cd code/coverage_FROM_CODEX
python run_task_e_parallel.py                    # E: n_T = 60, 80
python run_task_e2_parallel.py                   # E2: n_T = 686
```
Expected: two-sided = 0.859 / 0.847 / 0.943; outer-envelope = 0.861 / 0.847 / 0.946.

## Notes
- Seeds fixed per runner (42/52/62/72); B = 199 bootstrap; 1000 MC replications per
  setting (see each runner's config / `CONFIG_EXECUTED.json`).
- Target outcomes are used only to define the synthetic simulation truth.
- The generic NSCLC downstream runner is public, but the exact clinical/radiomic
  inputs are restricted and are not redistributed. The committed NSCLC tables
  are derived, configuration-level results rather than patient records.
- Horizon and diagnostic scripts require restricted inputs; their committed
  summaries can be inspected without those inputs.
- Figures are not part of this package; they are regenerable from the result files above.
