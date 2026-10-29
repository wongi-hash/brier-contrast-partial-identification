# Reproduction order and traceability

| Manuscript output | Committed result | Generator |
|---|---|---|
| Controlled-geometry illustration | `results/simulations/controlled_geometry.csv` | `code/simulations/generate_controlled_geometry.py` |
| Application-scale comparator | `results/simulations/application_scale_comparator.csv` | `code/simulations/application_scale_comparator/run_all.py` |
| Decision-regime study | `results/simulations/decision_regime/*.json` | `code/simulations/decision_regime/run_decision_regime.py` |
| Primary-implementation coverage | `results/simulations/primary_implementation_coverage.csv` | `code/simulations/primary_implementation_coverage/` |
| NSCLC horizon sensitivity | `results/nsclc/horizon_sensitivity.csv` | `code/nsclc/run_horizon_sensitivity.py` (restricted inputs required) |
| NSCLC outcome-free evaluation | `results/nsclc/outcome_free_evaluation.csv` | `code/nsclc/run_outcome_free_evaluation.py` (restricted inputs required; synthetic fixture public) |
| Rotterdam–GBSG supporting analysis | `results/applications/direct_vs_separate.csv` | `code/rotterdam_gbsg/run_public_pipeline.py` |
| One-step simulation (Web Table 16) | `results/onestep/simulation_summary.csv`, `results/onestep/simulation_replicates/`, `results/onestep/population_limit_bias.json` | `code/onestep/run_onestep_simulation.py`, `code/onestep/population_limit_bias.py` |
| One-step applications (Web Table 17) | `results/onestep/applications.csv` | `code/onestep/run_onestep_applications.py` (NSCLC rows: restricted inputs) |
| Scalar diagnostic κ and robustness values (Web Table 1) | `results/scalar_diagnostics/web_table_1.csv` | `code/scalar_diagnostics/run_scalar_kappa_robustness.py` (NSCLC rows: restricted inputs) |
| Biometrics Web Tables 1, 16, 17 | `tables/biometrics/` | `build_biometrics_tables.py` |
| Biometrics Figure 4 and Web Figure 2 | `figures/biometrics/` | `build_biometrics_figures.py` |

## Portable smoke checks

```bash
python code/simulations/generate_controlled_geometry.py
python code/nsclc/run_outcome_free_evaluation.py --synthetic --bootstrap 9 --out _check/nsclc
python code/simulations/application_scale_comparator/run_all.py --smoke --output _check/comparator
python build_tables.py
python build_figures.py
python validate_release.py
```

Regenerated figures are written to `build/figures/`; committed reference
figures under `figures/` are not overwritten. `validate_release.py` verifies
the complete `SHA256SUMS.txt` ledger before checking numerical claims.

## One-step estimator and scalar diagnostics (v1.2.0)

```bash
python code/onestep/run_onestep_simulation.py run --law A --scale nsclc   # repeat for B and for --scale rott
python code/onestep/run_onestep_simulation.py aggregate
python code/onestep/population_limit_bias.py
python code/onestep/run_onestep_applications.py --rotterdam                # after run_public_pipeline.py
python code/onestep/run_onestep_applications.py --nsclc-input-root <root>  # restricted inputs
python code/scalar_diagnostics/run_scalar_kappa_robustness.py --rotterdam --nsclc-input-root <root>
python build_biometrics_tables.py
python build_biometrics_figures.py
```

## Full application-scale simulation

```bash
python code/simulations/application_scale_comparator/run_all.py --full --output _full/comparator
```

## Rotterdam–GBSG public-data analysis

```bash
cd code/rotterdam_gbsg
python run_public_pipeline.py --bootstrap 999
```

The outcome-free stage rejects any target outcome column. Use
`--open-outcomes` only for the separate retrospective outcome comparison.

## Restricted NSCLC inputs

The institutional analysis requires the schema in `code/nsclc/input_schema.md`.
No institutional patient-level input is included in this public release. The
synthetic fixture verifies the same loader, firewall, censoring model, source
outcome regression, joint pairs bootstrap, and decision code without
reproducing manuscript data.

Supplementary S7, S8, S11, and S14 are restricted-input analyses. They are
auditable only through committed aggregate summaries and are not represented
as fully rerunnable public-data analyses.
