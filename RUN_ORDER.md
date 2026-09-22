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

## Portable smoke checks

```bash
python code/simulations/generate_controlled_geometry.py
python code/nsclc/run_outcome_free_evaluation.py --synthetic --bootstrap 9 --out _check/nsclc
python code/simulations/application_scale_comparator/run_all.py --smoke --output _check/comparator
python build_tables.py
python build_figures.py
python validate_release.py
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
