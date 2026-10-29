# One-step endpoint estimator (Web Appendix A.7; Web Tables 16–17)

The protocol is fixed in `PROTOCOL.md`. `onestep_core.py` implements the plug-in
and one-step endpoint estimators. Their nuisance fits reproduce the primary
implementation in `code/inference/direct_vs_separate.py`, and the plug-in
endpoints agree to floating-point precision.

## Simulation (fully public)

```bash
python code/onestep/run_onestep_simulation.py smoke                  # quick code-path check
python code/onestep/run_onestep_simulation.py run --law A --scale nsclc
python code/onestep/run_onestep_simulation.py run --law B --scale nsclc
python code/onestep/run_onestep_simulation.py run --law A --scale rott
python code/onestep/run_onestep_simulation.py run --law B --scale rott
python code/onestep/run_onestep_simulation.py aggregate              # -> results/onestep/simulation_summary.csv
python code/onestep/population_limit_bias.py                         # -> results/onestep/population_limit_bias.json
```

Runs are resumable. Add `--time-limit SECONDS` to split a long run into several
sessions.

The committed replicate files are in `results/onestep/simulation_replicates/`.
Rerunning them reproduces every record exactly.

## Applications (secondary analysis)

```bash
# Rotterdam–GBSG (public data): first run code/rotterdam_gbsg/run_public_pipeline.py
python code/onestep/run_onestep_applications.py --rotterdam
# NSCLC (restricted inputs following code/nsclc/input_schema.md)
python code/onestep/run_onestep_applications.py --nsclc-input-root <root>
# public code-path check with the synthetic fixture
python code/onestep/run_onestep_applications.py --nsclc-synthetic --bootstrap 9
```

The applications reuse the archived bootstrap seeds
(`code/inference/archived_bootstrap_seeds.json`) and the reported critical-value
grid. The plug-in columns of `results/onestep/applications.csv` therefore
reproduce `results/nsclc/outcome_free_evaluation.csv` and
`results/applications/direct_vs_separate.csv` exactly.

## Scalar diagnostics (Web Table 1)

```bash
python code/scalar_diagnostics/run_scalar_kappa_robustness.py --rotterdam
python code/scalar_diagnostics/run_scalar_kappa_robustness.py --nsclc-input-root <root>
```
