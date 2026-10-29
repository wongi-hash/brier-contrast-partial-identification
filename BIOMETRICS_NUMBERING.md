# Biometrics manuscript numbering

The v1.1.0 file names follow an earlier manuscript numbering. The table below maps
every figure and table of the Biometrics manuscript, *Sharp Partial Identification
for Survival Model Comparison Without Target Outcomes*, to its committed file and
generator. The `figures/biometrics/` directory holds the exact artwork used in the
Biometrics submission.

Decision labels in the Biometrics manuscript are CANDIDATE SUPPORTED, REFERENCE
SUPPORTED, and DEFER. Committed CSV and JSON files keep their internal codes:

| Manuscript label | Code in `action_*` columns | Code in decision-regime JSON |
|---|---|---|
| CANDIDATE SUPPORTED | `ADOPT` | `retain` |
| REFERENCE SUPPORTED | `KEEP_REFERENCE` | `reference` |
| DEFER | `DEFER` | `defer` |

Tables under `tables/biometrics/` use CAND, REF, and DEFER.

## Main text

| Biometrics item | Committed file | Generator / source |
|---|---|---|
| Figure 1 (workflow) | `figures/biometrics/Figure_1.*` | schematic artwork (no data) |
| Figure 2 (Proposition 2 geometry) | `figures/biometrics/Figure_2.*` (= `figures/Figure_2.*`) | `build_figures.py`, Figure 2 section |
| Figure 3 (NSCLC, seed 42, γ = 0.2) | `figures/biometrics/Figure_3.*` (= `figures/Figure_4.*`) | `build_figures.py`, Figure 4 section; `results/nsclc/outcome_free_evaluation.csv` |
| Figure 4 (Rotterdam–GBSG, seed 62, γ = 0.2) | `figures/biometrics/Figure_4.*` | `build_biometrics_figures.py`; `results/applications/direct_vs_separate.csv` |
| Table 1 (primary-implementation coverage) | `results/simulations/primary_implementation_coverage.csv` | `code/simulations/primary_implementation_coverage/` |
| Table 2 (Rotterdam–GBSG decisions) | `tables/Table_S4_Rotterdam_GBSG_decisions.csv` | `code/rotterdam_gbsg/run_public_pipeline.py` |

## Web Appendix

| Biometrics item | Committed file | Generator / source |
|---|---|---|
| Web Table 1 (scalar vs profile-varying, κ, robustness values) | `tables/biometrics/Web_Table_01_scalar_vs_profile_varying.csv`, `results/scalar_diagnostics/web_table_1.csv` | `code/scalar_diagnostics/run_scalar_kappa_robustness.py` (NSCLC rows: restricted inputs) |
| Web Tables 3–8 (simulations) | `results/simulations/` | `code/simulations/` |
| Web Tables 10–14 (NSCLC details) | `results/nsclc/` | restricted-input analyses; aggregate outputs only |
| Web Table 15 (Rotterdam–GBSG, all cells) | `tables/Table_S4_Rotterdam_GBSG_decisions.csv` | `code/rotterdam_gbsg/run_public_pipeline.py` |
| Web Table 16 (plug-in vs one-step, simulation) | `tables/biometrics/Web_Table_16_onestep_simulation.csv`, `results/onestep/simulation_summary.csv`, `results/onestep/simulation_replicates/`, `results/onestep/population_limit_bias.json` | `code/onestep/run_onestep_simulation.py`, `code/onestep/population_limit_bias.py` |
| Web Table 17 (plug-in vs one-step, applications) | `tables/biometrics/Web_Table_17_onestep_applications.csv`, `results/onestep/applications.csv` | `code/onestep/run_onestep_applications.py` (NSCLC rows: restricted inputs) |
| Web Figure 1 (simulation summary) | `figures/biometrics/Web_Figure_1.*` | same aggregates as `build_figures.py`, Figure 3 section (typography of the submitted artwork differs slightly) |
| Web Figure 2 (decision-regime frequencies) | `figures/biometrics/Web_Figure_2.*` | `build_biometrics_figures.py`; `results/simulations/decision_regime/` |
| Web Appendix A.7 (one-step estimator) | — | `code/onestep/`, `code/onestep/PROTOCOL.md` |
