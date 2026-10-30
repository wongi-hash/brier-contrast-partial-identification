# Changelog

## v1.2.3 — 2026-10-30

- added the finite-sample coverage diagnosis of Web Appendix C.4:
  `code/simulations/coverage_factorial/diagnose_coverage.py`, replicate
  records for the factorial cells n_S = 114 with 65% censoring and
  n_S = 516 with 50% censoring under
  `results/simulations/coverage_factorial/diagnosis/`, and
  `diagnosis_summary.json`. The diagnosis is exploratory, was added after the
  factorial results were available, reuses the factorial seeds, and reproduces
  the stored factorial envelopes and critical values exactly;
- `CITATION.cff` and the archive metadata record the exact v1.2.3 DOI
  (`10.5281/zenodo.23794202`); no previously committed result changed.
- There is no v1.2.2 release.

## v1.2.1 — 2026-10-30

- added the source-size by censoring coverage factorial (Web Table 8):
  `code/simulations/coverage_factorial/run_coverage_factorial.py`, its
  protocol, and replicate records for all 15 cells under
  `results/simulations/coverage_factorial/`; the original n_S = 114 and 516
  cells are now archived, and n_S = 200, 300, 400 were added under the
  protocol in `code/simulations/coverage_factorial/PROTOCOL.md`; Web Figure 3
  (`plot_web_figure_3.py`, `figures/biometrics/Web_Figure_3.*`) plots all 15 cells;
- added `code/rotterdam_gbsg/05_robustness_values.py` and
  `results/applications/robustness_direct_vs_separate.csv` (robustness values
  under direct identification and separate subtraction; Web Table 15 Notes);
- `CITATION.cff` and the archive metadata now record the exact v1.2.1 DOI
  (`10.5281/zenodo.23789365`); no previously committed result changed.

## v1.2.0 — 2026-10-29

- added the one-step endpoint estimator (Web Appendix A.7): core routines,
  prespecified protocol, simulation (Web Table 16; 4 conditions × 1,000
  replications, committed replicate records), population-limit bias, and the
  secondary application analysis (Web Table 17) computed from the archived
  bootstrap draws;
- added the scalar constant-shift versus profile-varying diagnostic κ and
  robustness values (Web Table 1);
- added `BIOMETRICS_NUMBERING.md`, `figures/biometrics/`, `tables/biometrics/`,
  `build_biometrics_figures.py`, and `build_biometrics_tables.py` for the
  Biometrics manuscript numbering and decision labels (CANDIDATE SUPPORTED,
  REFERENCE SUPPORTED, DEFER); committed CSV/JSON decision codes are unchanged;
- fixed exact reproduction of the committed NSCLC and Rotterdam–GBSG bootstrap
  envelopes: `code/inference/direct_vs_separate.py` now uses the archived
  integer seeds in `code/inference/archived_bootstrap_seeds.json` (B = 999);
  v1.1.0 derived seeds from the renamed public labels, which reproduced the
  identified sets and decisions but not every envelope digit;
- updated `validate_release.py` with checks for the new results.

## v1.1.0 — 2026-09-29

- aligned every public path and analysis label with manuscript terminology;
- removed local filesystem paths and execution-machine snapshots;
- added executable controlled-geometry, application-scale comparator,
  primary-implementation coverage, and Rotterdam–GBSG pipelines;
- added the NSCLC restricted-input schema, synthetic fixture, and explicit
  outcome-column firewall;
- added manuscript figure and table builders;
- added fail-closed privacy, credential, path, terminology, and numerical
  validation.
- aligned the controlled-geometry generator with its documented result path;
- changed ordinary figure rebuilds to `build/figures/` so reference artwork is
  not overwritten by environment-dependent font or metadata differences;
- made `validate_release.py` verify the complete SHA-256 ledger;
- documented the R 4.x / `survival` 3.8-3 environment and `Rscript` entry point;
- clarified that supplementary S7, S8, S11, and S14 are restricted-input
  NSCLC analyses auditable from aggregate summaries, not fully rerunnable from
  the public release;
- prepared repository metadata for a new GitHub/Zenodo release. The
  version-specific Zenodo DOI is inserted only after minting.

## Version 1.0.0

- initial reproducibility archive.
