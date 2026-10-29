# Changelog

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
