# Coverage factorial extension: protocol fixed 2026-10-29, before any extension run

Purpose: describe how simultaneous outer-envelope coverage of the primary conditional
Cox–Breslow procedure changes with source calibration information, filling the gap
between the original n_S = 114 and n_S = 516 cells (Web Table 8).

- Code: `run_coverage_factorial.py`; the replicate function, data-generating mechanism,
  estimator, joint pairs bootstrap, and sup-t critical value are unchanged from the
  original six cells.
- Added cells: n_S in {200, 300, 400} x source censoring in {0.65, 0.50, 0.38};
  n_T = 80; gamma grid {0, 0.1, 0.2, 0.3, 0.5}.
- 1,000 replications per cell, B = 199; replication seed
  `stable_seed("factorial", n_S, round(100 * censoring), rep)`, as in the original cells.
- Reported for all 15 cells: outer coverage with Wilson 95% interval, mean realized
  censoring, mean source events n_S x (1 - mean censoring), mean critical value.
- All cells are reported regardless of result. No single sample-size threshold is claimed.
