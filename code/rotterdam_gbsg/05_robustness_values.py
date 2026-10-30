"""Estimated population robustness values (Proposition 3) for candidate superiority in
Rotterdam -> GBSG under direct identification and under separate-risk subtraction
(Web Table 15, Notes).  Uses the locked analysis inputs written by run_public_pipeline.py
and the primary plug-in nuisance fit; target outcomes are never read.

The robustness value is the largest gamma on the grid {0, 0.005, ..., 1} at which the
upper endpoint of the identified set for R_T(q1) - R_T(q0) is still negative.  The upper
endpoint is nondecreasing in gamma (Theorem 1), so the scan stops at the first failure.

Usage
  python code/rotterdam_gbsg/05_robustness_values.py [--input-root code/rotterdam_gbsg/analysis_inputs]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]
sys.path.insert(0, str(REPOSITORY / "code" / "onestep"))
import onestep_core as C  # noqa: E402

SEEDS = (42, 52, 62, 72)
GRID = np.round(np.arange(0, 1.0001, 0.005), 3)
TAU = 365.25 * 5


def upper_endpoints(q0, q1, m, g):
    a = q1 * q1 - q0 * q0; b = q1 - q0; lp = C.logit(m); pm, pp = C.expit(lp - g), C.expit(lp + g)
    direct = np.mean(a - 2 * b * np.where(b >= 0, pm, pp))

    def single(q):
        c = 1 - 2 * q
        return np.mean(q * q + c * np.where(c > 0, pm, pp)), np.mean(q * q + c * np.where(c > 0, pp, pm))

    l1, u1 = single(q1); l0, u0 = single(q0)
    return direct, u1 - l0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-root", type=Path, default=HERE / "analysis_inputs")
    ap.add_argument("--out", type=Path, default=REPOSITORY / "results" / "applications" / "robustness_direct_vs_separate.csv")
    args = ap.parse_args()
    if not args.input_root.is_dir():
        raise FileNotFoundError("Run code/rotterdam_gbsg/run_public_pipeline.py first to create analysis_inputs/.")
    rows = []
    for s in SEEDS:
        d = args.input_root / f"seed{s}"
        cal = pd.read_csv(d / "calibration_GBSG.csv"); tg = pd.read_csv(d / "target_GBSG.csv")
        if {"time", "event"} & set(tg.columns):
            raise AssertionError("target outcome firewall violation")
        zs = cal[["q0_surv5y", "q1_surv5y"]].to_numpy(float); zt = tg[["q0_surv5y", "q1_surv5y"]].to_numpy(float)
        m = C.fit_nuisances(cal.time.to_numpy(float), cal.event.to_numpy(int), zs, zt, TAU)["mt"]
        rd = rs = 0.0; d_open = s_open = True
        for g in GRID:
            ud, us = upper_endpoints(zt[:, 0], zt[:, 1], m, g)
            if d_open and ud < 0: rd = float(g)
            else: d_open = False
            if s_open and us < 0: rs = float(g)
            else: s_open = False
        rows.append(dict(seed=s, robustness_direct=rd, robustness_separate=rs, ratio=round(rd / rs, 3)))
    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False, lineterminator="\n")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
