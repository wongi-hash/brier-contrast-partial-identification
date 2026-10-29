"""Scalar constant-shift versus profile-varying sensitivity models (Web Table 1; Web Appendix A.6).

For each application configuration, using the primary IPCW-weighted logistic estimate of m_S:
  * profile-varying (direct) identified-set width at gamma = 0.2;
  * scalar constant-shift range width at gamma = 0.2 over a 401-point eta grid on [-gamma, gamma];
  * width ratio and the first-order diagnostic kappa = |E_T(b w)| / E_T(|b| w), w = m_S(1 - m_S);
  * P_T(q1 > q0), the direction favored at gamma = 0, and population robustness values
    (grid step 0.005 up to 5) under both models.
All quantities use source calibration data and unlabeled target predictions only.

Usage
  python code/scalar_diagnostics/run_scalar_kappa_robustness.py --rotterdam
  python code/scalar_diagnostics/run_scalar_kappa_robustness.py --nsclc-input-root <root>   # restricted inputs
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]
sys.path.insert(0, str(REPOSITORY / "code" / "onestep"))
import onestep_core as C  # noqa: E402

_spec = importlib.util.spec_from_file_location("onestep_apps", REPOSITORY / "code" / "onestep" / "run_onestep_applications.py")
apps = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(apps)

GAMMA = 0.2
RV_GRID = np.round(np.arange(0, 5.0001, 0.005), 3)


def direct_set(q0, q1, m, g):
    a = q1 * q1 - q0 * q0; b = q1 - q0; lp = C.logit(m); pm, pp = C.expit(lp - g), C.expit(lp + g)
    return np.mean(a - 2 * b * np.where(b >= 0, pp, pm)), np.mean(a - 2 * b * np.where(b >= 0, pm, pp))


def scalar_range(q0, q1, m, g, n=401):
    a = q1 * q1 - q0 * q0; b = q1 - q0; lp = C.logit(m)
    values = [np.mean(a - 2 * b * C.expit(lp + eta)) for eta in np.linspace(-g, g, n)]
    return min(values), max(values)


def robustness(fn, favored):
    value = 0.0
    for g in RV_GRID:
        lo, hi = fn(g)
        if (lo > 0) if favored == "Reference" else (hi < 0):
            value = float(g)
        else:
            break
    return value


def row(dataset, seed, cal, target, tau):
    zs = cal[["q0", "q1"]].to_numpy(float); zt = target[["q0", "q1"]].to_numpy(float)
    nu = C.fit_nuisances(cal.time.to_numpy(float), cal.event.to_numpy(int), zs, zt, tau)
    m = nu["mt"]; q0, q1 = zt[:, 0], zt[:, 1]; b = q1 - q0; w = m * (1 - m)
    L, U = direct_set(q0, q1, m, GAMMA); sL, sU = scalar_range(q0, q1, m, GAMMA)
    favored = "Reference" if direct_set(q0, q1, m, 0.0)[0] > 0 else "Candidate"
    return {"dataset": dataset, "seed": seed, "n_target": len(zt), "P_T_q1_gt_q0": float(np.mean(b > 0)),
            "profile_varying_width": U - L, "scalar_width": sU - sL, "width_ratio": (sU - sL) / (U - L),
            "kappa": abs(np.mean(b * w)) / np.mean(np.abs(b) * w), "favored_at_gamma0": favored,
            "robustness_profile_varying": robustness(lambda g: direct_set(q0, q1, m, g), favored),
            "robustness_scalar": robustness(lambda g: scalar_range(q0, q1, m, g, 201), favored)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nsclc-input-root", type=Path)
    ap.add_argument("--nsclc-synthetic", action="store_true")
    ap.add_argument("--rotterdam", action="store_true")
    ap.add_argument("--rotterdam-input-root", type=Path, default=REPOSITORY / "code" / "rotterdam_gbsg" / "analysis_inputs")
    ap.add_argument("--out", type=Path, default=REPOSITORY / "_generated" / "scalar_diagnostics" / "web_table_1.csv")
    a = ap.parse_args()
    jobs = []
    if a.nsclc_synthetic:
        jobs += list(apps.nsclc_jobs(REPOSITORY / "synthetic_fixture" / "nsclc"))
    elif a.nsclc_input_root:
        jobs += list(apps.nsclc_jobs(a.nsclc_input_root))
    if a.rotterdam:
        jobs += list(apps.rotterdam_jobs(a.rotterdam_input_root))
    if not jobs:
        ap.error("select --rotterdam, --nsclc-input-root, or --nsclc-synthetic")
    rows = [row(ds, seed, cal, target, tau) for ds, _, seed, cal, target, tau in jobs]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(a.out, index=False)
    print(pd.DataFrame(rows).round(5).to_string(index=False))


if __name__ == "__main__":
    main()
