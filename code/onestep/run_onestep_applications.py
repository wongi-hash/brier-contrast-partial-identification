"""Secondary one-step analysis of the applications (Web Table 17; Web Appendix A.7).

Plug-in and one-step endpoints are computed from the SAME bootstrap draws (B = 999) and the
same critical-value grid {0, 0.01, ..., 1.00} as the reported plug-in analysis; the plug-in
columns therefore reproduce results/nsclc/outcome_free_evaluation.csv and
results/applications/direct_vs_separate.csv.  Target outcomes are never read.

Inputs
  NSCLC (restricted, not distributed): --nsclc-input-root following code/nsclc/input_schema.md
  Rotterdam-GBSG (public): code/rotterdam_gbsg/analysis_inputs/ written by run_public_pipeline.py

Usage
  python code/onestep/run_onestep_applications.py --rotterdam
  python code/onestep/run_onestep_applications.py --nsclc-input-root <root>
  python code/onestep/run_onestep_applications.py --nsclc-synthetic --bootstrap 9   # code-path check
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
sys.path.insert(0, str(HERE))
import onestep_core as C  # noqa: E402

_spec = importlib.util.spec_from_file_location("direct_vs_separate", REPOSITORY / "code" / "inference" / "direct_vs_separate.py")
core = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(core)

FORBIDDEN = {"time", "event", "status", "death", "recur", "rfstime", "dtime", "rtime"}
GRID = core.GAMMA_GRID
REPORT = [int(np.where(np.isclose(GRID, g))[0][0]) for g in core.REPORT_GAMMAS]


def bootstrap_seed(dataset, contrast, seed, boot):
    archived = core.archived_bootstrap_seed(dataset, contrast, seed, boot)
    return archived if archived is not None else core.stable_seed("harmonize-c4", dataset, contrast, seed, boot)


def analyze(dataset, contrast, seed, cal, target, tau, boot):
    if FORBIDDEN & set(target.columns):
        raise RuntimeError("target outcome firewall")
    zs = cal[["q0", "q1"]].to_numpy(float); zt = target[["q0", "q1"]].to_numpy(float)
    t = cal.time.to_numpy(float); e = cal.event.to_numpy(int)
    nu = C.fit_nuisances(t, e, zs, zt, tau); P, O = C.endpoints_both(zs, zt, nu, gammas=GRID)
    rng = np.random.default_rng(bootstrap_seed(dataset, contrast, seed, boot))
    dP = np.empty((boot, len(GRID), 2)); dO = np.empty_like(dP); v = rcap = floor = 0
    for _ in range(boot):
        si = rng.integers(0, len(t), len(t)); ti = rng.integers(0, len(zt), len(zt))
        try:
            nb = C.fit_nuisances(t[si], e[si], zs[si], zt[ti], tau)
            pb, ob = C.endpoints_both(zs[si], zt[ti], nb, gammas=GRID)
            dP[v] = pb; dO[v] = ob; v += 1; rcap += int(nb["rcap"] > 0); floor += int(nb["floor"] > 0)
        except Exception:
            pass
    env = {name: C.sup_t_envelope(p, d[:v]) for name, p, d in (("plugin", P, dP), ("onestep", O, dO))}
    rows = []
    for g, j in zip(core.REPORT_GAMMAS, REPORT):
        row = {"dataset": dataset, "contrast": contrast, "seed": seed, "gamma": float(g), "n_cal": len(t),
               "n_target": len(zt), "bootstrap_valid": v, "density_ratio_cap_point": int(nu["rcap"]),
               "density_ratio_cap_bootstrap_draws": rcap, "censoring_floor_point": int(nu["floor"]),
               "censoring_floor_bootstrap_draws": floor, "density_ratio_max_point": float(nu["r"].max())}
        for name, point in (("plugin", P), ("onestep", O)):
            lo, up, crit, _ = env[name]
            row.update({f"{name}_L": point[j, 0], f"{name}_U": point[j, 1], f"{name}_env_lower": lo[j],
                        f"{name}_env_upper": up[j], f"{name}_critical_value": crit,
                        f"{name}_action": C.action(lo[j], up[j])})
        rows.append(row)
    return rows


def nsclc_jobs(root):
    for cp in sorted(root.glob("pet_radiomic_candidate/seed*/calibration_*.csv")):
        seed = int(cp.parent.name.replace("seed", "")); site = cp.stem.replace("calibration_", "")
        tp = cp.with_name(f"target_{site}.csv")
        cal = pd.read_csv(cp).rename(columns={"q0_surv2y": "q0", "q1_surv2y": "q1"})
        target = pd.read_csv(tp).rename(columns={"q0_surv2y": "q0", "q1_surv2y": "q1"})
        yield f"NSCLC-{site}", "pet_radiomic_candidate", seed, cal, target, 730.0


def rotterdam_jobs(root):
    for seed in core.SEEDS:
        d = root / f"seed{seed}"
        cal = pd.read_csv(d / "calibration_GBSG.csv").rename(columns={"q0_surv5y": "q0", "q1_surv5y": "q1"})
        target = pd.read_csv(d / "target_GBSG.csv").rename(columns={"q0_surv5y": "q0", "q1_surv5y": "q1"})
        yield "Rotterdam-GBSG", "supporting_clinical_candidate", seed, cal, target, 365.25 * 5


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nsclc-input-root", type=Path)
    ap.add_argument("--nsclc-synthetic", action="store_true")
    ap.add_argument("--rotterdam", action="store_true")
    ap.add_argument("--rotterdam-input-root", type=Path, default=REPOSITORY / "code" / "rotterdam_gbsg" / "analysis_inputs")
    ap.add_argument("--bootstrap", type=int, default=999)
    ap.add_argument("--out", type=Path, default=REPOSITORY / "_generated" / "onestep" / "applications.csv")
    a = ap.parse_args()
    jobs = []
    if a.nsclc_synthetic:
        jobs += list(nsclc_jobs(REPOSITORY / "synthetic_fixture" / "nsclc"))
    elif a.nsclc_input_root:
        jobs += list(nsclc_jobs(a.nsclc_input_root))
    if a.rotterdam:
        if not a.rotterdam_input_root.is_dir():
            raise FileNotFoundError("Run code/rotterdam_gbsg/run_public_pipeline.py first to create analysis_inputs/.")
        jobs += list(rotterdam_jobs(a.rotterdam_input_root))
    if not jobs:
        ap.error("select --rotterdam, --nsclc-input-root, or --nsclc-synthetic")
    rows = []
    for job in jobs:
        rows += analyze(*job, a.bootstrap)
        print("done", job[0], job[2], flush=True)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(a.out, index=False)
    print("wrote", a.out, len(rows), "rows")


if __name__ == "__main__":
    main()
