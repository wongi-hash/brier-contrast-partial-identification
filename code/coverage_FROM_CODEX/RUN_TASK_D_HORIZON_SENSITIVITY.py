"""Task D: prespecified 1-year and 3-year NSCLC horizon sensitivity.

For each horizon, the source-only inductive upstream pipeline rebuilds the
locked survival probabilities while retaining the prespecified model family,
feature-selection procedure, seeds, and hyperparameter search.  The downstream
analysis uses the publication primary conditional Cox--Breslow censoring model,
IPCW logistic m_S, joint pairs bootstrap, and a simultaneous gamma-grid band.
Target outcomes are never read.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import platform
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
WORK = Path(os.environ.get("TASK_D_WORK_DIR", HERE / "restricted_inputs"))
PYTHON = Path(os.environ.get("TASK_D_PYTHON", sys.executable))
UPSTREAM = WORK / "HONEST_SPLIT_INDUCTIVE_UPSTREAM_20260901.py"
MATRIX = WORK / "ctdate_60_80_preprocessing" / "analysis_matrix_primary_outcome_eligible_ctdate_60_80.csv"
C4_RUNNER = Path(os.environ.get("TASK_D_C4_RUNNER", HERE.parent / "RUN_DIRECT_SEPARATE_HARMONIZATION_C4_20260902.py"))
OUT = Path(os.environ.get("TASK_D_OUT", HERE / "task_d_horizon_sensitivity"))
F_MANIFEST = HERE / "task_f_diagnostics" / "MANIFEST.json"
HORIZONS = (365.0, 1095.0)
MODELS = ("InductiveSparsePCA", "InductiveSparsePCA_SRDO")
SEEDS = (42, 52, 62, 72)
SITES = ("Stanford", "VA")
BOOT = int(os.environ.get("TASK_D_BOOT", "999"))
WORKERS = int(os.environ.get("TASK_D_WORKERS", "4"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_c4():
    spec = importlib.util.spec_from_file_location("locked_c4", C4_RUNNER)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def wait_for_task_f():
    if os.environ.get("TASK_D_SKIP_F_WAIT") == "1": return
    while True:
        try:
            m = json.loads(F_MANIFEST.read_text(encoding="utf-8"))
            if m.get("status") == "PASS_TASK_F" and m.get("bootstrap_per_configuration") == 999:
                return
        except (FileNotFoundError, json.JSONDecodeError):
            pass
        time.sleep(20)


def build_inputs(tau: float) -> Path:
    tag = str(int(tau)); out = OUT / f"tau{tag}" / "inputs"
    manifest = out / "MANIFEST.json"
    if manifest.is_file():
        current = json.loads(manifest.read_text(encoding="utf-8"))
        if float(current.get("tau", -1)) == tau and current.get("targets") == {"Stanford": 60, "VA": 80}:
            return out
    out.mkdir(parents=True, exist_ok=True)
    cmd = [str(PYTHON), "-u", str(UPSTREAM), "--matrix", str(MATRIX),
           "--output-root", str(out), "--tau", str(tau)]
    log = OUT / f"tau{tag}" / "upstream.log"
    with log.open("w", encoding="utf-8") as stream:
        code = subprocess.run(cmd, cwd=str(WORK), stdout=stream, stderr=subprocess.STDOUT,
                              creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0)).returncode
    if code != 0: raise RuntimeError(f"upstream failed tau={tau}; see {log}")
    return out


def run_one(payload):
    tau, model, seed, site, root = payload
    mod = load_c4(); d = Path(root) / model / f"seed{seed}"
    cal = pd.read_csv(d / f"calibration_{site}.csv").rename(
        columns={"q0_surv2y": "q0", "q1_surv2y": "q1"})
    target_path = d / f"target_{site}.csv"; target = pd.read_csv(target_path).rename(
        columns={"q0_surv2y": "q0", "q1_surv2y": "q1"})
    if {"time", "event"} & set(target.columns): raise AssertionError("target outcome firewall")
    rows, meta = mod.analyze(f"NSCLC-{site}-tau{int(tau)}", model, seed, cal, target, tau, BOOT)
    gmin = float(meta["G"]["G_tau_min"])
    unstable = bool(gmin <= mod.G_MIN + 1e-12 or meta["bootstrap_valid"] < BOOT or meta["bootstrap_failures"])
    for row in rows:
        row["tau_days"] = tau; row["G_tau_min"] = gmin; row["unstable_cell"] = unstable
    meta.update({"tau_days": tau, "G_tau_min": gmin, "unstable_cell": unstable,
                 "calibration_sha256": sha256(d / f"calibration_{site}.csv"),
                 "target_sha256": sha256(target_path)})
    return rows, meta


def main():
    OUT.mkdir(exist_ok=True); wait_for_task_f()
    roots = {tau: build_inputs(tau) for tau in HORIZONS}

    # Verify that changing tau did not silently retune the locked construction.
    summaries = {tau: pd.read_csv(root / "honest_split_inductive_upstream_summary.csv") for tau, root in roots.items()}
    lock_cols = [c for c in summaries[HORIZONS[0]].columns if c.startswith("lock_") or c in {"model", "seed", "site", "base_penalty"}]
    a = summaries[HORIZONS[0]][lock_cols].sort_values(["model", "seed", "site"]).reset_index(drop=True)
    b = summaries[HORIZONS[1]][lock_cols].sort_values(["model", "seed", "site"]).reset_index(drop=True)
    if not a.equals(b): raise AssertionError("model/hyperparameter locks changed across horizons")

    jobs = [(tau, model, seed, site, str(roots[tau]))
            for tau in HORIZONS for model in MODELS for seed in SEEDS for site in SITES]
    rows, metas = [], []
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futs = {pool.submit(run_one, job): job[:4] for job in jobs}
        for fut in as_completed(futs):
            rr, mm = fut.result(); rows.extend(rr); metas.append(mm)
            print("[DONE]", futs[fut], "unstable", mm["unstable_cell"], flush=True)
    frame = pd.DataFrame(rows).sort_values(["tau_days", "dataset", "contrast", "seed", "gamma"])
    result_path = OUT / "horizon_sensitivity.csv"; frame.to_csv(result_path, index=False)
    primary = frame[(frame.seed == 42) & frame.gamma.isin([0.0, 0.1, 0.2, 0.3, 0.5])]
    primary_path = OUT / "horizon_primary_seed42.csv"; primary.to_csv(primary_path, index=False)
    metadata_path = OUT / "job_metadata.json"
    metadata_path.write_text(json.dumps(metas, indent=2), encoding="utf-8")
    manifest = {
        "status": "PASS_TASK_D",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "target_outcomes_read": False,
        "horizons_days": list(HORIZONS), "bootstrap": BOOT,
        "gamma_grid": load_c4().GAMMA_GRID.tolist(),
        "report_gammas": load_c4().REPORT_GAMMAS.tolist(),
        "censoring_primary": "conditional Cox PH G(t|Z), Breslow; g_min=0.02; refit per pairs-bootstrap draw",
        "model_lock_equal_across_horizons": True,
        "unstable_cells": int(frame[["tau_days", "dataset", "contrast", "seed", "unstable_cell"]].drop_duplicates().unstable_cell.sum()),
        "runner_sha256": sha256(Path(__file__)),
        "upstream_runner_sha256": sha256(UPSTREAM), "c4_runner_sha256": sha256(C4_RUNNER),
        "input_manifests": {str(int(t)): sha256(r / "MANIFEST.json") for t, r in roots.items()},
        "outputs": {p.name: sha256(p) for p in (result_path, primary_path, metadata_path)},
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
    }
    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(primary[["tau_days", "dataset", "contrast", "gamma", "dir_sim_lower", "dir_sim_upper", "action_dir", "G_tau_min", "unstable_cell"]].to_string(index=False))
    print("[PASS]", OUT)


if __name__ == "__main__":
    main()
