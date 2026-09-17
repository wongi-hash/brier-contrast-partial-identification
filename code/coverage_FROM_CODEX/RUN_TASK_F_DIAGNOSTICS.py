"""Task F: real-data censoring-floor and nuisance-fit diagnostics.

This replays the locked pairs-bootstrap nuisance fits used by the conditional
Cox--Breslow primary implementation.  It never reads target outcomes.  The
reported point-fit floor frequency fills the manuscript placeholder; the
bootstrap diagnostics quantify exceptions and explicit failure to meet the
Newton convergence tolerance.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import platform
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
NSCLC = Path(os.environ.get("TASK_F_NSCLC_ROOT", HERE / "restricted_inputs" / "honest_split_inputs_inductive_ctdate_60_80"))
ROT = Path(os.environ.get("TASK_F_ROTTERDAM_ROOT", HERE / "restricted_inputs" / "rotterdam_gbsg_preregistered"))
OUT = Path(os.environ.get("TASK_F_OUT", HERE / "task_f_diagnostics"))
MODELS = ("InductiveSparsePCA", "InductiveSparsePCA_SRDO")
SITES = ("Stanford", "VA")
SEEDS = (42, 52, 62, 72)
G_MIN = 0.02
BOOT = int(os.environ.get("TASK_F_BOOT", "999"))
WORKERS = int(os.environ.get("TASK_F_WORKERS", "8"))


def stable_seed(*parts) -> int:
    return int(hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()[:8], 16)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def expit(x):
    return 1.0 / (1.0 + np.exp(-np.clip(np.asarray(x, float), -35, 35)))


def cox_components(time, censor_event, z):
    time = np.asarray(time, float)
    censor_event = np.asarray(censor_event, int)
    z = np.asarray(z, float)
    mu = z.mean(0); sd = z.std(0); sd[sd < 1e-8] = 1.0
    x = (z - mu) / sd
    event_times = np.unique(time[censor_event == 1])
    if len(event_times) == 0:
        return mu, sd, np.zeros(x.shape[1]), event_times, np.array([]), True, 0

    def score_info(beta):
        xb = np.clip(x @ beta, -25, 25); risk = np.exp(xb)
        score = np.zeros(x.shape[1]); info = np.zeros((x.shape[1], x.shape[1]))
        for t in event_times:
            ev = (time == t) & (censor_event == 1); d = int(ev.sum()); at = time >= t
            rw = risk[at]; xx = x[at]; s0 = rw.sum()
            mean = (rw[:, None] * xx).sum(0) / s0
            second = np.einsum("i,ij,ik->jk", rw, xx, xx) / s0
            score += x[ev].sum(0) - d * mean
            info += d * (second - np.outer(mean, mean))
        return score, info

    beta = np.zeros(x.shape[1]); converged = False; iterations = 0
    for iterations in range(1, 41):
        score, info = score_info(beta)
        step = np.linalg.pinv(info, rcond=1e-10) @ score
        if not np.all(np.isfinite(step)):
            raise FloatingPointError("nonfinite Cox step")
        mx = np.max(np.abs(step))
        if mx > 2.0: step *= 2.0 / mx
        beta += step
        if np.max(np.abs(step)) < 1e-8:
            converged = True; break
    xb = np.clip(x @ beta, -25, 25); risk = np.exp(xb)
    increments = []
    for t in event_times:
        d = int(((time == t) & (censor_event == 1)).sum())
        increments.append(d / risk[time >= t].sum())
    return mu, sd, beta, event_times, np.asarray(increments), converged, iterations


def conditional_ipcw_diag(time, event, z, tau):
    time = np.asarray(time, float); event = np.asarray(event, int); z = np.asarray(z, float)
    mu, sd, beta, et, dh, cox_ok, cox_iter = cox_components(time, 1 - event, z)
    ph = np.exp(np.clip(((z - mu) / sd) @ beta, -25, 25)); cum = np.cumsum(dh)

    def h0(query, left):
        idx = np.searchsorted(et, np.asarray(query, float), side="left" if left else "right") - 1
        out = np.zeros(len(np.atleast_1d(query))); ok = idx >= 0
        if len(cum): out[ok] = cum[np.minimum(idx[ok], len(cum) - 1)]
        return out

    raw_left = np.exp(-h0(time, True) * ph)
    raw_tau = np.exp(-h0(np.repeat(tau, len(time)), False) * ph)
    gtleft = np.clip(raw_left, G_MIN, 1.0); gtau = np.clip(raw_tau, G_MIN, 1.0)
    y = (time > tau).astype(float); w = np.zeros(len(time))
    survived = time > tau; failed = (time <= tau) & (event == 1); used = survived | failed
    w[survived] = 1.0 / gtau[survived]; w[failed] = 1.0 / gtleft[failed]
    raw_used = np.r_[raw_tau[survived], raw_left[failed]]
    diag = {
        "cox_converged": bool(cox_ok), "cox_iterations": int(cox_iter),
        "n_cal": int(len(time)), "n_ipcw_used": int(used.sum()),
        "floor_tau_n": int((raw_tau < G_MIN).sum()),
        "floor_left_n": int((raw_left < G_MIN).sum()),
        "floor_used_n": int((raw_used < G_MIN).sum()),
        "G_tau_min_raw": float(raw_tau.min()),
        "G_used_min_raw": float(raw_used.min()) if len(raw_used) else math.nan,
    }
    return y, w, diag


def weighted_logistic_diag(x, y, w):
    x = np.asarray(x, float); xd = np.column_stack([np.ones(len(x)), x]); beta = np.zeros(xd.shape[1])
    converged = False; iterations = 0
    for iterations in range(1, 61):
        p = expit(xd @ beta); ww = np.clip(w * p * (1 - p), 1e-10, None)
        info = (xd * ww[:, None]).T @ xd; score = xd.T @ (w * (y - p))
        step = np.linalg.pinv(info, rcond=1e-10) @ score
        if not np.all(np.isfinite(step)):
            raise FloatingPointError("nonfinite logistic step")
        mx = np.max(np.abs(step))
        if mx > 3.0: step *= 3.0 / mx
        beta += step
        if np.max(np.abs(step)) < 1e-8:
            converged = True; break
    return beta, converged, iterations


def fit_once(cal, tau):
    z = cal[["q0", "q1"]].to_numpy(float)
    y, w, gd = conditional_ipcw_diag(cal.time, cal.event, z, tau)
    mu = z.mean(0); sd = z.std(0); sd[sd < 1e-8] = 1.0
    _, log_ok, log_iter = weighted_logistic_diag((z - mu) / sd, y, w)
    gd.update({"logistic_converged": bool(log_ok), "logistic_iterations": int(log_iter)})
    return gd


def one_job(job):
    dataset, contrast, seed, cal_path, target_path, tau = job
    cal = pd.read_csv(cal_path).rename(columns={"q0_surv2y": "q0", "q1_surv2y": "q1",
                                                "q0_surv5y": "q0", "q1_surv5y": "q1"})
    target = pd.read_csv(target_path)
    if {"time", "event"} & set(target.columns):
        raise AssertionError("target-outcome firewall violation")
    point = fit_once(cal, tau)
    rng = np.random.default_rng(stable_seed("harmonize-c4", dataset, contrast, seed, BOOT))
    failures = 0; nonconv = 0; floor_draws = 0; floor_used_n = 0; used_n = 0
    error_types: dict[str, int] = {}
    for _ in range(BOOT):
        ci = rng.integers(0, len(cal), len(cal))
        rng.integers(0, len(target), len(target))  # preserve the locked RNG stream
        try:
            d = fit_once(cal.iloc[ci].reset_index(drop=True), tau)
            if not (d["cox_converged"] and d["logistic_converged"]): nonconv += 1
            if d["floor_used_n"]: floor_draws += 1
            floor_used_n += d["floor_used_n"]; used_n += d["n_ipcw_used"]
        except Exception as exc:
            failures += 1; key = type(exc).__name__; error_types[key] = error_types.get(key, 0) + 1
    return {
        "dataset": dataset, "contrast": contrast, "seed": seed,
        "tau": tau, "n_cal": len(cal), "n_target": len(target),
        "g_min": G_MIN,
        **{f"point_{k}": v for k, v in point.items()},
        "bootstrap_requested": BOOT,
        "bootstrap_fit_failures": failures,
        "bootstrap_fit_failure_fraction": failures / BOOT,
        "bootstrap_nonconverged": nonconv,
        "bootstrap_nonconvergence_fraction": nonconv / BOOT,
        "bootstrap_draws_with_floor": floor_draws,
        "bootstrap_draw_floor_fraction": floor_draws / max(1, BOOT - failures),
        "bootstrap_floor_used_n": floor_used_n,
        "bootstrap_ipcw_used_n": used_n,
        "bootstrap_floor_used_fraction": floor_used_n / used_n if used_n else math.nan,
        "error_types": json.dumps(error_types, sort_keys=True),
        "calibration_sha256": sha256(Path(cal_path)),
        "target_sha256": sha256(Path(target_path)),
    }


def jobs():
    out = []
    for model in MODELS:
        for seed in SEEDS:
            for site in SITES:
                d = NSCLC / model / f"seed{seed}"
                out.append((f"NSCLC-{site}", model, seed,
                            d / f"calibration_{site}.csv", d / f"target_{site}.csv", 730.0))
    for contrast, base in (
        ("ER-only", ROT),
        ("positive-control", ROT / "rotterdam_gbsg_positive_control"),
    ):
        for seed in SEEDS:
            d = base / "analysis_inputs" / f"seed{seed}"
            out.append(("Rotterdam-GBSG", contrast, seed,
                        d / "calibration_GBSG.csv", d / "target_GBSG.csv", 365.25 * 5))
    return out


def main():
    OUT.mkdir(exist_ok=True)
    all_jobs = jobs()
    for j in all_jobs:
        if not Path(j[3]).is_file() or not Path(j[4]).is_file(): raise FileNotFoundError(j)
    rows = []
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futs = {pool.submit(one_job, j): j[:3] for j in all_jobs}
        for fut in as_completed(futs):
            row = fut.result(); rows.append(row)
            print("[DONE]", futs[fut], "fail", row["bootstrap_fit_failures"],
                  "nonconv", row["bootstrap_nonconverged"], flush=True)
    frame = pd.DataFrame(rows).sort_values(["dataset", "contrast", "seed"]).reset_index(drop=True)
    details = OUT / "task_f_realdata_diagnostics.csv"; frame.to_csv(details, index=False)
    summary = frame.groupby(["dataset", "contrast"], as_index=False).agg(
        configurations=("seed", "size"),
        point_floor_used_n=("point_floor_used_n", "sum"),
        point_ipcw_used_n=("point_n_ipcw_used", "sum"),
        bootstrap_requested=("bootstrap_requested", "sum"),
        bootstrap_fit_failures=("bootstrap_fit_failures", "sum"),
        bootstrap_nonconverged=("bootstrap_nonconverged", "sum"),
        bootstrap_draws_with_floor=("bootstrap_draws_with_floor", "sum"),
        bootstrap_floor_used_n=("bootstrap_floor_used_n", "sum"),
        bootstrap_ipcw_used_n=("bootstrap_ipcw_used_n", "sum"),
    )
    summary["point_floor_used_fraction"] = summary.point_floor_used_n / summary.point_ipcw_used_n
    summary["bootstrap_fit_failure_fraction"] = summary.bootstrap_fit_failures / summary.bootstrap_requested
    summary["bootstrap_nonconvergence_fraction"] = summary.bootstrap_nonconverged / summary.bootstrap_requested
    summary["bootstrap_floor_used_fraction"] = summary.bootstrap_floor_used_n / summary.bootstrap_ipcw_used_n
    summary_path = OUT / "task_f_summary.csv"; summary.to_csv(summary_path, index=False)
    manifest = {
        "status": "PASS_TASK_F",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "target_outcomes_read": False,
        "g_min": G_MIN, "bootstrap_per_configuration": BOOT,
        "workers": WORKERS, "runner_sha256": sha256(Path(__file__)),
        "outputs": {p.name: sha256(p) for p in (details, summary_path)},
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
        "definitions": {
            "point_floor_used": "raw G below 0.02 among IPCW-contributing observations",
            "nonconverged": "Cox or weighted-logistic Newton loop did not meet 1e-8 tolerance",
        },
    }
    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(summary.to_string(index=False)); print("[PASS]", OUT)


if __name__ == "__main__":
    main()
