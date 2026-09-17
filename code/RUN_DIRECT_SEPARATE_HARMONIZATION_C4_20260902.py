"""Publication reporting harmonization for direct versus separate Brier bounds.

This script does not refit or tune q0/q1.  It reads the already locked
source-calibration and outcome-free target predictions, estimates p_S using a
conditional Cox censoring model G(t | Z) with Z=(q0,q1), and computes direct
and separate-risk endpoints with the same calibration/target bootstrap draws.

Covered here:
  * MAIN NSCLC inductive full-target analysis (Stanford 59, VA 80)
  * Rotterdam->GBSG ER-only analysis
  * Rotterdam->GBSG post-hoc positive-control analysis

Camelyon17 needs the per-slide locked q0/q1/pS artifact, which was not included
in the recovered compact archive; the script records that dependency rather
than reconstructing it from opened target outcomes.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import io
import json
import platform
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "direct_separate_harmonized_20260902"
NSCLC_ZIP = ROOT / "INDUCTIVE_FULLTARGET_GATE_PACKAGE_20260901.zip"
ROT = ROOT / "rotterdam_gbsg_preregistered"
REPORT_GAMMAS = np.array([0.0, 0.1, 0.2, 0.3, 0.5])
GAMMA_GRID = np.round(np.arange(0.0, 1.0001, 0.01), 2)
MODELS = ("InductiveSparsePCA", "InductiveSparsePCA_SRDO")
SEEDS = (42, 52, 62, 72)
SITES = ("Stanford", "VA")
G_MIN = 0.02


def expit(x):
    return 1.0 / (1.0 + np.exp(-np.clip(np.asarray(x, float), -35, 35)))


def logit(p):
    p = np.clip(np.asarray(p, float), 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


def stable_seed(*parts):
    return int(hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()[:8], 16)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _cox_components(time, censor_event, z):
    """Unpenalized Cox PH fit for censoring and Breslow baseline increments."""
    time = np.asarray(time, float)
    censor_event = np.asarray(censor_event, int)
    z = np.asarray(z, float)
    mu = z.mean(0)
    sd = z.std(0)
    sd[sd < 1e-8] = 1.0
    x = (z - mu) / sd
    event_times = np.unique(time[censor_event == 1])
    if len(event_times) == 0:
        return mu, sd, np.zeros(x.shape[1]), np.array([], float), np.array([], float)

    def score_info(beta):
        xb = np.clip(x @ beta, -25, 25)
        risk = np.exp(xb)
        score = np.zeros(x.shape[1])
        info = np.zeros((x.shape[1], x.shape[1]))
        for t in event_times:
            ev = (time == t) & (censor_event == 1)
            d = int(ev.sum())
            at = time >= t
            rw = risk[at]
            xx = x[at]
            s0 = rw.sum()
            mean = (rw[:, None] * xx).sum(0) / s0
            second = np.einsum("i,ij,ik->jk", rw, xx, xx) / s0
            score += x[ev].sum(0) - d * mean
            info += d * (second - np.outer(mean, mean))
        return score, info

    beta = np.zeros(x.shape[1])
    for _ in range(40):
        score, info = score_info(beta)
        step = np.linalg.pinv(info, rcond=1e-10) @ score
        if not np.all(np.isfinite(step)):
            raise FloatingPointError("nonfinite Cox step")
        # A bounded Newton step avoids rare resample separation explosions.
        mx = np.max(np.abs(step))
        if mx > 2.0:
            step *= 2.0 / mx
        beta += step
        if np.max(np.abs(step)) < 1e-8:
            break

    xb = np.clip(x @ beta, -25, 25)
    risk = np.exp(xb)
    increments = []
    for t in event_times:
        d = int(((time == t) & (censor_event == 1)).sum())
        increments.append(d / risk[time >= t].sum())
    return mu, sd, beta, event_times, np.asarray(increments)


def conditional_ipcw(time, event, z, tau):
    """Fixed-tau IPCW under T independent C conditional on H=Z=(q0,q1)."""
    time = np.asarray(time, float)
    event = np.asarray(event, int)
    z = np.asarray(z, float)
    mu, sd, beta, et, dh = _cox_components(time, 1 - event, z)
    ph = np.exp(np.clip(((z - mu) / sd) @ beta, -25, 25))
    cum = np.cumsum(dh)

    def h0(query, left):
        side = "left" if left else "right"
        idx = np.searchsorted(et, np.asarray(query, float), side=side) - 1
        out = np.zeros(len(np.atleast_1d(query)))
        ok = idx >= 0
        out[ok] = cum[np.minimum(idx[ok], len(cum) - 1)]
        return out

    gtleft = np.clip(np.exp(-h0(time, True) * ph), G_MIN, 1.0)
    gtau = np.clip(np.exp(-h0(np.repeat(tau, len(time)), False) * ph), G_MIN, 1.0)
    y = (time > tau).astype(float)
    w = np.zeros(len(time))
    w[time > tau] = 1.0 / gtau[time > tau]
    fail = (time <= tau) & (event == 1)
    w[fail] = 1.0 / gtleft[fail]
    diag = {
        "censoring_events": int((1 - event).sum()),
        "G_tau_min": float(gtau.min()),
        "G_tau_median": float(np.median(gtau)),
        "H": "Z=(q0,q1); no additional W_C was present in the locked package",
    }
    return y, w, diag


def weighted_logistic(x, y, w):
    x = np.asarray(x, float)
    xd = np.column_stack([np.ones(len(x)), x])
    beta = np.zeros(xd.shape[1])
    for _ in range(60):
        p = expit(xd @ beta)
        ww = np.clip(w * p * (1 - p), 1e-10, None)
        info = (xd * ww[:, None]).T @ xd
        score = xd.T @ (w * (y - p))
        step = np.linalg.pinv(info, rcond=1e-10) @ score
        mx = np.max(np.abs(step))
        if mx > 3.0:
            step *= 3.0 / mx
        beta += step
        if np.max(np.abs(step)) < 1e-8:
            break
    return beta


def fit_predict_ps(cal, target, tau):
    z = cal[["q0", "q1"]].to_numpy(float)
    zt = target[["q0", "q1"]].to_numpy(float)
    y, w, diag = conditional_ipcw(cal.time, cal.event, z, tau)
    mu, sd = z.mean(0), z.std(0)
    sd[sd < 1e-8] = 1.0
    beta = weighted_logistic((z - mu) / sd, y, w)
    ps = expit(beta[0] + ((zt - mu) / sd) @ beta[1:])
    return np.clip(ps, 1e-4, 1 - 1e-4), diag


def all_endpoints(q0, q1, ps, grid=GAMMA_GRID):
    """Direct and separate endpoints, sharing exactly the same p_S and draws."""
    q0, q1, ps = map(lambda v: np.asarray(v, float), (q0, q1, ps))
    a = q1**2 - q0**2
    b = q1 - q0
    lp = logit(ps)
    out = {k: np.empty((len(grid), 2)) for k in ("direct", "separate")}
    for j, gamma in enumerate(grid):
        pm, pp = expit(lp - gamma), expit(lp + gamma)
        out["direct"][j, 0] = np.mean(a - 2 * b * np.where(b >= 0, pp, pm))
        out["direct"][j, 1] = np.mean(a - 2 * b * np.where(b >= 0, pm, pp))

        def single(q):
            coef = 1 - 2 * q
            lo = np.mean(q**2 + coef * np.where(coef > 0, pm, pp))
            hi = np.mean(q**2 + coef * np.where(coef > 0, pp, pm))
            return lo, hi

        l1, u1 = single(q1)
        l0, u0 = single(q0)
        out["separate"][j] = (l1 - u0, u1 - l0)
    if np.any(out["direct"][:, 0] < out["separate"][:, 0] - 1e-10) or np.any(
        out["direct"][:, 1] > out["separate"][:, 1] + 1e-10
    ):
        raise AssertionError("direct set must be contained in separate set")
    return out


def action(lower, upper):
    if upper < 0:
        return "ADOPT"
    if lower > 0:
        return "KEEP_REFERENCE"
    return "DEFER"


def analyze(dataset, contrast, seed, cal, target, tau, boot):
    ps, gdiag = fit_predict_ps(cal, target, tau)
    q0 = target.q0.to_numpy(float)
    q1 = target.q1.to_numpy(float)
    points = all_endpoints(q0, q1, ps)
    rng = np.random.default_rng(stable_seed("harmonize-c4", dataset, contrast, seed, boot))
    draws = {m: np.empty((boot, len(GAMMA_GRID), 2)) for m in points}
    valid = 0
    failures = []
    for bi in range(boot):
        ci = rng.integers(0, len(cal), len(cal))
        ti = rng.integers(0, len(target), len(target))
        try:
            calb = cal.iloc[ci].reset_index(drop=True)
            psb, _ = fit_predict_ps(calb, target, tau)
            ep = all_endpoints(q0[ti], q1[ti], psb[ti])
            for m in draws:
                draws[m][valid] = ep[m]
            valid += 1
        except Exception as exc:
            failures.append(f"{type(exc).__name__}:{str(exc)[:120]}")
        if (bi + 1) % 200 == 0:
            print(f"[{dataset}/{contrast}/seed{seed}] {bi+1}/{boot} valid={valid}", flush=True)
    minimum_valid = max(2, int(0.9 * boot)) if boot < 200 else max(200, int(0.9 * boot))
    if valid < minimum_valid:
        raise RuntimeError(f"insufficient valid bootstrap draws: {valid}/{boot}")

    rows = []
    bands = {}
    for method in draws:
        d = draws[method][:valid]
        se = np.maximum(d.std(0, ddof=1), 1e-10)
        t = np.abs((d - points[method][None, :, :]) / se[None, :, :])
        crit = float(np.quantile(t.max(axis=(1, 2)), 0.95))
        point_crit = float(norm.ppf(0.975))
        bands[method] = {
            "se": se,
            "crit": crit,
            "point_lower": points[method][:, 0] - point_crit * se[:, 0],
            "point_upper": points[method][:, 1] + point_crit * se[:, 1],
            "sim_lower": points[method][:, 0] - crit * se[:, 0],
            "sim_upper": points[method][:, 1] + crit * se[:, 1],
        }

    for gamma in REPORT_GAMMAS:
        j = int(np.where(np.isclose(GAMMA_GRID, gamma))[0][0])
        ld, ud = points["direct"][j]
        ls, us = points["separate"][j]
        bd, bs = bands["direct"], bands["separate"]
        wd, ws = ud - ld, us - ls
        bwd = bd["sim_upper"][j] - bd["sim_lower"][j]
        bws = bs["sim_upper"][j] - bs["sim_lower"][j]
        rows.append({
            "dataset": dataset,
            "contrast": contrast,
            "seed": int(seed),
            "gamma": float(gamma),
            "n_cal": int(len(cal)),
            "n_target": int(len(target)),
            "censoring_model": "unpenalized Cox PH for censoring, H=Z=(q0,q1), Breslow baseline",
            "L_dir": ld, "U_dir": ud, "L_sep": ls, "U_sep": us,
            "W_dir": wd, "W_sep": ws,
            "ID_width_ratio": ws / wd if wd > 1e-12 else np.inf,
            "dir_pointwise_lower": bd["point_lower"][j],
            "dir_pointwise_upper": bd["point_upper"][j],
            "sep_pointwise_lower": bs["point_lower"][j],
            "sep_pointwise_upper": bs["point_upper"][j],
            "dir_sim_lower": bd["sim_lower"][j],
            "dir_sim_upper": bd["sim_upper"][j],
            "sep_sim_lower": bs["sim_lower"][j],
            "sep_sim_upper": bs["sim_upper"][j],
            "BandWidth_dir": bwd, "BandWidth_sep": bws,
            "band_ratio": bws / bwd if bwd > 1e-12 else np.inf,
            "action_dir": action(bd["sim_lower"][j], bd["sim_upper"][j]),
            "action_sep": action(bs["sim_lower"][j], bs["sim_upper"][j]),
            "direct_sim_critical": bd["crit"],
            "separate_sim_critical": bs["crit"],
            "bootstrap_requested": int(boot),
            "bootstrap_valid": int(valid),
        })
    meta = {"dataset": dataset, "contrast": contrast, "seed": seed, "G": gdiag,
            "bootstrap_valid": valid, "bootstrap_failures": failures[:20]}
    return rows, meta


def run_job(payload):
    dataset, contrast, seed, cal, target, tau, boot = payload
    return analyze(dataset, contrast, seed, cal, target, tau, boot)


def read_nsclc():
    if not NSCLC_ZIP.exists():
        raise FileNotFoundError(NSCLC_ZIP)
    with zipfile.ZipFile(NSCLC_ZIP) as zf:
        for model in MODELS:
            for seed in SEEDS:
                for site in SITES:
                    base = f"honest_split_inputs_inductive/{model}/seed{seed}"
                    cal = pd.read_csv(io.BytesIO(zf.read(f"{base}/calibration_{site}.csv"))).rename(
                        columns={"q0_surv2y": "q0", "q1_surv2y": "q1"})
                    target = pd.read_csv(io.BytesIO(zf.read(f"{base}/target_{site}.csv"))).rename(
                        columns={"q0_surv2y": "q0", "q1_surv2y": "q1"})
                    if {"time", "event"} & set(target.columns):
                        raise AssertionError("NSCLC target outcome firewall violation")
                    yield f"NSCLC-{site}", model, seed, cal, target, 730.0


def read_rotterdam(which):
    base = ROT if which == "ER-only" else ROT / "rotterdam_gbsg_positive_control"
    for seed in SEEDS:
        folder = base / "analysis_inputs" / f"seed{seed}"
        cal = pd.read_csv(folder / "calibration_GBSG.csv").rename(
            columns={"q0_surv5y": "q0", "q1_surv5y": "q1"})
        target = pd.read_csv(folder / "target_GBSG.csv").rename(
            columns={"q0_surv5y": "q0", "q1_surv5y": "q1"})
        if {"time", "event"} & set(target.columns):
            raise AssertionError("Rotterdam target outcome firewall violation")
        yield "Rotterdam-GBSG", which, seed, cal, target, 365.25 * 5


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=999)
    ap.add_argument("--only", choices=["all", "nsclc", "rotterdam"], default="all")
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    rows, metas = [], []
    jobs = []
    if args.only in ("all", "nsclc"):
        jobs.extend(read_nsclc())
    if args.only in ("all", "rotterdam"):
        jobs.extend(read_rotterdam("ER-only"))
        jobs.extend(read_rotterdam("positive-control"))
    payloads = []
    for dataset, contrast, seed, cal, target, tau in jobs:
        print(f"[QUEUE] {dataset} {contrast} seed={seed} n={len(cal)}/{len(target)}", flush=True)
        payloads.append((dataset, contrast, seed, cal, target, tau, args.bootstrap))
    if args.workers <= 1:
        results = map(run_job, payloads)
        for rr, mm in results:
            rows.extend(rr); metas.append(mm)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(run_job, p): p[:3] for p in payloads}
            for fut in as_completed(futs):
                rr, mm = fut.result()
                rows.extend(rr); metas.append(mm)
                print(f"[DONE] {futs[fut]}", flush=True)
    frame = pd.DataFrame(rows)
    out_csv = OUT / "harmonized_direct_separate.csv"
    frame.to_csv(out_csv, index=False)

    # Main Table 2 candidates: seed 42, gamma .2, plus positive-control transition.
    table = frame[(frame.seed == 42) & np.isclose(frame.gamma, 0.2)].copy()
    extra = frame[(frame.dataset == "Rotterdam-GBSG") &
                  (frame.contrast == "positive-control") &
                  (frame.seed == 42) & frame.gamma.isin([0.0, 0.1, 0.3])]
    table = pd.concat([table, extra], ignore_index=True).sort_values(["dataset", "contrast", "gamma"])
    table.to_csv(OUT / "TABLE2_candidates.csv", index=False)

    missing_camelyon = {
        "status": "BLOCKED_MISSING_LOCKED_PER_SLIDE_INPUTS",
        "required": [
            "target_scores_FOR_FINAL_EVALUATION_NO_LABELS.csv (slide/patient IDs plus q0/q1)",
            "source calibration outcome-reference predictions or per-slide p_S",
            "cluster identifier used for inference",
        ],
        "reason": "The recovered compact archive contains only aggregate endpoints/metrics and manifests. "
                  "Separate endpoints cannot be recovered from aggregate direct endpoints without q0/q1.",
        "prohibited_workaround": "Do not reconstruct q0/q1 from opened target outcomes.",
    }
    (OUT / "CAMELYON_HARMONIZATION_BLOCKER.json").write_text(
        json.dumps(missing_camelyon, indent=2), encoding="utf-8")

    outputs = [out_csv, OUT / "TABLE2_candidates.csv", OUT / "CAMELYON_HARMONIZATION_BLOCKER.json"]
    manifest = {
        "status": "PASS_NSCLC_ROTTERDAM_C4_DIRECT_SEPARATE__CAMELYON_INPUT_BLOCKED",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "bootstrap": args.bootstrap,
        "gamma_grid": GAMMA_GRID.tolist(),
        "report_gammas": REPORT_GAMMAS.tolist(),
        "target_outcomes_read": False,
        "censoring_primary": "Cox PH G(t|Z), Z=(q0,q1), unpenalized; Breslow baseline; refit per pairs-bootstrap draw",
        "same_draw_contract": True,
        "method_specific_sup_t": True,
        "runner_sha256": sha256(__file__),
        "input_sha256": {
            "NSCLC_zip": sha256(NSCLC_ZIP),
            **{
                str(p.relative_to(ROOT)): sha256(p)
                for p in sorted((ROT / "analysis_inputs").glob("seed*/*.csv"))
            },
            **{
                str(p.relative_to(ROOT)): sha256(p)
                for p in sorted(
                    (ROT / "rotterdam_gbsg_positive_control" / "analysis_inputs").glob("seed*/*.csv")
                )
            },
        },
        "outputs": {p.name: sha256(p) for p in outputs},
        "job_metadata": metas,
        "python": platform.python_version(),
        "numpy": np.__version__,
    }
    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[PASS] rows={len(frame)} -> {out_csv}")
    print(table[["dataset", "contrast", "gamma", "ID_width_ratio", "band_ratio", "action_dir", "action_sep"]].to_string(index=False))


if __name__ == "__main__":
    main()
