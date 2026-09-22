"""Targeted primary-implementation coverage study at NSCLC sample sizes.

Shard with COVERAGE_TARGET_SIZE and COVERAGE_SEED. Aggregate with
COVERAGE_AGGREGATE_ONLY=1.
The DGM has T independent of C conditional on Z=(q0,q1), a correctly specified
IPCW logistic source outcome regression, and an x-varying target shift bounded
by gamma=0.2.  Each bootstrap draw refits G and m_S and resamples target X.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import diagnostics as diag


HERE = Path(__file__).resolve().parent
OUT = HERE / "nsclc_scale_results"
SEEDS = (42, 52, 62, 72)
N_TARGETS = (60, 80)
N_SOURCE = 114
REPS_TOTAL = int(os.environ.get("COVERAGE_REPLICATES", "1000"))
BOOT = int(os.environ.get("COVERAGE_BOOTSTRAP", "199"))
ONLY_SEED = int(os.environ["COVERAGE_SEED"]) if "COVERAGE_SEED" in os.environ else None
ONLY_NT = int(os.environ["COVERAGE_TARGET_SIZE"]) if "COVERAGE_TARGET_SIZE" in os.environ else None
AGG_ONLY = os.environ.get("COVERAGE_AGGREGATE_ONLY") == "1"
GAMMAS = np.array([0.0, 0.1, 0.2, 0.3, 0.5])
TAU = 1.0
WEIBULL_SHAPE = 1.5
TRUE_SHIFT_BOUND = 0.2
TARGET_X_SHIFT = 0.5
ALPHA_C = np.array([1.0, -0.75])
TARGET_CENSORING = 0.65
G_MIN = 0.02


def expit(x): return 1.0 / (1.0 + np.exp(-np.clip(np.asarray(x, float), -35, 35)))
def logit(p):
    p = np.clip(np.asarray(p, float), 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))
def q_values(x):
    x = np.asarray(x, float)
    return expit(-0.15 + 0.45*x), expit(0.10 + 0.80*x)
def p_source_from_q(q0, q1): return expit(-0.35 + 1.25*q0 - 0.70*q1)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""): h.update(block)
    return h.hexdigest()


def calibrate_censoring_baseline():
    rng = np.random.default_rng(20260907)
    x = rng.normal(0, 1, 500_000); q0, q1 = q_values(x); ps = p_source_from_q(q0, q1)
    lam_t = -np.log(np.clip(ps, 1e-8, 1)) / TAU**WEIBULL_SHAPE
    t = (-np.log(rng.random(len(x))) / np.clip(lam_t, 1e-10, None))**(1/WEIBULL_SHAPE)
    z = np.column_stack([q0, q1]); z = (z-z.mean(0))/z.std(0)
    u = -np.log(rng.random(len(x))); mult = np.exp(z @ ALPHA_C)
    lo, hi = 1e-4, 20.0
    for _ in range(60):
        mid = (lo+hi)/2; c = u/(mid*mult); rate = np.mean(c < t)
        if rate < TARGET_CENSORING: lo = mid
        else: hi = mid
    base = (lo+hi)/2
    return float(base), float(np.mean(u/(base*mult) < t))


LAMBDA_C0, CALIBRATED_CENSORING = calibrate_censoring_baseline()


def endpoints(q0, q1, ps, gammas=GAMMAS):
    a = q1*q1-q0*q0; b = q1-q0; lp = logit(ps); out = np.empty((len(gammas), 2))
    for j, g in enumerate(gammas):
        pm, pp = expit(lp-g), expit(lp+g)
        out[j, 0] = np.mean(a-2*b*np.where(b >= 0, pp, pm))
        out[j, 1] = np.mean(a-2*b*np.where(b >= 0, pm, pp))
    return out


def population_truth():
    rng = np.random.default_rng(990011)
    x = rng.normal(TARGET_X_SHIFT, 1, 2_000_000); q0, q1 = q_values(x)
    ps = p_source_from_q(q0, q1)
    bounds = endpoints(q0, q1, ps)
    pt = expit(logit(ps) + TRUE_SHIFT_BOUND*np.tanh(x))
    dr = float(np.mean(q1*q1-q0*q0-2*(q1-q0)*pt))
    return bounds, dr


TRUE_ENDPOINTS, TRUE_DR = population_truth()


def simulate(rep_seed, n_target):
    rng = np.random.default_rng(rep_seed)
    xs = rng.normal(0, 1, N_SOURCE); xt = rng.normal(TARGET_X_SHIFT, 1, n_target)
    q0s, q1s = q_values(xs); q0t, q1t = q_values(xt); ps = p_source_from_q(q0s, q1s)
    lam_t = -np.log(np.clip(ps, 1e-8, 1)) / TAU**WEIBULL_SHAPE
    failure = (-np.log(rng.random(N_SOURCE))/np.clip(lam_t, 1e-10, None))**(1/WEIBULL_SHAPE)
    z = np.column_stack([q0s, q1s]); zstd = (z-z.mean(0))/np.where(z.std(0)<1e-8,1,z.std(0))
    censor = -np.log(rng.random(N_SOURCE))/(LAMBDA_C0*np.exp(zstd@ALPHA_C))
    time = np.minimum(failure, censor); event = (failure <= censor).astype(int)
    cal = pd.DataFrame({"time":time,"event":event,"q0":q0s,"q1":q1s})
    target = np.column_stack([q0t, q1t])
    return cal, target, float(1-event.mean())


def fit_predict(cal, target):
    z = cal[["q0","q1"]].to_numpy(float)
    y, w, gd = diag.conditional_ipcw_diag(cal.time, cal.event, z, TAU)
    mu = z.mean(0); sd = z.std(0); sd[sd<1e-8] = 1.0
    beta, log_ok, log_iter = diag.weighted_logistic_diag((z-mu)/sd, y, w)
    ps = expit(beta[0] + ((target-mu)/sd) @ beta[1:])
    gd.update({"logistic_converged":bool(log_ok),"logistic_iterations":int(log_iter)})
    return np.clip(ps, 1e-4, 1-1e-4), gd


def one_rep(rep_seed, n_target):
    cal, target, censoring = simulate(rep_seed, n_target)
    failures = 0; nonconv = 0; floor_draws = 0
    try:
        ps, point_diag = fit_predict(cal, target)
    except Exception as exc:
        return {"fatal":1,"error":type(exc).__name__,"censoring":censoring}
    point = endpoints(target[:,0], target[:,1], ps)
    draws = np.empty((BOOT, len(GAMMAS), 2)); valid = 0
    rng = np.random.default_rng(rep_seed + 800_000_003)
    for _ in range(BOOT):
        si = rng.integers(0, N_SOURCE, N_SOURCE); ti = rng.integers(0, n_target, n_target)
        try:
            psb, bd = fit_predict(cal.iloc[si].reset_index(drop=True), target)
            if not (bd["cox_converged"] and bd["logistic_converged"]): nonconv += 1
            if bd["floor_used_n"]: floor_draws += 1
            draws[valid] = endpoints(target[ti,0], target[ti,1], psb[ti]); valid += 1
        except Exception:
            failures += 1
    if valid < int(0.9*BOOT):
        return {"fatal":1,"error":"insufficient_bootstrap","censoring":censoring,
                "bootstrap_valid":valid,"bootstrap_failures":failures}
    draws = draws[:valid]; se = np.maximum(draws.std(0, ddof=1), 1e-10)
    stat = np.abs((draws-point[None,:,:])/se[None,:,:]).max(axis=(1,2)); crit = float(np.quantile(stat,.95))
    lower = point[:,0]-crit*se[:,0]; upper = point[:,1]+crit*se[:,1]
    simultaneous = int(np.all((lower <= TRUE_ENDPOINTS[:,0]) & (TRUE_ENDPOINTS[:,0] <= point[:,0]+crit*se[:,0]) &
                              (point[:,1]-crit*se[:,1] <= TRUE_ENDPOINTS[:,1]) & (TRUE_ENDPOINTS[:,1] <= upper)))
    per_gamma = ((lower <= TRUE_ENDPOINTS[:,0]) & (TRUE_ENDPOINTS[:,0] <= point[:,0]+crit*se[:,0]) &
                 (point[:,1]-crit*se[:,1] <= TRUE_ENDPOINTS[:,1]) & (TRUE_ENDPOINTS[:,1] <= upper)).astype(int)
    return {
        "fatal":0,"simultaneous_coverage":simultaneous,
        "outer_coverage":int(np.all((lower <= TRUE_ENDPOINTS[:,0]) & (TRUE_ENDPOINTS[:,1] <= upper))),
        **{f"coverage_gamma_{g:.1f}":int(v) for g,v in zip(GAMMAS,per_gamma)},
        "inclass_simultaneous_coverage":int(np.all(per_gamma[GAMMAS>=TRUE_SHIFT_BOUND-1e-12])),
        "censoring":censoring,"point_cox_converged":int(point_diag["cox_converged"]),
        "point_logistic_converged":int(point_diag["logistic_converged"]),
        "point_floor_used_n":point_diag["floor_used_n"],"point_ipcw_used_n":point_diag["n_ipcw_used"],
        "bootstrap_valid":valid,"bootstrap_failures":failures,"bootstrap_nonconverged":nonconv,
        "bootstrap_floor_draws":floor_draws,"critical_value":crit,
    }


def wilson(k,n,z=1.959963984540054):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return c-h,c+h


def checkpoint_path(nt,seed): return OUT/f"checkpoint_nt{nt}_seed{seed}.json"


def run_shard(nt, seed):
    per = REPS_TOTAL//len(SEEDS); records=[]
    for i in range(per):
        records.append(one_rep(seed*1_000_000+nt*10_000+i,nt))
        if (i+1)%10==0: print(f"[nt={nt} seed={seed}] {i+1}/{per}",flush=True)
    payload={"meta":{"n_target":nt,"seed":seed,"reps":per,"boot":BOOT},"records":records}
    checkpoint_path(nt,seed).write_text(json.dumps(payload),encoding="utf-8")


def aggregate():
    rows=[]
    for nt in N_TARGETS:
        records=[]
        for seed in SEEDS:
            p=checkpoint_path(nt,seed)
            if not p.is_file(): raise FileNotFoundError(p)
            v=json.loads(p.read_text()); meta=v["meta"]
            if meta!={"n_target":nt,"seed":seed,"reps":REPS_TOTAL//4,"boot":BOOT}: raise RuntimeError(f"invalid {p}")
            records.extend(v["records"])
        good=[r for r in records if not r["fatal"]]; n=len(good); k=sum(r["simultaneous_coverage"] for r in good)
        ko=sum(r["outer_coverage"] for r in good); olo,ohi=wilson(ko,n)
        ki=sum(r["inclass_simultaneous_coverage"] for r in good); lo,hi=wilson(k,n); ilo,ihi=wilson(ki,n)
        rows.append({"n_source":N_SOURCE,"n_target":nt,"reps_requested":len(records),"reps_valid":n,
                     "fatal_fit_failures":len(records)-n,
                     "simultaneous_endpoint_coverage":k/n,"coverage_ci_low":lo,"coverage_ci_high":hi,
                     "outer_envelope_coverage":ko/n,"outer_ci_low":olo,"outer_ci_high":ohi,
                     "inclass_gamma_coverage":ki/n,"inclass_ci_low":ilo,"inclass_ci_high":ihi,
                     "mean_censoring":np.mean([r["censoring"] for r in good]),
                     "point_nonconvergence":sum((not r["point_cox_converged"]) or (not r["point_logistic_converged"]) for r in good),
                     "bootstrap_requested":n*BOOT,"bootstrap_fit_failures":sum(r["bootstrap_failures"] for r in good),
                     "bootstrap_nonconverged":sum(r["bootstrap_nonconverged"] for r in good),
                     "bootstrap_floor_draws":sum(r["bootstrap_floor_draws"] for r in good),
                     "mean_critical_value":np.mean([r["critical_value"] for r in good])})
    frame=pd.DataFrame(rows); result=OUT/"coverage_summary.csv"; frame.to_csv(result,index=False)
    manifest={"status":"PASS","study":"targeted primary-implementation stress study","created_utc":datetime.now(timezone.utc).isoformat(),
              "n_source":N_SOURCE,"n_targets":list(N_TARGETS),"reps":REPS_TOTAL,"boot":BOOT,
              "gammas":GAMMAS.tolist(),"true_shift":"0.2*tanh(x)","inclass_gammas":[.2,.3,.5],
              "censoring":{"form":"exponential PH given Z=(q0,q1)","alpha":ALPHA_C.tolist(),
                            "lambda0":LAMBDA_C0,"calibrated_rate":CALIBRATED_CENSORING,"g_min":G_MIN},
              "target_outcomes_used":"simulation truth only; no real target outcomes",
              "runner_sha256":sha256(Path(__file__)),"output_sha256":sha256(result),
              "python":platform.python_version(),"numpy":np.__version__,"pandas":pd.__version__}
    (OUT/"MANIFEST.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    print(frame.to_string(index=False)); print("[PASS]",OUT)


def main():
    OUT.mkdir(exist_ok=True)
    if AGG_ONLY: aggregate(); return
    nts=(ONLY_NT,) if ONLY_NT is not None else N_TARGETS
    seeds=(ONLY_SEED,) if ONLY_SEED is not None else SEEDS
    for nt in nts:
        for seed in seeds: run_shard(nt,seed)
    if ONLY_NT is None and ONLY_SEED is None: aggregate()


if __name__=="__main__": main()
