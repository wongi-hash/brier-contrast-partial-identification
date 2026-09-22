"""Locked helpers for the preregistered Rotterdam -> GBSG experiment."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict

TAU = 365.25 * 5
CAP = 365.25 * 7
G_MIN = 0.02
SEEDS = [42, 52, 62, 72]
REPORT_GAMMAS = np.array([0.0, 0.1, 0.2, 0.3, 0.5])
GAMMA_GRID = np.round(np.arange(0.0, 1.0001, 0.01), 2)
Q0_COLS = ["age3", "age3log", "meno", "size_mid", "size_large", "nodes_inv_sqrt", "hormon"]
Q1_COLS = Q0_COLS + ["er"]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump_json(obj, path):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -35, 35)))


def logit(p):
    p = np.clip(np.asarray(p, dtype=float), 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


def source_outcome(raw):
    time = np.where(raw["recur"].to_numpy() == 1, raw["rtime"].to_numpy(), raw["dtime"].to_numpy()).astype(float)
    event = ((raw["recur"].to_numpy() == 1) | (raw["death"].to_numpy() == 1)).astype(int)
    event = np.where(time <= CAP, event, 0)
    time = np.minimum(time, CAP)
    return time, event


def target_outcome(raw):
    time = np.minimum(raw["rfstime"].to_numpy(dtype=float), CAP)
    event = np.where(raw["rfstime"].to_numpy(dtype=float) <= CAP, raw["status"].to_numpy(dtype=int), 0)
    return time, event


def feature_frame(raw, cohort):
    age10 = raw["age"].to_numpy(dtype=float) / 10.0
    if cohort == "rotterdam":
        size = raw["size"].astype(str)
        size_mid = size.eq("20-50").astype(float).to_numpy()
        size_large = size.eq(">50").astype(float).to_numpy()
    elif cohort == "gbsg":
        size = raw["size"].to_numpy(dtype=float)
        size_mid = ((size > 20) & (size <= 50)).astype(float)
        size_large = (size > 50).astype(float)
    else:
        raise ValueError(cohort)
    nodes = raw["nodes"].to_numpy(dtype=float)
    if np.any(nodes <= 0):
        raise ValueError("feature_frame requires node-positive rows")
    return pd.DataFrame({
        "age3": age10 ** 3,
        "age3log": (age10 ** 3) * np.log(age10),
        "meno": raw["meno"].to_numpy(dtype=float),
        "size_mid": size_mid,
        "size_large": size_large,
        "nodes_inv_sqrt": 1.0 / np.sqrt(nodes),
        "hormon": raw["hormon"].to_numpy(dtype=float),
        "er": raw["er"].to_numpy(dtype=float),
    })


def fit_scaler(frame, cols):
    mu = frame[cols].mean(axis=0)
    sd = frame[cols].std(axis=0, ddof=0).replace(0, 1.0)
    return mu, sd


def apply_scaler(frame, cols, mu, sd):
    return (frame[cols] - mu) / sd


def fit_predictor(train_features, train_time, train_event, cols, penalizer=0.01):
    mu, sd = fit_scaler(train_features, cols)
    x = apply_scaler(train_features, cols, mu, sd).copy()
    x["time"] = train_time
    x["event"] = train_event
    model = CoxPHFitter(penalizer=penalizer)
    model.fit(x, duration_col="time", event_col="event", show_progress=False)
    return {"model": model, "cols": list(cols), "mu": mu, "sd": sd}


def predict_survival(lock, features, tau=TAU):
    x = apply_scaler(features, lock["cols"], lock["mu"], lock["sd"])
    pred = lock["model"].predict_survival_function(x, times=[tau]).iloc[0].to_numpy(dtype=float)
    return np.clip(pred, 1e-4, 1 - 1e-4)


def serialize_lock(lock):
    model = lock["model"]
    return {
        "columns": lock["cols"],
        "mean": {k: float(v) for k, v in lock["mu"].items()},
        "sd": {k: float(v) for k, v in lock["sd"].items()},
        "coef": {k: float(v) for k, v in model.params_.items()},
        "penalizer": float(model.penalizer),
        "baseline_survival_tau": float(model.baseline_survival_.reindex(model.baseline_survival_.index.union([TAU])).sort_index().ffill().loc[TAU].iloc[0]),
    }


def _baseline_hazard_at(model, query, left=False):
    base = model.baseline_cumulative_hazard_.iloc[:, 0]
    times = base.index.to_numpy(dtype=float)
    vals = base.to_numpy(dtype=float)
    q = np.asarray(query, dtype=float)
    side = "left" if left else "right"
    idx = np.searchsorted(times, q, side=side) - 1
    out = np.zeros(len(q), dtype=float)
    ok = idx >= 0
    out[ok] = vals[np.minimum(idx[ok], len(vals) - 1)]
    return out


def conditional_ipcw(time, event, z, penalizer=0.05):
    """Returns fixed-tau binary outcome and IPCW under T independent C | Z."""
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    z = np.asarray(z, dtype=float)
    mu = z.mean(axis=0)
    sd = z.std(axis=0)
    sd[sd < 1e-8] = 1.0
    zs = (z - mu) / sd
    frame = pd.DataFrame(zs, columns=["z0", "z1"])
    frame["time"] = time
    frame["censor_event"] = 1 - event
    cox = CoxPHFitter(penalizer=penalizer)
    cox.fit(frame, duration_col="time", event_col="censor_event", show_progress=False)
    ph = cox.predict_partial_hazard(frame[["z0", "z1"]]).to_numpy(dtype=float)
    h_tleft = _baseline_hazard_at(cox, time, left=True)
    h_tau = _baseline_hazard_at(cox, np.repeat(TAU, len(time)), left=False)
    g_tleft = np.clip(np.exp(-h_tleft * ph), G_MIN, 1.0)
    g_tau = np.clip(np.exp(-h_tau * ph), G_MIN, 1.0)
    y = (time > TAU).astype(int)
    w = np.zeros(len(time), dtype=float)
    w[time > TAU] = 1.0 / g_tau[time > TAU]
    failed = (time <= TAU) & (event == 1)
    w[failed] = 1.0 / g_tleft[failed]
    return y, w, {"g_tau_min": float(g_tau.min()), "g_tau_median": float(np.median(g_tau))}


def fit_ps(time, event, cal_q0, cal_q1, target_q0, target_q1):
    z = np.column_stack([cal_q0, cal_q1])
    y, w, gdiag = conditional_ipcw(time, event, z)
    mu = z.mean(axis=0)
    sd = z.std(axis=0)
    sd[sd < 1e-8] = 1.0
    model = LogisticRegression(C=10.0, penalty="l2", solver="lbfgs", max_iter=1000, random_state=0)
    model.fit((z - mu) / sd, y, sample_weight=w)
    zt = np.column_stack([target_q0, target_q1])
    ps = model.predict_proba((zt - mu) / sd)[:, 1]
    return np.clip(ps, 1e-4, 1 - 1e-4), gdiag


def sharp_endpoints(q0, q1, ps, gamma_grid=GAMMA_GRID):
    q0 = np.asarray(q0, dtype=float)
    q1 = np.asarray(q1, dtype=float)
    a = q1 ** 2 - q0 ** 2
    b = q1 - q0
    lp = logit(ps)
    lower = np.empty(len(gamma_grid))
    upper = np.empty(len(gamma_grid))
    for j, gamma in enumerate(gamma_grid):
        pm = sigmoid(lp - gamma)
        pp = sigmoid(lp + gamma)
        p_for_lower = np.where(b >= 0, pp, pm)
        p_for_upper = np.where(b >= 0, pm, pp)
        lower[j] = np.mean(a - 2 * b * p_for_lower)
        upper[j] = np.mean(a - 2 * b * p_for_upper)
    return lower, upper


def stable_seed(*parts):
    text = "|".join(str(x) for x in parts)
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:8], 16)


def gate_one(cal, target, seed, boot=999, progress=True):
    ct = cal["time"].to_numpy(dtype=float)
    ce = cal["event"].to_numpy(dtype=int)
    cq0 = cal["q0_surv5y"].to_numpy(dtype=float)
    cq1 = cal["q1_surv5y"].to_numpy(dtype=float)
    tq0 = target["q0_surv5y"].to_numpy(dtype=float)
    tq1 = target["q1_surv5y"].to_numpy(dtype=float)
    ps, gdiag = fit_ps(ct, ce, cq0, cq1, tq0, tq1)
    lower, upper = sharp_endpoints(tq0, tq1, ps)
    rng = np.random.default_rng(stable_seed("rotterdam_gbsg", seed, boot))
    lb = []
    ub = []
    failures = []
    for b in range(boot):
        ci = rng.integers(0, len(cal), len(cal))
        ti = rng.integers(0, len(target), len(target))
        try:
            ps_b, _ = fit_ps(ct[ci], ce[ci], cq0[ci], cq1[ci], tq0, tq1)
            l_b, u_b = sharp_endpoints(tq0[ti], tq1[ti], ps_b[ti])
            if np.all(np.isfinite(l_b)) and np.all(np.isfinite(u_b)):
                lb.append(l_b)
                ub.append(u_b)
            else:
                failures.append({"rep": b, "reason": "nonfinite"})
        except Exception as exc:
            failures.append({"rep": b, "reason": type(exc).__name__ + ":" + str(exc)[:160]})
        if progress and (b + 1) % 100 == 0:
            print("  seed {} bootstrap {}/{} valid={}".format(seed, b + 1, boot, len(lb)), flush=True)
    # Small B is allowed only for public code-path smoke tests; manuscript runs
    # use B=999 and therefore still require at least max(200, 90%) valid draws.
    minimum_valid = max(2, int(0.9 * boot)) if boot < 200 else max(200, int(0.9 * boot))
    if len(lb) < minimum_valid:
        raise RuntimeError("Too many bootstrap failures: {}/{} valid".format(len(lb), boot))
    lb = np.asarray(lb)
    ub = np.asarray(ub)
    se_l = np.maximum(lb.std(axis=0, ddof=1), 1e-10)
    se_u = np.maximum(ub.std(axis=0, ddof=1), 1e-10)
    max_t = np.maximum(np.abs(lb - lower) / se_l, np.abs(ub - upper) / se_u).max(axis=1)
    critical = float(np.quantile(max_t, 0.95))
    lower_ci = lower - critical * se_l
    upper_ci = upper + critical * se_u
    rows = []
    for gamma in REPORT_GAMMAS:
        j = int(np.where(np.isclose(GAMMA_GRID, gamma))[0][0])
        decision = "ADOPT" if upper_ci[j] < 0 else ("KEEP_REFERENCE" if lower_ci[j] > 0 else "DEFER")
        rows.append({
            "seed": int(seed), "gamma": float(gamma),
            "L_point": float(lower[j]), "U_point": float(upper[j]),
            "L_simCI": float(lower_ci[j]), "U_simCI": float(upper_ci[j]),
            "decision": decision,
        })
    cand = GAMMA_GRID[upper < 0]
    ref = GAMMA_GRID[lower > 0]
    summary = {
        "seed": int(seed), "n_cal": int(len(cal)), "n_target": int(len(target)),
        "bootstrap_requested": int(boot), "bootstrap_valid": int(len(lb)),
        "bootstrap_failures": failures[:20], "simultaneous_critical": critical,
        "gamma_star_candidate_point": float(cand.max()) if len(cand) else 0.0,
        "gamma_star_reference_point": float(ref.max()) if len(ref) else 0.0,
        "mean_abs_q1_minus_q0": float(np.mean(np.abs(tq1 - tq0))),
        "censoring_diagnostics": gdiag,
    }
    dense = pd.DataFrame({
        "seed": seed, "gamma": GAMMA_GRID, "L_point": lower, "U_point": upper,
        "L_simCI": lower_ci, "U_simCI": upper_ci,
    })
    return pd.DataFrame(rows), dense, summary


def overlap_diagnostics(cal, target, seed):
    zs = cal[["q0_surv5y", "q1_surv5y"]].to_numpy(dtype=float)
    zt = target[["q0_surv5y", "q1_surv5y"]].to_numpy(dtype=float)
    z = np.vstack([zs, zt])
    y = np.r_[np.zeros(len(zs), dtype=int), np.ones(len(zt), dtype=int)]
    mu = z.mean(axis=0)
    sd = z.std(axis=0)
    sd[sd < 1e-8] = 1.0
    zz = (z - mu) / sd
    cv = StratifiedKFold(5, shuffle=True, random_state=seed)
    clf = LogisticRegression(C=1.0, solver="lbfgs", max_iter=1000, random_state=seed)
    prop = cross_val_predict(clf, zz, y, cv=cv, method="predict_proba")[:, 1]
    auc = roc_auc_score(y, prop)
    p_source = np.clip(prop[:len(zs)], 1e-5, 1 - 1e-5)
    odds = p_source / (1 - p_source) * (len(zs) / len(zt))
    ess = odds.sum() ** 2 / np.sum(odds ** 2)
    outside = np.mean(
        (zt[:, 0] < zs[:, 0].min()) | (zt[:, 0] > zs[:, 0].max()) |
        (zt[:, 1] < zs[:, 1].min()) | (zt[:, 1] > zs[:, 1].max())
    )
    zss = (zs - zs.mean(axis=0)) / np.where(zs.std(axis=0) < 1e-8, 1, zs.std(axis=0))
    ztt = (zt - zs.mean(axis=0)) / np.where(zs.std(axis=0) < 1e-8, 1, zs.std(axis=0))
    source_nn = np.sqrt(((zss[:, None, :] - zss[None, :, :]) ** 2).sum(axis=2) + np.eye(len(zss)) * 1e12).min(axis=1)
    target_nn = np.sqrt(((ztt[:, None, :] - zss[None, :, :]) ** 2).sum(axis=2)).min(axis=1)
    return {
        "seed": int(seed), "domain_auc_orientation_invariant": float(max(auc, 1 - auc)),
        "target_ess": float(ess), "target_ess_fraction": float(ess / len(zs)),
        "target_outside_source_rectangle_fraction": float(outside),
        "target_knn_above_source_p95_fraction": float(np.mean(target_nn > np.quantile(source_nn, 0.95))),
        "source_n": int(len(zs)), "target_n": int(len(zt)),
    }


def oracle_metrics(time, event, q0, q1):
    _, w, _ = conditional_ipcw(time, event, np.column_stack([q0, q1]))
    y = (np.asarray(time) > TAU).astype(float)
    b0 = float(np.mean(w * (y - q0) ** 2))
    b1 = float(np.mean(w * (y - q1) ** 2))
    c0 = float(concordance_index(time, q0, event))
    c1 = float(concordance_index(time, q1, event))
    return {"brier_q0": b0, "brier_q1": b1, "delta_brier_q1_minus_q0": b1 - b0,
            "cindex_q0": c0, "cindex_q1": c1, "delta_cindex_q1_minus_q0": c1 - c0}
