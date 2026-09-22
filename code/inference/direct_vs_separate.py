"""Shared direct-versus-separate Brier-risk contrast computations.

The module estimates the source outcome regression using a conditional Cox
censoring model with Breslow baseline and computes direct and separate-risk
endpoints from the same calibration and outcome-free target bootstrap draws.
It contains no application-specific data paths.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm


REPORT_GAMMAS = np.array([0.0, 0.1, 0.2, 0.3, 0.5])
GAMMA_GRID = np.round(np.arange(0.0, 1.0001, 0.01), 2)
SEEDS = (42, 52, 62, 72)
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

    raw_left = np.exp(-h0(time, True) * ph)
    raw_tau = np.exp(-h0(np.repeat(tau, len(time)), False) * ph)
    gtleft = np.clip(raw_left, G_MIN, 1.0)
    gtau = np.clip(raw_tau, G_MIN, 1.0)
    y = (time > tau).astype(float)
    w = np.zeros(len(time))
    w[time > tau] = 1.0 / gtau[time > tau]
    fail = (time <= tau) & (event == 1)
    w[fail] = 1.0 / gtleft[fail]
    diag = {
        "censoring_events": int((1 - event).sum()),
        "G_tau_min": float(gtau.min()),
        "G_tau_median": float(np.median(gtau)),
        "G_floor": G_MIN,
        "floor_tau_n": int((raw_tau < G_MIN).sum()),
        "floor_left_n": int((raw_left < G_MIN).sum()),
        "floor_used_n": int((raw_tau[time > tau] < G_MIN).sum() + (raw_left[fail] < G_MIN).sum()),
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
