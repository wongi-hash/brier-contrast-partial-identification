"""Publication-grade application-matched comparator simulation.

This supersedes RUN_COMPARATOR_APPMATCHED_SIM.py.  The v1 program parsed
FB_WITH_CI but never executed a bootstrap CI, used process-randomized Python
hashes in DGP seeds, had a no-op pS_misspec branch, and changed the true target
law with the analysis gamma.  V2 fixes all four issues.

The target law is fixed within each scenario.  At every Monte Carlo replicate,
the program estimates source p_S from right-censored data, computes four
gamma-indexed identified sets, and (when FB_WITH_CI=1) obtains method-specific
gamma-simultaneous endpoint bands using a joint source/target pairs bootstrap.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPS = int(os.environ.get("FB_REPS", "1000"))
WITH_CI = os.environ.get("FB_WITH_CI", "1") == "1"
BOOT = int(os.environ.get("FB_BOOT", "199" if WITH_CI else "0"))
ONLY_SCENARIO = os.environ.get("FB_SCENARIO", "").strip()
BASE_SEED = int(os.environ.get("FB_SEED", "42"))
DEFER_COST = float(os.environ.get("FB_DEFER_COST", "0.003"))

N_SOURCE = int(os.environ.get("FB_N_SOURCE", "114"))
N_TARGET = int(os.environ.get("FB_N_TARGET", "30"))
TAU = 1.0
CENSOR_SOURCE_DESIGN = float(os.environ.get("FB_CENSOR_SOURCE", "0.69"))
CENSOR_SCALE = float(os.environ.get("FB_CENSOR_SCALE", "0.70"))
GAMMAS = np.asarray([0.1, 0.2, 0.3, 0.5], float)
METHODS = ("transport_assume", "scalar_tilt", "separate_risk", "proposed_direct")
SCENARIOS = (
    "no_shift",
    "constant_shift",
    "xvary_in_class",
    "sign_change",
    "out_of_class",
    "low_overlap",
    "pS_misspec",
)
if ONLY_SCENARIO:
    if ONLY_SCENARIO not in SCENARIOS:
        raise ValueError(ONLY_SCENARIO)
    ACTIVE_SCENARIOS = (ONLY_SCENARIO,)
else:
    ACTIVE_SCENARIOS = SCENARIOS


def expit(z):
    return 1.0 / (1.0 + np.exp(-np.clip(z, -35.0, 35.0)))


def logit(p):
    p = np.clip(np.asarray(p, float), 1e-7, 1 - 1e-7)
    return np.log(p / (1 - p))


def q0f(x):
    return np.clip(expit(0.10 + 0.50 * x[:, 0]), 1e-4, 1 - 1e-4)


def q1f(x):
    return np.clip(expit(0.35 + 0.90 * x[:, 0]), 1e-4, 1 - 1e-4)


def ps_linear_predictor(x):
    q0, q1 = q0f(x), q1f(x)
    return -0.25 + 1.40 * q0 - 0.80 * q1


def ps_true(x, scenario):
    eta = ps_linear_predictor(x)
    if scenario == "pS_misspec":
        # The fitted nuisance remains linear in (q0,q1), while the true law
        # contains an omitted smooth nonlinear term.
        eta = eta + 2.5 * (q1f(x) - 0.55) ** 2
    return np.clip(expit(eta), 1e-5, 1 - 1e-5)


def delta_true(x, scenario):
    z = x[:, 0]
    if scenario == "no_shift":
        return np.zeros(len(x))
    if scenario == "constant_shift":
        return np.full(len(x), 0.20)
    if scenario in {"xvary_in_class", "low_overlap", "pS_misspec"}:
        return 0.20 * np.tanh(z)
    if scenario == "sign_change":
        return 0.40 * np.sin(2.0 * z)
    if scenario == "out_of_class":
        return 0.80 * np.tanh(z)
    raise ValueError(scenario)


def target_shift(scenario):
    return 1.20 if scenario == "low_overlap" else 0.0


def class_contains(method, scenario, gamma):
    gamma = float(gamma)
    if method == "transport_assume":
        return scenario == "no_shift"
    if method == "scalar_tilt":
        if scenario == "no_shift":
            return True
        return scenario == "constant_shift" and gamma >= 0.20 - 1e-12
    if method in {"separate_risk", "proposed_direct"}:
        needed = {
            "no_shift": 0.0,
            "constant_shift": 0.20,
            "xvary_in_class": 0.20,
            "sign_change": 0.40,
            "out_of_class": 0.80,
            "low_overlap": 0.20,
            "pS_misspec": 0.20,
        }[scenario]
        return gamma >= needed - 1e-12
    raise ValueError(method)


def population_truth(scenario):
    # Deterministic Gauss-Hermite integration under X_T ~ N(mu_T,1).
    nodes, weights = np.polynomial.hermite.hermgauss(120)
    x = (target_shift(scenario) + np.sqrt(2.0) * nodes)[:, None]
    w = weights / np.sqrt(np.pi)
    q0, q1 = q0f(x), q1f(x)
    pt = expit(logit(ps_true(x, scenario)) + delta_true(x, scenario))
    a, b = q1**2 - q0**2, q1 - q0
    return float(np.sum(w * (a - 2.0 * b * pt)))


def generate(rep_seed, scenario):
    rng = np.random.default_rng(int(rep_seed))
    xs = rng.normal(0.0, 1.0, (N_SOURCE, 1))
    xt = rng.normal(target_shift(scenario), 1.0, (N_TARGET, 1))
    ps = ps_true(xs, scenario)
    shape = 1.5
    rate = -np.log(ps) / (TAU**shape)
    failure = (-np.log(rng.random(N_SOURCE)) / np.clip(rate, 1e-8, None)) ** (1 / shape)
    # Fixed censoring mechanism independent of failure conditional on X.  The
    # scale was locked before the full run; observed censoring is reported.
    censor = rng.exponential(CENSOR_SCALE, N_SOURCE)
    observed = np.minimum(failure, censor)
    event = (failure <= censor).astype(int)
    return xs, observed, event, xt


def censoring_km(time, censored):
    order = np.argsort(time, kind="mergesort")
    t, c = np.asarray(time)[order], np.asarray(censored)[order]
    unique, first = np.unique(t, return_index=True)
    at_risk = len(t) - first
    deaths = np.add.reduceat(c, first)
    return unique, np.cumprod(1.0 - deaths / at_risk)


def km_left(unique, survival, query):
    query = np.atleast_1d(query)
    index = np.searchsorted(unique, query, side="left") - 1
    out = np.ones(len(query), float)
    good = index >= 0
    out[good] = survival[np.clip(index[good], 0, len(survival) - 1)]
    return out


def fit_ps_model(xs, time, event):
    unique, ghat = censoring_km(time, 1 - event)
    gt = np.clip(km_left(unique, ghat, time), 1e-3, None)
    gtau = max(float(km_left(unique, ghat, [TAU])[0]), 1e-3)
    outcome = (time > TAU).astype(float)
    weight = np.zeros(len(time), float)
    weight[time > TAU] = 1.0 / gtau
    observed_failure = (time <= TAU) & (event == 1)
    weight[observed_failure] = 1.0 / gt[observed_failure]
    feature = np.column_stack([q0f(xs), q1f(xs)])
    center = feature.mean(axis=0)
    scale = feature.std(axis=0, ddof=0)
    scale[scale < 1e-8] = 1.0
    design = np.column_stack([np.ones(len(xs)), (feature - center) / scale])
    beta = np.zeros(design.shape[1], float)
    for _ in range(40):
        fitted = expit(design @ beta)
        work_weight = np.clip(weight * fitted * (1 - fitted), 1e-10, None)
        hessian = (design * work_weight[:, None]).T @ design + 1e-6 * np.eye(design.shape[1])
        step = np.linalg.solve(hessian, design.T @ (weight * (outcome - fitted)))
        beta += step
        if float(np.max(np.abs(step))) < 1e-8:
            break
    return beta, center, scale


def predict_ps(model, x):
    beta, center, scale = model
    feature = np.column_stack([q0f(x), q1f(x)])
    design = np.column_stack([np.ones(len(x)), (feature - center) / scale])
    return np.clip(expit(design @ beta), 1e-4, 1 - 1e-4)


def endpoints_all(q0, q1, ps, gamma):
    gamma = float(gamma)
    a, b = q1**2 - q0**2, q1 - q0
    lp = logit(ps)
    lower_p, upper_p = expit(lp - gamma), expit(lp + gamma)

    p_for_lower = np.where(b >= 0, upper_p, lower_p)
    p_for_upper = np.where(b >= 0, lower_p, upper_p)
    direct = (
        float(np.mean(a - 2 * b * p_for_lower)),
        float(np.mean(a - 2 * b * p_for_upper)),
    )

    def single_risk(q):
        coefficient = 1 - 2 * q
        p_min = np.where(coefficient >= 0, lower_p, upper_p)
        p_max = np.where(coefficient >= 0, upper_p, lower_p)
        return (
            float(np.mean(q**2 + coefficient * p_min)),
            float(np.mean(q**2 + coefficient * p_max)),
        )

    l1, u1 = single_risk(q1)
    l0, u0 = single_risk(q0)
    separate = (l1 - u0, u1 - l0)

    scalar_values = [
        float(np.mean(a - 2 * b * expit(lp + eta)))
        for eta in np.linspace(-gamma, gamma, 101)
    ]
    scalar = (min(scalar_values), max(scalar_values))
    transported = float(np.mean(a - 2 * b * ps))
    return {
        "transport_assume": (transported, transported),
        "scalar_tilt": scalar,
        "separate_risk": separate,
        "proposed_direct": direct,
    }


def decision(lower, upper):
    if upper < 0:
        return "retain"
    if lower > 0:
        return "reference"
    return "defer"


def regret(action, truth):
    if action == "retain":
        return max(float(truth), 0.0)
    if action == "reference":
        return max(-float(truth), 0.0)
    return DEFER_COST


def one_rep(rep, scenario, truth):
    scenario_index = SCENARIOS.index(scenario)
    rep_seed = BASE_SEED + 1000003 * scenario_index + 7919 * int(rep)
    xs, time, event, xt = generate(rep_seed, scenario)
    model = fit_ps_model(xs, time, event)
    ps = predict_ps(model, xt)
    q0, q1 = q0f(xt), q1f(xt)
    points = {method: np.empty((len(GAMMAS), 2), float) for method in METHODS}
    for gi, gamma in enumerate(GAMMAS):
        current = endpoints_all(q0, q1, ps, gamma)
        for method in METHODS:
            points[method][gi] = current[method]

    bands = {}
    if WITH_CI:
        rng = np.random.default_rng(rep_seed + 500000003)
        bootstrap = {method: np.empty((BOOT, len(GAMMAS), 2), float) for method in METHODS}
        for bi in range(BOOT):
            source_index = rng.integers(0, N_SOURCE, N_SOURCE)
            target_index = rng.integers(0, N_TARGET, N_TARGET)
            model_b = fit_ps_model(xs[source_index], time[source_index], event[source_index])
            xt_b = xt[target_index]
            ps_b = predict_ps(model_b, xt_b)
            q0_b, q1_b = q0f(xt_b), q1f(xt_b)
            for gi, gamma in enumerate(GAMMAS):
                current = endpoints_all(q0_b, q1_b, ps_b, gamma)
                for method in METHODS:
                    bootstrap[method][bi, gi] = current[method]
        for method in METHODS:
            se = bootstrap[method].std(axis=0, ddof=1)
            se = np.maximum(se, 1e-10)
            studentized = np.abs((bootstrap[method] - points[method][None, :, :]) / se[None, :, :])
            critical = float(np.quantile(np.max(studentized, axis=(1, 2)), 0.95))
            bands[method] = np.column_stack(
                [points[method][:, 0] - critical * se[:, 0], points[method][:, 1] + critical * se[:, 1]]
            )
    else:
        bands = {method: points[method].copy() for method in METHODS}

    rows = []
    for method in METHODS:
        for gi, gamma in enumerate(GAMMAS):
            point_l, point_u = points[method][gi]
            ci_l, ci_u = bands[method][gi]
            action = decision(ci_l, ci_u)
            rows.append(
                {
                    "scenario": scenario,
                    "rep": int(rep),
                    "method": method,
                    "gamma": float(gamma),
                    "class_contains_true": bool(class_contains(method, scenario, gamma)),
                    "truth": float(truth),
                    "point_lower": float(point_l),
                    "point_upper": float(point_u),
                    "point_coverage": bool(point_l <= truth <= point_u),
                    "point_width": float(point_u - point_l),
                    "ci_lower": float(ci_l),
                    "ci_upper": float(ci_u),
                    "ci_coverage": bool(ci_l <= truth <= ci_u),
                    "ci_width": float(ci_u - ci_l),
                    "decision": action,
                    "regret": regret(action, truth),
                    "observed_censoring": float(1 - event.mean()),
                }
            )
    return rows


def summarize(frame):
    return (
        frame.groupby(["scenario", "method", "gamma", "class_contains_true"], as_index=False)
        .agg(
            reps=("rep", "nunique"),
            truth=("truth", "mean"),
            point_coverage=("point_coverage", "mean"),
            point_width=("point_width", "mean"),
            ci_coverage=("ci_coverage", "mean"),
            ci_width=("ci_width", "mean"),
            retain=("decision", lambda x: float(np.mean(np.asarray(x) == "retain"))),
            reference=("decision", lambda x: float(np.mean(np.asarray(x) == "reference"))),
            defer=("decision", lambda x: float(np.mean(np.asarray(x) == "defer"))),
            mean_regret=("regret", "mean"),
            observed_censoring=("observed_censoring", "mean"),
        )
    )


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    config = {
        "version": "application_matched_comparator_v2",
        "reps": REPS,
        "with_ci": WITH_CI,
        "bootstrap": BOOT,
        "base_seed": BASE_SEED,
        "n_source": N_SOURCE,
        "n_target": N_TARGET,
        "source_censoring_design": CENSOR_SOURCE_DESIGN,
        "censor_scale": CENSOR_SCALE,
        "gammas": GAMMAS.tolist(),
        "scenarios": list(ACTIVE_SCENARIOS),
        "methods": list(METHODS),
        "defer_cost": DEFER_COST,
        "truth": "fixed scenario-specific target population contrast via Gauss-Hermite quadrature",
    }
    code_hash = sha256(Path(__file__))
    config_hash = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()[:12]
    out = HERE / "comparator_out_v2" / f"run_{config_hash}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "CONFIG_EXECUTED.json").write_text(
        json.dumps({**config, "code_sha256": code_hash}, indent=2), encoding="utf-8"
    )
    all_rows = []
    for scenario in ACTIVE_SCENARIOS:
        truth = population_truth(scenario)
        checkpoint = out / f"replicates_{scenario}.csv"
        completed = set()
        existing = None
        if checkpoint.is_file():
            existing = pd.read_csv(checkpoint)
            completed = set(existing.rep.astype(int).unique())
            all_rows.append(existing)
        for rep in range(REPS):
            if rep in completed:
                continue
            rows = one_rep(rep, scenario, truth)
            block = pd.DataFrame(rows)
            block.to_csv(checkpoint, mode="a", index=False, header=not checkpoint.exists())
            all_rows.append(block)
            if (rep + 1) % max(1, REPS // 10) == 0:
                print(f"[{scenario}] {rep + 1}/{REPS}", flush=True)
    frame = pd.concat(all_rows, ignore_index=True)
    frame = frame.drop_duplicates(["scenario", "rep", "method", "gamma"], keep="last")
    summary = summarize(frame)
    replicate_path = out / "all_replicates.csv"
    summary_path = out / "summary.csv"
    frame.to_csv(replicate_path, index=False)
    summary.to_csv(summary_path, index=False)
    manifest = {
        "status": "PASS" if len(frame.scenario.unique()) == len(ACTIVE_SCENARIOS) else "PARTIAL",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config": config,
        "code_sha256": code_hash,
        "python": sys.version,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "platform": platform.platform(),
        "outputs": {
            str(replicate_path): sha256(replicate_path),
            str(summary_path): sha256(summary_path),
        },
    }
    (out / "MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(summary.to_string(index=False))
    print("[PASS]", out)


if __name__ == "__main__":
    main()
