"""Plug-in versus one-step endpoint estimators in simulation (Web Table 16; Web Appendix C.2).

Two source outcome laws (A: correctly specified; B: omitted term 2.5 (q1 - 0.55)^2 in the
source logit) at two information scales (NSCLC: n_S=114, n_T=80, ~65% censoring;
Rotterdam-GBSG: n_S=516, n_T=686, ~38% censoring).  Plug-in and one-step endpoints are
computed from the same data and the same bootstrap draws (B=199).

Usage
  python code/onestep/run_onestep_simulation.py run --law A --scale nsclc [--reps 1000] [--time-limit SEC]
  python code/onestep/run_onestep_simulation.py aggregate
  python code/onestep/run_onestep_simulation.py smoke          # 3 replications, B=19, temporary output

Runs are resumable: completed replications are skipped.  Replication seeds are fixed by
stable_seed("onestep-sim-20261028", law, scale_key, rep) so committed results are reproduced.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]
sys.path.insert(0, str(HERE))
import onestep_core as C  # noqa: E402

RESULTS = REPOSITORY / "results" / "onestep"
TAU = 1.0
WEIBULL_SHAPE = 1.5
TARGET_X_SHIFT = 0.5
ALPHA_C = np.array([1.0, -0.75])
# scale key (used in seeds) -> (n_S, n_T, target censoring, file stem)
SCALES = {"nsclc": (114, 80, 0.65, "nsclc_scale"), "rott": (516, 686, 0.38, "rotterdam_gbsg_scale")}
GAMMAS = C.GAMMAS


def q_values(x):
    x = np.asarray(x, float)
    return C.expit(-0.15 + 0.45 * x), C.expit(0.10 + 0.80 * x)


def source_law(q0, q1, law):
    eta = -0.35 + 1.25 * q0 - 0.70 * q1
    if law == "B":
        eta = eta + 2.5 * (q1 - 0.55) ** 2
    return C.expit(eta)


def calibrate_censoring(law, target_censoring):
    rng = np.random.default_rng(20260907)
    x = rng.normal(0, 1, 500_000); q0, q1 = q_values(x); ps = source_law(q0, q1, law)
    lam = -np.log(np.clip(ps, 1e-8, 1)) / TAU ** WEIBULL_SHAPE
    t = (-np.log(rng.random(len(x))) / np.clip(lam, 1e-10, None)) ** (1 / WEIBULL_SHAPE)
    z = np.column_stack([q0, q1]); z = (z - z.mean(0)) / z.std(0)
    u = -np.log(rng.random(len(x))); mult = np.exp(z @ ALPHA_C)
    lo, hi = 1e-4, 20.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if np.mean(u / (mid * mult) < t) < target_censoring:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def population_truth(law):
    rng = np.random.default_rng(990011)
    x = rng.normal(TARGET_X_SHIFT, 1, 2_000_000); q0, q1 = q_values(x)
    ps = source_law(q0, q1, law); out = np.empty((len(GAMMAS), 2))
    for j, g in enumerate(GAMMAS):
        for k, endpoint in enumerate(("L", "U")):
            f, _ = C._f_and_df(q0, q1, ps, g, endpoint)
            out[j, k] = f.mean()
    return out


def simulate(seed, law, n_s, n_t, lam_c0):
    rng = np.random.default_rng(seed)
    xs = rng.normal(0, 1, n_s); xt = rng.normal(TARGET_X_SHIFT, 1, n_t)
    q0s, q1s = q_values(xs); q0t, q1t = q_values(xt); ps = source_law(q0s, q1s, law)
    lam = -np.log(np.clip(ps, 1e-8, 1)) / TAU ** WEIBULL_SHAPE
    failure = (-np.log(rng.random(n_s)) / np.clip(lam, 1e-10, None)) ** (1 / WEIBULL_SHAPE)
    z = np.column_stack([q0s, q1s]); zs = (z - z.mean(0)) / np.where(z.std(0) < 1e-8, 1, z.std(0))
    censor = -np.log(rng.random(n_s)) / (lam_c0 * np.exp(zs @ ALPHA_C))
    time_ = np.minimum(failure, censor); event = (failure <= censor).astype(int)
    return time_, event, z, np.column_stack([q0t, q1t])


def one_replication(seed, law, scale, lam_c0, truth, boot):
    n_s, n_t, _, _ = SCALES[scale]
    t, e, zs, zt = simulate(seed, law, n_s, n_t, lam_c0)
    rec = {"seed": seed, "cens": float(1 - e.mean())}
    try:
        nu = C.fit_nuisances(t, e, zs, zt, TAU); P, O = C.endpoints_both(zs, zt, nu)
    except Exception as exc:  # pragma: no cover - recorded, not hidden
        rec.update({"fatal": 1, "err": type(exc).__name__}); return rec
    rb = np.random.default_rng(seed + 800_000_003)
    dP = np.empty((boot, len(GAMMAS), 2)); dO = np.empty_like(dP); v = rcap = floor = 0
    for _ in range(boot):
        si = rb.integers(0, n_s, n_s); ti = rb.integers(0, n_t, n_t)
        try:
            nb = C.fit_nuisances(t[si], e[si], zs[si], zt[ti], TAU)
            pb, ob = C.endpoints_both(zs[si], zt[ti], nb)
            if not (np.all(np.isfinite(pb)) and np.all(np.isfinite(ob))):
                raise FloatingPointError
            dP[v] = pb; dO[v] = ob; v += 1; rcap += int(nb["rcap"] > 0); floor += int(nb["floor"] > 0)
        except Exception:
            pass
    if v < int(0.9 * boot):
        rec.update({"fatal": 1, "err": "boot"}); return rec
    rec.update({"fatal": 0, "boot_valid": v, "rcap_point": int(nu["rcap"]), "rcap_boot": rcap,
                "floor_point": int(nu["floor"]), "floor_boot": floor, "P": P.tolist(), "O": O.tolist()})
    for name, point, draws in (("P", P, dP[:v]), ("O", O, dO[:v])):
        lo, up, crit, _ = C.sup_t_envelope(point, draws)
        rec[name + "_outer"] = int(np.all((lo <= truth[:, 0]) & (truth[:, 1] <= up)))
        rec[name + "_width"] = (up - lo).tolist(); rec[name + "_crit"] = crit
    return rec


def run(law, scale, reps, boot, out_dir, time_limit):
    n_s, n_t, cens, stem = SCALES[scale]
    lam = calibrate_censoring(law, cens); truth = population_truth(law)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"law{law}_{stem}.jsonl"; done = set()
    if path.exists():
        for line in path.read_text().splitlines():
            try:
                done.add(json.loads(line)["rep"])
            except Exception:
                pass
    t0 = time.time()
    with path.open("a") as handle:
        for rep in range(reps):
            if rep in done:
                continue
            if time.time() - t0 > time_limit:
                break
            rec = one_replication(C.stable_seed("onestep-sim-20261028", law, scale, rep), law, scale, lam, truth, boot)
            rec["rep"] = rep
            handle.write(json.dumps(rec) + "\n"); handle.flush()
    print(f"law {law}, scale {scale}: {len(done)} previously completed; output {path}")


def wilson(k, n, z=1.959963984540054):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def aggregate(rep_dir, out_csv):
    rows = []
    for law in ("A", "B"):
        truth = population_truth(law)
        for scale, (_, _, _, stem) in SCALES.items():
            path = rep_dir / f"law{law}_{stem}.jsonl"
            if not path.exists():
                continue
            recs = [json.loads(l) for l in path.read_text().splitlines() if l.strip()]
            good = [r for r in recs if not r.get("fatal")]; n = len(good)
            for est, label in (("P", "plug-in"), ("O", "one-step")):
                E = np.array([r[est] for r in good]); err = E - truth[None]
                k = sum(r[est + "_outer"] for r in good); lo, hi = wilson(k, n)
                W = np.array([r[est + "_width"] for r in good])
                for j, g in enumerate(GAMMAS):
                    rows.append({"law": law, "scale": stem, "estimator": label, "gamma": g,
                                 "n_valid": n, "n_fatal": len(recs) - n,
                                 "true_L": truth[j, 0], "true_U": truth[j, 1],
                                 "bias_L": err[:, j, 0].mean(), "bias_U": err[:, j, 1].mean(),
                                 "sd_L": E[:, j, 0].std(ddof=1), "sd_U": E[:, j, 1].std(ddof=1),
                                 "rmse_L": float(np.sqrt((err[:, j, 0] ** 2).mean())),
                                 "rmse_U": float(np.sqrt((err[:, j, 1] ** 2).mean())),
                                 "outer_cov": k / n, "wilson_lo": lo, "wilson_hi": hi,
                                 "mean_env_width": W[:, j].mean(),
                                 "mean_cens": float(np.mean([r["cens"] for r in good])),
                                 "rcap_point_reps": int(sum(r["rcap_point"] > 0 for r in good)),
                                 "floor_point_reps": int(sum(r["floor_point"] > 0 for r in good))})
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    print("wrote", out_csv, len(rows), "rows")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=["run", "aggregate", "smoke"])
    ap.add_argument("--law", choices=["A", "B"])
    ap.add_argument("--scale", choices=list(SCALES))
    ap.add_argument("--reps", type=int, default=1000)
    ap.add_argument("--boot", type=int, default=199)
    ap.add_argument("--time-limit", type=float, default=float("inf"))
    ap.add_argument("--replicate-dir", type=Path, default=RESULTS / "simulation_replicates")
    ap.add_argument("--summary", type=Path, default=RESULTS / "simulation_summary.csv")
    a = ap.parse_args()
    if a.command == "run":
        if not (a.law and a.scale):
            ap.error("run requires --law and --scale")
        run(a.law, a.scale, a.reps, a.boot, a.replicate_dir, a.time_limit)
    elif a.command == "aggregate":
        aggregate(a.replicate_dir, a.summary)
    else:
        tmp = REPOSITORY / "_check" / "onestep_smoke"
        for law in ("A", "B"):
            for scale in SCALES:
                run(law, scale, 3, 19, tmp, float("inf"))
        aggregate(tmp, tmp / "simulation_summary.csv")


if __name__ == "__main__":
    main()
