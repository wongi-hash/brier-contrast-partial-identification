"""Finite-sample coverage diagnosis (Web Appendix C.4).

Exploratory analysis added after the factorial results were available. For a factorial
cell it re-runs the replicate function of run_coverage_factorial.py with the same seeds and
stores, per replication, the plug-in endpoints, bootstrap standard errors, sup-t critical
value, and the target-sample endpoints computed with the TRUE source survival probability.
The recomputed envelopes reproduce the stored factorial records exactly (checked in
`summarize`).

Usage
  python code/simulations/coverage_factorial/diagnose_coverage.py run 114 65 [--time-limit SEC]
  python code/simulations/coverage_factorial/diagnose_coverage.py run 516 50 [--time-limit SEC]
  python code/simulations/coverage_factorial/diagnose_coverage.py summarize
Runs are resumable: completed replications in the output file are skipped.
"""
import argparse, json, math, runpy, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RES = ROOT / "results" / "simulations" / "coverage_factorial"
OUT = RES / "diagnosis"
F = runpy.run_path(str(HERE / "run_coverage_factorial.py"), run_name="factorial")
CELLS = ((114, 65), (516, 50))


def replicate(seed, nS, nT, lam_c0):
    """Identical random-number stream to run_coverage_factorial.one_rep."""
    g = F; GAM = g["GAMMAS"]; BOOT = g["BOOT"]
    rng = np.random.default_rng(seed); xs = rng.normal(0, 1, nS); xt = rng.normal(g["TXSHIFT"], 1, nT)
    q0s, q1s = g["qv"](xs); q0t, q1t = g["qv"](xt); ps = g["psrc"](q0s, q1s)
    lam = -np.log(np.clip(ps, 1e-8, 1)) / g["TAU"] ** g["WK"]
    fail = (-np.log(rng.random(nS)) / np.clip(lam, 1e-10, None)) ** (1 / g["WK"])
    z = np.column_stack([q0s, q1s]); zs = (z - z.mean(0)) / np.where(z.std(0) < 1e-8, 1, z.std(0))
    cens = -np.log(rng.random(nS)) / (lam_c0 * np.exp(zs @ g["ALPHA_C"]))
    time_ = np.minimum(fail, cens); ev = (fail <= cens).astype(int)
    psh = g["fitps"](time_, ev, q0s, q1s, q0t, q1t); point = g["endpoints"](q0t, q1t, psh)
    oracle_m = g["endpoints"](q0t, q1t, g["psrc"](q0t, q1t))
    draws = np.empty((BOOT, len(GAM), 2)); v = 0; rb = np.random.default_rng(seed + 800_000_003)
    for _ in range(BOOT):
        si = rb.integers(0, nS, nS); ti = rb.integers(0, nT, nT)
        try:
            psb = g["fitps"](time_[si], ev[si], q0s[si], q1s[si], q0t, q1t)
            draws[v] = g["endpoints"](q0t[ti], q1t[ti], psb[ti]); v += 1
        except Exception:
            pass
    d = draws[:v]; se = np.maximum(d.std(0, ddof=1), 1e-10)
    crit = float(np.quantile(np.abs((d - point[None]) / se[None]).max(axis=(1, 2)), 0.95))
    return dict(point=point.tolist(), se=se.tolist(), crit=crit, oracle_m=oracle_m.tolist(), events=int(ev.sum()))


def run(nS, pct, limit):
    OUT.mkdir(parents=True, exist_ok=True); path = OUT / f"replicates_{nS}_{pct}.jsonl"
    done = {json.loads(l)["rep"] for l in path.read_text().splitlines()} if path.exists() else set()
    lam = F["calib_lambda"](pct / 100); t0 = time.time()
    with path.open("a", encoding="utf-8", newline="\n") as fh:
        for i in range(F["REPS"]):
            if i in done: continue
            if limit and time.time() - t0 > limit: break
            r = replicate(F["stable_seed"]("factorial", nS, pct, i), nS, F["N_T"], lam); r["rep"] = i
            fh.write(json.dumps(r) + "\n"); fh.flush()
    print(nS, pct, "replications:", sum(1 for _ in path.open()))


def summarize():
    T = F["TRUE_EP"]; GAM = F["GAMMAS"]; out = {}
    for nS, pct in CELLS:
        R = sorted((json.loads(l) for l in (OUT / f"replicates_{nS}_{pct}.jsonl").read_text().splitlines()), key=lambda r: r["rep"])
        P = np.array([r["point"] for r in R]); S = np.array([r["se"] for r in R]); C = np.array([r["crit"] for r in R])
        O = np.array([r["oracle_m"] for r in R])
        lo = P[:, :, 0] - C[:, None] * S[:, :, 0]; up = P[:, :, 1] + C[:, None] * S[:, :, 1]
        lf = (lo > T[None, :, 0]).any(1); uf = (up < T[None, :, 1]).any(1); miss = lf | uf
        stored = json.loads((RES / f"cell_{nS}_{pct}.json").read_text())["records"]
        reproduces = bool(all(int(not m) == s["outer"] for m, s in zip(miss, stored)) and np.allclose(C, [s["crit"] for s in stored]))
        bias = P.mean(0) - T; sd = P.std(0, ddof=1); se = S.mean(0); obias = O.mean(0) - T

        def cov(l, u): return float(1 - ((l > T[None, :, 0]) | (u < T[None, :, 1])).any(1).mean())
        # (a) Monte Carlo SD in place of each replication's bootstrap SEs; each replication keeps its own bootstrap critical value
        cov_sd = cov(P[:, :, 0] - C[:, None] * sd[None, :, 0], P[:, :, 1] + C[:, None] * sd[None, :, 1])
        # (b) Monte Carlo mean bias (per gamma and endpoint) subtracted from both limits; bootstrap SEs and critical values unchanged
        cov_bias = cov(lo - bias[None, :, 0], up - bias[None, :, 1])
        j = list(GAM).index(0.2)
        out[f"{nS}_{pct}"] = {
            "n_rep": len(R), "mean_source_events": float(np.mean([r["events"] for r in R])),
            "outer_coverage": float(1 - miss.mean()), "reproduces_stored_factorial_records": reproduces,
            "misses_lower_only": int((lf & ~uf).sum()), "misses_upper_only": int((uf & ~lf).sum()), "misses_both": int((lf & uf).sum()),
            "gamma": GAM.tolist(),
            "bias_lower": bias[:, 0].tolist(), "bias_upper": bias[:, 1].tolist(),
            "bias_lower_over_sd": (bias[:, 0] / sd[:, 0]).tolist(), "bias_upper_over_sd": (bias[:, 1] / sd[:, 1]).tolist(),
            "bias_with_true_mS_lower": obias[:, 0].tolist(), "bias_with_true_mS_upper": obias[:, 1].tolist(),
            "mc_sd_lower": sd[:, 0].tolist(), "mc_sd_upper": sd[:, 1].tolist(),
            "se_over_sd_lower": (se[:, 0] / sd[:, 0]).tolist(), "se_over_sd_upper": (se[:, 1] / sd[:, 1]).tolist(),
            "corr_signed_error_bootstrap_se_lower_gamma0.2": float(np.corrcoef(P[:, j, 0] - T[j, 0], S[:, j, 0])[0, 1]),
            "coverage_mc_sd_in_place_of_bootstrap_se": cov_sd, "coverage_mc_bias_removed": cov_bias,
        }
    # Path.write_text() has no ``newline`` argument on the Python 3.8/3.9
    # versions supported by this release.  Use an explicit text handle so the
    # archived summary is still written with deterministic LF line endings.
    with (OUT / "diagnosis_summary.json").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write(json.dumps(out, indent=2))
    for k, v in out.items():
        print(k, {a: v[a] for a in ("outer_coverage", "reproduces_stored_factorial_records", "misses_lower_only", "misses_upper_only",
                                     "misses_both", "coverage_mc_sd_in_place_of_bootstrap_se", "coverage_mc_bias_removed")})


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("action", choices=["run", "summarize"])
    ap.add_argument("nS", nargs="?", type=int); ap.add_argument("pct", nargs="?", type=int); ap.add_argument("--time-limit", type=float, default=0)
    a = ap.parse_args()
    if a.action == "run": run(a.nS, a.pct, a.time_limit)
    else: summarize()
