"""
Factorial coverage study (Web Table 8): outer-envelope coverage of the primary
conditional Cox-Breslow implementation as source calibration size (n_S) and source
censoring vary, holding n_T = 80 fixed. Same data-generating mechanism as
code/simulations/primary_implementation_coverage/ (conditional Cox-Breslow IPCW +
weighted logistic + joint pairs-bootstrap sup-t). Deterministic seeds; Wilson CIs.

Cells: n_S in {114, 200, 300, 400, 516} x source censoring in {0.65, 0.50, 0.38}.
The n_S = 114 and 516 cells were run first; the 200/300/400 cells were added under
PROTOCOL.md (fixed before they were run). The code path is identical for all cells.

Full run (results/simulations/coverage_factorial/):
    python code/simulations/coverage_factorial/run_coverage_factorial.py
One cell, then aggregate:
    FACT_ONLY_CELL=300_50 python code/simulations/coverage_factorial/run_coverage_factorial.py
    FACT_AGG_ONLY=1 python code/simulations/coverage_factorial/run_coverage_factorial.py
Smoke (writes to _check/, ignored by Git):
    FACT_REPS=20 FACT_BOOT=29 FACT_OUT=_check/factorial python code/simulations/coverage_factorial/run_coverage_factorial.py
"""
import os, json, math, time, hashlib, platform
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = Path(os.environ["FACT_OUT"]) if os.environ.get("FACT_OUT") else ROOT / "results" / "simulations" / "coverage_factorial"
OUT.mkdir(parents=True, exist_ok=True)
REPS = int(os.environ.get("FACT_REPS", "1000"))
BOOT = int(os.environ.get("FACT_BOOT", "199"))
ONLY_CELL = os.environ.get("FACT_ONLY_CELL")           # e.g. "114_65"
AGG_ONLY = os.environ.get("FACT_AGG_ONLY") == "1"

GAMMAS = np.array([0.0, 0.1, 0.2, 0.3, 0.5]); TAU = 1.0; WK = 1.5
TRUE_SHIFT = 0.2; TXSHIFT = 0.5; ALPHA_C = np.array([1.0, -0.75]); G_MIN = 0.02
N_T = 80
CELLS = [(n, c) for n in (114, 200, 300, 400, 516) for c in (0.65, 0.50, 0.38)]

def expit(x): return 1.0/(1.0+np.exp(-np.clip(np.asarray(x, float), -35, 35)))
def logit(p): p = np.clip(np.asarray(p, float), 1e-6, 1-1e-6); return np.log(p/(1-p))
def qv(x): x = np.asarray(x, float); return expit(-0.15+0.45*x), expit(0.10+0.80*x)
def psrc(q0, q1): return expit(-0.35+1.25*q0-0.70*q1)
def stable_seed(*parts): return int(hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()[:8], 16)

def endpoints(q0, q1, ps):
    a = q1*q1-q0*q0; b = q1-q0; lp = logit(ps); out = np.empty((len(GAMMAS), 2))
    for j, g in enumerate(GAMMAS):
        pm, pp = expit(lp-g), expit(lp+g)
        out[j, 0] = np.mean(a-2*b*np.where(b >= 0, pp, pm))
        out[j, 1] = np.mean(a-2*b*np.where(b >= 0, pm, pp))
    return out

def calib_lambda(target_cens):
    rng = np.random.default_rng(20260907); x = rng.normal(0, 1, 500_000); q0, q1 = qv(x); ps = psrc(q0, q1)
    lam = -np.log(np.clip(ps, 1e-8, 1))/TAU**WK; t = (-np.log(rng.random(len(x)))/np.clip(lam, 1e-10, None))**(1/WK)
    z = np.column_stack([q0, q1]); z = (z-z.mean(0))/z.std(0); u = -np.log(rng.random(len(x))); mult = np.exp(z@ALPHA_C)
    lo, hi = 1e-4, 20.0
    for _ in range(60):
        mid = (lo+hi)/2; c = u/(mid*mult)
        if np.mean(c < t) < target_cens: lo = mid
        else: hi = mid
    return (lo+hi)/2

def population_truth():
    rng = np.random.default_rng(990011); x = rng.normal(TXSHIFT, 1, 2_000_000); q0, q1 = qv(x); ps = psrc(q0, q1)
    return endpoints(q0, q1, ps)
TRUE_EP = population_truth()

def _cox(time, ce, z):
    time = np.asarray(time, float); ce = np.asarray(ce, int); z = np.asarray(z, float)
    mu = z.mean(0); sd = z.std(0); sd[sd < 1e-8] = 1.0; x = (z-mu)/sd; et = np.unique(time[ce == 1])
    if len(et) == 0: return mu, sd, np.zeros(x.shape[1]), np.array([]), np.array([])
    asc = np.argsort(time); ta = time[asc]; xa = x[asc]; i0 = np.searchsorted(ta, et, 'left')
    d = np.array([int(((time == t) & (ce == 1)).sum()) for t in et], float)
    sumx = np.array([x[(time == t) & (ce == 1)].sum(0) for t in et]); xaxa = xa[:, :, None]*xa[:, None, :]
    beta = np.zeros(x.shape[1])
    for _ in range(40):
        ra = np.exp(np.clip(xa@beta, -25, 25)); r0 = np.cumsum(ra[::-1])[::-1]
        r1 = np.cumsum((ra[:, None]*xa)[::-1], 0)[::-1]; r2 = np.cumsum((ra[:, None, None]*xaxa)[::-1], 0)[::-1]
        S0 = r0[i0]; m = r1[i0]/S0[:, None]
        info = np.sum(d[:, None, None]*(r2[i0]/S0[:, None, None]-m[:, :, None]*m[:, None, :]), 0)
        sc = sumx.sum(0)-np.sum(d[:, None]*m, 0); step = np.linalg.pinv(info, rcond=1e-10)@sc
        mx = np.max(np.abs(step))
        if mx > 2: step *= 2/mx
        beta += step
        if np.max(np.abs(step)) < 1e-8: break
    ra = np.exp(np.clip(xa@beta, -25, 25)); r0 = np.cumsum(ra[::-1])[::-1]; inc = d/r0[i0]
    return mu, sd, beta, et, inc

def cipcw(time, event, z, tau):
    time = np.asarray(time, float); event = np.asarray(event, int); z = np.asarray(z, float)
    mu, sd, beta, et, dh = _cox(time, 1-event, z); ph = np.exp(np.clip(((z-mu)/sd)@beta, -25, 25)); cum = np.cumsum(dh)
    def h0(q, left):
        idx = np.searchsorted(et, np.asarray(q, float), side=('left' if left else 'right'))-1
        o = np.zeros(len(np.atleast_1d(q))); ok = idx >= 0; o[ok] = cum[np.minimum(idx[ok], len(cum)-1)]; return o
    gtl = np.clip(np.exp(-h0(time, True)*ph), G_MIN, 1.0); gtau = np.clip(np.exp(-h0(np.repeat(tau, len(time)), False)*ph), G_MIN, 1.0)
    y = (time > tau).astype(float); w = np.zeros(len(time)); w[time > tau] = 1/gtau[time > tau]
    f = (time <= tau) & (event == 1); w[f] = 1/gtl[f]
    return y, w

def wlog(x, y, w):
    xd = np.column_stack([np.ones(len(x)), x]); b = np.zeros(xd.shape[1])
    for _ in range(60):
        p = expit(xd@b); ww = np.clip(w*p*(1-p), 1e-10, None)
        step = np.linalg.pinv((xd*ww[:, None]).T@xd, rcond=1e-10)@(xd.T@(w*(y-p)))
        mx = np.max(np.abs(step))
        if mx > 3: step *= 3/mx
        b += step
        if np.max(np.abs(step)) < 1e-8: break
    return b

def fitps(ct, ce, cq0, cq1, tq0, tq1):
    z = np.column_stack([cq0, cq1]); y, w = cipcw(ct, ce, z, TAU); mu = z.mean(0); sd = z.std(0); sd[sd < 1e-8] = 1.0
    beta = wlog((z-mu)/sd, y, w); zt = np.column_stack([tq0, tq1])
    return np.clip(expit(beta[0]+((zt-mu)/sd)@beta[1:]), 1e-4, 1-1e-4)

def one_rep(seed, nS, nT, lam_c0):
    rng = np.random.default_rng(seed); xs = rng.normal(0, 1, nS); xt = rng.normal(TXSHIFT, 1, nT)
    q0s, q1s = qv(xs); q0t, q1t = qv(xt); ps = psrc(q0s, q1s)
    lam = -np.log(np.clip(ps, 1e-8, 1))/TAU**WK; fail = (-np.log(rng.random(nS))/np.clip(lam, 1e-10, None))**(1/WK)
    z = np.column_stack([q0s, q1s]); zs = (z-z.mean(0))/np.where(z.std(0) < 1e-8, 1, z.std(0))
    cens = -np.log(rng.random(nS))/(lam_c0*np.exp(zs@ALPHA_C)); time = np.minimum(fail, cens); ev = (fail <= cens).astype(int)
    crate = float(1-ev.mean())
    try: psh = fitps(time, ev, q0s, q1s, q0t, q1t)
    except Exception: return {"fatal": 1, "cens": crate}
    point = endpoints(q0t, q1t, psh); draws = np.empty((BOOT, len(GAMMAS), 2)); v = 0
    rb = np.random.default_rng(seed + 800_000_003)
    for _ in range(BOOT):
        si = rb.integers(0, nS, nS); ti = rb.integers(0, nT, nT)
        try: psb = fitps(time[si], ev[si], q0s[si], q1s[si], q0t, q1t); draws[v] = endpoints(q0t[ti], q1t[ti], psb[ti]); v += 1
        except Exception: pass
    if v < int(0.9*BOOT): return {"fatal": 1, "cens": crate}
    d = draws[:v]; se = np.maximum(d.std(0, ddof=1), 1e-10)
    crit = float(np.quantile(np.abs((d-point[None])/se[None]).max(axis=(1, 2)), 0.95))
    lo = point[:, 0]-crit*se[:, 0]; up = point[:, 1]+crit*se[:, 1]
    outer = int(np.all((lo <= TRUE_EP[:, 0]) & (TRUE_EP[:, 1] <= up)))
    return {"fatal": 0, "outer": outer, "cens": crate, "crit": crit}

def wilson(k, n, z=1.959963984540054):
    if n == 0: return (0.0, 1.0)
    p = k/n; d = 1+z*z/n; c = (p+z*z/(2*n))/d; h = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return c-h, c+h

def cell_tag(nS, cens): return f"{nS}_{int(round(cens*100))}"

def run_cell(nS, cens):
    lam = calib_lambda(cens); recs = []
    for i in range(REPS):
        recs.append(one_rep(stable_seed("factorial", nS, int(round(cens*100)), i), nS, N_T, lam))
        if (i+1) % 100 == 0: print(f"  [{cell_tag(nS,cens)}] {i+1}/{REPS}", flush=True)
    payload = {"meta": {"nS": nS, "target_cens": cens, "nT": N_T, "reps": REPS, "boot": BOOT, "lambda_c0": lam}, "records": recs}
    (OUT / f"cell_{cell_tag(nS,cens)}.json").write_text(json.dumps(payload), encoding="utf-8")
    print(f"  saved cell_{cell_tag(nS,cens)}.json")

def aggregate():
    import csv
    rows = []
    for nS, cens in CELLS:
        p = OUT / f"cell_{cell_tag(nS,cens)}.json"
        if not p.is_file(): raise FileNotFoundError(p)
        v = json.loads(p.read_text()); good = [r for r in v["records"] if not r["fatal"]]
        n = len(good); k = sum(r["outer"] for r in good); lo, hi = wilson(k, n)
        rows.append({"nS": nS, "target_cens": cens, "nT": N_T, "reps_valid": n,
                     "fatal": len(v["records"])-n, "mean_censoring": round(float(np.mean([r["cens"] for r in good])), 3),
                     "mean_source_events": round(nS*(1-float(np.mean([r["cens"] for r in good]))), 1),
                     "outer_envelope_coverage": round(k/n, 4), "wilson_lo": round(lo, 4), "wilson_hi": round(hi, 4),
                     "mean_critical_value": round(float(np.mean([r["crit"] for r in good])), 3)})
    with open(OUT / "factorial_coverage_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); [w.writerow(r) for r in rows]
    manifest = {"status": "PASS_FACTORIAL", "reps": REPS, "boot": BOOT, "n_T": N_T, "cells": [list(c) for c in CELLS],
                "gammas": GAMMAS.tolist(), "true_shift": "0.2*tanh(x)", "dgm": "identical to code/simulations/primary_implementation_coverage",
                "python": platform.python_version(), "numpy": np.__version__}
    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("\n n_S  cens  outer_cov   Wilson CI        mean_cens")
    for r in rows:
        print(f" {r['nS']:>4} {r['target_cens']:.2f}  {r['outer_envelope_coverage']:.4f}   [{r['wilson_lo']:.4f},{r['wilson_hi']:.4f}]  {r['mean_censoring']:.3f}")
    print("[PASS]", OUT / "factorial_coverage_summary.csv")

def main():
    t0 = time.time()
    if AGG_ONLY: aggregate(); return
    cells = CELLS
    if ONLY_CELL:
        nS, c = ONLY_CELL.split("_"); cells = [(int(nS), int(c)/100.0)]
    for nS, cens in cells:
        print(f"[cell n_S={nS} cens={cens}] reps={REPS} boot={BOOT}", flush=True); run_cell(nS, cens)
    if not ONLY_CELL and not AGG_ONLY: aggregate()
    print(f"[elapsed {time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()
