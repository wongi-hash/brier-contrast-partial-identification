"""Core routines: primary plug-in and one-step (orthogonalized) endpoint estimators.

Nuisance fits replicate the primary implementation (conditional Cox-Breslow censoring
on H=Z, G floor 0.02, IPCW-weighted logistic on standardized (q0,q1)).
See code/onestep/PROTOCOL.md.
"""
import hashlib
import numpy as np

G_MIN = 0.02
R_CAP = 20.0
GAMMAS = np.array([0.0, 0.1, 0.2, 0.3, 0.5])


def expit(x):
    return 1.0 / (1.0 + np.exp(-np.clip(np.asarray(x, float), -35, 35)))


def logit(p):
    p = np.clip(np.asarray(p, float), 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


def stable_seed(*parts):
    return int(hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()[:8], 16)


def _cox(time, ce, z):
    """Unpenalized Cox PH (Breslow ties) for censoring; vectorized risk sets."""
    time = np.asarray(time, float); ce = np.asarray(ce, int); z = np.asarray(z, float)
    mu = z.mean(0); sd = z.std(0); sd[sd < 1e-8] = 1.0; x = (z - mu) / sd
    et = np.unique(time[ce == 1])
    if len(et) == 0:
        return mu, sd, np.zeros(x.shape[1]), np.array([]), np.array([])
    asc = np.argsort(time, kind="stable"); ta = time[asc]; xa = x[asc]
    i0 = np.searchsorted(ta, et, "left")
    d = np.array([int(((time == t) & (ce == 1)).sum()) for t in et], float)
    sumx = np.array([x[(time == t) & (ce == 1)].sum(0) for t in et])
    xaxa = xa[:, :, None] * xa[:, None, :]
    beta = np.zeros(x.shape[1])
    for _ in range(40):
        ra = np.exp(np.clip(xa @ beta, -25, 25)); r0 = np.cumsum(ra[::-1])[::-1]
        r1 = np.cumsum((ra[:, None] * xa)[::-1], 0)[::-1]
        r2 = np.cumsum((ra[:, None, None] * xaxa)[::-1], 0)[::-1]
        S0 = r0[i0]; m = r1[i0] / S0[:, None]
        info = np.sum(d[:, None, None] * (r2[i0] / S0[:, None, None] - m[:, :, None] * m[:, None, :]), 0)
        sc = sumx.sum(0) - np.sum(d[:, None] * m, 0)
        step = np.linalg.pinv(info, rcond=1e-10) @ sc
        if not np.all(np.isfinite(step)):
            raise FloatingPointError("nonfinite Cox step")
        mx = np.max(np.abs(step))
        if mx > 2:
            step *= 2 / mx
        beta += step
        if np.max(np.abs(step)) < 1e-8:
            break
    ra = np.exp(np.clip(xa @ beta, -25, 25)); r0 = np.cumsum(ra[::-1])[::-1]
    inc = d / r0[i0]
    return mu, sd, beta, et, inc


def cipcw(time, event, z, tau):
    time = np.asarray(time, float); event = np.asarray(event, int); z = np.asarray(z, float)
    mu, sd, beta, et, dh = _cox(time, 1 - event, z)
    ph = np.exp(np.clip(((z - mu) / sd) @ beta, -25, 25)); cum = np.cumsum(dh)

    def h0(q, left):
        idx = np.searchsorted(et, np.asarray(q, float), side=("left" if left else "right")) - 1
        o = np.zeros(len(np.atleast_1d(q))); ok = idx >= 0
        o[ok] = cum[np.minimum(idx[ok], len(cum) - 1)]
        return o

    graw_l = np.exp(-h0(time, True) * ph); graw_t = np.exp(-h0(np.repeat(tau, len(time)), False) * ph)
    gtl = np.clip(graw_l, G_MIN, 1.0); gtau = np.clip(graw_t, G_MIN, 1.0)
    y = (time > tau).astype(float); w = np.zeros(len(time))
    w[time > tau] = 1 / gtau[time > tau]
    f = (time <= tau) & (event == 1); w[f] = 1 / gtl[f]
    floor = int(np.sum((time > tau) & (graw_t < G_MIN)) + np.sum(f & (graw_l < G_MIN)))
    return y, w, floor


def logistic_newton(xd, y, w):
    b = np.zeros(xd.shape[1])
    for _ in range(60):
        p = expit(xd @ b); ww = np.clip(w * p * (1 - p), 1e-10, None)
        step = np.linalg.pinv((xd * ww[:, None]).T @ xd, rcond=1e-10) @ (xd.T @ (w * (y - p)))
        mx = np.max(np.abs(step))
        if mx > 3:
            step *= 3 / mx
        b += step
        if np.max(np.abs(step)) < 1e-8:
            break
    return b


def fit_nuisances(time, event, zs, zt, tau):
    """Returns m_hat on source and target, IPCW weights/outcomes, normalized r_hat on source."""
    zs = np.asarray(zs, float); zt = np.asarray(zt, float)
    y, w, floor = cipcw(time, event, zs, tau)
    mu = zs.mean(0); sd = zs.std(0); sd[sd < 1e-8] = 1.0
    xs = np.column_stack([np.ones(len(zs)), (zs - mu) / sd]); xt = np.column_stack([np.ones(len(zt)), (zt - mu) / sd])
    beta = logistic_newton(xs, y, w)
    ms = np.clip(expit(xs @ beta), 1e-4, 1 - 1e-4); mt = np.clip(expit(xt @ beta), 1e-4, 1 - 1e-4)
    # domain classifier on standardized (logit q0, logit q1), pooled
    lz = np.vstack([logit(zs), logit(zt)]); lmu = lz.mean(0); lsd = lz.std(0); lsd[lsd < 1e-8] = 1.0
    xd = np.column_stack([np.ones(len(lz)), (lz - lmu) / lsd])
    dlab = np.r_[np.zeros(len(zs)), np.ones(len(zt))]
    a = logistic_newton(xd, dlab, np.ones(len(lz)))
    eta_s = xd[: len(zs)] @ a
    r = np.exp(np.clip(eta_s, -30, 30)) * (len(zs) / len(zt))
    r = r / r.mean()
    cap = int(np.sum(r > R_CAP)); r = np.minimum(r, R_CAP)
    return {"ms": ms, "mt": mt, "y": y, "w": w, "r": r, "floor": floor, "rcap": cap}


def _f_and_df(q0, q1, m, gamma, which):
    a = q1 * q1 - q0 * q0; b = q1 - q0
    sL = np.where(b >= 0, 1.0, -1.0); s = sL if which == "L" else -sL
    lin = logit(m) + s * gamma; g = expit(lin)
    f = a - 2 * b * g
    df = -2 * b * g * (1 - g) / (m * (1 - m))
    return f, df


def endpoints_both(src_q, tgt_q, nu, gammas=GAMMAS, tgt_idx=None, src_idx=None):
    """Return (plugin[K,2], onestep[K,2]).  Optional index arrays reuse fitted nuisances on resampled rows."""
    q0t, q1t = tgt_q[:, 0], tgt_q[:, 1]; q0s, q1s = src_q[:, 0], src_q[:, 1]
    mt, ms, r, w, y = nu["mt"], nu["ms"], nu["r"], nu["w"], nu["y"]
    if tgt_idx is not None:
        q0t, q1t, mt = q0t[tgt_idx], q1t[tgt_idx], mt[tgt_idx]
    plug = np.empty((len(gammas), 2)); one = np.empty((len(gammas), 2))
    for j, g in enumerate(gammas):
        for k, which in enumerate(("L", "U")):
            ft, _ = _f_and_df(q0t, q1t, mt, g, which)
            _, dfs = _f_and_df(q0s, q1s, ms, g, which)
            plug[j, k] = ft.mean()
            one[j, k] = ft.mean() + np.mean(r * dfs * w * (y - ms))
    return plug, one


def sup_t_envelope(point, draws, level=0.95):
    se = np.maximum(draws.std(0, ddof=1), 1e-10)
    crit = float(np.quantile(np.abs((draws - point[None]) / se[None]).max(axis=(1, 2)), level))
    return point[:, 0] - crit * se[:, 0], point[:, 1] + crit * se[:, 1], crit, se


def action(lo, up):
    return "ADOPT" if up < 0 else ("KEEP_REFERENCE" if lo > 0 else "DEFER")
