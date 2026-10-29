"""Population limits of the plug-in and one-step biases (Web Table 16 notes; Web Appendix A.7).

The working logistic model is fitted to the true source law on 2,000,000 source draws; the
plug-in limit is E_T f(Z, m~) - E_T f(Z, m_S) and the one-step limit (true density ratio) is
E_T[f(Z, m~) + f'(Z, m~)(m_S - m~)] - E_T f(Z, m_S), over the grid and both endpoints.

Usage: python code/onestep/population_limit_bias.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import onestep_core as C  # noqa: E402
import run_onestep_simulation as S  # noqa: E402

OUT = HERE.parents[1] / "results" / "onestep" / "population_limit_bias.json"


def main():
    rng = np.random.default_rng(7)
    xs = rng.normal(0, 1, 2_000_000); xt = rng.normal(S.TARGET_X_SHIFT, 1, 2_000_000)
    out = {}
    for law in ("A", "B"):
        q0s, q1s = S.q_values(xs); q0t, q1t = S.q_values(xt)
        ms = S.source_law(q0s, q1s, law); mt = S.source_law(q0t, q1t, law)
        z = np.column_stack([q0s, q1s]); mu, sd = z.mean(0), z.std(0)
        beta = C.logistic_newton(np.column_stack([np.ones(len(z)), (z - mu) / sd]), ms, np.ones(len(z)))
        mtil = C.expit(np.column_stack([np.ones(len(xt)), (np.column_stack([q0t, q1t]) - mu) / sd]) @ beta)
        plug, one = [], []
        for g in C.GAMMAS:
            for endpoint in ("L", "U"):
                f_true, _ = C._f_and_df(q0t, q1t, mt, g, endpoint)
                f_til, df = C._f_and_df(q0t, q1t, mtil, g, endpoint)
                plug.append(f_til.mean() - f_true.mean())
                one.append(np.mean(f_til + df * (mt - mtil)) - f_true.mean())
        out[law] = {"plugin_bias_x1e3_range": [min(plug) * 1e3, max(plug) * 1e3],
                    "onestep_bias_x1e3_range": [min(one) * 1e3, max(one) * 1e3]}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
