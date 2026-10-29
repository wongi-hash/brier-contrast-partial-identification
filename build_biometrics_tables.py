"""Build the v1.2.0 Biometrics Web Tables from committed results.

  tables/biometrics/Web_Table_01_scalar_vs_profile_varying.csv  <- results/scalar_diagnostics/web_table_1.csv
  tables/biometrics/Web_Table_16_onestep_simulation.csv         <- results/onestep/simulation_summary.csv
                                                                    results/onestep/population_limit_bias.json
  tables/biometrics/Web_Table_17_onestep_applications.csv       <- results/onestep/applications.csv
Decision codes ADOPT / KEEP_REFERENCE / DEFER are written as CAND / REF / DEFER
(CANDIDATE SUPPORTED / REFERENCE SUPPORTED / DEFER in the manuscript).
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "tables" / "biometrics"
CODE = {"ADOPT": "CAND", "KEEP_REFERENCE": "REF", "DEFER": "DEFER"}


def web_table_1():
    d = pd.read_csv(ROOT / "results" / "scalar_diagnostics" / "web_table_1.csv")
    order = {"NSCLC-Stanford": 0, "NSCLC-VA": 1, "Rotterdam-GBSG": 2}
    d = d.sort_values(["dataset", "seed"], key=lambda s: s.map(order) if s.name == "dataset" else s)
    out = pd.DataFrame({
        "Setting": d.dataset, "Seed": d.seed, "P_T(q1>q0)": d.P_T_q1_gt_q0.round(2),
        "Profile-varying width": d.profile_varying_width.round(5), "Scalar width": d.scalar_width.round(5),
        "Width ratio": d.width_ratio.round(3), "kappa": d.kappa.round(3), "Favored at gamma=0": d.favored_at_gamma0,
        "Robustness value profile-varying": d.robustness_profile_varying.round(3),
        "Robustness value scalar": d.robustness_scalar.round(3)})
    out.to_csv(OUT / "Web_Table_01_scalar_vs_profile_varying.csv", index=False)


def web_table_16():
    s = pd.read_csv(ROOT / "results" / "onestep" / "simulation_summary.csv")
    rows = []
    for (law, scale, est), g in s.groupby(["law", "scale", "estimator"], sort=False):
        bias = pd.concat([g.bias_L, g.bias_U]) * 1e3; sd = pd.concat([g.sd_L, g.sd_U]) * 1e3
        rmse = pd.concat([g.rmse_L, g.rmse_U]) * 1e3; first = g.iloc[0]
        rows.append({"Law": law, "Scale": scale, "Estimator": est,
                     "Bias x1e3 min": round(bias.min(), 1), "Bias x1e3 max": round(bias.max(), 1),
                     "SD x1e3 min": round(sd.min(), 1), "SD x1e3 max": round(sd.max(), 1),
                     "RMSE x1e3 min": round(rmse.min(), 1), "RMSE x1e3 max": round(rmse.max(), 1),
                     "Outer coverage": round(first.outer_cov, 3), "Wilson lower": round(first.wilson_lo, 3),
                     "Wilson upper": round(first.wilson_hi, 3),
                     "Mean envelope width at gamma=0.2 x1e3": round(g[g.gamma == 0.2].mean_env_width.iloc[0] * 1e3, 1),
                     "Replications": int(first.n_valid)})
    out = pd.DataFrame(rows)
    order = {"nsclc_scale": 0, "rotterdam_gbsg_scale": 1}
    out = out.sort_values(["Scale", "Law", "Estimator"], key=lambda c: c.map(order) if c.name == "Scale" else c,
                          ascending=[True, True, False])
    out.to_csv(OUT / "Web_Table_16_onestep_simulation.csv", index=False)
    limits = json.loads((ROOT / "results" / "onestep" / "population_limit_bias.json").read_text())
    (OUT / "Web_Table_16_population_limit_bias.json").write_text(json.dumps(limits, indent=2) + "\n")


def web_table_17():
    d = pd.read_csv(ROOT / "results" / "onestep" / "applications.csv")
    d = d[d.gamma.round(6) == 0.2]
    out = pd.DataFrame({
        "Setting": d.dataset, "Seed": d.seed,
        "Plug-in L": d.plugin_L.round(4), "Plug-in U": d.plugin_U.round(4),
        "One-step L": d.onestep_L.round(4), "One-step U": d.onestep_U.round(4),
        "Plug-in envelope lower": d.plugin_env_lower.round(4), "Plug-in envelope upper": d.plugin_env_upper.round(4),
        "One-step envelope lower": d.onestep_env_lower.round(4), "One-step envelope upper": d.onestep_env_upper.round(4),
        "Plug-in decision": d.plugin_action.map(CODE), "One-step decision": d.onestep_action.map(CODE)})
    out.to_csv(OUT / "Web_Table_17_onestep_applications.csv", index=False)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    web_table_1(); web_table_16(); web_table_17()
    print("wrote", sorted(p.name for p in OUT.iterdir()))
