"""Retrospective outcome comparison for the Rotterdam–GBSG supporting analysis."""
import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import SEEDS, dump_json, oracle_metrics, sha256, stable_seed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=1999)
    ap.add_argument("--confirm-open-target-outcomes", action="store_true")
    args = ap.parse_args()
    if not args.confirm_open_target_outcomes:
        raise RuntimeError("use --confirm-open-target-outcomes")
    gate_path = HERE / "outcome_free_results" / "GATE_MANIFEST.json"
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    if gate.get("status") != "PASS" or gate.get("bootstrap") != 999:
        raise RuntimeError("complete B=999 outcome-free analysis required")
    oracle_path = HERE / "sealed_oracle" / "target_GBSG_outcomes.csv"
    oracle = pd.read_csv(oracle_path)
    rows = []
    for seed in SEEDS:
        target = pd.read_csv(HERE / "analysis_inputs" / ("seed{}".format(seed)) / "target_GBSG.csv")
        data = target.merge(oracle, on="patient_id", validate="one_to_one")
        t, e = data.time.to_numpy(float), data.event.to_numpy(int)
        q0, q1 = data.q0_surv5y.to_numpy(float), data.q1_surv5y.to_numpy(float)
        point = oracle_metrics(t, e, q0, q1)
        rng = np.random.default_rng(stable_seed("retrospective-outcome-comparison", seed, args.bootstrap))
        boots, failures = [], 0
        for b in range(args.bootstrap):
            ix = rng.integers(0, len(data), len(data))
            try: boots.append(oracle_metrics(t[ix], e[ix], q0[ix], q1[ix]))
            except Exception: failures += 1
            if (b + 1) % 200 == 0:
                print("seed {} oracle {}/{} valid={}".format(seed, b + 1, args.bootstrap, len(boots)), flush=True)
        if len(boots) < int(0.9 * args.bootstrap):
            raise RuntimeError("oracle bootstrap failure")
        for key, value in point.items():
            vals = np.array([x[key] for x in boots])
            rows.append({"seed": seed, "metric": key, "estimate": value,
                         "CI_low": float(np.quantile(vals, .025)), "CI_high": float(np.quantile(vals, .975)),
                         "bootstrap_valid": len(vals), "bootstrap_failures": failures})
    result = pd.DataFrame(rows)
    out = HERE / "retrospective_outcome_results"; out.mkdir(exist_ok=True)
    result_path = out / "outcome_metrics.csv"; result.to_csv(result_path, index=False)
    dump_json({"status": "PASS",
               "created_utc": datetime.now(timezone.utc).isoformat(),
               "target_outcomes_read": True, "target_outcomes_had_been_accessed_in_antecedent_exploratory_work": True,
               "gate_manifest_sha256": sha256(gate_path), "oracle_input_sha256": sha256(oracle_path),
               "runner_sha256": sha256(__file__), "common_sha256": sha256(HERE / "common.py"),
               "result_sha256": sha256(result_path), "bootstrap": args.bootstrap,
               "interpretation": "Supporting external demonstration; not prospectively blinded clinical validation.",
               "python": sys.version, "platform": platform.platform()}, out / "OUTCOME_COMPARISON_MANIFEST.json")
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
