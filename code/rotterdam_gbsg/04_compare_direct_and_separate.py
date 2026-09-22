"""Compute direct and separate-risk bounds for the Rotterdam–GBSG analysis."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
CORE = HERE.parent / "inference" / "direct_vs_separate.py"
spec = importlib.util.spec_from_file_location("direct_vs_separate", CORE)
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bootstrap", type=int, default=999)
    parser.add_argument("--out", default="direct_vs_separate_results")
    args = parser.parse_args()
    rows, metadata = [], []
    for seed in core.SEEDS:
        folder = HERE / "analysis_inputs" / f"seed{seed}"
        calibration_path = folder / "calibration_GBSG.csv"
        target_path = folder / "target_GBSG.csv"
        calibration = pd.read_csv(calibration_path).rename(
            columns={"q0_surv5y": "q0", "q1_surv5y": "q1"}
        )
        target = pd.read_csv(target_path).rename(
            columns={"q0_surv5y": "q0", "q1_surv5y": "q1"}
        )
        if {"time", "event"} & set(target):
            raise RuntimeError("target outcome firewall violation")
        result, meta = core.analyze(
            "Rotterdam-GBSG", "supporting_clinical_candidate", seed,
            calibration, target, 365.25 * 5, args.bootstrap
        )
        rows.extend(result)
        metadata.append(meta)
    output_directory = HERE / args.out
    output_directory.mkdir(exist_ok=True)
    output = output_directory / "direct_vs_separate.csv"
    pd.DataFrame(rows).to_csv(output, index=False)
    manifest = {
        "status": "PASS",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "target_outcomes_read": False,
        "bootstrap_draws": args.bootstrap,
        "core_sha256": sha(CORE),
        "output_sha256": sha(output),
        "jobs": metadata,
    }
    (output_directory / "MANIFEST.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print("[PASS]", len(rows), "rows ->", output)


if __name__ == "__main__":
    main()
