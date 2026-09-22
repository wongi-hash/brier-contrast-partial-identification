"""Run the outcome-free Rotterdam–GBSG supporting analysis."""
import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import SEEDS, dump_json, gate_one, overlap_diagnostics, sha256

FORBIDDEN = {"time", "event", "status", "rfstime", "death", "recur", "dtime", "rtime"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=999)
    ap.add_argument("--out", default="outcome_free_results")
    args = ap.parse_args()
    packaging_path = HERE / "PACKAGING_MANIFEST.json"
    packaging = json.loads(packaging_path.read_text(encoding="utf-8"))
    if packaging["target_outcomes_exported_to_analysis"] is not False:
        raise RuntimeError("packaging firewall failed")
    if packaging["specification_sha256"] != sha256(HERE / "supporting_analysis_specification.md"):
        raise RuntimeError("supporting-analysis specification changed after packaging")
    out = HERE / args.out
    out.mkdir(exist_ok=True)
    reports, dense_frames, summaries, overlaps, inputs = [], [], [], [], []
    for seed in SEEDS:
        d = HERE / "analysis_inputs" / ("seed{}".format(seed))
        cp, tp = d / "calibration_GBSG.csv", d / "target_GBSG.csv"
        cal, target = pd.read_csv(cp), pd.read_csv(tp)
        if FORBIDDEN.intersection(target.columns):
            raise RuntimeError("target outcome leakage")
        if set(target.columns) != {"patient_id", "q0_surv5y", "q1_surv5y"}:
            raise RuntimeError("unexpected target schema")
        print("[Rotterdam-GBSG seed {}] n_cal={} n_target={} B={}".format(seed, len(cal), len(target), args.bootstrap), flush=True)
        report, dense, summary = gate_one(cal, target, seed, boot=args.bootstrap, progress=True)
        reports.append(report); dense_frames.append(dense); summaries.append(summary)
        overlaps.append(overlap_diagnostics(cal, target, seed))
        inputs += [{"path": str(cp.relative_to(HERE)), "sha256": sha256(cp)}, {"path": str(tp.relative_to(HERE)), "sha256": sha256(tp)}]
    report = pd.concat(reports, ignore_index=True)
    dense = pd.concat(dense_frames, ignore_index=True)
    overlap = pd.DataFrame(overlaps)
    rp, dp, op, sp = out / "gamma_report.csv", out / "gamma_dense.csv", out / "overlap_diagnostics.csv", out / "seed_summary.json"
    report.to_csv(rp, index=False); dense.to_csv(dp, index=False); overlap.to_csv(op, index=False); dump_json(summaries, sp)
    manifest = {
        "status": "PASS", "created_utc": datetime.now(timezone.utc).isoformat(),
        "study_role": "supporting Rotterdam–GBSG survival analysis",
        "target_outcomes_read_by_this_script": False,
        "target_outcomes_previously_opened_in_prior_analysis": True,
        "bootstrap": args.bootstrap,
        "specification_sha256": sha256(HERE / "supporting_analysis_specification.md"),
        "packaging_manifest_sha256": sha256(packaging_path),
        "runner_sha256": sha256(__file__), "common_sha256": sha256(HERE / "common.py"),
        "inputs": inputs,
        "outputs": [{"path": str(p.relative_to(HERE)), "sha256": sha256(p)} for p in [rp, dp, op, sp]],
        "python": sys.version, "platform": platform.platform(),
    }
    dump_json(manifest, out / "GATE_MANIFEST.json")
    print("\n[OUTCOME-FREE EVALUATION]"); print(report.to_string(index=False))
    print("\n[OVERLAP]"); print(overlap.to_string(index=False))


if __name__ == "__main__":
    main()
