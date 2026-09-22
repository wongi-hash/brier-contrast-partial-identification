"""Run the public Rotterdam–GBSG supporting-analysis pipeline."""
from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent


def run(command):
    print("[RUN]", " ".join(map(str, command)), flush=True)
    subprocess.run(command, cwd=HERE, check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bootstrap", type=int, default=999)
    parser.add_argument("--skip-export", action="store_true")
    parser.add_argument("--open-outcomes", action="store_true")
    args = parser.parse_args()
    if not args.skip_export:
        run(["Rscript", str(HERE / "00_export_survival_data.R")])
        for name in ("rotterdam.csv", "gbsg.csv"):
            path = HERE / name
            pd.read_csv(path).to_csv(path, index=False, lineterminator="\n")
    run([sys.executable, "01_prepare_locked_predictors.py"])
    run([sys.executable, "02_run_outcome_free_evaluation.py", "--bootstrap", str(args.bootstrap)])
    run([sys.executable, "04_compare_direct_and_separate.py", "--bootstrap", str(args.bootstrap)])
    if args.open_outcomes:
        run([
            sys.executable,
            "03_run_retrospective_outcome_comparison.py",
            "--confirm-open-target-outcomes",
        ])
    print("[PASS] Rotterdam–GBSG supporting analysis")


if __name__ == "__main__":
    main()
