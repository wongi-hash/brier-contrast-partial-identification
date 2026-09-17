"""Task E2: Rotterdam-scale coverage arm.

This is intentionally a thin, provenance-preserving wrapper around the locked
Task E implementation.  Only the sample sizes, output directory, and target
censoring rate are changed:

  * n_source = 516 (exact Rotterdam calibration split size)
  * n_target = 686 (GBSG)
  * censoring = 194 / 516, observed in the locked Rotterdam calibration data

All DGM structure, gamma grid, conditional Cox--Breslow implementation,
IPCW-weighted logistic outcome regression, pairs bootstrap, and simultaneous
endpoint envelope are inherited unchanged from Task E.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import RUN_TASK_E_PRIMARY_COVERAGE as task_e


HERE = Path(__file__).resolve().parent
OUT = HERE / "task_e2_rotterdam_scale_coverage"
N_SOURCE = 516
N_TARGET = 686
OBSERVED_CENSORING = 194 / 516


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def configure() -> None:
    task_e.OUT = OUT
    task_e.N_SOURCE = N_SOURCE
    task_e.N_TARGETS = (N_TARGET,)
    task_e.TARGET_CENSORING = OBSERVED_CENSORING
    # The baseline calibration depends on the requested censoring fraction.
    task_e.LAMBDA_C0, task_e.CALIBRATED_CENSORING = (
        task_e.calibrate_censoring_baseline()
    )


def finalize_manifest() -> None:
    manifest_path = OUT / "MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.update(
        {
            "status": "PASS_TASK_E2",
            "task": "Rotterdam-scale conditional Cox--Breslow coverage",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "n_source": N_SOURCE,
            "n_targets": [N_TARGET],
            "sample_size_basis": {
                "source": "exact locked Rotterdam calibration split",
                "target": "GBSG cohort",
            },
            "censoring_rate_basis": {
                "observed_events": 322,
                "observed_censored": 194,
                "observed_n": 516,
                "target_rate": OBSERVED_CENSORING,
                "calibrated_rate": task_e.CALIBRATED_CENSORING,
            },
            "wrapper_sha256": sha256(Path(__file__)),
            "parent_task_e_runner_sha256": sha256(
                HERE / "RUN_TASK_E_PRIMARY_COVERAGE.py"
            ),
            "target_outcomes_used": "simulation truth only; no real target outcomes",
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> None:
    configure()
    task_e.main()
    if task_e.AGG_ONLY:
        finalize_manifest()
        print("[PASS TASK E2]", OUT, flush=True)


if __name__ == "__main__":
    main()
