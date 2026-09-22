"""Targeted primary-implementation coverage study at Rotterdam–GBSG scale.

This is intentionally a thin, provenance-preserving wrapper around the locked
NSCLC-scale implementation. Only the sample sizes, output directory, and target
censoring rate are changed:

  * n_source = 516 (exact Rotterdam calibration split size)
  * n_target = 686 (GBSG)
  * censoring = 194 / 516, observed in the locked Rotterdam calibration data

All DGM structure, sensitivity grid, conditional Cox--Breslow implementation,
IPCW-weighted logistic outcome regression, pairs bootstrap, and simultaneous
endpoint envelope are inherited unchanged.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import run_nsclc_scale as coverage


HERE = Path(__file__).resolve().parent
OUT = HERE / "rotterdam_gbsg_scale_results"
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
    coverage.OUT = OUT
    coverage.N_SOURCE = N_SOURCE
    coverage.N_TARGETS = (N_TARGET,)
    coverage.TARGET_CENSORING = OBSERVED_CENSORING
    # The baseline calibration depends on the requested censoring fraction.
    coverage.LAMBDA_C0, coverage.CALIBRATED_CENSORING = (
        coverage.calibrate_censoring_baseline()
    )


def finalize_manifest() -> None:
    manifest_path = OUT / "MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.update(
        {
            "status": "PASS",
            "study": "targeted primary-implementation stress study at Rotterdam–GBSG scale",
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
                "calibrated_rate": coverage.CALIBRATED_CENSORING,
            },
            "wrapper_sha256": sha256(Path(__file__)),
            "parent_runner_sha256": sha256(
                HERE / "run_nsclc_scale.py"
            ),
            "target_outcomes_used": "simulation truth only; no real target outcomes",
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> None:
    configure()
    coverage.main()
    if coverage.AGG_ONLY:
        finalize_manifest()
        print("[PASS]", OUT, flush=True)


if __name__ == "__main__":
    main()
