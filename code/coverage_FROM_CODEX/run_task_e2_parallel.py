"""Run the four locked Task E2 seed shards, then aggregate them."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
RUNNER = HERE / "RUN_TASK_E2_ROTTERDAM_SCALE_COVERAGE.py"
LOG_DIR = HERE / "task_e2_logs"
STATUS = HERE / "TASK_E2_STATUS.json"
SEEDS = (42, 52, 62, 72)


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_status(status: str, completed: int, failed: int, **extra) -> None:
    payload = {
        "task": "E2_rotterdam_scale_coverage",
        "status": status,
        "completed": completed,
        "total": len(SEEDS),
        "failed": failed,
        "updated_utc": utcnow(),
        **extra,
    }
    STATUS.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main() -> None:
    LOG_DIR.mkdir(exist_ok=True)
    write_status("RUNNING", 0, 0, started_utc=utcnow())
    jobs = []
    for seed in SEEDS:
        env = os.environ.copy()
        env.update(
            {
                "TASK_E_NT": "686",
                "TASK_E_SEED": str(seed),
                "TASK_E_REPS": "1000",
                "TASK_E_BOOT": "199",
            }
        )
        stdout = (LOG_DIR / f"seed{seed}.stdout.log").open("w", encoding="utf-8")
        stderr = (LOG_DIR / f"seed{seed}.stderr.log").open("w", encoding="utf-8")
        proc = subprocess.Popen(
            [sys.executable, "-u", str(RUNNER)],
            cwd=HERE,
            env=env,
            stdout=stdout,
            stderr=stderr,
        )
        jobs.append((seed, proc, stdout, stderr))

    completed = 0
    failed = 0
    pending = list(jobs)
    while pending:
        next_pending = []
        for seed, proc, stdout, stderr in pending:
            code = proc.poll()
            if code is None:
                next_pending.append((seed, proc, stdout, stderr))
                continue
            stdout.close()
            stderr.close()
            completed += 1
            failed += int(code != 0)
            write_status("RUNNING", completed, failed, last_seed=seed, returncode=code)
        pending = next_pending
        if pending:
            time.sleep(15)

    if failed:
        write_status("FAILED", completed, failed)
        raise SystemExit(1)

    env = os.environ.copy()
    env.update(
        {
            "TASK_E_AGG_ONLY": "1",
            "TASK_E_REPS": "1000",
            "TASK_E_BOOT": "199",
        }
    )
    aggregate = subprocess.run(
        [sys.executable, "-u", str(RUNNER)], cwd=HERE, env=env, check=False
    )
    final_status = "PASS" if aggregate.returncode == 0 else "FAILED_AGGREGATION"
    write_status(
        final_status,
        completed,
        failed,
        aggregation_returncode=aggregate.returncode,
        finished_utc=utcnow(),
    )
    raise SystemExit(aggregate.returncode)


if __name__ == "__main__":
    main()
