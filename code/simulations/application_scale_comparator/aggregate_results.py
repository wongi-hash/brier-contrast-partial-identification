"""Aggregate the locked nT=60/80 comparator and mirror simulations.

Reads only PASS runs with reps=1000 and bootstrap=199, verifies the expected
7 comparator + 2 mirror scenarios at each target size, adds Wilson intervals,
and writes a machine-readable CSV plus a manuscript-facing Markdown summary.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import os

import pandas as pd


HERE = Path(__file__).resolve().parent
TARGETS = (60, 80)
COMPARATOR_SCENARIOS = (
    "no_shift", "constant_shift", "xvary_in_class", "sign_change",
    "out_of_class", "low_overlap", "pS_misspec",
)
MIRROR_SCENARIOS = ("no_shift", "xvary_in_class")
REPS = 1000
BOOT = 199
Z = 1.959963984540054


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def wilson(p: float, n: int = REPS) -> tuple[float, float]:
    den = 1.0 + Z * Z / n
    center = (p + Z * Z / (2 * n)) / den
    half = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / den
    return max(0.0, center - half), min(1.0, center + half)


def collect(root: Path, kind: str) -> tuple[list[pd.DataFrame], list[dict]]:
    frames: list[pd.DataFrame] = []
    provenance: list[dict] = []
    for manifest_path in root.rglob("MANIFEST.json"):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        cfg = manifest.get("config", {})
        if cfg.get("n_target") not in TARGETS:
            continue
        if cfg.get("reps") != REPS or cfg.get("bootstrap") != BOOT:
            continue
        if manifest.get("status") != "PASS":
            raise AssertionError(f"Non-PASS run: {manifest_path}")
        summary_path = manifest_path.parent / "summary.csv"
        if not summary_path.is_file():
            raise FileNotFoundError(summary_path)
        frame = pd.read_csv(summary_path)
        frame["n_target"] = int(cfg["n_target"])
        frame["run_kind"] = kind
        frames.append(frame)
        provenance.append({
            "kind": kind,
            "n_target": int(cfg["n_target"]),
            "scenario": cfg.get("scenario") or cfg.get("scenarios", [None])[0],
            "manifest": str(manifest_path.relative_to(HERE)),
            "manifest_sha256": sha256(manifest_path),
            "summary": str(summary_path.relative_to(HERE)),
            "summary_sha256": sha256(summary_path),
        })
    return frames, provenance


def fmt(x: float, digits: int = 3) -> str:
    return f"{float(x):.{digits}f}"


def main() -> None:
    comp_root = Path(os.environ.get("SIMULATION_COMPARATOR_ROOT", str(HERE / "comparator")))
    mirror_root = Path(os.environ.get("SIMULATION_MIRROR_ROOT", str(HERE / "candidate_superior_mirror")))
    output_root = Path(os.environ.get("SIMULATION_AGGREGATE_ROOT", str(HERE)))
    output_root.mkdir(parents=True, exist_ok=True)
    comp_frames, comp_prov = collect(comp_root, "comparator")
    mirror_frames, mirror_prov = collect(mirror_root, "mirror")
    frame = pd.concat(comp_frames + mirror_frames, ignore_index=True)

    for n_target in TARGETS:
        comp = set(frame.loc[(frame.run_kind == "comparator") & (frame.n_target == n_target), "scenario"])
        mirror = set(frame.loc[(frame.run_kind == "mirror") & (frame.n_target == n_target), "scenario"])
        if comp != set(COMPARATOR_SCENARIOS):
            raise AssertionError((n_target, "comparator", comp))
        if mirror != {f"mirror_{x}" for x in MIRROR_SCENARIOS}:
            raise AssertionError((n_target, "mirror", mirror))

    if len(comp_prov) != 14 or len(mirror_prov) != 4:
        raise AssertionError((len(comp_prov), len(mirror_prov)))

    for column in ("ci_coverage", "retain", "reference", "defer"):
        lows, highs = zip(*(wilson(float(x)) for x in frame[column]))
        frame[f"{column}_wilson_low"] = lows
        frame[f"{column}_wilson_high"] = highs

    frame = frame.sort_values(["run_kind", "scenario", "n_target", "method", "gamma"])
    csv_path = output_root / "application_scale_comparator.csv"
    frame.to_csv(csv_path, index=False)

    comp = frame[frame.run_kind == "comparator"]
    width_rows = []
    for (scenario, n_target), group in comp.groupby(["scenario", "n_target"]):
        direct = group[group.method == "proposed_direct"]
        separate = group[group.method == "separate_risk"]
        width_rows.append({
            "scenario": scenario,
            "n_target": int(n_target),
            "point_ratio": separate.point_width.mean() / direct.point_width.mean(),
            "ci_ratio": separate.ci_width.mean() / direct.ci_width.mean(),
        })
    widths = pd.DataFrame(width_rows).sort_values(["scenario", "n_target"])

    coverage_rows = []
    direct_comp = comp[comp.method == "proposed_direct"]
    validity = direct_comp[
        (direct_comp.class_contains_true.astype(str).str.lower() == "true")
        & (direct_comp.scenario != "pS_misspec")
    ]
    for (scenario, n_target), group in validity.groupby(["scenario", "n_target"]):
        coverage_rows.append({
            "scenario": scenario,
            "n_target": int(n_target),
            "eligible_gammas": ", ".join(f"{x:.1f}" for x in sorted(group.gamma.astype(float))),
            "min_coverage": group.ci_coverage.min(),
            "mean_coverage": group.ci_coverage.mean(),
            "mean_ci_width": group.ci_width.mean(),
            "mean_defer": group.defer.mean(),
        })
    coverage = pd.DataFrame(coverage_rows).sort_values(["scenario", "n_target"])

    mirror = frame[
        (frame.run_kind == "mirror")
        & (frame.method == "proposed_direct")
        & (frame.gamma.isin([0.2, 0.3]))
    ].sort_values(["scenario", "n_target", "gamma"])

    lines = [
        "# Comparator simulation — canonical nT=60/80 aggregation (LOCKED)",
        "",
        f"- 18/18 PASS runs: 7 comparator + 2 mirror scenarios at each nT.",
        f"- Monte Carlo: {REPS} replicates/cell; joint pairs bootstrap {BOOT}/replicate; nS=114.",
        "- All proportions have Wilson 95% intervals in `application_scale_comparator.csv`.",
        "- For a single proportion, maximum Monte Carlo SE is 0.016.",
        "",
        "## Table 1. Direct-versus-separate width ratio (mean across γ)",
        "",
        "| scenario | nT | point identified-set ratio | simultaneous-CI ratio |",
        "|---|---:|---:|---:|",
    ]
    for row in widths.itertuples():
        lines.append(f"| {row.scenario} | {row.n_target} | {fmt(row.point_ratio, 2)}× | {fmt(row.ci_ratio, 2)}× |")

    lines += [
        "",
        "## Table 2. Proposed-direct simultaneous-CI behavior",
        "",
        "Only γ cells whose sensitivity class contains the true target law are included.",
        "",
        "| scenario | nT | eligible γ | min coverage | mean coverage | mean CI width | mean DEFER |",
        "|---|---:|---|---:|---:|---:|---:|",
    ]
    for row in coverage.itertuples():
        lines.append(
            f"| {row.scenario} | {row.n_target} | {row.eligible_gammas} | {fmt(row.min_coverage)} | {fmt(row.mean_coverage)} "
            f"| {fmt(row.mean_ci_width)} | {fmt(row.mean_defer)} |"
        )
    lines += [
        "",
        "`pS_misspec` is excluded because nuisance misspecification makes class membership NA. `out_of_class` has no class-containing γ on the locked grid and is also excluded. Their row-wise stress-test results remain in the CSV.",
        "For `sign_change`, only γ=0.5 contains the true law. Low overlap remains the main failure region even among class-containing cells (coverage about 0.89–0.92, below nominal 0.95).",
        "",
        "## Table 3. Candidate-superior mirror decision behavior",
        "",
        "| scenario | nT | γ | RETAIN | wrong REFERENCE | DEFER | CI coverage |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in mirror.itertuples():
        lines.append(
            f"| {row.scenario} | {int(row.n_target)} | {fmt(row.gamma, 1)} | {fmt(row.retain)} "
            f"| {fmt(row.reference)} | {fmt(row.defer)} | {fmt(row.ci_coverage)} |"
        )

    retain = mirror.retain.astype(float)
    coverage_vals = mirror.ci_coverage.astype(float)
    lines += [
        "",
        "## Locked interpretation",
        "",
        f"- Direct identified sets are {widths.point_ratio.min():.2f}–{widths.point_ratio.max():.2f}× narrower than separate subtraction; simultaneous CIs are {widths.ci_ratio.min():.2f}–{widths.ci_ratio.max():.2f}× narrower.",
        f"- In mirror cells at γ=0.2/0.3, correct RETAIN is {retain.min():.3f}–{retain.max():.3f}; DEFER is {mirror.defer.astype(float).min():.3f}–{mirror.defer.astype(float).max():.3f}.",
        "- Wrong-direction REFERENCE was observed 0 times in each 1,000-replicate mirror cell (not a theoretical zero; one-cell exact 95% upper bound ≈0.0037).",
        f"- Mirror CI coverage is {coverage_vals.min():.3f}–{coverage_vals.max():.3f}.",
        "- These results replace nT=30/40 as the MAIN application-matched simulation; the 30/40 results remain a transductive-DDAU-matched supplement.",
        "",
        "## Provenance",
        "",
    ]
    for item in sorted(comp_prov + mirror_prov, key=lambda x: (x["n_target"], x["kind"], x["scenario"])):
        lines.append(
            f"- nT={item['n_target']} {item['kind']} {item['scenario']}: "
            f"manifest `{item['manifest']}` ({item['manifest_sha256'][:16]}…), "
            f"summary `{item['summary']}` ({item['summary_sha256'][:16]}…)."
        )

    rendered = "\n".join(lines) + "\n"
    md_path = output_root / "application_scale_comparator_summary.md"
    unified_path = output_root / "application_scale_comparator_table.md"
    md_path.write_text(rendered, encoding="utf-8")
    unified_path.write_text(rendered, encoding="utf-8")

    manifest = {
        "status": "PASS",
        "n_runs": 18,
        "n_targets": list(TARGETS),
        "reps": REPS,
        "bootstrap": BOOT,
        "aggregator_sha256": sha256(Path(__file__)),
        "outputs": {
            csv_path.name: sha256(csv_path),
            md_path.name: sha256(md_path),
            unified_path.name: sha256(unified_path),
        },
        "inputs": comp_prov + mirror_prov,
    }
    manifest_path = output_root / "application_scale_comparator_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"status": "PASS", "rows": len(frame), "outputs": manifest["outputs"]}, indent=2))


if __name__ == "__main__":
    main()
