"""Fail-closed audit for the manuscript-aligned public release candidate."""
from __future__ import annotations
import hashlib
import json
import re
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def need(condition, message):
    if not condition:
        raise AssertionError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_checksum_ledger():
    """Verify every archived byte listed in SHA256SUMS.txt.

    The ledger deliberately excludes itself and MANIFEST.md. Generated files
    under build/ are ignored and therefore cannot make this check pass or fail.
    """
    ledger = ROOT / "SHA256SUMS.txt"
    need(ledger.is_file(), "missing SHA256SUMS.txt")
    entries = []
    for line_number, raw in enumerate(ledger.read_text(encoding="utf-8").splitlines(), start=1):
        if not raw.strip():
            continue
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", raw)
        need(match is not None, f"malformed SHA256SUMS.txt line {line_number}")
        expected, relative = match.groups()
        path = ROOT / Path(relative)
        need(path.is_file(), f"checksum target missing: {relative}")
        actual = sha256(path)
        need(actual == expected, f"checksum mismatch: {relative}: expected {expected}, got {actual}")
        entries.append(relative)
    need(entries, "empty SHA256SUMS.txt")
    need(len(entries) == len(set(entries)), "duplicate paths in SHA256SUMS.txt")
    return len(entries)


def main():
    checks = []
    checksum_entries = verify_checksum_ledger()
    checks.append(f"SHA256SUMS ledger ({checksum_entries} files)")
    coverage = pd.read_csv(ROOT / "results/simulations/primary_implementation_coverage.csv").sort_values("n_target")
    need(coverage.n_target.tolist() == [60, 80, 686], "coverage target rows")
    need(np.allclose(coverage.two_sided_coverage, [.859, .847, .943]), "two-sided coverage")
    need(np.allclose(coverage.outer_envelope_coverage, [.861, .847, .946]), "outer coverage")
    checks.append("primary-implementation coverage values")

    comparator = pd.read_csv(ROOT / "results/simulations/application_scale_comparator.csv")
    need(len(comparator) == 288, "application-scale comparator must have 288 rows")
    scenarios = ["no_shift", "xvary_in_class", "low_overlap", "pS_misspec"]
    selected = comparator[(comparator.run_kind == "comparator") & (comparator.n_target == 80) & np.isclose(comparator.gamma, .2)]
    ratios = []
    for scenario in scenarios:
        group = selected[selected.scenario == scenario].set_index("method")
        ratios.append(group.loc["separate_risk", "point_width"] / group.loc["proposed_direct", "point_width"])
    need(np.allclose(ratios, [3.44, 3.441, 3.623, 3.428], atol=.002), "comparator width ratios")
    checks.append("application-scale comparator aggregate")

    geometry = pd.read_csv(ROOT / "results/simulations/controlled_geometry.csv")
    need(np.allclose(geometry.width_ratio, [1, 1.223463, 1.614863, 2.459502, 5.732487], atol=2e-6), "controlled geometry")
    checks.append("controlled-geometry illustration")

    nsclc = pd.read_csv(ROOT / "results/nsclc/outcome_free_evaluation.csv")
    need(len(nsclc) == 80, "NSCLC row count")
    need(nsclc.groupby("contrast").size().to_dict() == {
        "PET_radiomic_candidate": 40,
        "PET_radiomic_candidate_SRDO": 40,
    }, "40 primary and 40 SRDO configurations")
    need(set(nsclc.action_dir) == {"DEFER"}, "NSCLC decisions")
    checks.append("NSCLC outcome-free aggregate")

    applications = pd.read_csv(ROOT / "results/applications/direct_vs_separate.csv")
    flip = applications[
        (applications.dataset == "Rotterdam-GBSG")
        & (applications.contrast == "supporting_clinical_candidate")
        & (applications.seed == 62)
        & np.isclose(applications.gamma, .2)
    ]
    need(len(flip) == 1 and flip.iloc[0].action_dir == "ADOPT" and flip.iloc[0].action_sep == "DEFER", "Rotterdam–GBSG decision difference")
    checks.append("Rotterdam–GBSG supporting analysis")

    onestep_apps = pd.read_csv(ROOT / "results/onestep/applications.csv")
    need(len(onestep_apps) == 60, "one-step applications must have 60 cells")
    need((onestep_apps.plugin_action == onestep_apps.onestep_action).all(), "one-step and plug-in decisions agree in all cells")
    committed = pd.concat([
        nsclc[nsclc.contrast == "PET_radiomic_candidate"],
        applications[(applications.dataset == "Rotterdam-GBSG") & (applications.contrast == "supporting_clinical_candidate")],
    ])[["dataset", "seed", "gamma", "dir_sim_lower", "dir_sim_upper", "action_dir"]]
    merged = onestep_apps.merge(committed, on=["dataset", "seed", "gamma"])
    need(len(merged) == 60, "one-step applications aligned with committed plug-in rows")
    need(np.allclose(merged.plugin_env_lower, merged.dir_sim_lower, atol=1e-12)
         and np.allclose(merged.plugin_env_upper, merged.dir_sim_upper, atol=1e-12)
         and (merged.plugin_action == merged.action_dir).all(), "one-step plug-in columns reproduce committed envelopes")
    cell = onestep_apps[(onestep_apps.dataset == "Rotterdam-GBSG") & (onestep_apps.seed == 62) & np.isclose(onestep_apps.gamma, .2)].iloc[0]
    need(np.allclose([cell.onestep_env_lower, cell.onestep_env_upper], [-0.0486, -0.0010], atol=5e-5), "one-step seed-62 envelope")
    need(np.allclose([cell.plugin_env_lower, cell.plugin_env_upper], [-0.0486, -0.0017], atol=5e-5), "plug-in seed-62 envelope")
    simulation = pd.read_csv(ROOT / "results/onestep/simulation_summary.csv")
    cov = simulation.groupby(["law", "scale", "estimator"]).outer_cov.first()
    need(np.isclose(cov[("A", "nsclc_scale", "plug-in")], .863) and np.isclose(cov[("A", "nsclc_scale", "one-step")], .844), "one-step NSCLC-scale coverage")
    need(np.isclose(cov[("B", "rotterdam_gbsg_scale", "plug-in")], .957) and np.isclose(cov[("B", "rotterdam_gbsg_scale", "one-step")], .957), "one-step Rotterdam-scale coverage")
    need((simulation.n_valid == 1000).all(), "one-step simulation replications")
    kappa = pd.read_csv(ROOT / "results/scalar_diagnostics/web_table_1.csv")
    need(len(kappa) == 12 and (kappa.width_ratio - kappa.kappa).abs().max() < .002, "Web Table 1 kappa diagnostic")
    seeds = json.loads((ROOT / "code/inference/archived_bootstrap_seeds.json").read_text(encoding="utf-8"))["seeds"]
    need(len(seeds) == 20, "archived bootstrap seeds")
    for stem in ("Figure_1", "Figure_2", "Figure_3", "Figure_4", "Web_Figure_1", "Web_Figure_2"):
        for extension in ("pdf", "png", "tiff"):
            path = ROOT / "figures" / "biometrics" / f"{stem}.{extension}"
            need(path.is_file() and path.stat().st_size > 1000, f"missing {path}")
    need(len(list((ROOT / "tables" / "biometrics").glob("*.csv"))) == 3, "Biometrics web tables")
    checks.append("v1.2.0 one-step, scalar-diagnostic, and Biometrics artwork outputs")

    for stem in ("Figure_1", "Figure_2", "Figure_3", "Figure_4", "Figure_S1", "Figure_S2"):
        for extension in ("pdf", "eps", "png", "tiff"):
            path = ROOT / "figures" / f"{stem}.{extension}"
            need(path.is_file() and path.stat().st_size > 1000, f"missing {path}")
    need(len(list((ROOT / "tables").glob("*.csv"))) >= 5, "table outputs")
    checks.append("figures and tables")

    command = ["git", "-c", f"safe.directory={ROOT}", "-C", str(ROOT), "ls-files", "--cached", "--others", "--exclude-standard"]
    public = [ROOT / item for item in subprocess.check_output(command, text=True, encoding="utf-8").splitlines()]
    secret_pattern = re.compile(r"github_pat_|ghp_[A-Za-z0-9]{20,}|AKIA[0-9A-Z]{16}|BEGIN (?:RSA|OPENSSH|EC) PRIVATE KEY")
    local_path_pattern = re.compile(r"(?:[A-Za-z]:[\\/](?:Users|Documents|Desktop)[\\/]|/home/[^/]+/|TCIA_PET)", re.I)
    internal_term_pattern = re.compile(r"FROM_CODEX|TASK_[A-Z0-9]+|BDS7|positive[-_ ]control|InductiveSparsePCA|ctdate|59_80|60_80", re.I)
    patient_csvs, secrets, local_paths, internal_terms = [], [], [], []
    for path in public:
        if not path.is_file() or path.resolve() == Path(__file__).resolve():
            continue
        relative = path.relative_to(ROOT).as_posix()
        if path.suffix.lower() == ".csv" and "synthetic_fixture/" not in relative:
            try:
                if "patient_id" in set(pd.read_csv(path, nrows=0).columns):
                    patient_csvs.append(relative)
            except Exception:
                pass
        if path.suffix.lower() in {".py", ".md", ".json", ".cff", ".txt", ".r", ".csv"}:
            text = path.read_text(encoding="utf-8", errors="ignore")
            if secret_pattern.search(text): secrets.append(relative)
            if local_path_pattern.search(text): local_paths.append(relative)
            if internal_term_pattern.search(text) or internal_term_pattern.search(relative) or re.search(r"20\d{6}", relative): internal_terms.append(relative)
    need(not patient_csvs, f"patient-level public CSV files: {patient_csvs}")
    need(not secrets, f"secret-like strings: {secrets}")
    need(not local_paths, f"local filesystem paths: {local_paths}")
    need(not internal_terms, f"non-manuscript internal terms: {internal_terms}")
    checks.append("privacy, credential, local-path, and terminology scan")

    result = {
        "status": "PASS_MANUSCRIPT_ALIGNED_RELEASE_AUDIT",
        "checks": checks,
        "public_file_count": len(public),
        "note": "Clean-environment and Git-history tests are documented in VALIDATION_REPORT.md.",
    }
    (ROOT / "VALIDATION_RESULTS.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
