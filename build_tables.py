"""Build manuscript-facing tables from committed aggregate results."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "tables"


def main():
    OUT.mkdir(exist_ok=True)
    comparator = pd.read_csv(ROOT / "results/simulations/application_scale_comparator.csv")
    table_1 = comparator[
        (comparator.run_kind == "comparator")
        & (comparator.n_target == 80)
        & np.isclose(comparator.gamma, 0.2)
        & comparator.scenario.isin(["no_shift", "xvary_in_class", "low_overlap", "pS_misspec"])
    ].copy()
    table_1.to_csv(OUT / "Table_1_application_scale_comparator.csv", index=False)

    nsclc = pd.read_csv(ROOT / "results/nsclc/outcome_free_evaluation.csv")
    nsclc[(nsclc.seed == 42) & np.isclose(nsclc.gamma, 0.2)].to_csv(
        OUT / "Table_3_NSCLC_outcome_free_evaluation.csv", index=False
    )

    supporting = pd.read_csv(ROOT / "results/applications/direct_vs_separate.csv")
    supporting[supporting.dataset == "Rotterdam-GBSG"].to_csv(
        OUT / "Table_S4_Rotterdam_GBSG_decisions.csv", index=False
    )

    coverage = pd.read_csv(ROOT / "results/simulations/primary_implementation_coverage.csv")
    coverage.to_csv(OUT / "Table_S9_primary_implementation_coverage.csv", index=False)

    regimes = []
    for path in sorted((ROOT / "results/simulations/decision_regime").glob("*.json")):
        regimes.append({"result": path.stem, **json.loads(path.read_text(encoding="utf-8"))})
    pd.json_normalize(regimes).to_csv(OUT / "Table_S5_decision_regimes.csv", index=False)
    print("[PASS] manuscript tables ->", OUT)


if __name__ == "__main__":
    main()
