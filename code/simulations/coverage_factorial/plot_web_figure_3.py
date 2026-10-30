"""Web Figure 3: outer-envelope coverage in the 15-cell coverage factorial (Web Table 8).
Reads results/simulations/coverage_factorial/factorial_coverage_summary.csv and writes
figures/biometrics/Web_Figure_3.{pdf,png,tiff} (or to $FIGDIR if set)."""
import os
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
FIGDIR = Path(os.environ.get("FIGDIR", ROOT / "figures" / "biometrics")); FIGDIR.mkdir(parents=True, exist_ok=True)
FONT = "DejaVu Sans"
for f in ("Arial", "Helvetica", "Liberation Sans"):
    try:
        font_manager.fontManager.addfont(font_manager.findfont(f, fallback_to_default=False)); FONT = f; break
    except Exception:
        pass
plt.rcParams.update({"font.size": 8.5, "font.family": FONT, "axes.linewidth": 0.8, "pdf.fonttype": 42,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8})
d = pd.read_csv(ROOT / "results/simulations/coverage_factorial/factorial_coverage_summary.csv")
STYLE = {0.65: ("#1f3b73", "o", "65% censoring", -6), 0.50: ("#2e8b57", "s", "50% censoring", 0),
         0.38: ("#c0504d", "^", "38% censoring", 6)}
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), sharey=True)
for cens, (col, mk, lab, dx) in STYLE.items():
    s = d[abs(d.target_cens - cens) < 1e-9].sort_values("nS")
    err = [s.outer_envelope_coverage - s.wilson_lo, s.wilson_hi - s.outer_envelope_coverage]
    axes[0].errorbar(s.nS + dx, s.outer_envelope_coverage, yerr=err, color=col, marker=mk, ms=4, lw=1.1, capsize=2, label=lab)
    axes[1].errorbar(s.mean_source_events, s.outer_envelope_coverage, yerr=err, color=col, marker=mk, ms=4, lw=1.1, capsize=2, label=lab)
for ax in axes:
    ax.axhline(0.95, ls="--", lw=0.9, color="0.35"); ax.set_ylim(0.82, 0.975)
    ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="y", lw=0.4, color="0.88")
axes[0].set_xticks([114, 200, 300, 400, 516]); axes[0].set_xlabel("Source calibration size $n_S$")
axes[1].set_xlabel("Mean number of source events"); axes[0].set_ylabel("Outer-envelope coverage")
axes[0].set_title("(A)", loc="left", fontsize=9, fontweight="bold"); axes[1].set_title("(B)", loc="left", fontsize=9, fontweight="bold")
axes[0].text(118, 0.953, "nominal 0.95", fontsize=7, color="0.35")
axes[0].legend(frameon=False, loc="lower right")
fig.tight_layout()
for ext in ("pdf", "png", "tiff"):
    kw = {"dpi": 600} if ext == "png" else ({"dpi": 800, "pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {})
    fig.savefig(FIGDIR / f"Web_Figure_3.{ext}", bbox_inches="tight", facecolor="white", **kw)
print("wrote", FIGDIR / "Web_Figure_3.*")
