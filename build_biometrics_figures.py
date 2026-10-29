"""Regenerate the Biometrics-numbered figures that changed in v1.2.0 from committed aggregates.

  Figure 4      Rotterdam-GBSG seed 62, gamma = 0.2: direct and separate identified sets and
                simultaneous envelopes  <- results/applications/direct_vs_separate.csv
  Web Figure 2  decision frequencies across prespecified simulation regimes
                <- results/simulations/decision_regime/*.json

Decision labels follow the Biometrics manuscript: CANDIDATE SUPPORTED, REFERENCE SUPPORTED,
DEFER (the committed CSV/JSON files keep the internal codes ADOPT/KEEP/DEFER and
retain/reference/defer).  Output: build/figures_biometrics/ (override with BIOMETRICS_FIGDIR).
Figures 1-3 and Web Figure 1 are listed in BIOMETRICS_NUMBERING.md.
"""
import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402

BASE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.environ.get("BIOMETRICS_FIGDIR", os.path.join(BASE, "build", "figures_biometrics"))
os.makedirs(FIGDIR, exist_ok=True)
DOUBLE = 180 / 25.4

for _font in ("Arial", "Helvetica", "Liberation Sans"):
    try:
        _path = font_manager.findfont(_font, fallback_to_default=False)
    except ValueError:
        continue
    font_manager.fontManager.addfont(_path); FONT = _font; break
else:
    raise RuntimeError("Install Arial, Helvetica, or Liberation Sans.")

plt.rcParams.update({"font.size": 8.5, "font.family": FONT, "font.sans-serif": [FONT],
                     "mathtext.fontset": "stixsans", "axes.titlesize": 9, "axes.labelsize": 8.5,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
                     "axes.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42,
                     "ps.useafm": True, "svg.fonttype": "none"})
LABEL = {"ADOPT": "CANDIDATE SUPPORTED", "KEEP_REFERENCE": "REFERENCE SUPPORTED",
         "KEEP": "REFERENCE SUPPORTED", "DEFER": "DEFER"}


def save(fig, name):
    for ext in ("pdf", "eps"):
        fig.savefig(os.path.join(FIGDIR, f"{name}.{ext}"), bbox_inches="tight")
    fig.savefig(os.path.join(FIGDIR, f"{name}.tiff"), dpi=800, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(os.path.join(FIGDIR, f"{name}.png"), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig); print("saved", name)


def figure_4():
    src = os.path.join(BASE, "results", "applications", "direct_vs_separate.csv")
    row = next(r for r in csv.DictReader(open(src)) if r["dataset"] == "Rotterdam-GBSG"
               and r["contrast"] == "supporting_clinical_candidate" and int(r["seed"]) == 62
               and abs(float(r["gamma"]) - 0.2) < 1e-9)
    g = lambda c: float(row[c])
    idr, ids = (g("L_dir"), g("U_dir")), (g("L_sep"), g("U_sep"))
    denv, senv = (g("dir_sim_lower"), g("dir_sim_upper")), (g("sep_sim_lower"), g("sep_sim_upper"))
    c_dir, c_sep = "#1f4e79", "#c0504d"
    fig, ax = plt.subplots(1, 1, figsize=(DOUBLE, 2.25)); yD, yS = 0.68, 0.32
    for env, iset, y, col in ((senv, ids, yS, c_sep), (denv, idr, yD, c_dir)):
        ax.plot(env, [y, y], color=col, lw=2.2, solid_capstyle="butt", zorder=2)
        for x in env:
            ax.plot([x], [y], "|", color=col, ms=10, mew=1.8, zorder=3)
        ax.plot(iset, [y, y], color=col, lw=7, solid_capstyle="round", zorder=4)
    ax.axvline(0, color="0.35", ls="--", lw=1.1, zorder=1)
    ax.set_xlim(-0.065, 0.015); ax.set_ylim(0, 1)
    ax.set_yticks([yS, yD]); ax.set_yticklabels(["separate", "direct"], fontsize=8)
    ax.tick_params(axis="x", labelsize=10.5)
    ax.set_title(f"Rotterdam → GBSG, seed 62, $\\gamma$ = 0.2    |    ID width ratio {g('ID_width_ratio'):.2f}×    |    "
                 f"direct {LABEL[row['action_dir'].strip()]} / separate {LABEL[row['action_sep'].strip()]}",
                 fontsize=8.5, loc="left")
    ax.tick_params(axis="y", length=0); ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_xlabel(r"Target Brier-risk contrast  $\Delta R_T$  (candidate $-$ reference)", fontsize=9)
    fig.tight_layout(); save(fig, "Figure_4")


def web_figure_2():
    bdir = os.path.join(BASE, "results", "simulations", "decision_regime")
    regimes = ["candidate_superior", "finite_crossing", "no_crossing_within_range", "reference_superior", "ambiguous"]
    labels = ["Candidate superior", "Finite crossing", "No crossing within range", "Reference superior", "Ambiguous"]
    rates = lambda r, c: [json.load(open(os.path.join(bdir, f"{r}_{c}.json")))["decision_rates"][k]
                          for k in ("retain", "reference", "defer")]
    cats = ["CANDIDATE SUPPORTED", "REFERENCE SUPPORTED", "DEFER"]; colors = ["#2e7d32", "#2d6cb5", "#b8bdc6"]
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 3.9), sharey=True); x = np.arange(len(regimes))
    for ax, cens, title in zip(axes, ["marginal", "conditional"],
                               ["(A) Independent censoring", "(B) Covariate-dependent censoring"]):
        M = np.array([rates(r, cens) for r in regimes]); bottom = np.zeros(len(regimes))
        for k in range(3):
            ax.bar(x, M[:, k], bottom=bottom, color=colors[k], width=0.62, edgecolor="white", linewidth=0.6, label=cats[k])
            for xi, (v, b) in enumerate(zip(M[:, k], bottom)):
                if v >= 0.06:
                    ax.text(xi, b + v / 2, f"{v:.2f}", ha="center", va="center", fontsize=7,
                            color="white" if k != 2 else "#333333")
            bottom += M[:, k]
        ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7, rotation=30, ha="right")
        ax.tick_params(axis="y", labelsize=8); ax.set_ylim(0, 1.0)
        ax.set_title(title, fontsize=9, color="#1f3b73", fontweight="bold")
        ax.grid(axis="y", color="0.9", linewidth=0.8); ax.set_axisbelow(True)
    axes[0].set_ylabel("Monte Carlo decision frequency", fontsize=8.5)
    handles, legend_labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.015), frameon=False, fontsize=7.5)
    fig.tight_layout(rect=[0, 0.16, 1, 1]); save(fig, "Web_Figure_2")


if __name__ == "__main__":
    figure_4()
    web_figure_2()
