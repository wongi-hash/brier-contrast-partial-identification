# %% [markdown]
# # All manuscript figures — single editable notebook (SiM submission format)
#
# Percent-format notebook. In VSCode or Jupyter each `# %%` block is a cell.
# (To get a classic .ipynb:  `jupytext --to notebook ALL_FIGURES.py`)
#
# **Format (per SiM / Wiley artwork guidelines):**
# - Each figure saved as its own **vector PDF + EPS**, **800-dpi LZW TIFF**
#   (submission), and **600-dpi PNG** (preview).
# - Figures are designed directly at the Wiley 180-mm full-page width rather than
#   being drawn oversized and reduced later.
# - Arial/Helvetica/Liberation Sans is required explicitly (no silent fallback).
# - Panel labels unified to (A), (B), (C). No pale background tints ("Tints are not
#   acceptable"): white fills + colored outlines only.
# - direct = blue, separate = red throughout; distinctions also carried by row/label
#   so they survive greyscale / color-vision deficiency.
#
# **Manuscript ↔ code mapping / data sources**
# | Manuscript | Content | Data source |
# |---|---|---|
# | Figure 1 | workflow schematic | none (schematic) |
# | Figure 2 | Prop 2 geometry | computed from closed forms |
# | Figure 3 | simulation summary | controlled geometry and application-scale comparator aggregates |
# | Figure 4 | NSCLC direct versus separate | final outcome-free evaluation aggregate (seed 42, γ=0.2) |
# | Figure S1 | Rotterdam→GBSG decisions | supporting-analysis aggregate |
# | Figure S2 | decision-regime frequencies | prespecified regime aggregates |

# %%
# ---- shared setup (run this cell first) ----
import os, csv, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless, reproducible build on CI/servers
import matplotlib.pyplot as plt
from matplotlib import font_manager

try:
    BASE = os.path.dirname(os.path.abspath(__file__))
except NameError:
    BASE = os.environ.get("SIM_REPRO_ROOT", os.getcwd())

FIGDIR = os.path.join(BASE, "figures")
os.makedirs(FIGDIR, exist_ok=True)

MM = 1 / 25.4
DOUBLE = 180 * MM

_FONT_CANDIDATES = ("Arial", "Helvetica", "Liberation Sans")
SELECTED_FONT = None
SELECTED_FONT_PATH = None
for _font in _FONT_CANDIDATES:
    try:
        _font_path = font_manager.findfont(_font, fallback_to_default=False)
    except ValueError:
        continue
    if os.path.isfile(_font_path):
        font_manager.fontManager.addfont(_font_path)
        SELECTED_FONT, SELECTED_FONT_PATH = _font, _font_path
        break
if SELECTED_FONT is None:
    raise RuntimeError(
        "No approved figure font found. Install Arial, Helvetica, or Liberation Sans; "
        "silent DejaVu fallback is disabled."
    )

plt.rcParams.update({
    "font.size": 9.5,
    "font.family": SELECTED_FONT,
    "font.sans-serif": [SELECTED_FONT],
    "mathtext.fontset": "stixsans",
    "axes.titlesize": 10,
    "axes.labelsize": 9.5,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,   # embed TrueType (editable text) in PDF
    "ps.fonttype": 42,
    # Matplotlib on this Windows host cannot embed Arial as PostScript Type 42.
    # EPS therefore uses the standard Helvetica AFM face explicitly; TIFF/PDF
    # retain Arial. This is an explicit compatibility choice, not silent fallback.
    "ps.useafm": True,
    "svg.fonttype": "none",
})


def expit(x):
    return 1.0 / (1.0 + np.exp(-x))


def logit(p):
    return np.log(p / (1.0 - p))


def wilson(p, n=1000, z=1.959963984540054):
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


def save(fig, name, preview_dpi=600, submission_dpi=800):
    """Vector masters, 800-dpi LZW TIFF submission file, and PNG preview."""
    for ext in ("pdf", "eps"):
        fig.savefig(os.path.join(FIGDIR, f"{name}.{ext}"), bbox_inches="tight")
    fig.savefig(
        os.path.join(FIGDIR, f"{name}.tiff"), dpi=submission_dpi,
        bbox_inches="tight", facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    fig.savefig(
        os.path.join(FIGDIR, f"{name}.png"), dpi=preview_dpi,
        bbox_inches="tight", facecolor="white",
    )
    print("saved", name, "(pdf/eps/tiff/png)")


print("BASE =", BASE)
print("FIGURE FONT =", SELECTED_FONT, "|", SELECTED_FONT_PATH)

# %% [markdown]
# ## Figure 1 — Overview (4-step, illustrated)
# Step 1 uses hospital icons: a source site's locked models q0, q1 are compared at
# target sites with no observed outcomes. (1) problem, (2) sensitivity assumption,
# (3) direct sharp identification (core; right-censored source handled as a note),
# (4) simultaneous inference → decision (on the confidence envelope).

# %%
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle

NAVY = "#1f3b6e"; EDGE = "#1f3b6e"; FILL = "#ffffff"; GREY = "#5a6270"; TXT = "#111111"
GREEN = "#2e9e4f"; BLUE = "#2d6cb5"; RED = "#c0392b"; HOSP = "#ffffff"

fig, ax = plt.subplots(figsize=(DOUBLE, 9.2))
ax.set_xlim(0, 10); ax.set_ylim(0, 18); ax.axis("off")


def box(x, y, w, h, ec=EDGE, lw=2.2, dashed=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.22",
                 linewidth=lw, edgecolor=ec, facecolor=FILL, linestyle=("--" if dashed else "-"), joinstyle="round"))


def badge(x, y, n):
    ax.add_patch(Circle((x, y), 0.44, facecolor=NAVY, edgecolor="white", linewidth=1.6, zorder=8))
    ax.text(x, y, str(n), ha="center", va="center", color="white", fontsize=10, fontweight="bold", zorder=9)


def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=20, lw=2.0,
                 color=NAVY, shrinkA=1, shrinkB=1, zorder=1))


def icon(x, y, color, glyph):
    ax.add_patch(Circle((x, y), 0.26, facecolor=color, edgecolor="none", zorder=6))
    ax.text(x, y, glyph, ha="center", va="center", color="white", fontsize=8.5, fontweight="bold", zorder=7)


def hospital(cx, cy, s, label=None):
    ax.add_patch(FancyBboxPatch((cx - 1.0 * s, cy - 0.75 * s), 2.0 * s, 1.5 * s,
                 boxstyle="round,pad=0.01,rounding_size=0.06", fc=HOSP, ec=NAVY, lw=1.6, zorder=3))
    ax.add_patch(Rectangle((cx - 0.11 * s, cy + 0.02 * s), 0.22 * s, 0.55 * s, color=RED, zorder=5))
    ax.add_patch(Rectangle((cx - 0.28 * s, cy + 0.20 * s), 0.56 * s, 0.19 * s, color=RED, zorder=5))
    ax.add_patch(Rectangle((cx - 0.16 * s, cy - 0.75 * s), 0.32 * s, 0.45 * s, color="#8aa6cc", zorder=5))
    if label:
        ax.text(cx, cy - 0.98 * s, label, ha="center", va="top", fontsize=8.5, color=TXT, fontweight="bold")


def chip(x, y, txt):
    ax.add_patch(FancyBboxPatch((x - 0.42, y - 0.26), 0.84, 0.52, boxstyle="round,pad=0.02,rounding_size=0.12",
                 fc=FILL, ec=NAVY, lw=1.4, zorder=6))
    ax.text(x, y, txt, ha="center", va="center", fontsize=9, color=NAVY, fontweight="bold", zorder=7)


# 1 Problem setup (illustrated)
box(0.5, 13.3, 9.0, 4.3)
ax.text(5.0, 17.42, "Problem setup", ha="center", va="top", fontsize=10, fontweight="bold", color=NAVY)
hospital(2.25, 15.75, 0.92, "Source hospital A")
chip(1.6, 14.3, r"$q_0$"); chip(2.9, 14.3, r"$q_1$")
arrow(3.7, 15.6, 5.75, 15.6)
ax.text(4.72, 16.02, "deploy", ha="center", va="center", fontsize=8.5, color=GREY, style="italic")
hospital(6.95, 15.95, 0.62); hospital(8.05, 15.2, 0.62)
ax.text(7.5, 16.85, "?", ha="center", va="center", fontsize=14, color=NAVY, fontweight="bold")
ax.text(7.5, 14.5, "target hospitals", ha="center", va="top", fontsize=8.5, color=TXT, fontweight="bold")
ax.text(7.5, 14.12, r"($Y_\tau$ not observed)", ha="center", va="top", fontsize=8.3, color=RED)
ax.text(5.0, 13.5, r"Which of $q_0, q_1$ is better?      $\Delta R_T = R_T(q_1) - R_T(q_0)$", ha="center", va="center", fontsize=9, color=TXT)
badge(0.85, 17.5, 1)

# 2 Sensitivity assumption
box(1.5, 10.8, 7.0, 1.85)
ax.text(5.0, 12.42, "Sensitivity assumption", ha="center", va="top", fontsize=9.5, fontweight="bold", color=NAVY)
ax.text(5.0, 11.35, "between-site conditional shift bounded:   |δ(z)| ≤ γ", ha="center", va="center", fontsize=9, color=TXT)
badge(1.85, 12.55, 2)

# 3 Direct sharp identification (center, emphasized)
box(0.8, 4.7, 8.4, 4.6, lw=3.2)
ax.text(5.0, 9.05, "Direct sharp identification", ha="center", va="top", fontsize=10, fontweight="bold", color=NAVY)
ax.text(5.0, 8.55, "preserving the shared target law", ha="center", va="top", fontsize=9, fontweight="bold", color=NAVY)
ax.text(5.0, 7.85, r"identified set for $\Delta R_T$    ·    right-censored source (IPCW + Cox)", ha="center", va="center", fontsize=8.3, color=GREY, style="italic")
sep_y, dir_y = 7.0, 6.4
ax.plot([3.4, 7.9], [sep_y, sep_y], color="#7d9bc4", lw=3.4, solid_capstyle="butt", zorder=5)
for xx in (3.4, 7.9): ax.plot([xx, xx], [sep_y - 0.11, sep_y + 0.11], color="#7d9bc4", lw=3.4, zorder=5)
ax.text(3.15, sep_y, "Separate", ha="right", va="center", fontsize=8.5, color="#7d9bc4", fontweight="bold")
ax.plot([4.7, 6.6], [dir_y, dir_y], color=NAVY, lw=3.4, solid_capstyle="butt", zorder=5)
for xx in (4.7, 6.6): ax.plot([xx, xx], [dir_y - 0.11, dir_y + 0.11], color=NAVY, lw=3.4, zorder=5)
ax.text(3.15, dir_y, "Direct", ha="right", va="center", fontsize=8.5, color=NAVY, fontweight="bold")
ax.text(5.65, 5.72, r"$W_{\mathrm{dir}}$  ≤  $W_{\mathrm{sep}}$", ha="center", va="center", fontsize=9, color=TXT)
ax.text(5.0, 5.22, r"$\rightarrow$ sharper identification of the better model", ha="center", va="center", fontsize=8.5, color=NAVY)
badge(1.15, 9.2, 3)

# 4 Simultaneous inference -> decision
box(1.1, 0.4, 7.8, 3.6)
ax.text(5.0, 3.82, r"Simultaneous inference $\rightarrow$ decision", ha="center", va="top", fontsize=9.5, fontweight="bold", color=NAVY)
ax.text(5.0, 3.18, r"pairs bootstrap envelope $[\,L_\gamma,\ \bar{U}_\gamma\,]$", ha="center", va="center", fontsize=9, color=TXT)
icon(2.2, 2.5, GREEN, "A"); ax.text(2.62, 2.5, "upper envelope < 0", ha="left", va="center", fontsize=8.5, color=TXT); ax.text(5.6, 2.5, r"$\rightarrow$ ADOPT", ha="left", va="center", fontsize=8.5, color=GREEN, fontweight="bold")
icon(2.2, 1.9, BLUE, "K"); ax.text(2.62, 1.9, "lower envelope > 0", ha="left", va="center", fontsize=8.5, color=TXT); ax.text(5.6, 1.9, r"$\rightarrow$ KEEP", ha="left", va="center", fontsize=8.5, color=BLUE, fontweight="bold")
icon(2.2, 1.3, RED, "D"); ax.text(2.62, 1.3, "otherwise", ha="left", va="center", fontsize=8.5, color=TXT); ax.text(5.6, 1.3, r"$\rightarrow$ DEFER", ha="left", va="center", fontsize=8.5, color=RED, fontweight="bold")
ax.text(5.0, 0.72, r"report robustness value $\gamma^{*}$", ha="center", va="center", fontsize=8.3, color=GREY, style="italic")
badge(1.45, 3.9, 4)

arrow(5.0, 13.3, 5.0, 12.65)
arrow(5.0, 10.8, 5.0, 9.3)
arrow(5.0, 4.7, 5.0, 4.0)

plt.tight_layout()
save(fig, "Figure_1")
plt.show()

# %% [markdown]
# ## Figure 2 — Direct vs. separate identified sets (Prop 2 geometry)
# Section 2.4. Three point-mass profiles at γ=0.2 (computed from closed forms).
# Direct = navy, separate = red; dashed line = 0. Panels (A)(B)(C).

# %%
GAMMA = 0.2
panels = [("(A) Same side of 1/2", 0.60, 0.75, 0.65),
          ("(B) Straddling 1/2", 0.40, 0.70, 0.55),
          ("(C) Identical predictors", 0.70, 0.70, 0.65)]


def _compute(q0, q1, mS, g=GAMMA):
    mp = expit(logit(mS) + g); mm = expit(logit(mS) - g)
    a = q1**2 - q0**2; b = q1 - q0
    d1 = a - 2 * b * mm; d2 = a - 2 * b * mp
    Ldir, Udir = min(d1, d2), max(d1, d2)

    def Rr(qj):
        r1 = qj**2 + (1 - 2 * qj) * mm; r2 = qj**2 + (1 - 2 * qj) * mp
        return min(r1, r2), max(r1, r2)

    L1, U1 = Rr(q1); L0, U0 = Rr(q0)
    Lsep, Usep = L1 - U0, U1 - L0
    return dict(Ldir=Ldir, Udir=Udir, Lsep=Lsep, Usep=Usep, Wdir=Udir - Ldir, Wsep=Usep - Lsep)


res = [_compute(q0, q1, mS) for _, q0, q1, mS in panels]
NAVY2 = "#1f4e79"; RED = "#c0504d"
fig, axes = plt.subplots(3, 1, figsize=(DOUBLE, 6.3))
for ax, (t, q0, q1, mS), r in zip(axes, panels, res):
    ax.plot([r["Ldir"], r["Udir"]], [0.70, 0.70], color=NAVY2, lw=7, solid_capstyle="round", zorder=3)
    ax.plot([r["Lsep"], r["Usep"]], [0.30, 0.30], color=RED, lw=7, solid_capstyle="round", zorder=3)
    for x in (r["Ldir"], r["Udir"]): ax.plot([x], [0.70], "|", color=NAVY2, ms=13, mew=2.2, zorder=4)
    for x in (r["Lsep"], r["Usep"]): ax.plot([x], [0.30], "|", color=RED, ms=13, mew=2.2, zorder=4)
    if r["Wdir"] < 1e-12: ax.plot([r["Ldir"]], [0.70], "o", color=NAVY2, ms=9, zorder=5)
    ax.axvline(0.0, color="0.35", ls="--", lw=1.1, zorder=1)
    lo = min(r["Ldir"], r["Lsep"]); hi = max(r["Udir"], r["Usep"]); pad = 0.06 * max(hi - lo, 0.02)
    ax.set_xlim(lo - pad - 0.01, hi + pad + 0.01); ax.set_ylim(0, 1)
    ax.set_yticks([0.30, 0.70]); ax.set_yticklabels(["separate", "direct"], fontsize=10)
    ax.tick_params(axis="x", labelsize=10)
    rtxt = r"width ratio undefined ($W_{\mathrm{dir}}=0$)" if r["Wdir"] < 1e-12 else f"$W_{{\\mathrm{{sep}}}}/W_{{\\mathrm{{dir}}}}$ = {r['Wsep']/r['Wdir']:.2f}"
    ax.set_title(f"{t}    |    $q_0$={q0}, $q_1$={q1}    |    {rtxt}", fontsize=10, loc="left")
    ax.spines[["top", "right", "left"]].set_visible(False); ax.tick_params(axis="y", length=0)
axes[-1].set_xlabel(r"Target Brier-risk contrast  $\Delta R_T$", fontsize=10)
fig.tight_layout()
save(fig, "Figure_2")
plt.show()

# %% [markdown]
# ## Figure 3 — Simulation summary
# Section 3.5. (A) analytic geometry sweep (γ=0.3, stated constants). (B),(C) read
# from COMPARATOR_6080_AGGREGATED.csv at n_T=80, γ=0.2. (C) = per-γ true-contrast
# containment (NOT endpoint coverage; no nominal-calibration reference line is shown).

# %%
COMP_CSV = os.path.join(BASE, "results", "simulations", "application_scale_comparator.csv")
SCN = ["no_shift", "xvary_in_class", "low_overlap", "pS_misspec"]
SCN_LABEL = ["no\nshift", "Z-varying", "low\noverlap", r"$m_S$" + "\nmisspec."]
GAMMA_FIG3 = 0.2

_rows = [r for r in csv.DictReader(open(COMP_CSV))
         if r["run_kind"] == "comparator" and int(float(r["n_target"])) == 80
         and abs(float(r["gamma"]) - GAMMA_FIG3) < 1e-9 and r["scenario"] in SCN]


def _cell(scn, method, col):
    for r in _rows:
        if r["scenario"] == scn and r["method"] == method:
            return float(r[col])
    raise ValueError(f"Fig3 missing {scn}/{method}/{col} at gamma={GAMMA_FIG3}, n_T=80")


idr = [_cell(s, "separate_risk", "point_width") / _cell(s, "proposed_direct", "point_width") for s in SCN]
cir = [_cell(s, "separate_risk", "ci_width") / _cell(s, "proposed_direct", "ci_width") for s in SCN]
cov = [_cell(s, "proposed_direct", "ci_coverage") for s in SCN]
ci_lo = [_cell(s, "proposed_direct", "ci_coverage_wilson_low") for s in SCN]
ci_hi = [_cell(s, "proposed_direct", "ci_coverage_wilson_high") for s in SCN]
print("Fig3 B/C read from:", COMP_CSV, "| n_T=80, gamma=", GAMMA_FIG3)
print("  idr:", [round(v, 3) for v in idr], "| cir:", [round(v, 3) for v in cir], "| containment:", cov)

NAVY3 = "#1f4e79"; STEEL3 = "#8fb3d9"; RED3 = "#c0504d"
# Preserve one compound manuscript figure while separating the conceptual and
# empirical layers: panel A spans the first row; B and C share the second row.
fig = plt.figure(figsize=(DOUBLE, 6.1))
gs = fig.add_gridspec(2, 2, height_ratios=[0.92, 1.08], hspace=0.62, wspace=0.38)
axes = [fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
TT = 10

# (A) controlled geometry sweep (gamma=0.3), generated by code/GENERATE_CONTROLLED_GEOMETRY.py
GEOM_CSV = os.path.join(BASE, "results", "simulations", "controlled_geometry.csv")
with open(GEOM_CSV, newline="", encoding="utf-8-sig") as _fh:
    _geom = list(csv.DictReader(_fh))
f = [float(r["same_side_fraction"]) for r in _geom]
ratio = [float(r["width_ratio"]) for r in _geom]
print("Fig3 A read from:", GEOM_CSV, "| ratios:", [round(v, 2) for v in ratio])
axes[0].plot(f, ratio, "o-", color=NAVY3, lw=2.2, ms=7); axes[0].axhline(1.0, color="0.6", ls="--", lw=1)
axes[0].set_xlabel("Proportion on same side of 1/2", fontsize=9.5); axes[0].set_ylabel(r"$W_{\mathrm{sep}}/W_{\mathrm{dir}}$", fontsize=9.5)
axes[0].set_title(r"(A) Geometry ($\gamma=0.3$)", loc="left", fontsize=TT)
axes[0].set_xlim(-0.05, 1.05); axes[0].set_ylim(0.75, 6.15)
for xi, yi in zip(f, ratio):
    axes[0].annotate(f"{yi:.2f}", (xi, yi), textcoords="offset points",
                     xytext=(0, 10), ha="center", va="bottom", fontsize=8.5,
                     bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none"},
                     zorder=5, clip_on=False)
axes[0].tick_params(labelsize=9); axes[0].spines[["top", "right"]].set_visible(False)

# (B) width ratios (n_T=80, γ=0.2)
x = np.arange(len(SCN)); w = 0.38
axes[1].bar(x - w / 2, idr, w, color=NAVY3, label="Identified set")
axes[1].bar(x + w / 2, cir, w, color=STEEL3, hatch="//", edgecolor=NAVY3, label="Confidence envelope")
axes[1].axhline(1.0, color="0.6", ls="--", lw=1)
axes[1].set_xticks(x); axes[1].set_xticklabels(SCN_LABEL, fontsize=9)
axes[1].set_ylabel("Separate / direct width ratio", fontsize=9.5); axes[1].tick_params(axis="y", labelsize=9)
axes[1].set_title("(B) Width reduction", loc="left", fontsize=TT)
axes[1].set_ylim(0, 4.4)
axes[1].legend(fontsize=8.5, frameon=False, loc="upper center", ncol=2,
               bbox_to_anchor=(0.5, 0.99), columnspacing=1.0, handlelength=2.0)
axes[1].spines[["top", "right"]].set_visible(False)

# (C) per-γ true-contrast containment (n_T=80, γ=0.2); error bars = stored MC Wilson CI
lo = [c - l for c, l in zip(cov, ci_lo)]; hi = [h - c for c, h in zip(cov, ci_hi)]
colors = [NAVY3] * len(cov)
axes[2].bar(x, cov, 0.6, color=colors, yerr=[lo, hi], capsize=5, ecolor="0.3", error_kw={"lw": 1.1})
axes[2].set_ylim(0.80, 1.015)
axes[2].set_xticks(x); axes[2].set_xticklabels(SCN_LABEL, fontsize=9); axes[2].tick_params(axis="y", labelsize=9)
axes[2].set_ylabel("True-contrast containment", fontsize=9.5)
axes[2].set_title("(C) Contrast containment", loc="left", fontsize=TT)
for xi, yi, y_top in zip(x, cov, ci_hi):
    axes[2].annotate(f"{yi:.3f}", (xi, y_top), textcoords="offset points",
                     xytext=(0, 6), ha="center", va="bottom", fontsize=8.5,
                     bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none"},
                     zorder=5, clip_on=False)
axes[2].spines[["top", "right"]].set_visible(False)
fig.subplots_adjust(left=0.09, right=0.985, bottom=0.09, top=0.96)
save(fig, "Figure_3")
plt.show()

# %% [markdown]
# ## Figure 4 — NSCLC outcome-free evaluation
# Section 4.6. Direct vs separate sharp identified sets (thick bars) + simultaneous
# confidence envelopes (thin capped), seed 42, γ=0.2. Panels
# (A) Stanford, (B) VA. Read from the final manuscript aggregate.

# %%
FIG4_CSV = os.path.join(BASE, "results", "nsclc", "outcome_free_evaluation.csv")


def _f4(ds_key):
    for r in csv.DictReader(open(FIG4_CSV)):
        if (r["dataset"] == ds_key and r["contrast"] == "PET_radiomic_candidate"
                and int(r["seed"]) == 42 and abs(float(r["gamma"]) - 0.2) < 1e-9):
            def g(c):
                return float(r[c])
            idr = (g("L_dir"), g("U_dir")); ids = (g("L_sep"), g("U_sep"))
            denv = (g("dir_sim_lower"), g("dir_sim_upper"))
            senv = (g("sep_sim_lower"), g("sep_sim_upper"))
            ratio = g("ID_width_ratio"); dec = r["action_dir"].strip()
            samp = (not (idr[0] <= 0 <= idr[1])) and (denv[0] < 0 < denv[1])
            kind = "sampling-limited" if samp else dec.lower()
            return idr, ids, denv, senv, kind, ratio
    raise ValueError(f"Fig4 row not found: {ds_key}, seed 42, gamma 0.2")


cfg = [("(A) Stanford",) + _f4("NSCLC-Stanford"), ("(B) VA",) + _f4("NSCLC-VA")]
print("Fig4 read from:", FIG4_CSV)
c_dir = "#1f4e79"; c_sep = "#c0504d"
fig, axes = plt.subplots(2, 1, figsize=(DOUBLE, 4.2))
for ax, (title, idr, ids, denv, senv, kind, ratio_tbl) in zip(axes, cfg):
    yD, yS = 0.68, 0.32
    ax.plot(senv, [yS, yS], color=c_sep, lw=2.2, solid_capstyle="butt", zorder=2)
    for x in senv: ax.plot([x], [yS], "|", color=c_sep, ms=10, mew=1.8, zorder=3)
    ax.plot(ids, [yS, yS], color=c_sep, lw=7, solid_capstyle="round", zorder=4)
    ax.plot(denv, [yD, yD], color=c_dir, lw=2.2, solid_capstyle="butt", zorder=2)
    for x in denv: ax.plot([x], [yD], "|", color=c_dir, ms=10, mew=1.8, zorder=3)
    ax.plot(idr, [yD, yD], color=c_dir, lw=7, solid_capstyle="round", zorder=4)
    ax.axvline(0, color="0.35", ls="--", lw=1.1, zorder=1)
    ax.set_xlim(-0.06, 0.14); ax.set_ylim(0, 1)
    ax.set_yticks([yS, yD]); ax.set_yticklabels(["separate", "direct"], fontsize=9)
    ax.tick_params(axis="x", labelsize=10.5)
    ax.set_title(f"{title}    |    ID width ratio {ratio_tbl:.2f}×    |    DEFER ({kind})", fontsize=9.5, loc="left")
    ax.tick_params(axis="y", length=0); ax.spines[["top", "right", "left"]].set_visible(False)
axes[-1].set_xlabel(r"Target Brier-risk contrast  $\Delta R_T$  (candidate $-$ reference)", fontsize=9.5)
fig.tight_layout()
save(fig, "Figure_4")
plt.show()

# %% [markdown]
# ## Figure S1 — Rotterdam→GBSG decisions across locked configs and γ
# (A) direct identification, (B) separate-risk subtraction. seed × γ grid; read from
# supporting-analysis CSV. Starred/bordered cell (seed 62, γ=0.2): direct = ADOPT, separate = DEFER.

# %%
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, Rectangle

SRC = os.path.join(BASE, "results", "applications", "direct_vs_separate.csv")
DATASET = "Rotterdam-GBSG"; CAND_TAG = "supporting_clinical_candidate"
seeds = [42, 52, 62, 72]; gammas = [0.0, 0.1, 0.2, 0.3, 0.5]

rows = [r for r in csv.DictReader(open(SRC)) if r["dataset"] == DATASET and r["contrast"] == CAND_TAG]
assert len(rows) == len(seeds) * len(gammas), f"expected {len(seeds)*len(gammas)} rows, got {len(rows)}"
code = {"ADOPT": 1, "DEFER": 0, "KEEP": 2}
Adir = np.full((len(seeds), len(gammas)), -1, int)
Asep = np.full((len(seeds), len(gammas)), -1, int)
for r in rows:
    i = seeds.index(int(r["seed"])); j = gammas.index(round(float(r["gamma"]), 3))
    Adir[i, j] = code.get(r["action_dir"].strip().upper().split()[0], 0)
    Asep[i, j] = code.get(r["action_sep"].strip().upper().split()[0], 0)
assert (Adir >= 0).all() and (Asep >= 0).all(), "missing cell in source CSV"

GREEN = "#2e7d32"; GREY = "#c9ccd1"; NAVYS = "#1f3b73"; STAR = "#c0392b"
cmap = ListedColormap([GREY, GREEN])
fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 3.2))


def _draw(ax, M, title):
    ax.imshow(M, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(gammas))); ax.set_xticklabels([f"{g:g}" for g in gammas], fontsize=9)
    ax.set_yticks(range(len(seeds))); ax.set_yticklabels([f"seed {s}" for s in seeds], fontsize=9)
    ax.set_xlabel(r"Sensitivity value $\gamma$", fontsize=9.5)
    ax.set_title(title, fontsize=10, color=NAVYS, fontweight="bold")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, "ADOPT" if M[i, j] == 1 else "DEFER", ha="center", va="center",
                fontsize=8, color="white" if M[i, j] == 1 else "#333333")
    ax.set_xticks(np.arange(-.5, len(gammas), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(seeds), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5); ax.tick_params(which="minor", length=0)


_draw(axes[0], Adir, "(A) Direct identification")
_draw(axes[1], Asep, "(B) Separate-risk subtraction")
si, gj = seeds.index(62), gammas.index(0.2)
for ax in axes:
    ax.add_patch(Rectangle((gj - .5, si - .5), 1, 1, fill=False, edgecolor=STAR, linewidth=2.8, zorder=5))
    ax.text(gj + .34, si - .34, "*", ha="center", va="center", fontsize=12, color=STAR, fontweight="bold", zorder=6)
fig.legend(handles=[Patch(facecolor=GREEN, label="ADOPT CANDIDATE"), Patch(facecolor=GREY, label="DEFER")],
           loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.02), fontsize=9)
fig.tight_layout(rect=[0, 0.05, 1, 1])
save(fig, "Figure_S1")
print("S1 rows used:", len(rows), "| SRC:", SRC)
plt.show()

# %% [markdown]
# ## Figure S2 — Decision frequencies across prespecified simulation regimes
# (A) independent (marginal) censoring, (B) covariate-dependent (conditional).
# Frequencies read from locked aggregate JSON. retain=ADOPT, reference=KEEP, defer=DEFER.

# %%
BDIR = os.path.join(BASE, "results", "simulations", "decision_regime")
regimes = ["candidate_superior", "finite_crossing", "no_crossing_within_range", "reference_superior", "ambiguous"]
labels = ["Candidate\nsuperior", "Finite\ncrossing", "No\ncrossing", "Reference\nsuperior", "Ambiguous"]


def _rates(regime, cens):
    d = json.load(open(os.path.join(BDIR, f"{regime}_{cens}.json")))["decision_rates"]
    return [d["retain"], d["reference"], d["defer"]]   # ADOPT, KEEP, DEFER


GREEN2 = "#2e7d32"; BLUE2 = "#2d6cb5"; GREY2 = "#b8bdc6"; NAVY4 = "#1f3b73"
cats = ["ADOPT CANDIDATE", "KEEP REFERENCE", "DEFER"]; colors = [GREEN2, BLUE2, GREY2]
fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 3.4), sharey=True)
x = np.arange(len(regimes))
for ax, cens, title in zip(axes, ["marginal", "conditional"],
                           ["(A) Independent censoring", "(B) Covariate-dependent censoring"]):
    M = np.array([_rates(r, cens) for r in regimes])
    bottom = np.zeros(len(regimes))
    for k in range(3):
        ax.bar(x, M[:, k], bottom=bottom, color=colors[k], width=0.62, edgecolor="white", linewidth=0.6, label=cats[k])
        for xi, (v, b) in enumerate(zip(M[:, k], bottom)):
            if v >= 0.06:
                ax.text(xi, b + v / 2, f"{v:.2f}", ha="center", va="center", fontsize=8,
                        color="white" if k != 2 else "#333333")
        bottom += M[:, k]
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8.5)
    ax.tick_params(axis="y", labelsize=8)
    ax.set_ylim(0, 1.0); ax.set_title(title, fontsize=10, color=NAVY4, fontweight="bold")
    ax.grid(axis="y", color="0.9", linewidth=0.8); ax.set_axisbelow(True)
axes[0].set_ylabel("Monte Carlo decision frequency", fontsize=9.5)
handles, legend_labels = axes[1].get_legend_handles_labels()
fig.legend(handles, legend_labels, loc="lower center", ncol=3,
           bbox_to_anchor=(0.5, -0.015), frameon=False, fontsize=8.5)
fig.tight_layout(rect=[0, 0.12, 1, 1])
save(fig, "Figure_S2")
print("S2 JSON dir:", BDIR)
plt.show()
