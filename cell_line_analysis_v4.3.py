"""
Cell Line Immortalization Analysis
===================================
Analyzes association between Drosophila species/fly line and:
  1. Primary line formation rate
  2. Immortalization rate (conditional on primary line formation)
  3. Time to immortalization (mel vs sim; Mann-Whitney U)

Data notes:
  - Accepts CSV (comma-separated) input files.
  - Contaminated flasks (Contamination? = 'Y') and maybe-contaminated ('M') are
    excluded from all analyses.
  - "ongoing" in Immortalized? means outcome not yet known; excluded from the
    immortalization denominator entirely.
  - "harvested" in Immortalized? means the primary culture was deliberately
    collected before immortalization; plotted as a separate outcome category.
  - Immortalization rate is calculated only among flasks that formed primary cultures
    AND have a resolved immortalization outcome (Y or N).

Statistical tests:
  - Between-species: Chi-square test of independence (2×4 contingency tables)
  - Riv84 vs Mer23 (wRi strains): Fisher's Exact Test (2×2 tables, small sample sizes)
  - Within-species comparisons: assessed but not formally tested where sample sizes
    are too small or variance is too low for meaningful inference.
  - Time to immortalization (mel vs sim): Mann-Whitney U (non-parametric).

Outputs (all PDF):
  01_overview.pdf
  02_between_species.pdf
  03_within_species.pdf
  04_wRi_Riv84_vs_Mer23.pdf
  07_outcomes_with_harvested.pdf
  08_time_to_immortalization.pdf
  09_mel_vs_sim_tti_test.pdf
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy import stats
import warnings
import sys
warnings.filterwarnings("ignore")

# ── Styling ──────────────────────────────────────────────────────────────────
SPECIES_ORDER   = ["melanogaster", "simulans", "willistoni", "F1hybrid:melxsim"]
SPECIES_LABELS  = ["D. melanogaster", "D. simulans", "D. willistoni", "F1 hybrid\n(mel × sim)"]
SPECIES_COLORS  = ["#0077BB", "#BB2200", "#228833", "#7B3FA0"]
SPECIES_PALETTE = dict(zip(SPECIES_ORDER, SPECIES_COLORS))

WRI_COLORS = {"D._simulans_Riv84": "#BB2200", "D._simulans_Mer23_#5": "#EE6644"}
WRI_LABELS = {"D._simulans_Riv84": "Riv84 (wRi)", "D._simulans_Mer23_#5": "Mer23 (wRi)"}

LINE_SHORT = {
    "Dsim_w[501]":                                   "w[501]",
    "D._simulans_Riv84":                             "Riv84",
    "D._simulans_Mer23_#5":                          "Mer23",
    "D._melanogaster_ub:Khc":                        "ub:Khc",
    "D._melogaster_PlasMem:GFP":                     "PlasMem:GFP",
    "D._melanogaster_mei-P26[1]":                    "mei-P26[1]",
    "D._melanogaster_Jupiter:GFP":                   "Jupiter:GFP",
    "D._willistoni_(genome_line)":                   "genome line",
    "D._melanogaster_virginator_x_D._simulans_w[501]_F1_hybrid": "F1 hybrid",
}

plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "axes.spines.top":  False,
    "axes.spines.right":False,
    "axes.titlesize":   13,
    "axes.labelsize":   11,
    "xtick.labelsize":  10,
    "ytick.labelsize":  10,
    "figure.dpi":       150,
})


# ── Data loading and cleaning ────────────────────────────────────────────────
def load_and_clean(path: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(path)   # CSV (comma-separated)

    # Normalise Immortalized? to lowercase for consistent matching
    df["Immortalized?"] = df["Immortalized?"].astype(str).str.strip().str.lower()

    # Exclude contaminated (Y) and maybe-contaminated (M) flasks
    df_clean = df[df["Contamination?"] == "N"].copy()

    # Binary flag for primary culture formation
    df_clean["primary"] = (df_clean["Primary_Culture?"] == "Y").astype(int)

    # Time-to-immortalization (numeric days; "NY" and non-numeric → NaN)
    df["tti"] = pd.to_numeric(df["Time to immortalization"], errors="coerce")

    # Subset that formed primary cultures AND have a resolved immortalization outcome.
    # "ongoing" and "harvested" are excluded from the immortalization rate denominator:
    #   ongoing  → outcome not yet known
    #   harvested → deliberately collected; not a failed or successful immortalization
    df_primary = df_clean[
        (df_clean["Primary_Culture?"] == "Y") &
        (~df_clean["Immortalized?"].isin(["ongoing", "harvested"]))
    ].copy()
    df_primary["immortalized_flag"] = (df_primary["Immortalized?"] == "y").astype(int)

    return df, df_clean, df_primary


# ── Helper: error bar from proportions ───────────────────────────────────────
def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score confidence interval for a proportion."""
    if n == 0:
        return 0.0, 0.0
    p_hat = k / n
    denom = 1 + z**2 / n
    centre = (p_hat + z**2 / (2 * n)) / denom
    margin = (z * np.sqrt(p_hat * (1 - p_hat) / n + z**2 / (4 * n**2))) / denom
    lo = max(0.0, centre - margin)
    hi = min(1.0, centre + margin)
    return lo, hi


def report_categorical_test(ct: list, label: str, mode: str = "both") -> None:
    """
    Run and print chi-square and/or Fisher's exact test on a contingency table.

    Parameters
    ----------
    ct   : list of [successes, failures] per group (rows of the contingency table)
    label: short description printed before the results
    mode : "chi2"   — chi-square only
           "fisher" — Fisher's exact only
           "both"   — both (default)

    Automatically checks expected cell counts for chi-square and prints a warning
    when any expected count < 5, recommending Fisher's exact in that case.
    """
    ct_arr = np.array(ct)
    is_2x2 = ct_arr.shape == (2, 2)

    # ── Chi-square ────────────────────────────────────────────────────────────
    if mode in ("chi2", "both"):
        chi2, p_chi, dof, expected = stats.chi2_contingency(ct_arr)
        min_exp = expected.min()
        low_exp = min_exp < 5
        sig = ("***" if p_chi < 0.001 else "**" if p_chi < 0.01
               else "*" if p_chi < 0.05 else "ns")
        p_str = f"p = {p_chi:.2e}" if p_chi < 0.001 else f"p = {p_chi:.4f}"
        print(f"    → Chi-square: χ²({dof}) = {chi2:.2f}, {p_str}  ({sig})")
        if low_exp:
            print(f"       ⚠  Min expected cell count = {min_exp:.1f} (< 5); "
                  f"chi-square approximation may be unreliable — see Fisher's exact below")

    # ── Fisher's exact ────────────────────────────────────────────────────────
    if mode in ("fisher", "both"):
        if is_2x2:
            _, p_fish = stats.fisher_exact(ct_arr)
            sig = ("***" if p_fish < 0.001 else "**" if p_fish < 0.01
                   else "*" if p_fish < 0.05 else "ns")
            p_str = f"p = {p_fish:.2e}" if p_fish < 0.001 else f"p = {p_fish:.4f}"
            print(f"    → Fisher's exact: {p_str}  ({sig})")
        else:
            # scipy's fisher_exact only handles 2×2; use barnard or report chi2 only
            # For larger tables use the chi2 result with an advisory note
            if mode == "fisher":
                chi2, p_chi, dof, expected = stats.chi2_contingency(ct_arr)
                min_exp = expected.min()
                sig = ("***" if p_chi < 0.001 else "**" if p_chi < 0.01
                       else "*" if p_chi < 0.05 else "ns")
                p_str = f"p = {p_chi:.2e}" if p_chi < 0.001 else f"p = {p_chi:.4f}"
                print(f"    → Fisher's exact: not available for {ct_arr.shape[0]}×2 tables; "
                      f"reporting chi-square: χ²({dof}) = {chi2:.2f}, {p_str}  ({sig})")
                if min_exp < 5:
                    print(f"       ⚠  Min expected cell count = {min_exp:.1f} (< 5)")
            else:
                print(f"    → Fisher's exact: not available for {ct_arr.shape[0]}×2 tables "
                      f"(scipy limitation); chi-square result shown above")


def rate_stats(df: pd.DataFrame, group_col: str, flag_col: str) -> pd.DataFrame:
    rows = []
    for grp, sub in df.groupby(group_col):
        k = sub[flag_col].sum()
        n = len(sub)
        p = k / n if n else np.nan
        lo, hi = wilson_ci(int(k), n)
        rows.append({"group": grp, "successes": k, "trials": n,
                     "rate": p, "ci_lo": lo, "ci_hi": hi})
    return pd.DataFrame(rows).set_index("group")


# ── PLOT 1: Overview / Exploratory ───────────────────────────────────────────
def plot_overview(df: pd.DataFrame, df_clean: pd.DataFrame,
                  df_primary: pd.DataFrame, save_path: str):

    fig = plt.figure(figsize=(16, 12))
    fig.suptitle("Overview: Drosophila Cell Line Immortalization Dataset",
                 fontsize=15, fontweight="bold", y=0.98)

    gs = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.38)
    ax1 = fig.add_subplot(gs[0, 0])   # flask counts by species (all data)
    ax2 = fig.add_subplot(gs[0, 1])   # contamination breakdown by species
    ax3 = fig.add_subplot(gs[0, 2])   # primary + immortalization funnel
    ax4 = fig.add_subplot(gs[1, 0])   # primary rate by species
    ax5 = fig.add_subplot(gs[1, 1])   # immortalization rate by species
    ax6 = fig.add_subplot(gs[1, 2])   # stacked outcome by species

    # ── ax1: total flasks by species ────────────────────────────────────────
    counts = df.groupby("Drosophila_species").size().reindex(SPECIES_ORDER)
    bars = ax1.bar(SPECIES_LABELS, counts.values,
                   color=[SPECIES_PALETTE[s] for s in SPECIES_ORDER], edgecolor="white")
    for bar, v in zip(bars, counts.values):
        ax1.text(bar.get_x() + bar.get_width()/2, v + 0.5, str(v),
                 ha="center", va="bottom", fontsize=9)
    ax1.set_title("Total Flasks Attempted")
    ax1.set_ylabel("Number of flasks")
    ax1.set_ylim(0, max(counts.values) * 1.18)

    # ── ax2: contamination by species ───────────────────────────────────────
    cont_map = {"N": "Clean", "Y": "Contaminated", "M": "Maybe contam."}
    cont_colors = {"Clean": "#2ECC71", "Contaminated": "#E74C3C", "Maybe contam.": "#F39C12"}
    cont_df = df.copy()
    cont_df["Contam_label"] = cont_df["Contamination?"].map(cont_map)
    ct = (cont_df.groupby(["Drosophila_species", "Contam_label"])
                 .size().unstack(fill_value=0)
                 .reindex(SPECIES_ORDER))
    ct = ct.reindex(columns=["Clean", "Contaminated", "Maybe contam."], fill_value=0)
    bottoms = np.zeros(len(SPECIES_ORDER))
    for status, color in cont_colors.items():
        if status in ct.columns:
            vals = ct[status].values
            ax2.bar(SPECIES_LABELS, vals, bottom=bottoms,
                    color=color, label=status, edgecolor="white")
            bottoms += vals
    ax2.set_title("Flask Outcomes by Contamination")
    ax2.set_ylabel("Number of flasks")
    ax2.legend(fontsize=8, loc="upper right")

    # ── ax3: outcome funnel (all clean flasks) ───────────────────────────────
    total_clean = len(df_clean)
    total_primary = (df_clean["primary"] == 1).sum()
    total_imm = int(df_primary["immortalized_flag"].sum())

    funnel_vals  = [total_clean, total_primary, total_imm]
    funnel_labels = [f"Clean flasks\n(n={total_clean}, "
                     f"{total_clean/len(df)*100:.0f}% of total)",
                     f"Formed primary\nculture\n(n={total_primary}, "
                     f"{total_primary/total_clean*100:.0f}% of uncontaminated)",
                     f"Immortalized\n(n={total_imm}, "
                     f"{total_imm/total_primary*100:.0f}% of primary)"]
    funnel_colors = ["#3498DB", "#2ECC71", "#9B59B6"]
    widths = [v / total_clean for v in funnel_vals]
    for i, (w, lab, col) in enumerate(zip(widths, funnel_labels, funnel_colors)):
        ax3.barh(2 - i, w, color=col, alpha=0.85, edgecolor="white", height=0.7)
        ax3.text(w + 0.01, 2 - i, lab, va="center", fontsize=8.5)
    ax3.set_xlim(0, 1.65)
    ax3.set_yticks([])
    ax3.set_xlabel("Proportion of clean flasks")
    ax3.set_title("Overall Outcome Funnel\n(non-contaminated flasks)")
    ax3.spines["left"].set_visible(False)

    # ── ax4: primary culture rate by species ─────────────────────────────────
    pr = rate_stats(df_clean, "Drosophila_species", "primary").reindex(SPECIES_ORDER)
    bars = ax4.bar(SPECIES_LABELS, pr["rate"].values,
                   color=[SPECIES_PALETTE[s] for s in SPECIES_ORDER], edgecolor="white")
    for i, (bar, row) in enumerate(zip(bars, pr.itertuples())):
        ax4.errorbar(bar.get_x() + bar.get_width()/2, row.rate,
                     yerr=[[row.rate - row.ci_lo], [row.ci_hi - row.rate]],
                     fmt="none", color="black", capsize=4, linewidth=1.5)
        ax4.text(bar.get_x() + bar.get_width()/2,
                 row.ci_hi + 0.03, f"{int(row.successes)}/{int(row.trials)}",
                 ha="center", va="bottom", fontsize=8.5)
    ax4.set_ylim(0, 1.25)
    ax4.set_ylabel("Primary culture rate")
    ax4.set_title("Primary Culture Formation Rate\nby Species (clean flasks)")
    ax4.axhline(0.5, ls="--", color="grey", alpha=0.4, lw=1)

    # ── ax5: immortalization rate by species ──────────────────────────────────
    ir = rate_stats(df_primary, "Drosophila_species", "immortalized_flag").reindex(SPECIES_ORDER)
    bars = ax5.bar(SPECIES_LABELS, ir["rate"].values,
                   color=[SPECIES_PALETTE[s] for s in SPECIES_ORDER], edgecolor="white")
    for bar, row in zip(bars, ir.itertuples()):
        ax5.errorbar(bar.get_x() + bar.get_width()/2, row.rate,
                     yerr=[[row.rate - row.ci_lo], [row.ci_hi - row.rate]],
                     fmt="none", color="black", capsize=4, linewidth=1.5)
        ax5.text(bar.get_x() + bar.get_width()/2,
                 row.ci_hi + 0.03, f"{int(row.successes)}/{int(row.trials)}",
                 ha="center", va="bottom", fontsize=8.5)
    ax5.set_ylim(0, 1.25)
    ax5.set_ylabel("Immortalization rate")
    ax5.set_title("Immortalization Rate by Species\n(of primary cultures formed)")
    ax5.axhline(0.5, ls="--", color="grey", alpha=0.4, lw=1)

    # ── ax6: stacked bar of all outcomes ─────────────────────────────────────
    outcome_colors = {"No primary culture": "#E74C3C",
                      "Primary only": "#F39C12",
                      "Immortalized": "#2ECC71"}
    stacks = []
    for sp in SPECIES_ORDER:
        sub = df_clean[df_clean["Drosophila_species"] == sp]
        no_prim = (sub["Primary_Culture?"] == "N").sum()
        prim_only = ((sub["Primary_Culture?"] == "Y") &
                     (~sub["Immortalized?"].isin(["y", "harvested"]))).sum()
        imm = (sub["Immortalized?"] == "y").sum()
        stacks.append([no_prim, prim_only, imm])
    stacks = np.array(stacks)  # (4 species × 3 outcomes)
    bottoms = np.zeros(4)
    for j, (label, color) in enumerate(outcome_colors.items()):
        ax6.bar(SPECIES_LABELS, stacks[:, j], bottom=bottoms,
                color=color, label=label, edgecolor="white")
        bottoms += stacks[:, j]
    ax6.set_title("Outcome Stacks by Species\n(clean flasks)")
    ax6.set_ylabel("Number of flasks")
    ax6.legend(fontsize=8, loc="upper right")

    plt.savefig(save_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {save_path}")


# ── PLOT 2: Between-species statistical tests ────────────────────────────────
def plot_between_species(df_clean: pd.DataFrame, df_primary: pd.DataFrame,
                         save_path: str):

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    fig.suptitle("Between-Species Comparisons", fontsize=14, fontweight="bold")

    # ── Chi-square for primary culture rate ──────────────────────────────────
    ax = axes[0]
    pr = rate_stats(df_clean, "Drosophila_species", "primary").reindex(SPECIES_ORDER)
    # Chi-square
    ct = [[r.successes, r.trials - r.successes] for r in pr.itertuples()]
    chi2, p_prim, dof, _ = stats.chi2_contingency(ct)
    bars = ax.bar(SPECIES_LABELS, pr["rate"].values,
                  color=[SPECIES_PALETTE[s] for s in SPECIES_ORDER],
                  edgecolor="white", width=0.6)
    for bar, row in zip(bars, pr.itertuples()):
        ax.errorbar(bar.get_x() + bar.get_width()/2, row.rate,
                    yerr=[[row.rate - row.ci_lo], [row.ci_hi - row.rate]],
                    fmt="none", color="black", capsize=5, linewidth=2)
        ax.text(bar.get_x() + bar.get_width()/2,
                min(row.ci_hi + 0.04, 1.08),
                f"{int(row.successes)}/{int(row.trials)}\n({row.rate*100:.0f}%)",
                ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, 1.35)
    ax.set_ylabel("Primary culture formation rate", fontsize=11)
    ax.set_title("Primary Culture Formation Rate\nby Drosophila Species")
    p_str = f"p = {p_prim:.2e}" if p_prim < 0.001 else f"p = {p_prim:.4f}"
    ax.text(0.5, 1.28, f"χ²({dof}) = {chi2:.2f}, {p_str}",
            transform=ax.transAxes, ha="center", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", fc="#EBF5FB", ec="#85C1E9"))
    ax.set_xlabel("Error bars = 95% Wilson CI", fontsize=9)

    # ── Chi-square for immortalization rate ──────────────────────────────────
    ax = axes[1]
    ir = rate_stats(df_primary, "Drosophila_species", "immortalized_flag").reindex(SPECIES_ORDER)
    ct2 = [[r.successes, r.trials - r.successes] for r in ir.itertuples()]
    chi2b, p_imm, dof2, _ = stats.chi2_contingency(ct2)
    bars = ax.bar(SPECIES_LABELS, ir["rate"].values,
                  color=[SPECIES_PALETTE[s] for s in SPECIES_ORDER],
                  edgecolor="white", width=0.6)
    for bar, row in zip(bars, ir.itertuples()):
        ax.errorbar(bar.get_x() + bar.get_width()/2, row.rate,
                    yerr=[[row.rate - row.ci_lo], [row.ci_hi - row.rate]],
                    fmt="none", color="black", capsize=5, linewidth=2)
        ax.text(bar.get_x() + bar.get_width()/2,
                min(row.ci_hi + 0.04, 1.08),
                f"{int(row.successes)}/{int(row.trials)}\n({row.rate*100:.0f}%)",
                ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, 1.35)
    ax.set_ylabel("Immortalization rate (of primary cultures)", fontsize=11)
    ax.set_title("Immortalization Rate\nby Drosophila Species")
    p_str2 = f"p = {p_imm:.2e}" if p_imm < 0.001 else f"p = {p_imm:.4f}"
    ax.text(0.5, 1.28, f"χ²({dof2}) = {chi2b:.2f}, {p_str2}",
            transform=ax.transAxes, ha="center", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", fc="#EBF5FB", ec="#85C1E9"))
    ax.set_xlabel("Error bars = 95% Wilson CI", fontsize=9)

    plt.tight_layout()
    plt.savefig(save_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {save_path}")


# ── PLOT 3: Within-species sample size assessment ────────────────────────────
def plot_within_species(df_clean: pd.DataFrame, df_primary: pd.DataFrame,
                        save_path: str):

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Within-Species Fly Line Comparisons\n"
                 "(formal testing only where sample size permits)",
                 fontsize=13, fontweight="bold")

    for ax_idx, (sp, ax_pc, ax_im) in enumerate([
        ("melanogaster",    axes[0, 0], axes[1, 0]),
        ("simulans",        axes[0, 1], axes[1, 1]),
    ]):
        sub_clean = df_clean[df_clean["Drosophila_species"] == sp]
        sub_prim  = df_primary[df_primary["Drosophila_species"] == sp]
        lines     = sub_clean["Fly_Line"].unique()
        short_labels = [LINE_SHORT.get(l, l) for l in lines]
        colors = [WRI_COLORS.get(l, SPECIES_PALETTE[sp]) for l in lines]

        # Primary culture rate per line
        pc_vals, pc_lo, pc_hi, pc_counts = [], [], [], []
        for l in lines:
            sub = sub_clean[sub_clean["Fly_Line"] == l]
            k, n = sub["primary"].sum(), len(sub)
            pc_vals.append(k / n if n else 0)
            lo, hi = wilson_ci(int(k), n)
            pc_lo.append(lo); pc_hi.append(hi)
            pc_counts.append((int(k), n))

        bars = ax_pc.bar(short_labels, pc_vals, color=colors, edgecolor="white", width=0.6)
        for bar, rate, lo, hi, (k, n) in zip(bars, pc_vals, pc_lo, pc_hi, pc_counts):
            ax_pc.errorbar(bar.get_x() + bar.get_width()/2, rate,
                           yerr=[[rate - lo], [hi - rate]],
                           fmt="none", color="black", capsize=4, lw=1.5)
            ax_pc.text(bar.get_x() + bar.get_width()/2,
                       min(hi + 0.03, 1.08), f"{k}/{n}",
                       ha="center", va="bottom", fontsize=8.5)
        ax_pc.set_ylim(0, 1.35)
        ax_pc.set_ylabel("Primary culture rate")
        ax_pc.set_title(f"D. {sp} — Primary Culture Formation Rate")
        ax_pc.tick_params(axis="x", rotation=20)

        # Immortalization rate per line
        im_vals, im_lo, im_hi, im_counts = [], [], [], []
        for l in lines:
            sub = sub_prim[sub_prim["Fly_Line"] == l]
            k, n = sub["immortalized_flag"].sum(), len(sub)
            im_vals.append(k / n if n else 0)
            lo, hi = wilson_ci(int(k), n)
            im_lo.append(lo); im_hi.append(hi)
            im_counts.append((int(k), n))

        bars = ax_im.bar(short_labels, im_vals, color=colors, edgecolor="white", width=0.6)
        for bar, rate, lo, hi, (k, n) in zip(bars, im_vals, im_lo, im_hi, im_counts):
            if n > 0:
                ax_im.errorbar(bar.get_x() + bar.get_width()/2, rate,
                               yerr=[[rate - lo], [hi - rate]],
                               fmt="none", color="black", capsize=4, lw=1.5)
            ax_im.text(bar.get_x() + bar.get_width()/2,
                       min(hi + 0.03, 1.08), f"{k}/{n}",
                       ha="center", va="bottom", fontsize=8.5)
        ax_im.set_ylim(0, 1.35)
        ax_im.set_ylabel("Immortalization rate")
        ax_im.set_title(f"D. {sp} — Immortalization Rate (of primary cultures)")
        ax_im.tick_params(axis="x", rotation=20)

        # Annotation for insufficient data
        if sp == "melanogaster":
            ax_pc.text(0.5, 1.20,
                       "⚠ Sample sizes too small / variance too low\nfor formal within-species test",
                       transform=ax_pc.transAxes, ha="center", fontsize=8.5,
                       color="#7D6608",
                       bbox=dict(boxstyle="round,pad=0.3", fc="#FEF9E7", ec="#F0B27A"))
            ax_im.text(0.5, 1.20,
                       "⚠ ub:Khc only line with immortalized cells\n(n too small for formal test)",
                       transform=ax_im.transAxes, ha="center", fontsize=8.5,
                       color="#7D6608",
                       bbox=dict(boxstyle="round,pad=0.3", fc="#FEF9E7", ec="#F0B27A"))

    plt.tight_layout()
    plt.savefig(save_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {save_path}")


# ── PLOT 4: Riv84 vs Mer23 (wRi strains) ────────────────────────────────────
def plot_wri_comparison(df_clean: pd.DataFrame, df_primary: pd.DataFrame,
                        save_path: str):

    df_sim   = df_clean[df_clean["Drosophila_species"] == "simulans"]
    df_sim_p = df_primary[df_primary["Drosophila_species"] == "simulans"]

    wri_lines = ["D._simulans_Riv84", "D._simulans_Mer23_#5"]
    df_wri    = df_sim[df_sim["Fly_Line"].isin(wri_lines)]
    df_wri_p  = df_sim_p[df_sim_p["Fly_Line"].isin(wri_lines)]

    # Fisher's exact — primary culture
    riv  = df_wri[df_wri["Fly_Line"] == "D._simulans_Riv84"]
    mer  = df_wri[df_wri["Fly_Line"] == "D._simulans_Mer23_#5"]
    riv_k, riv_n = int(riv["primary"].sum()), len(riv)
    mer_k, mer_n = int(mer["primary"].sum()), len(mer)
    _, p_pc = stats.fisher_exact([[riv_k, riv_n - riv_k],
                                   [mer_k, mer_n - mer_k]])

    # Fisher's exact — immortalization
    riv_p = df_wri_p[df_wri_p["Fly_Line"] == "D._simulans_Riv84"]
    mer_p = df_wri_p[df_wri_p["Fly_Line"] == "D._simulans_Mer23_#5"]
    riv_ik, riv_in = int(riv_p["immortalized_flag"].sum()), len(riv_p)
    mer_ik, mer_in = int(mer_p["immortalized_flag"].sum()), len(mer_p)
    _, p_im = stats.fisher_exact([[riv_ik, riv_in - riv_ik],
                                   [mer_ik, mer_in - mer_ik]])

    fig, axes = plt.subplots(1, 2, figsize=(11, 6))
    fig.suptitle("D. simulans wRi Strains: Riv84 vs Mer23\n(Fisher's Exact Test)",
                 fontsize=13, fontweight="bold")

    for ax, (k_vals, n_vals, lines, p_val, ylabel, title) in zip(axes, [
        (
            [riv_k / riv_n, mer_k / mer_n],
            [(riv_k, riv_n), (mer_k, mer_n)],
            wri_lines,
            p_pc,
            "Primary culture formation rate",
            "Primary Culture Formation Rate",
        ),
        (
            [riv_ik / riv_in if riv_in else 0, mer_ik / mer_in if mer_in else 0],
            [(riv_ik, riv_in), (mer_ik, mer_in)],
            wri_lines,
            p_im,
            "Immortalization rate (of primary cultures)",
            "Immortalization Rate",
        ),
    ]):
        bar_colors = [WRI_COLORS[l] for l in lines]
        bar_labels = [WRI_LABELS[l] for l in lines]
        bars = ax.bar(bar_labels, k_vals, color=bar_colors, edgecolor="white",
                      width=0.45)
        for bar, rate, (k, n) in zip(bars, k_vals, n_vals):
            lo, hi = wilson_ci(k, n)
            if n > 0:
                ax.errorbar(bar.get_x() + bar.get_width()/2, rate,
                            yerr=[[rate - lo], [hi - rate]],
                            fmt="none", color="black", capsize=6, lw=2)
            label = f"{k}/{n}\n({rate*100:.0f}%)" if n else "0/0"
            ax.text(bar.get_x() + bar.get_width()/2,
                    min((hi if n else 0) + 0.04, 1.07), label,
                    ha="center", va="bottom", fontsize=11, fontweight="bold")

        ax.set_ylim(0, 1.35)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title)

        # Significance annotation
        p_fmt = f"p = {p_val:.2e}" if p_val < 0.001 else f"p = {p_val:.4f}"
        sig = "***" if p_val < 0.001 else ("**" if p_val < 0.01 else
               ("*" if p_val < 0.05 else "ns"))
        color = "#E74C3C" if sig not in ("ns",) else "#7F8C8D"
        ax.text(0.5, 1.20, f"Fisher's Exact Test\n{p_fmt}  ({sig})",
                transform=ax.transAxes, ha="center", fontsize=10, color=color,
                bbox=dict(boxstyle="round,pad=0.35", fc="#FDEDEC" if sig != "ns" else "#F2F3F4",
                           ec=color))
        ax.set_xlabel("Error bars = 95% Wilson CI", fontsize=9)

    plt.tight_layout()
    plt.savefig(save_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {save_path}")


# ── Summary printout ─────────────────────────────────────────────────────────
def print_summary(df: pd.DataFrame, df_clean: pd.DataFrame,
                  df_primary: pd.DataFrame, stats_mode: str = "both"):
    print("\n" + "="*62)
    print("  ANALYSIS SUMMARY")
    print("="*62)
    print(f"  Total flasks:              {len(df)}")
    cont_y = len(df[df['Contamination?'] == 'Y'])
    cont_m = len(df[df['Contamination?'] == 'M'])
    print(f"  Contaminated excl. (Y):    {cont_y}")
    print(f"  Maybe-contam. excl. (M):   {cont_m}")
    og_n      = (df['Immortalized?'] == 'ongoing').sum()
    harv_n    = (df['Immortalized?'] == 'harvested').sum()
    print(f"  Ongoing/unresolved:        {og_n}  (excluded from immortalization denom.)")
    print(f"  Harvested:                 {harv_n}  (excluded from immortalization denom.)")
    print(f"  Clean flasks analyzed:     {len(df_clean)}")
    print()

    print("  PRIMARY CULTURE RATE (clean flasks)")
    pr = rate_stats(df_clean, "Drosophila_species", "primary").reindex(SPECIES_ORDER)
    for sp, row in zip(SPECIES_ORDER, pr.itertuples()):
        print(f"    {sp:22s}  {row.successes:.0f}/{row.trials:.0f}  "
              f"({row.rate*100:.1f}%)  95% CI [{row.ci_lo*100:.1f}–{row.ci_hi*100:.1f}%]")
    ct = [[r.successes, r.trials - r.successes] for r in pr.itertuples()]
    report_categorical_test(ct, "primary culture rate by species", stats_mode)
    print()

    print("  IMMORTALIZATION RATE (of primary cultures)")
    ir = rate_stats(df_primary, "Drosophila_species", "immortalized_flag").reindex(SPECIES_ORDER)
    for sp, row in zip(SPECIES_ORDER, ir.itertuples()):
        print(f"    {sp:22s}  {row.successes:.0f}/{row.trials:.0f}  "
              f"({row.rate*100:.1f}%)  95% CI [{row.ci_lo*100:.1f}–{row.ci_hi*100:.1f}%]")
    ct2 = [[r.successes, r.trials - r.successes] for r in ir.itertuples()]
    report_categorical_test(ct2, "immortalization rate by species", stats_mode)
    print()

    print("  WITHIN-SPECIES ASSESSMENT")
    print("    D. melanogaster: Lines have very high + uniform primary rates;")
    print("      only ub:Khc immortalizes (7/16). No formal test possible —")
    print("      insufficient variance and n.")
    print("    D. simulans: Three lines (w[501], Riv84, Mer23) with variable")
    print("      rates; formal testing feasible for wRi pair (Riv84 vs Mer23).")
    print("    D. willistoni: Single line only — no within-species comparison.")
    print("      (2 OG/ongoing flasks excluded from immortalization denominator.)")
    print("    F1 hybrid: Single line only — no within-species comparison.")
    print()

    print("  wRi STRAIN COMPARISON (Riv84 vs Mer23, D. simulans)")
    df_sim   = df_clean[df_clean["Drosophila_species"] == "simulans"]
    df_sim_p = df_primary[df_primary["Drosophila_species"] == "simulans"]
    riv   = df_sim[df_sim["Fly_Line"] == "D._simulans_Riv84"]
    mer   = df_sim[df_sim["Fly_Line"].str.contains("Mer23")]
    riv_p = df_sim_p[df_sim_p["Fly_Line"] == "D._simulans_Riv84"]
    mer_p = df_sim_p[df_sim_p["Fly_Line"].str.contains("Mer23")]

    riv_k,  riv_n  = int(riv["primary"].sum()), len(riv)
    mer_k,  mer_n  = int(mer["primary"].sum()), len(mer)
    riv_ik, riv_in = int(riv_p["immortalized_flag"].sum()), len(riv_p)
    mer_ik, mer_in = int(mer_p["immortalized_flag"].sum()), len(mer_p)

    print(f"    Primary culture: Riv84 {riv_k}/{riv_n}, Mer23 {mer_k}/{mer_n}")
    report_categorical_test([[riv_k, riv_n - riv_k], [mer_k, mer_n - mer_k]],
                            "primary culture Riv84 vs Mer23", stats_mode)
    print(f"    Immortalization: Riv84 {riv_ik}/{riv_in}, Mer23 {mer_ik}/{mer_in}")
    report_categorical_test([[riv_ik, riv_in - riv_ik], [mer_ik, mer_in - mer_ik]],
                            "immortalization Riv84 vs Mer23", stats_mode)
    print()

    # ── Time-to-immortalization ───────────────────────────────────────────────
    tti = df[df["tti"].notna()].copy()
    print("  TIME TO IMMORTALIZATION (days: primary culture → soft splits)")
    print(f"  {'Group':<30}  {'n':>3}  {'Mean':>7}  {'Median':>7}  {'SD':>7}  Range")
    print(f"  {'-'*30}  {'---':>3}  {'-------':>7}  {'-------':>7}  {'-------':>7}  -----")

    def _desc(vals):
        return (len(vals), vals.mean(), vals.median(),
                vals.std(ddof=1) if len(vals) > 1 else float("nan"),
                vals.min(), vals.max())

    # By species
    print()
    print("  By species:")
    for sp, grp in tti.groupby("Drosophila_species"):
        n, mn, med, sd, lo, hi = _desc(grp["tti"])
        print(f"    {'D. ' + sp:<28}  {n:>3}  {mn:>7.1f}  {med:>7.1f}  "
              f"{sd:>7.1f}  {lo:.0f}–{hi:.0f}")
    sp_groups = [g["tti"].values for _, g in tti.groupby("Drosophila_species")]
    if len(sp_groups) == 2:
        mel_v, sim_v = sp_groups
        U, p_mw = stats.mannwhitneyu(mel_v, sim_v, alternative="two-sided")
        t,  p_t  = stats.ttest_ind(mel_v, sim_v, equal_var=False)
        sig = "***" if p_mw < 0.001 else ("**" if p_mw < 0.01 else
              ("*"   if p_mw < 0.05  else "ns"))
        print(f"    → Mann-Whitney U = {U:.0f}, p = {p_mw:.4f}  ({sig})")
        print(f"    → Welch t-test   t = {t:.3f}, p = {p_t:.4f}  (secondary)")
        if p_mw >= 0.05:
            print("    → No significant difference detected (small n; interpret with caution)")

    # By fly line
    print()
    print("  By fly line:")
    for line, grp in tti.groupby("Fly_Line"):
        n, mn, med, sd, lo, hi = _desc(grp["tti"])
        label = LINE_SHORT.get(line, line)
        print(f"    {label:<28}  {n:>3}  {mn:>7.1f}  {med:>7.1f}  "
              f"{'  n/a' if n < 2 else f'{sd:>7.1f}'}  {lo:.0f}–{hi:.0f}")
    # Within simulans: Mer23 vs w[501] — report with caveat
    mer23 = tti[tti["Fly_Line"].str.contains("Mer23")]["tti"].values
    w501  = tti[tti["Fly_Line"] == "Dsim_w[501]"]["tti"].values
    if len(mer23) >= 2 and len(w501) >= 2:
        U2, p2 = stats.mannwhitneyu(mer23, w501, alternative="two-sided")
        sig2 = "***" if p2 < 0.001 else ("**" if p2 < 0.01 else
               ("*"   if p2 < 0.05  else "ns"))
        print(f"    → Mer23 vs w[501] (both D. simulans): Mann-Whitney U = {U2:.0f}, "
              f"p = {p2:.4f}  ({sig2})")
        print(f"       ⚠  Very small n (Mer23 n={len(mer23)}, w[501] n={len(w501)}); "
              f"interpret with caution")
    else:
        print("    → Formal line comparison not possible: n < 2 in one or more groups")

    # By Wolbachia — descriptive only (mirrors fly line 1:1 here)
    print()
    print("  By Wolbachia infection:")
    for wolb, grp in tti.groupby("Wolbachia infection"):
        n, mn, med, sd, lo, hi = _desc(grp["tti"])
        print(f"    {wolb:<28}  {n:>3}  {mn:>7.1f}  {med:>7.1f}  "
              f"{'  n/a' if n < 2 else f'{sd:>7.1f}'}  {lo:.0f}–{hi:.0f}")
    print("    → Wolbachia groups mirror fly lines 1:1 here; no additional test performed")

    print("="*62 + "\n")


# ── PLOT 7: Outcome overview with Harvested ───────────────────────────────────
OUT_PAL = {
    "Immortalized": "#009988",
    "Primary Only": "#EE7733",
    "Harvested":    "#DDAA33",
    "Ongoing":      "#AA4499",
    "No Primary":   "#BBBBBB",
    "Contaminated": "#332288",
}
OUT_ORDER = ["Immortalized", "Harvested", "Ongoing", "Primary Only", "No Primary", "Contaminated"]

LINE_SP = {
    "ub:Khc": "melanogaster", "mei-P26[1]": "melanogaster",
    "PlasMem:GFP": "melanogaster", "Jupiter:GFP": "melanogaster",
    "Riv84": "simulans", "Mer23": "simulans", "w[501]": "simulans",
    "genome line": "willistoni", "F1 hybrid": "F1hybrid:melxsim",
}

def _flask_outcome(row):
    imm  = str(row["Immortalized?"]).strip().lower()
    cont = str(row["Contamination?"]).strip()
    prim = str(row["Primary_Culture?"]).strip()
    if cont in ("Y", "M"):    return "Contaminated"
    if imm == "harvested":    return "Harvested"
    if imm == "ongoing":      return "Ongoing"
    if prim == "N":           return "No Primary"
    if imm == "y":            return "Immortalized"
    return "Primary Only"

def plot_outcomes(df: pd.DataFrame, save_path: str):
    df = df.copy()
    df["_outcome"] = df.apply(_flask_outcome, axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle("Flask Outcomes Including 'Harvested' Category",
                 fontsize=13, fontweight="bold")

    line_order = ["ub:Khc", "mei-P26[1]", "PlasMem:GFP", "Jupiter:GFP",
                  "Riv84", "Mer23", "w[501]", "genome line", "F1 hybrid"]

    for ax, (group_col, groups, labels, title) in zip(axes, [
        ("Drosophila_species", SPECIES_ORDER, SPECIES_LABELS, "By Species"),
        ("line_short",         line_order,    line_order,      "By Fly Line  (total in species colour)"),
    ]):
        for xi, (grp, lbl) in enumerate(zip(groups, labels)):
            sub = df[df[group_col] == grp]
            bottom = 0
            for out in OUT_ORDER:
                n = (sub["_outcome"] == out).sum()
                if n > 0:
                    ax.bar(xi, n, bottom=bottom, color=OUT_PAL[out],
                           edgecolor="white", linewidth=0.5,
                           label=out if xi == 0 else "_nolegend_")
                    if n >= 2:
                        ax.text(xi, bottom + n/2, str(n), ha="center", va="center",
                                fontsize=8,
                                color="white" if out in ("Immortalized","Contaminated","Ongoing") else "#333")
                    bottom += n
            sp_col = SPECIES_PALETTE.get(LINE_SP.get(grp, grp), "#555")
            ax.text(xi, bottom + 0.3, str(int(bottom)),
                    ha="center", va="bottom", fontsize=8.5, fontweight="bold",
                    color=sp_col if group_col == "line_short" else "#333")

        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(labels, rotation=25 if group_col == "line_short" else 0,
                           ha="right" if group_col == "line_short" else "center")
        ax.set_ylabel("Number of flasks")
        ax.set_title(title)

    axes[0].legend(title="Outcome", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    harv_n = (df["_outcome"] == "Harvested").sum()
    axes[1].text(0.99, 0.99, f"Harvested total: {harv_n}",
                 transform=axes[1].transAxes, ha="right", va="top", fontsize=8.5,
                 bbox=dict(boxstyle="round,pad=0.3", fc="#FEF3CD", ec="#DDAA33"))

    plt.tight_layout()
    plt.savefig(save_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {save_path}")


# ── PLOT 8: Time-to-immortalization by species, line, Wolbachia ───────────────
def plot_tti_overview(df: pd.DataFrame, save_path: str):
    tti = df[df["tti"].notna()].copy()
    jitter = 0.12

    fig, axes = plt.subplots(1, 3, figsize=(15, 6))
    fig.suptitle(
        "Time to Immortalization (days from primary culture to weekly/biweekly soft splits)",
        fontsize=12, fontweight="bold")

    groups = [
        ("Drosophila_species",
         [s for s in SPECIES_ORDER if s in tti["Drosophila_species"].values],
         SPECIES_PALETTE, "By Species",
         lambda s: "D. mel." if s=="melanogaster" else ("D. sim." if s=="simulans" else s)),
        ("line_short",
         [l for l in ["ub:Khc","mei-P26[1]","PlasMem:GFP","Jupiter:GFP",
                       "Riv84","Mer23","w[501]","genome line","F1 hybrid"]
          if l in tti["line_short"].values],
         LINE_SHORT, "By Fly Line", lambda l: l),
        ("Wolbachia infection",
         [w for w in ["wMel_(DB)","wRi","wRi_(Riv84)","wRi_(Mer23)","wWil","Unknown/mixed"]
          if w in tti["Wolbachia infection"].values],
         {}, "By Wolbachia Infection", lambda w: w),
    ]
    WOLB_PAL = {"wMel_(DB)":"#0077BB","wRi":"#BB2200","wRi_(Riv84)":"#BB2200",
                "wRi_(Mer23)":"#EE6644","wWil":"#228833","Unknown/mixed":"#999999"}

    for ax, (col, order, pal, title, labelfn) in zip(axes, groups):
        for xi, grp in enumerate(order):
            sub = tti[tti[col] == grp]["tti"]
            col_pal = WOLB_PAL if title.startswith("By Wolb") else pal
            color = col_pal.get(grp, SPECIES_PALETTE.get(grp, "#888"))
            jx = xi + np.random.uniform(-jitter, jitter, size=len(sub))
            ax.scatter(jx, sub, color=color, s=60, zorder=3,
                       edgecolors="white", linewidths=0.6)
            ax.hlines(sub.mean(), xi-0.25, xi+0.25, colors=color, linewidths=2.5, zorder=4)
            ax.text(xi, sub.max() + 8, f"n={len(sub)}\nμ={sub.mean():.0f}d",
                    ha="center", va="bottom", fontsize=8, color=color)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([labelfn(g) for g in order], rotation=20, ha="right")
        ax.set_ylabel("Days to immortalization")
        ax.set_title(title)
        ax.set_ylim(0, tti["tti"].max() + 60)

    plt.tight_layout()
    plt.savefig(save_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {save_path}")


# ── PLOT 9: mel vs sim time-to-immortalization statistical test ───────────────
def plot_tti_test(df: pd.DataFrame, save_path: str):
    tti = df[df["tti"].notna()].copy()
    mel = tti[tti["Drosophila_species"] == "melanogaster"]["tti"].values
    sim = tti[tti["Drosophila_species"] == "simulans"]["tti"].values

    U, p_mw = stats.mannwhitneyu(mel, sim, alternative="two-sided")
    t, p_t  = stats.ttest_ind(mel, sim, equal_var=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6),
                             gridspec_kw={"width_ratios": [1.6, 1]})
    fig.suptitle("D. melanogaster vs D. simulans: Time to Immortalization",
                 fontsize=13, fontweight="bold")

    np.random.seed(42)
    jitter = 0.10
    ax = axes[0]
    for xi, (vals, sp, col) in enumerate([
        (mel, "melanogaster", "#0077BB"),
        (sim, "simulans",     "#BB2200"),
    ]):
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        iqr = q3 - q1
        wlo = max(vals.min(), q1 - 1.5*iqr)
        whi = min(vals.max(), q3 + 1.5*iqr)
        ax.fill_between([xi-0.2, xi+0.2], q1, q3, color=col, alpha=0.18, zorder=1)
        ax.hlines(med, xi-0.2, xi+0.2, colors=col, linewidths=2.5, zorder=3)
        ax.vlines(xi, wlo, whi, colors=col, linewidths=1.2, zorder=2)
        ax.hlines([wlo, whi], xi-0.1, xi+0.1, colors=col, linewidths=1.2, zorder=2)
        jx = xi + np.random.uniform(-jitter, jitter, size=len(vals))
        ax.scatter(jx, vals, color=col, s=70, zorder=4, edgecolors="white", linewidths=0.8)
        ax.text(xi, vals.max() + 15,
                f"n = {len(vals)}\nmean = {vals.mean():.0f}d\nmedian = {np.median(vals):.0f}d",
                ha="center", va="bottom", fontsize=9, color=col)

    ymax = max(mel.max(), sim.max()) + 55
    ax.annotate("", xy=(1, ymax), xytext=(0, ymax),
                arrowprops=dict(arrowstyle="-", color="#333", lw=1.5))
    p_str = f"p = {p_mw:.4f}" if p_mw >= 0.001 else f"p = {p_mw:.2e}"
    sig   = "***" if p_mw < 0.001 else ("**" if p_mw < 0.01 else ("*" if p_mw < 0.05 else "ns"))
    ax.text(0.5, ymax + 5, f"Mann-Whitney U  {p_str}  ({sig})",
            ha="center", va="bottom", fontsize=9.5, color="#333")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["D. melanogaster\n(ub:Khc only)", "D. simulans\n(Mer23 + w[501])"],
                       fontsize=10)
    ax.set_ylabel("Days from primary culture\nto weekly/biweekly soft splits")
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0, ymax + 40)
    ax.set_title("Distribution & Test")

    ax2 = axes[1]
    ax2.axis("off")
    rows = [
        ["", "D. melanogaster", "D. simulans"],
        ["n",             str(len(mel)),               str(len(sim))],
        ["Mean (days)",   f"{mel.mean():.1f}",         f"{sim.mean():.1f}"],
        ["Median (days)", f"{np.median(mel):.1f}",     f"{np.median(sim):.1f}"],
        ["SD",            f"{mel.std(ddof=1):.1f}",    f"{sim.std(ddof=1):.1f}"],
        ["Range",         f"{mel.min()}–{mel.max()}",  f"{sim.min()}–{sim.max()}"],
        ["", "", ""],
        ["Test", "Mann-Whitney U", ""],
        ["U statistic", f"{U:.0f}", ""],
        ["p-value (2-sided)", p_str, ""],
        ["Result", sig, ""],
        ["", "", ""],
        ["Welch t-test", f"t = {t:.3f}", f"p = {p_t:.4f}"],
    ]
    col_x = [0.02, 0.38, 0.72]
    y = 0.97
    for ri, row in enumerate(rows):
        for ci, (cell, cx) in enumerate(zip(row, col_x)):
            weight = "bold" if ri == 0 or (ri > 0 and ci == 0 and cell) else "normal"
            color  = ("#0077BB" if ci == 1 and ri == 0
                      else "#BB2200" if ci == 2 and ri == 0 else "#222")
            ax2.text(cx, y, cell, transform=ax2.transAxes,
                     fontsize=9, va="top", ha="left", fontweight=weight, color=color)
        y -= 0.072

    ax2.text(0.02, y - 0.03,
             "⚠  Small sample sizes (n = 5, 6).\nMann-Whitney U is non-parametric and\npreferred; interpret with caution.",
             transform=ax2.transAxes, fontsize=8, va="top", color="#7D6608",
             bbox=dict(boxstyle="round,pad=0.4", fc="#FEF9E7", ec="#F0B27A"))
    ax2.set_title("Summary Statistics", fontsize=11)

    plt.tight_layout()
    plt.savefig(save_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {save_path}")


# ── mel vs sim subset comparison ─────────────────────────────────────────────
def print_mel_sim_comparison(df_clean: pd.DataFrame, df_primary: pd.DataFrame,
                             stats_mode: str = "both"):
    """
    Print primary culture and immortalization rate comparisons restricted to
    D. melanogaster and D. simulans only, for use in protocol paper contexts
    where willistoni and F1 hybrid data are not included.

    Dropping to a 2×2 table means Fisher's exact is available for both tests,
    and expected cell counts can be verified directly.
    """
    MEL_SIM = ["melanogaster", "simulans"]
    MEL_SIM_LABELS = ["D. melanogaster", "D. simulans"]

    dc = df_clean[df_clean["Drosophila_species"].isin(MEL_SIM)].copy()
    dp = df_primary[df_primary["Drosophila_species"].isin(MEL_SIM)].copy()

    print("\n" + "="*62)
    print("  MEL vs SIM COMPARISON (protocol paper subset)")
    print("="*62)

    # Primary culture rate
    print("  PRIMARY CULTURE RATE (clean flasks, mel & sim only)")
    pr = rate_stats(dc, "Drosophila_species", "primary").reindex(MEL_SIM)
    for sp, lbl, row in zip(MEL_SIM, MEL_SIM_LABELS, pr.itertuples()):
        print(f"    {lbl:22s}  {row.successes:.0f}/{row.trials:.0f}  "
              f"({row.rate*100:.1f}%)  95% CI [{row.ci_lo*100:.1f}–{row.ci_hi*100:.1f}%]")
    ct_pc = [[r.successes, r.trials - r.successes] for r in pr.itertuples()]
    report_categorical_test(ct_pc, "primary culture rate mel vs sim", stats_mode)
    print()

    # Immortalization rate
    print("  IMMORTALIZATION RATE (of primary cultures, mel & sim only)")
    print("  (harvested and ongoing excluded from denominator)")
    ir = rate_stats(dp, "Drosophila_species", "immortalized_flag").reindex(MEL_SIM)
    for sp, lbl, row in zip(MEL_SIM, MEL_SIM_LABELS, ir.itertuples()):
        # Also report how many were excluded from the denominator
        n_formed = (dc[dc["Drosophila_species"] == sp]["primary"] == 1).sum()
        n_excluded = int(n_formed - row.trials)
        excl_note = f"  [{n_excluded} harvested/ongoing excluded]" if n_excluded else ""
        print(f"    {lbl:22s}  {row.successes:.0f}/{row.trials:.0f}  "
              f"({row.rate*100:.1f}%)  95% CI [{row.ci_lo*100:.1f}–{row.ci_hi*100:.1f}%]"
              f"{excl_note}")
    ct_im = [[r.successes, r.trials - r.successes] for r in ir.itertuples()]
    report_categorical_test(ct_im, "immortalization rate mel vs sim", stats_mode)
    print("="*62 + "\n")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():

    # ── Argument parsing ──────────────────────────────────────────────────────
    # Usage: python cell_line_analysis_v4.py <input_csv> <output_dir> [--stats chi2|fisher|both]
    #   --stats chi2    report chi-square only
    #   --stats fisher  report Fisher's exact only (falls back to chi-square for >2×2 tables)
    #   --stats both    report both (default)
    if len(sys.argv) < 3:
        print("Error: Please provide a file path and output directory.")
        print("Usage: python cell_line_analysis_v4.py <input_csv> <output_dir> [--stats chi2|fisher|both]")
        sys.exit(1)

    data_path  = sys.argv[1]
    out_dir    = sys.argv[2].rstrip("/") + "/"
    stats_mode = "both"

    # Parse optional --stats flag
    if "--stats" in sys.argv:
        idx = sys.argv.index("--stats")
        if idx + 1 >= len(sys.argv):
            print("Error: --stats requires an argument: chi2, fisher, or both")
            sys.exit(1)
        stats_mode = sys.argv[idx + 1].lower()
        if stats_mode not in ("chi2", "fisher", "both"):
            print(f"Error: --stats must be chi2, fisher, or both (got '{stats_mode}')")
            sys.exit(1)

    try:
        df, df_clean, df_primary = load_and_clean(data_path)
    except FileNotFoundError:
        print(f"Error: The file '{data_path}' was not found.")
        sys.exit(1)

    # Add line_short column needed by plot_outcomes
    df["line_short"] = df["Fly_Line"].map(LINE_SHORT)

    print_summary(df, df_clean, df_primary, stats_mode)
    print_mel_sim_comparison(df_clean, df_primary, stats_mode)
    plot_overview(df, df_clean, df_primary, f"{out_dir}01_overview.pdf")
    plot_between_species(df_clean, df_primary, f"{out_dir}02_between_species.pdf")
    plot_within_species(df_clean, df_primary, f"{out_dir}03_within_species.pdf")
    plot_wri_comparison(df_clean, df_primary, f"{out_dir}04_wRi_Riv84_vs_Mer23.pdf")
    plot_outcomes(df,     f"{out_dir}07_outcomes_with_harvested.pdf")
    plot_tti_overview(df, f"{out_dir}08_time_to_immortalization.pdf")
    plot_tti_test(df,     f"{out_dir}09_mel_vs_sim_tti_test.pdf")
    print("All plots saved.")

if __name__ == "__main__":
    main()
