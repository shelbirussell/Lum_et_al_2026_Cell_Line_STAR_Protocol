#!/usr/bin/env python3
"""
02_plot_stable_titer.py
=======================
Plots Wolbachia titer for stable cell-line samples.

Reads stable_infections_23C.csv (or any stable CSV).
White-food rows (fly-derived samples) are pulled out and plotted
separately with a 'fly_whitefood_' filename prefix, then excluded
from the main cell-line plots.

Usage:
  python 02_plot_stable_titer.py \
      --input   ./output/stable_infections_23C.csv \
      --out_dir ./output/plots_stable_23C
"""

import argparse
import os
import re
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning)

# ── Media ──────────────────────────────────────────────────────────────────
MEDIA_LINESTYLE = {
    "Shields/Sang M3": "-",
    "Schneiders":      "--",
    "Transitioning":   (0, (3, 1, 1, 1)),
    "White food":      ":",
    "Unknown":         ":",
}
MEDIA_LABEL = {
    "Shields/Sang M3": "Shields & Sang M3",
    "Schneiders":      "Schneider's",
    "Transitioning":   "Transitioning (mixed ratio)",
    "White food":      "White food (fly diet)",
    "Unknown":         "Unknown",
}

def classify_media(media_str):
    if pd.isna(media_str) or str(media_str).strip().lower() in {"", "nan"}:
        return "Unknown"
    m = str(media_str).strip()
    if re.search(r"\d:\d", m) or re.search(r"M3.*Schneider|Schneider.*M3", m, re.IGNORECASE):
        return "Transitioning"
    ml = m.lower()
    if "white food" in ml:
        return "White food"
    if "shield" in ml or ml in {"m3"}:
        return "Shields/Sang M3"
    if "schneider" in ml:
        return "Schneiders"
    return "Unknown"

# ── Titer columns ──────────────────────────────────────────────────────────
STRAIN_TITER = {
    "wMel":          ("wmel_titer",   "wMel titer\n(Wolbachia / Drosophila genomes)", True),
    "wRi":           ("wri_titer",    "wRi titer\n(Wolbachia / Drosophila genomes)",  True),
    "wRi_Merrill23": ("wri_titer",    "wRi titer\n(Wolbachia / Drosophila genomes)",  True),
    "wWil":          ("wwil_titer",   "wWil titer\n(Wolbachia / Drosophila genomes)", True),
    # DOX/uninf: computed titer from Wolbachia/host mean depth;
    # NaN Wolbachia depth → titer = 0 (no Wolbachia detected above threshold).
    # Column "_uninf_titer" is added by add_host_depth_col(); linear axis.
    "DOX":           ("_uninf_titer", "Titer (Wolbachia / Drosophila genomes)\n[0 = below detection]", False),
    "uninf":         ("_uninf_titer", "Titer (Wolbachia / Drosophila genomes)\n[0 = below detection]", False),
    "uninf.":        ("_uninf_titer", "Titer (Wolbachia / Drosophila genomes)\n[0 = below detection]", False),
    "uninfected":    ("_uninf_titer", "Titer (Wolbachia / Drosophila genomes)\n[0 = below detection]", False),
}

def add_host_depth_col(df):
    df = df.copy()
    is_dsim = df["Cell Line"].fillna("").str.lower().str.startswith("dsim")
    df["_host_depth"] = np.where(is_dsim, df["dsim_mean_depth"], df["dmel_mean_depth"])

    # _uninf_titer: for DOX/uninf lines, compute Wolbachia/host ratio using
    # whichever Wolbachia mean_depth column has a value (any strain could appear
    # as low-level contamination). Where Wolbachia depth is NaN or zero,
    # titer = 0 (below detection). Linear axis is used for these.
    host = pd.to_numeric(df["_host_depth"], errors="coerce")
    wmel = pd.to_numeric(df.get("wmel_mean_depth", pd.Series(np.nan, index=df.index)), errors="coerce")
    wri  = pd.to_numeric(df.get("wri_mean_depth",  pd.Series(np.nan, index=df.index)), errors="coerce")
    wwil = pd.to_numeric(df.get("wwil_mean_depth", pd.Series(np.nan, index=df.index)), errors="coerce")
    # Take the maximum detected Wolbachia depth across all strains
    wol  = pd.concat([wmel, wri, wwil], axis=1).max(axis=1)
    # Where Wolbachia depth is NaN/0 → titer 0; otherwise divide by host depth
    titer = np.where(
        wol.isna() | (wol == 0) | host.isna() | (host == 0),
        0.0,
        wol / host
    )
    df["_uninf_titer"] = titer
    return df

def pick_titer_info(infection, df_cols):
    info = STRAIN_TITER.get(str(infection).strip())
    if info is None:
        return None
    col, ylabel, use_log = info
    if col in ("_host_depth", "_uninf_titer") or col in df_cols:
        return col, ylabel, use_log
    return None

# ── Replicate letter ───────────────────────────────────────────────────────
def extract_rep_letter(row):
    """
    Extract the replicate letter from Split (DDMMRep).
    If no letter is found (date-only splits, no split value, or single-flask
    maintainers like early SLR and JJ), return "A" — these are all the same
    continuous flask, just recorded without a letter designation.
    """
    split_str = row["Split (DDMMRep)"] if "Split (DDMMRep)" in row.index else ""

    if pd.isna(split_str) or str(split_str).strip().lower() in {"", "nan"}:
        return "A"
    s = str(split_str).strip().strip("\"'")
    m = re.search(r"[A-Z][A-Z'\"]*$", s, re.IGNORECASE)
    return m.group(0).upper() if m else "A"

# ── Colours ────────────────────────────────────────────────────────────────
REP_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    "#bcbd22", "#17becf", "#aec7e8", "#ffbb78",
]

# ── Drawing ────────────────────────────────────────────────────────────────
def plot_rep_line(ax, x, y, media_classes, color, use_log):
    if len(x) == 0:
        return set()
    media_seen = set()
    for i in range(len(x) - 1):
        mc = media_classes.iloc[i]
        ls = MEDIA_LINESTYLE.get(mc, ":")
        media_seen.add(mc)
        if use_log:
            # Skip segments involving zero values on log axis
            if y.iloc[i] > 0 and y.iloc[i+1] > 0:
                ax.semilogy([x.iloc[i], x.iloc[i+1]], [y.iloc[i], y.iloc[i+1]],
                            color=color, linestyle=ls, linewidth=1.6)
        else:
            ax.plot([x.iloc[i], x.iloc[i+1]], [y.iloc[i], y.iloc[i+1]],
                    color=color, linestyle=ls, linewidth=1.6)
    media_seen.add(media_classes.iloc[-1])
    # Markers: on log axis use downward triangle for zeros, circle otherwise
    if use_log:
        pos = y > 0
        if pos.any():
            ax.semilogy(x[pos], y[pos], color=color, linestyle="none",
                        marker="o", markersize=5)
        if (~pos).any():
            # Plot zeros at axis minimum as downward triangles
            ymin = ax.get_ylim()[0] if ax.get_ylim()[0] > 0 else 1e-4
            ax.semilogy(x[~pos], [ymin] * (~pos).sum(), color=color,
                        linestyle="none", marker="v", markersize=6, alpha=0.7)
    else:
        ax.plot(x, y, color=color, linestyle="none", marker="o", markersize=5)
        # Mark zeros with a downward triangle on top of the circle
        zeros = y == 0
        if zeros.any():
            ax.plot(x[zeros], y[zeros], color=color, linestyle="none",
                    marker="v", markersize=8, alpha=0.7)
    return media_seen

# ── Load ───────────────────────────────────────────────────────────────────
def load_data(path):
    df = pd.read_csv(path, low_memory=False)
    df["_date"] = pd.to_datetime(
        df["Date Collected (MM/DD/YYYY)"], errors="coerce", dayfirst=False)
    df.loc[df["_date"].dt.year.fillna(9999) < 2000, "_date"] = pd.NaT
    df["_media_class"] = df["Media Type"].apply(classify_media)
    df["_rep_letter"]  = df.apply(extract_rep_letter, axis=1)
    df = add_host_depth_col(df)
    return df

# ── Core plot function ─────────────────────────────────────────────────────
def _make_plots(df, out_dir, prefix, title_tag=""):
    """
    Generate one PNG per Cell Line × Infection from df.
    prefix  : filename prefix (e.g. 'stable', 'fly_whitefood')
    title_tag: appended to figure suptitle (e.g. ' — White food / fly')
    """
    cell_lines = sorted(df["Cell Line"].dropna().unique())
    infections = sorted(df["Infection"].dropna().unique())

    for cell_line in cell_lines:
        for infection in infections:
            titer_info = pick_titer_info(infection, df.columns)
            if titer_info is None:
                continue
            titer_col, ylabel, use_log = titer_info

            # For _uninf_titer, include zeros (they represent below-detection samples)
            if titer_col == "_uninf_titer":
                sub = df[
                    (df["Cell Line"].fillna("") == cell_line) &
                    (df["Infection"].fillna("") == infection) &
                    df[titer_col].notna()
                ].copy()
            else:
                sub = df[
                    (df["Cell Line"].fillna("") == cell_line) &
                    (df["Infection"].fillna("") == infection) &
                    df[titer_col].notna() &
                    (df[titer_col] > 0)
                ].copy()

            if sub.empty:
                continue

            persons = sorted(sub["Initials"].dropna().unique())
            if not persons:
                continue

            n_persons = len(persons)
            fig, axes = plt.subplots(
                1, n_persons,
                figsize=(max(5 * n_persons, 6), 5),
                sharey=True, squeeze=False
            )
            fig.suptitle(
                f"{cell_line}  ×  {infection}{title_tag}",
                fontsize=13, fontweight="bold", y=1.02
            )

            any_data   = False
            media_seen = set()

            for ax, person in zip(axes[0], persons):
                psub = sub[sub["Initials"] == person].sort_values("_date")
                rep_letters = sorted(psub["_rep_letter"].unique())
                color_map   = {r: REP_COLORS[i % len(REP_COLORS)]
                               for i, r in enumerate(rep_letters)}

                for rep in rep_letters:
                    rsub = psub[psub["_rep_letter"] == rep].dropna(
                        subset=["_date", titer_col]).sort_values("_date")
                    rsub = rsub[rsub[titer_col] > 0]
                    if rsub.empty:
                        continue

                    seen = plot_rep_line(ax,
                                        rsub["_date"], rsub[titer_col],
                                        rsub["_media_class"],
                                        color=color_map[rep],
                                        use_log=use_log)
                    media_seen.update(seen)
                    any_data = True

                ax.set_title(person, fontsize=11)
                ax.set_xlabel("Collection date", fontsize=9)
                ax.tick_params(axis="x", rotation=30, labelsize=8)
                ax.grid(True, which="both", linestyle=":", alpha=0.4)

                if rep_letters:
                    rep_handles = [
                        mlines.Line2D([], [], color=color_map[r],
                                      marker="o", markersize=5, linewidth=1.6,
                                      label=f"Rep {r}" if r else "Rep A")
                        for r in rep_letters
                    ]
                    ax.legend(handles=rep_handles, fontsize=7,
                              loc="upper left", title="Replicate", title_fontsize=7)

            if not any_data:
                plt.close(fig)
                continue

            axes[0][0].set_ylabel(ylabel, fontsize=10)
            if use_log:
                for ax in axes[0]:
                    ax.set_ylim(1e-3, 1e3)

            media_handles = [
                mlines.Line2D([], [], color="grey",
                              linestyle=MEDIA_LINESTYLE[m], linewidth=1.5,
                              label=MEDIA_LABEL.get(m, m))
                for m in sorted(media_seen) if m in MEDIA_LINESTYLE
            ]
            # tight_layout first, then add the figure legend so it isn't displaced
            plt.tight_layout()

            if media_handles:
                fig.legend(handles=media_handles,
                           loc="upper right",
                           fontsize=9, title="Media", title_fontsize=9,
                           framealpha=0.9)

            safe = re.sub(r"[^\w\-]", "_", f"{cell_line}_{infection}")
            out_path = os.path.join(out_dir, f"{prefix}_{safe}.pdf")
            plt.savefig(out_path, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved: {out_path}")


def plot_stable(df, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    # Exclude GP from plots
    df = df[df["Initials"].fillna("").str.upper() != "GP"].copy()

    # Premix / pre-sv / post-sv-stable rows are baseline samples, not continuous
    # replicates — exclude them so they don't create a spurious "" rep group.
    premix_types = {"premix", "post-sv-stable"}
    premix_mask = df["Type"].fillna("").str.lower().isin(premix_types)
    if premix_mask.any():
        print(f"  Excluding {premix_mask.sum()} premix/pre-sv rows from line plots")
        df = df[~premix_mask].copy()

    # White-food rows are fly-derived — plot separately, then exclude from main plots
    wf_mask = df["_media_class"] == "White food"
    if wf_mask.any():
        print(f"  Plotting {wf_mask.sum()} white-food (fly) rows separately…")
        _make_plots(df[wf_mask], out_dir,
                    prefix="fly_whitefood",
                    title_tag="  [White food / fly]")

    _make_plots(df[~wf_mask], out_dir, prefix="stable")


def main():
    parser = argparse.ArgumentParser(description="Plot stable-infection titer over time")
    parser.add_argument("--input",   default="./output/stable_infections_23C.csv")
    parser.add_argument("--out_dir", default="./output/plots_stable_23C")
    args = parser.parse_args()
    if not os.path.isfile(args.input):
        sys.exit(f"ERROR: {args.input} not found")
    print(f"Loading {args.input}…")
    df = load_data(args.input)
    print(f"  {len(df)} rows, {df['Cell Line'].nunique()} cell lines, "
          f"{df['Infection'].nunique()} infections")
    print("Generating plots…")
    plot_stable(df, args.out_dir)
    print("Done.")

if __name__ == "__main__":
    main()
