#!/usr/bin/env python3
"""
04_plot_freeze_thaw.py
======================
Plots freeze-thaw revival titer over time.

Layout:
  - One figure per infection state (DOX, wMel, etc.)
  - One subplot column per person (Initials)
  - Both cell lines (JW18, S2, …) are on the same axes:
      JW18  → solid line   (─)
      S2    → dashed line  (--)
      other → dotted line  (·)
  - Each replicate letter gets its own colour, shared across cell lines
    so Rep A JW18 and Rep A S2 are the same colour (different dash style)
  - Media type is indicated by marker shape:
      M3 / Shields & Sang M3  → circle  (o)
      Schneider's              → square  (s)
      Unknown                 → diamond (D)
  - Y-axis: titer (log10) for infected strains;
            host genome mean depth (linear) for DOX / uninf
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
M3_KEYWORDS     = {"m3", "shields", "s&sm3", "shields/sang"}
SCHNDR_KEYWORDS = {"schneider", "schneid"}

def classify_media(media_str):
    if pd.isna(media_str) or str(media_str).strip().lower() in {"", "nan"}:
        return "Unknown"
    m = str(media_str).strip().lower()
    if any(k in m for k in M3_KEYWORDS):
        return "Shields/Sang M3"
    if any(k in m for k in SCHNDR_KEYWORDS):
        return "Schneiders"
    return "Unknown"

MEDIA_MARKER = {"M3": "o", "Schneiders": "s", "Unknown": "D"}
MEDIA_LABEL  = {
    "M3":         "Shields & Sang M3",
    "Schneiders": "Schneider's",
    "Unknown":    "Unknown / mixed media",
}

# ── Cell-line line style ───────────────────────────────────────────────────
CELLLINE_LS    = {"JW18": "-", "S2": "--"}
CELLLINE_LABEL = {"JW18": "JW18", "S2": "S2"}

def cellline_ls(cl):
    return CELLLINE_LS.get(cl, ":")

# ── Titer columns ──────────────────────────────────────────────────────────
STRAIN_TITER = {
    "wMel":          ("wmel_titer",  "wMel titer\n(Wolbachia / Drosophila genomes)", True),
    "wRi":           ("wri_titer",   "wRi titer\n(Wolbachia / Drosophila genomes)",  True),
    "wRi_Merrill23": ("wri_titer",   "wRi titer\n(Wolbachia / Drosophila genomes)",  True),
    "wWil":          ("wwil_titer",  "wWil titer\n(Wolbachia / Drosophila genomes)", True),
    "DOX":           ("_host_depth", "Host genome mean depth\n(cell-health proxy)",  False),
    "uninf":         ("_host_depth", "Host genome mean depth\n(cell-health proxy)",  False),
    "uninf.":        ("_host_depth", "Host genome mean depth\n(cell-health proxy)",  False),
    "uninfected":    ("_host_depth", "Host genome mean depth\n(cell-health proxy)",  False),
}

def add_host_depth_col(df):
    df = df.copy()
    is_dsim = df["Cell Line"].fillna("").str.lower().str.startswith("dsim")
    df["_host_depth"] = np.where(is_dsim, df["dsim_mean_depth"], df["dmel_mean_depth"])
    return df

def pick_titer_info(infection, df_cols):
    info = STRAIN_TITER.get(str(infection).strip())
    if info is None:
        return None
    col, ylabel, use_log = info
    if col == "_host_depth" or col in df_cols:
        return col, ylabel, use_log
    return None

# ── Replicate letter ───────────────────────────────────────────────────────
def extract_rep_letter(split_str):
    if pd.isna(split_str) or str(split_str).strip().lower() in {"", "nan"}:
        return ""
    s = str(split_str).strip().strip("\"'")
    # JJ-style: letter immediately before '-Thaw'
    m = re.search(r"([A-Z][A-Z']*)-[Tt]haw$", s, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    # AL-style: trailing letter(s) with optional apostrophes
    m = re.search(r"([A-Z]['\u2019\u2018]*)$", s, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    return ""

# ── Colours (one per replicate letter) ────────────────────────────────────
REP_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    "#bcbd22", "#17becf",
]


def load_data(path):
    df = pd.read_csv(path, low_memory=False)
    df["_date"] = pd.to_datetime(
        df["Date Collected (MM/DD/YYYY)"], errors="coerce", dayfirst=False)
    df.loc[df["_date"].dt.year.fillna(9999) < 2000, "_date"] = pd.NaT
    df["_media_class"] = df["Media Type"].apply(classify_media)
    df["_rep_letter"]  = df["Split (DDMMRep)"].apply(extract_rep_letter)
    df = add_host_depth_col(df)
    return df


def plot_series(ax, x, y, media_classes, color, cell_line, use_log):
    """
    Draw one (replicate × cell-line) series.
    Line style encodes cell line; marker shape encodes media type.
    Segments connect consecutive points with the line style of the earlier point.
    """
    if len(x) == 0:
        return
    ls = cellline_ls(cell_line)
    plot_fn = ax.semilogy if use_log else ax.plot

    # Draw connecting segments
    for i in range(len(x) - 1):
        plot_fn([x.iloc[i], x.iloc[i+1]],
                [y.iloc[i], y.iloc[i+1]],
                color=color, linestyle=ls, linewidth=1.6)

    # Draw markers per-point (shape = media)
    for _, row in pd.DataFrame({"x": x, "y": y, "mc": media_classes}).iterrows():
        mk = MEDIA_MARKER.get(row["mc"], "D")
        plot_fn([row["x"]], [row["y"]],
                color=color, linestyle="none",
                marker=mk, markersize=7, markeredgecolor="white",
                markeredgewidth=0.5)


def plot_freeze_thaw(df, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    infections = sorted(df["Infection"].dropna().unique())

    for infection in infections:
        titer_info = pick_titer_info(infection, df.columns)
        if titer_info is None:
            continue
        titer_col, ylabel, use_log = titer_info

        sub = df[
            (df["Infection"].fillna("") == infection) &
            df[titer_col].notna() &
            (df[titer_col] > 0)
        ].copy()

        if sub.empty:
            continue

        persons     = sorted(sub["Initials"].fillna("Unknown").unique())
        n_persons   = len(persons)
        cell_lines  = sorted(sub["Cell Line"].dropna().unique())

        # Build a global replicate colour map across all cell lines / persons
        all_reps  = sorted(sub["_rep_letter"].unique())
        color_map = {r: REP_COLORS[i % len(REP_COLORS)]
                     for i, r in enumerate(all_reps)}

        fig, axes = plt.subplots(
            1, n_persons,
            figsize=(max(5 * n_persons, 6), 5),
            sharey=True, squeeze=False
        )
        fig.suptitle(
            f"Freeze-Thaw Revival  |  {infection}",
            fontsize=13, fontweight="bold", y=1.02
        )

        any_data   = False
        media_seen = set()

        for ax, person in zip(axes[0], persons):
            psub = sub[sub["Initials"].fillna("Unknown") == person]

            for cell_line in cell_lines:
                csub = psub[psub["Cell Line"] == cell_line].sort_values("_date")

                for rep in all_reps:
                    rsub = csub[csub["_rep_letter"] == rep].dropna(
                        subset=["_date", titer_col]).sort_values("_date")
                    rsub = rsub[rsub[titer_col] > 0]
                    if rsub.empty:
                        continue

                    color = color_map[rep]
                    plot_series(ax,
                                rsub["_date"], rsub[titer_col],
                                rsub["_media_class"],
                                color=color,
                                cell_line=cell_line,
                                use_log=use_log)
                    media_seen.update(rsub["_media_class"].unique())
                    any_data = True

            ax.set_title(person, fontsize=11)
            ax.set_xlabel("Collection date", fontsize=9)
            ax.tick_params(axis="x", rotation=30, labelsize=8)
            ax.grid(True, which="both", linestyle=":", alpha=0.4)

        if not any_data:
            plt.close(fig)
            continue

        axes[0][0].set_ylabel(ylabel, fontsize=10)

        # ── Figure-level legend ────────────────────────────────────────
        legend_handles = []

        # Replicate colours
        for rep in all_reps:
            legend_handles.append(
                mlines.Line2D([], [], color=color_map[rep],
                              linewidth=2,
                              label=f"Rep {rep}" if rep else "Rep ?")
            )

        # Spacer
        legend_handles.append(
            mlines.Line2D([], [], color="none", label="")
        )

        # Cell-line line styles
        for cl in cell_lines:
            legend_handles.append(
                mlines.Line2D([], [], color="grey",
                              linestyle=cellline_ls(cl),
                              linewidth=1.8,
                              label=CELLLINE_LABEL.get(cl, cl))
            )

        # Spacer
        legend_handles.append(
            mlines.Line2D([], [], color="none", label="")
        )

        # Media marker shapes
        for mc in sorted(media_seen):
            if mc in MEDIA_MARKER:
                legend_handles.append(
                    mlines.Line2D([], [], color="grey",
                                  linestyle="none",
                                  marker=MEDIA_MARKER[mc], markersize=7,
                                  label=MEDIA_LABEL.get(mc, mc))
                )

        fig.legend(handles=legend_handles,
                   loc="upper right", fontsize=8,
                   title="Legend", title_fontsize=9,
                   bbox_to_anchor=(1.0, 1.0))

        plt.tight_layout()
        safe = re.sub(r"[^\w\-]", "_", infection)
        out_path = os.path.join(out_dir, f"freeze_thaw_{safe}.pdf")
        plt.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {out_path}")

        # ── JJ-only wMel version saved as PDF ─────────────────────────
        if infection == "wMel":
            jj_sub = sub[sub["Initials"].fillna("Unknown") == "JJ"]
            if not jj_sub.empty:
                jj_cell_lines = sorted(jj_sub["Cell Line"].dropna().unique())
                jj_reps       = sorted(jj_sub["_rep_letter"].unique())
                jj_color_map  = {r: REP_COLORS[i % len(REP_COLORS)]
                                 for i, r in enumerate(jj_reps)}

                fig_jj, ax_jj = plt.subplots(1, 1, figsize=(6, 5))
                fig_jj.suptitle(
                    "Freeze-Thaw Revival  |  wMel  |  JJ",
                    fontsize=13, fontweight="bold", y=1.02
                )

                jj_media_seen = set()
                for cell_line in jj_cell_lines:
                    csub = jj_sub[jj_sub["Cell Line"] == cell_line].sort_values("_date")
                    for rep in jj_reps:
                        rsub = csub[csub["_rep_letter"] == rep].dropna(
                            subset=["_date", titer_col]).sort_values("_date")
                        rsub = rsub[rsub[titer_col] > 0]
                        if rsub.empty:
                            continue
                        plot_series(ax_jj,
                                    rsub["_date"], rsub[titer_col],
                                    rsub["_media_class"],
                                    color=jj_color_map[rep],
                                    cell_line=cell_line,
                                    use_log=use_log)
                        jj_media_seen.update(rsub["_media_class"].unique())

                ax_jj.set_title("JJ", fontsize=11)
                ax_jj.set_xlabel("Collection date", fontsize=9)
                ax_jj.set_ylabel(ylabel, fontsize=10)
                ax_jj.tick_params(axis="x", rotation=30, labelsize=8)
                ax_jj.grid(True, which="both", linestyle=":", alpha=0.4)

                # Legend
                jj_handles = []
                for rep in jj_reps:
                    jj_handles.append(
                        mlines.Line2D([], [], color=jj_color_map[rep],
                                      linewidth=2,
                                      label=f"Rep {rep}" if rep else "Rep ?")
                    )
                jj_handles.append(mlines.Line2D([], [], color="none", label=""))
                for cl in jj_cell_lines:
                    jj_handles.append(
                        mlines.Line2D([], [], color="grey",
                                      linestyle=cellline_ls(cl), linewidth=1.8,
                                      label=CELLLINE_LABEL.get(cl, cl))
                    )
                jj_handles.append(mlines.Line2D([], [], color="none", label=""))
                for mc in sorted(jj_media_seen):
                    if mc in MEDIA_MARKER:
                        jj_handles.append(
                            mlines.Line2D([], [], color="grey", linestyle="none",
                                          marker=MEDIA_MARKER[mc], markersize=7,
                                          label=MEDIA_LABEL.get(mc, mc))
                        )
                fig_jj.legend(handles=jj_handles, loc="upper right", fontsize=8,
                              title="Legend", title_fontsize=9,
                              bbox_to_anchor=(1.0, 1.0))

                plt.tight_layout()
                pdf_path = os.path.join(out_dir, "freeze_thaw_wMel_JJ.pdf")
                fig_jj.savefig(pdf_path, bbox_inches="tight")
                plt.close(fig_jj)
                print(f"  Saved: {pdf_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot freeze-thaw revival titer over time")
    parser.add_argument("--input",   default="./output/freeze_thaw.csv")
    parser.add_argument("--out_dir", default="./output/plots_freeze_thaw")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f"ERROR: input not found: {args.input}")

    print(f"Loading {args.input}…")
    df = load_data(args.input)
    print(f"  {len(df)} rows, infections: {sorted(df['Infection'].dropna().unique())}")

    print("Generating freeze-thaw plots…")
    plot_freeze_thaw(df, args.out_dir)
    print("Done.")


if __name__ == "__main__":
    main()
