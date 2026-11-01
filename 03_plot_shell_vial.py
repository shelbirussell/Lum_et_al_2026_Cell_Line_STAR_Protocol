#!/usr/bin/env python3
"""
03_plot_shell_vial.py
=====================
Plots titer over time for shell-vial (SV) experiments.

Split (DDMMRep) format for SV data:
  DDMM_N[_label]   e.g.  2010_1   2010_4   1209_1_ctl   0209_3-cipro
  or legacy letter  e.g.  0407A   2606A

The replicate identifier is:
  - For numeric splits: the number N (and optional treatment label) after
    the first underscore, e.g. '1', '2', '1_ctl', '3_cipro'.
    Each unique N (ignoring DDMM date prefix) is one continuous series
    tracked across all timepoints sorted by collection date.
  - For letter splits: the trailing letter(s), as before.

Layout:
  - One figure per Cell Line × Infection.
  - One subplot per person (Initials).
  - Each replicate (numeric or letter) is a separate coloured line.
  - Shields & Sang M3  → solid line
  - Schneider's        → dashed line
  - Unknown / mixed    → dotted line
  - Y-axis: titer (log10) for infected strains; host depth (linear) for uninf/DOX.
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

MEDIA_LINESTYLE = {"M3": "-", "Schneiders": "--", "Unknown": ":"}
MEDIA_LABEL     = {
    "M3":         "Shields & Sang M3",
    "Schneiders": "Schneider's",
    "Unknown":    "Unknown / mixed media",
}

# ── Titer columns ──────────────────────────────────────────────────────────
STRAIN_TITER = {
    "wMel":            ("wmel_titer",  "wMel titer\n(Wolbachia / Drosophila genomes)", True),
    "wMel(1/10)":      ("wmel_titer",  "wMel (1/10) titer\n(Wolbachia / Drosophila genomes)", True),
    "wMel SV-1":       ("wmel_titer",  "wMel titer\n(Wolbachia / Drosophila genomes)", True),
    "wMel SV-1'":      ("wmel_titer",  "wMel titer\n(Wolbachia / Drosophila genomes)", True),
    "wMel SV-2":       ("wmel_titer",  "wMel titer\n(Wolbachia / Drosophila genomes)", True),
    "wMel SV-2'":      ("wmel_titer",  "wMel titer\n(Wolbachia / Drosophila genomes)", True),
    "wRi":           ("wri_titer",   "wRi titer\n(Wolbachia / Drosophila genomes)",  True),
    "wRi_Merrill23": ("wri_titer",   "wRi titer\n(Wolbachia / Drosophila genomes)",  True),
    "wWil":          ("wwil_titer",  "wWil titer\n(Wolbachia / Drosophila genomes)", True),
    "wMel+wRi":      ("wmel_titer",  "wMel titer\n(Wolbachia / Drosophila genomes)", True),
    "wMel:wRi":      ("wmel_titer",  "wMel titer\n(Wolbachia / Drosophila genomes)", True),
    "DOX":           ("_host_depth", "Host genome mean depth\n(cell-health proxy)",  False),
    "uninf":         ("_host_depth", "Host genome mean depth\n(cell-health proxy)",  False),
    "uninfected":    ("_host_depth", "Host genome mean depth\n(cell-health proxy)",  False),
}

def add_host_depth_col(df):
    df = df.copy()
    is_dsim = df["Cell Line"].fillna("").str.lower().str.startswith("dsim")
    df["_host_depth"] = np.where(is_dsim, df["dsim_mean_depth"], df["dmel_mean_depth"])
    return df

def pick_titer_info(infection, df_cols):
    # Exact match only — substring fallback caused cross-experiment mixing
    # (e.g. wMel(1/10) being included in the wMel plot)
    info = STRAIN_TITER.get(str(infection).strip())
    if info is None:
        return None
    col, ylabel, use_log = info
    if col == "_host_depth" or col in df_cols:
        return col, ylabel, use_log
    return None

# ── Replicate ID extraction ────────────────────────────────────────────────
def extract_rep_id(split_str):
    """
    Extract the replicate identifier from a Split (DDMMRep) value.

    Numeric format  DDMM_N[_label|-label]  →  'N'  or  'N_label'
      e.g.  '2010_1'         → '1'
            '1209_3_cipro'   → '3_cipro'
            '0209_4-oflox'   → '4_oflox'   (normalise - to _)
            '0310_10'        → '10'

    Letter format   DDMM[A-Z]   →  trailing letter
      e.g.  '0407A'   → 'A'
            '2606A'   → 'A'
    """
    if pd.isna(split_str) or str(split_str).strip().lower() in {"", "nan"}:
        return ""
    s = str(split_str).strip().strip("\"'")

    # Numeric: DDMM_N or DDMM_N_label or DDMM_N-label
    m = re.match(r"^\d{4}_(\d+)(?:[_-](.+))?$", s)
    if m:
        num  = m.group(1)
        label = m.group(2)
        if label:
            # normalise separator to underscore
            return f"{num}_{label}"
        return num

    # Letter: trailing capital after digits
    m = re.search(r"[A-Z][A-Z']*$", s, re.IGNORECASE)
    if m:
        return m.group(0).upper()

    return s   # fallback: use as-is

# ── Colours ────────────────────────────────────────────────────────────────
REP_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    "#bcbd22", "#17becf", "#aec7e8", "#ffbb78",
]

def plot_rep_line(ax, x, y, media_classes, color, use_log):
    """Draw one replicate as a continuous coloured line with per-segment media style."""
    if len(x) == 0:
        return set()
    media_seen = set()
    for i in range(len(x) - 1):
        mc = media_classes.iloc[i]
        ls = MEDIA_LINESTYLE.get(mc, ":")
        media_seen.add(mc)
        pts_x = [x.iloc[i], x.iloc[i+1]]
        pts_y = [y.iloc[i], y.iloc[i+1]]
        if use_log:
            ax.semilogy(pts_x, pts_y, color=color, linestyle=ls, linewidth=1.6)
        else:
            ax.plot(pts_x, pts_y, color=color, linestyle=ls, linewidth=1.6)
    media_seen.add(media_classes.iloc[-1])
    plot_fn = ax.semilogy if use_log else ax.plot
    plot_fn(x, y, color=color, linestyle="none", marker="o", markersize=5)
    return media_seen


def load_data(path):
    df = pd.read_csv(path, low_memory=False)
    df["_date"] = pd.to_datetime(
        df["Date Collected (MM/DD/YYYY)"], errors="coerce", dayfirst=False)
    df.loc[df["_date"].dt.year.fillna(9999) < 2000, "_date"] = pd.NaT
    df["_media_class"] = df["Media Type"].apply(classify_media)
    df["_rep_id"]      = df["Split (DDMMRep)"].apply(extract_rep_id)
    df = add_host_depth_col(df)
    return df


PRESV_TITER_COLS = {
    "wMel": "wmel_titer",
    "wRi":  "wri_titer",
    "wWil": "wwil_titer",
}

def load_presv(path):
    """Load pre-sv source samples and return a lookup:
    (Cell Line, Initials) → list of (date, titer, infection) tuples.

    Titer is taken from the appropriate Wolbachia titer column for the
    source infection (e.g. wWil titer for a wWil source line).
    For uninf/DOX sources all Wolbachia titers are NaN → assigned 0,
    meaning zero Wolbachia genomes detected (correct baseline).
    These points will be plotted on every SV subplot for that cell line
    as a pre-treatment zero-titer reference.
    """
    if path is None or not os.path.isfile(path):
        return {}
    df = pd.read_csv(path, low_memory=False)
    df["_date"] = pd.to_datetime(
        df["Date Collected (MM/DD/YYYY)"], errors="coerce", dayfirst=False)
    df.loc[df["_date"].dt.year.fillna(9999) < 2000, "_date"] = pd.NaT
    df["Initials"] = df["Initials"].fillna("SLR")

    # Build SV first-date lookup per (cell line, target infection) so we match
    # the right SV experiment — e.g. JW18 DOX 2025 is a different experiment
    # from JW18 wMel 2023.
    sv_path = os.path.join(os.path.dirname(path), "shell_vial.csv")
    sv_first = {}   # (cell_line, infection) → earliest SV date
    sv_first_cl = {}  # cell_line → earliest SV date across all infections
    if os.path.isfile(sv_path):
        sv = pd.read_csv(sv_path, low_memory=False)
        sv["_date"] = pd.to_datetime(
            sv["Date Collected (MM/DD/YYYY)"], errors="coerce", dayfirst=False)
        sv.loc[sv["_date"].dt.year.fillna(9999) < 2000, "_date"] = pd.NaT
        for (cl, inf), grp in sv.groupby(["Cell Line", "Infection"]):
            first = grp["_date"].dropna().min()
            if pd.notna(first):
                sv_first[(cl, inf)] = first
        for cl, grp in sv.groupby("Cell Line"):
            first = grp["_date"].dropna().min()
            if pd.notna(first):
                sv_first_cl[cl] = first

    # For each pre-sv row, find the SV experiment it belongs to:
    # The relevant SV start is the earliest SV date for the same cell line
    # that comes AFTER the pre-sv date (i.e. this pre-sv precedes that experiment).
    # Among valid candidates per cell line, keep only the latest (closest to SV).
    from collections import defaultdict
    candidates = defaultdict(list)  # cell_line → [(date, titer, media, initials)]

    for _, row in df.iterrows():
        if pd.isna(row["_date"]):
            continue
        cl  = row["Cell Line"]
        inf = str(row.get("Infection", "") or "").strip()

        # Find the earliest SV date for this cell line that is after this pre-sv date
        # (i.e. which SV experiment does this pre-sv sample belong to?)
        relevant_sv_starts = [
            d for (c, i), d in sv_first.items()
            if c == cl and d > row["_date"]
        ]
        if not relevant_sv_starts:
            continue  # no SV experiment follows this pre-sv date
        nearest_sv_start = min(relevant_sv_starts)

        titer_col = PRESV_TITER_COLS.get(inf)
        if titer_col and titer_col in row.index and pd.notna(row[titer_col]):
            titer = float(row[titer_col])
        else:
            titer = 0.0
        candidates[cl].append(
            (row["_date"], titer, row.get("Media Type", ""),
             row["Initials"], nearest_sv_start)
        )

    # Per cell line, keep the single pre-sv sample closest to (but before)
    # its matched SV start date.
    lookup = {}
    for cl, pts in candidates.items():
        # Sort by gap to nearest SV start ascending, pick smallest gap
        pts.sort(key=lambda x: (x[4] - x[0]).days)
        best_date, best_titer, best_media, best_initials, _ = pts[0]
        key = (cl, best_initials)
        lookup[key] = [(best_date, best_titer, best_media)]

    return lookup


def plot_sv(df, out_dir, presv_lookup=None, out_suffix=""):
    os.makedirs(out_dir, exist_ok=True)

    cell_lines = sorted(df["Cell Line"].dropna().unique())
    infections = sorted(df["Infection"].dropna().unique())

    for cell_line in cell_lines:
        for infection in infections:
            titer_info = pick_titer_info(infection, df.columns)
            if titer_info is None:
                continue
            titer_col, ylabel, use_log = titer_info

            sub = df[
                (df["Cell Line"].fillna("") == cell_line) &
                (df["Infection"].fillna("") == infection) &
                df[titer_col].notna() &
                (df[titer_col] > 0)
            ].copy()

            if sub.empty:
                continue

            persons = sorted(sub["Initials"].fillna("Unknown").unique())
            n_persons = len(persons)

            fig, axes = plt.subplots(
                1, n_persons,
                figsize=(max(5 * n_persons, 6), 5),
                sharey=True, squeeze=False
            )
            fig.suptitle(
                f"Shell Vial  |  {cell_line}  ×  {infection}",
                fontsize=13, fontweight="bold", y=1.02
            )

            any_data   = False
            media_seen = set()

            for ax, person in zip(axes[0], persons):
                psub = sub[sub["Initials"].fillna("Unknown") == person]

                def _rep_sort_key(r):
                    r = str(r)
                    parts = r.split("_")
                    return (int(parts[0]), "_".join(parts[1:])) if parts[0].isdigit() else (999, r)
                rep_ids = sorted(psub["_rep_id"].unique(), key=_rep_sort_key)
                color_map = {r: REP_COLORS[i % len(REP_COLORS)]
                             for i, r in enumerate(rep_ids)}

                # Get zero anchor for this cell_line/infection/person if available
                anchor_pts = presv_lookup.get((cell_line, person, infection), []) \
                             if presv_lookup else []
                # Fall back to (cell_line, person) key for real pre-sv CSV data
                if not anchor_pts and presv_lookup:
                    anchor_pts = presv_lookup.get((cell_line, person), [])
                # Use the single earliest anchor point
                anchor = min(anchor_pts, key=lambda p: p[0]) \
                         if anchor_pts else None

                for rep in rep_ids:
                    rsub = psub[psub["_rep_id"] == rep].dropna(
                        subset=["_date", titer_col]).sort_values("_date")
                    rsub = rsub[rsub[titer_col] > 0]
                    if rsub.empty:
                        continue

                    # Prepend anchor point to connect the line from zero,
                    # but only if this rep's first sample is within 4 weeks
                    # of the anchor (i.e. it belongs to the same experiment)
                    rep_first_date = rsub["_date"].min()
                    anchor_relevant = (anchor is not None and
                        (rep_first_date - anchor[0]).days <= 28)
                    if anchor_relevant:
                        anchor_date, anchor_titer, anchor_media = anchor
                        if anchor_titer == 0 and use_log:
                            # For log axis: draw a dotted line from axis floor
                            # to the first real point, then hand off to plot_rep_line
                            yfloor = 1e-4
                            ax.semilogy(
                                [anchor_date, rsub["_date"].iloc[0]],
                                [yfloor, rsub[titer_col].iloc[0]],
                                color=color_map[rep], linestyle=":",
                                linewidth=1.4, zorder=3)
                            # Marker at anchor date at axis floor
                            ax.semilogy(anchor_date, yfloor,
                                        color=color_map[rep], linestyle="none",
                                        marker="v", markersize=8, zorder=4)
                        elif anchor_titer > 0:
                            # Real titer anchor — prepend as first data row
                            anchor_row = pd.DataFrame({
                                "_date": [anchor_date],
                                titer_col: [anchor_titer],
                                "_media_class": [anchor_media or rsub["_media_class"].iloc[0]]
                            })
                            rsub = pd.concat([anchor_row, rsub], ignore_index=True)

                    seen = plot_rep_line(
                        ax,
                        rsub["_date"], rsub[titer_col],
                        rsub["_media_class"],
                        color=color_map[rep],
                        use_log=use_log
                    )
                    media_seen.update(seen)
                    any_data = True

                presv_plotted = anchor is not None and anchor_pts

                ax.set_title(person, fontsize=11)
                ax.set_xlabel("Collection date", fontsize=9)
                ax.tick_params(axis="x", rotation=30, labelsize=8)
                ax.grid(True, which="major", linestyle=":", alpha=0.4)
                ax.minorticks_off()

                # Per-subplot replicate legend
                rep_handles = [
                    mlines.Line2D([], [], color=color_map[r],
                                  marker="o", markersize=5, linewidth=1.6,
                                  label=f"Rep {r}")
                    for r in rep_ids
                ]
                if presv_plotted:
                    rep_handles.append(
                        mlines.Line2D([], [], color="grey",
                                      marker="v", markersize=7,
                                      linestyle=":", linewidth=1.4,
                                      label="Start (titer = 0)")
                    )
                if rep_handles:
                    ax.legend(handles=rep_handles, fontsize=7,
                              loc="upper left",
                              title="Replicate", title_fontsize=7)

            if not any_data:
                plt.close(fig)
                continue

            for ax in axes[0]:
                if use_log:
                    ax.set_ylim(1e-4, 1e2)

            axes[0][0].set_ylabel(ylabel, fontsize=10)

            plt.tight_layout()

            # Figure-level media legend
            media_handles = [
                mlines.Line2D([], [], color="grey",
                              linestyle=MEDIA_LINESTYLE[m], linewidth=1.5,
                              label=MEDIA_LABEL.get(m, m))
                for m in sorted(media_seen) if m in MEDIA_LINESTYLE
            ]
            if media_handles:
                fig.legend(handles=media_handles, loc="upper right",
                           fontsize=9, title="Media", title_fontsize=9,
                           framealpha=0.9)

            safe = re.sub(r"[^\w\-]", "_", f"{cell_line}_{infection}")
            out_path = os.path.join(out_dir, f"sv_{safe}{out_suffix}.pdf")
            plt.savefig(out_path, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot shell-vial titer over time")
    parser.add_argument("--input",   default="./output/shell_vial.csv")
    parser.add_argument("--pre_sv",  default="./output/pre_sv.csv",
                        help="Pre-SV source sample CSV (pre_sv.csv from parser). "
                             "Pass --pre_sv '' to disable.")
    parser.add_argument("--out_dir", default="./output/plots_sv")
    parser.add_argument("--filter_reps", default=None,
                        help="Comma-separated rep ID prefixes to keep, e.g. '1,5' "
                             "keeps rep IDs starting with 1 or 5 (1_ctl, 5_ctl, etc.)")
    parser.add_argument("--filter_cell_line", default=None,
                        help="Only plot this cell line (e.g. JW18)")
    parser.add_argument("--filter_infection", default=None,
                        help="Only plot rows whose Infection contains this string")
    parser.add_argument("--out_suffix", default="",
                        help="Suffix appended to output filenames before .pdf")
    parser.add_argument("--zero_anchor", action="store_true",
                        help="Instead of real pre-sv data, place a titer=0 point "
                             "one week before the first sample for each "
                             "cell-line × infection × person group")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f"ERROR: input not found: {args.input}")

    print(f"Loading {args.input}…")
    df = load_data(args.input)
    print(f"  {len(df)} rows")
    print("  Replicate ID extraction sample:")
    sample = df[["Split (DDMMRep)", "_rep_id"]].drop_duplicates().sort_values("Split (DDMMRep)")
    print(sample.to_string(index=False))
    print()

    # Optional filters
    if args.filter_cell_line:
        df = df[df["Cell Line"].fillna("") == args.filter_cell_line]
        print(f"  Filtered to cell line: {args.filter_cell_line}")
    if args.filter_infection:
        df = df[df["Infection"].fillna("").str.contains(args.filter_infection, case=False, regex=False)]
        print(f"  Filtered to infections containing: {args.filter_infection}")
    if args.filter_reps:
        prefixes = [r.strip() for r in args.filter_reps.split(",")]
        # Match rep IDs that equal or start with any of the given prefixes
        # e.g. prefix "1" matches "1", "1_ctl"; prefix "5" matches "5", "5_ctl"
        def rep_matches(rep_id):
            r = str(rep_id)
            return any(r == p or r.startswith(p + "_") for p in prefixes)
        df = df[df["_rep_id"].apply(rep_matches)]
        print(f"  Filtered to rep ID prefixes: {prefixes} ({len(df)} rows remaining)")

    if args.zero_anchor:
        # Build synthetic zero-anchor lookup: titer=0 one week before
        # the earliest timepoint (T0 if present, else T1/min T) for each
        # (Cell Line, Infection, Initials) combination.
        import re as _re

        def _extract_T(sid):
            m = _re.search(r'_T(\d+)_', str(sid))
            return int(m.group(1)) if m else None

        df["_T_num"] = df["SampleID"].apply(_extract_T)

        presv_lookup = {}
        for (cl, inf, init), grp in df.groupby(["Cell Line", "Infection", "Initials"]):
            grp2 = grp.dropna(subset=["_date"])
            if grp2.empty:
                continue
            # Find earliest date at lowest T number (or min date if no T#)
            has_T = grp2["_T_num"].notna()
            if has_T.any():
                min_T  = grp2.loc[has_T, "_T_num"].min()
                first  = grp2.loc[grp2["_T_num"] == min_T, "_date"].min()
            else:
                first = grp2["_date"].min()
            zero_date = first - pd.Timedelta(days=7)
            # Key by (cell_line, infection, initials) to avoid mixing experiments
            key = (cl, init)
            # Use infection-specific sub-key to store separately, then merge
            # into the (cl, init) lookup only if this anchor is earlier than
            # any existing anchor for that cell line + person
            presv_lookup.setdefault(key, {})
            presv_lookup[key][inf] = (zero_date, 0.0, "")

        # Convert to the list format expected by the plotting code,
        # keeping one anchor per (cell_line, infection, person) triplet
        final_lookup = {}
        for (cl, init), inf_dict in presv_lookup.items():
            # Store all infection-specific anchors; the plotting code will
            # match them to the correct subplot via the cell_line + person key
            for inf, pt in inf_dict.items():
                final_key = (cl, init, inf)
                final_lookup[final_key] = [pt]
        presv_lookup = final_lookup

        print(f"  Zero anchors built for {len(presv_lookup)} cell-line/infection/person combinations")
        for (cl, init, inf), pts in sorted(presv_lookup.items()):
            for pt_date, _, _ in pts:
                print(f"    {cl} / {inf} / {init}: anchor date = {pt_date.date()}")
    else:
        presv_lookup = load_presv(args.pre_sv if args.pre_sv else None)
        if presv_lookup:
            print(f"  Pre-SV reference points loaded for "
                  f"{len(presv_lookup)} cell-line/person combinations")

    print("Generating shell-vial plots…")
    plot_sv(df, args.out_dir, presv_lookup=presv_lookup,
            out_suffix=args.out_suffix)
    print("Done.")


if __name__ == "__main__":
    main()
