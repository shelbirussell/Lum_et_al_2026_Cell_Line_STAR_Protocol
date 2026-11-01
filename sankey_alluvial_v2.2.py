"""
Sankey and Alluvial diagrams for Drosophila cell line immortalization data.

Sankey  : Species → Fly Line → Wolbachia infection → Outcome
          (flow width ∝ flask count; outcome is the final fate of each flask)

Alluvial: Species → Fly Line → Wolbachia infection
          (shows how experimental effort was distributed across categories)

Both are drawn with matplotlib using cubic-Bezier ribbons.
Three species subsets are produced for each diagram:
  all       — all four species  (05_sankey_all, 06_alluvial_all)
  noF1      — mel, sim, wil     (05_sankey_noF1, 06_alluvial_noF1)
  melSim    — mel and sim only  (05_sankey_melSim, 06_alluvial_melSim)

New outcome categories vs prior versions:
  "ongoing"   → lowercase replacement for "OG"; outcome not yet resolved
  "harvested" → primary culture deliberately collected before immortalization
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import matplotlib.patheffects as pe
import warnings
import sys
warnings.filterwarnings("ignore")

# ── Colour palettes ───────────────────────────────────────────────────────────
# Designed for colour-blind accessibility (Okabe-Ito / IBM / Paul Tol frameworks).
# Species families share hue with their associated fly lines and Wolbachia strains:
#   D. melanogaster / wMel  → blue family
#   D. simulans / wRi       → red family
#   D. willistoni / wWil    → green
#   F1 hybrid               → purple (perceptual midpoint of blue and red)
# Outcome colours are drawn from a neutral set that avoids red/green conflict.

SP_COL = {
    "D. mel.":   "#0077BB",   # IBM blue
    "D. sim.":   "#BB2200",   # deep brick red
    "D. wil.":   "#228833",   # Okabe teal-green
    "F1 hybrid": "#7B3FA0",   # purple
}
LINE_COL = {
    # mel lines — blue shades, dark → light
    "ub:Khc":        "#0077BB",  # anchor blue
    "mei-P26[1]":    "#4DAADD",  # mid sky blue
    "PlasMem:GFP":   "#88CCEE",  # pale blue
    "Jupiter:GFP":   "#C8E8F5",  # very pale blue
    # sim lines — red shades, dark → light
    "Riv84":         "#BB2200",  # anchor red
    "Mer23":         "#EE6644",  # coral-red
    "w[501]":        "#FFAA99",  # salmon
    # wil line — green
    "genome line":   "#228833",
    # F1 hybrid line — purple
    "F1 hybrid":     "#7B3FA0",
}
WOLB_COL = {
    "wMel_(DB)":     "#0077BB",  # blue = mel
    "wRi_(Riv84)":   "#BB2200",  # red = Riv84
    "wRi_(Mer23)":   "#EE6644",  # coral = Mer23
    "wWil":          "#228833",  # green = willistoni
}
WOLB_COL = {
    "wMel_(DB)":     "#0077BB",
    "wRi":           "#BB2200",
    "wRi_(Riv84)":   "#BB2200",
    "wRi_(Mer23)":   "#EE6644",
    "wWil":          "#228833",
    "Unknown/mixed": "#999999",
}
OUT_COL = {
    # Neutral set — avoids red/green conflict for protanopes/deuteranopes;
    # Contaminated uses dark grey so it doesn't clash with sim species colour.
    "Immortalized": "#009988",   # Okabe teal
    "Harvested":    "#DDAA33",   # amber-gold
    "Primary Only": "#EE7733",   # Okabe amber
    "Ongoing":      "#AA4499",   # Okabe purple
    "No Primary":   "#BBBBBB",   # mid grey
    "Contaminated": "#555555",   # dark grey
}

plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.spines.left":  False,
    "axes.spines.bottom":False,
    "figure.dpi":        180,
})


# ── Data preparation ──────────────────────────────────────────────────────────
def load_data(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)   # CSV (comma-separated)
    df["Immortalized?"] = df["Immortalized?"].astype(str).str.strip().str.lower()

    def outcome(row):
        imm  = row["Immortalized?"]
        cont = str(row["Contamination?"]).strip()
        prim = str(row["Primary_Culture?"]).strip()
        if cont in ("Y", "M"):   return "Contaminated"
        if imm == "harvested":   return "Harvested"
        if imm == "ongoing":     return "Ongoing"
        if prim == "N":          return "No Primary"
        if imm == "y":           return "Immortalized"
        return "Primary Only"

    sp_short = {
        "melanogaster":    "D. mel.",
        "simulans":        "D. sim.",
        "willistoni":      "D. wil.",
        "F1hybrid:melxsim":"F1 hybrid",
    }
    line_short = {
        "Dsim_w[501]":                                                   "w[501]",
        "D._simulans_Riv84":                                             "Riv84",
        "D._simulans_Mer23_#5":                                          "Mer23",
        "D._melanogaster_ub:Khc":                                        "ub:Khc",
        "D._melogaster_PlasMem:GFP":                                     "PlasMem:GFP",
        "D._melanogaster_mei-P26[1]":                                    "mei-P26[1]",
        "D._melanogaster_Jupiter:GFP":                                   "Jupiter:GFP",
        "D._willistoni_(genome_line)":                                   "genome line",
        "D._melanogaster_virginator_x_D._simulans_w[501]_F1_hybrid":    "F1 hybrid",
    }

    df["outcome"] = df.apply(outcome, axis=1)
    df["sp"]      = df["Drosophila_species"].map(sp_short)
    df["line"]    = df["Fly_Line"].map(line_short)
    df["wolb"]    = df["Wolbachia infection"].fillna("Unknown/mixed")
    return df


# ── Bezier ribbon helper ──────────────────────────────────────────────────────
def bezier_ribbon(ax, x0, y0_bot, y0_top, x1, y1_bot, y1_top,
                  color, alpha=0.55, zorder=1):
    """Draw a filled cubic-Bezier ribbon between two vertical stacks."""
    cx = (x0 + x1) / 2
    verts = [
        (x0,  y0_bot),
        (cx,  y0_bot),
        (cx,  y1_bot),
        (x1,  y1_bot),
        (x1,  y1_top),
        (cx,  y1_top),
        (cx,  y0_top),
        (x0,  y0_top),
        (x0,  y0_bot),
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    patch = mpatches.PathPatch(
        Path(verts, codes),
        facecolor=color, edgecolor="none", alpha=alpha, zorder=zorder
    )
    ax.add_patch(patch)


def draw_node(ax, x, y_bot, y_top, color, width=0.018, label=None,
              label_side="right", fontsize=8):
    """Draw a coloured node bar and optional label."""
    rect = mpatches.FancyBboxPatch(
        (x - width/2, y_bot), width, y_top - y_bot,
        boxstyle="square,pad=0",
        facecolor=color, edgecolor="white", linewidth=0.6, zorder=3
    )
    ax.add_patch(rect)
    if label:
        mid = (y_bot + y_top) / 2
        xoff = x + width/2 + 0.012 if label_side == "right" else x - width/2 - 0.012
        ha   = "left" if label_side == "right" else "right"
        n_txt = f"  {label}" if label_side == "right" else f"{label}  "
        ax.text(xoff, mid, n_txt, va="center", ha=ha, fontsize=fontsize,
                color="#222222", zorder=4)


# ── Node-stack builder ────────────────────────────────────────────────────────
def build_stacks(groups: list[str], counts: pd.Series,
                 palette: dict, gap: float = 0.01) -> dict:
    """
    Given an ordered list of group names and their counts,
    return {name: (y_bot, y_top, color)} with small gaps between nodes.
    Total height = sum(counts); groups stacked bottom-to-top.
    """
    total = counts.sum()
    stacks = {}
    cursor = 0.0
    for i, g in enumerate(groups):
        h = counts.get(g, 0) / total
        stacks[g] = (cursor, cursor + h, palette.get(g, "#AAAAAA"))
        cursor += h + gap
    # Normalize so total height = 1 (gaps included)
    total_h = cursor - gap
    stacks = {g: (bot/total_h, top/total_h, col)
              for g, (bot, top, col) in stacks.items()}
    return stacks


# ═══════════════════════════════════════════════════════════════════════════════
# SANKEY: Species → Line → Wolbachia → Outcome
# ═══════════════════════════════════════════════════════════════════════════════
def plot_sankey(df: pd.DataFrame, save_path: str):
    # Drop rows whose sp/line mapped to NaN (excluded species)
    df = df.dropna(subset=["sp", "line"])

    # ── aggregation ──────────────────────────────────────────────────────────
    grp = (df.groupby(["sp","line","wolb","outcome"])
             .size().reset_index(name="n"))
    total = grp["n"].sum()

    # Ordered node lists — filter to what is present in this subset
    ALL_SP   = ["D. mel.", "D. sim.", "D. wil.", "F1 hybrid"]
    ALL_LINE = ["ub:Khc","mei-P26[1]","PlasMem:GFP","Jupiter:GFP",
                "Riv84","Mer23","w[501]","genome line","F1 hybrid"]
    ALL_WOLB = ["wMel_(DB)","wRi","wRi_(Riv84)","wRi_(Mer23)","wWil","Unknown/mixed"]
    ALL_OUT  = ["Immortalized","Harvested","Primary Only","Ongoing","No Primary","Contaminated"]

    present = lambda lst, col: [x for x in lst if x in grp[col].values]
    SP_ORD   = present(ALL_SP,   "sp")
    LINE_ORD = present(ALL_LINE, "line")
    WOLB_ORD = present(ALL_WOLB, "wolb")
    OUT_ORD  = present(ALL_OUT,  "outcome")

    # Column x positions
    xs = [0.05, 0.33, 0.61, 0.90]

    # Build per-column stacks
    sp_counts   = grp.groupby("sp")["n"].sum()
    line_counts = grp.groupby("line")["n"].sum()
    wolb_counts = grp.groupby("wolb")["n"].sum()
    out_counts  = grp.groupby("outcome")["n"].sum()

    sp_stk   = build_stacks(SP_ORD,   sp_counts,   SP_COL,   gap=0.012)
    line_stk = build_stacks(LINE_ORD, line_counts, LINE_COL, gap=0.012)
    wolb_stk = build_stacks(WOLB_ORD, wolb_counts, WOLB_COL, gap=0.012)
    out_stk  = build_stacks(OUT_ORD,  out_counts,  OUT_COL,  gap=0.012)

    # ── figure ───────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(16, 9))
    ax.set_xlim(-0.02, 1.18)
    ax.set_ylim(-0.04, 1.10)
    ax.set_aspect("auto")
    ax.axis("off")
    fig.patch.set_facecolor("#FAFAFA")

    fig.text(0.5, 1.01, "Sankey Diagram: Flask Outcomes by Species, Line & Wolbachia Infection",
             ha="center", va="top", fontsize=14, fontweight="bold", transform=ax.transAxes)
    fig.text(0.5, 0.97, f"Width ∝ number of flasks  (total n = {total})",
             ha="center", va="top", fontsize=9.5, color="#555", transform=ax.transAxes)

    # Column headers
    for x, lbl in zip(xs, ["Species", "Fly Line", "Wolbachia", "Outcome"]):
        ax.text(x, 1.055, lbl, ha="center", va="bottom", fontsize=11,
                fontweight="bold", color="#333", transform=ax.transAxes)

    # ── draw nodes ────────────────────────────────────────────────────────────
    NODE_W = 0.016
    for g, (bot, top, col) in sp_stk.items():
        draw_node(ax, xs[0], bot, top, col, NODE_W,
                  label=f"{g}\n(n={int(sp_counts.get(g,0))})",
                  label_side="left", fontsize=8.5)
    for g, (bot, top, col) in line_stk.items():
        draw_node(ax, xs[1], bot, top, col, NODE_W, label=g,
                  label_side="left", fontsize=8)
    for g, (bot, top, col) in wolb_stk.items():
        draw_node(ax, xs[2], bot, top, col, NODE_W, label=g,
                  label_side="left", fontsize=8)
    for g, (bot, top, col) in out_stk.items():
        draw_node(ax, xs[3], bot, top, col, NODE_W,
                  label=f"{g}\n(n={int(out_counts.get(g,0))})",
                  label_side="right", fontsize=8.5)

    # ── draw ribbons ──────────────────────────────────────────────────────────
    # We need per-node "cursors" tracking how much of each node has been allocated
    sp_cur   = {g: bot for g, (bot, top, col) in sp_stk.items()}
    line_cur = {g: bot for g, (bot, top, col) in line_stk.items()}
    wolb_cur = {g: bot for g, (bot, top, col) in wolb_stk.items()}
    out_cur  = {g: bot for g, (bot, top, col) in out_stk.items()}

    def scale(n): return n / total  # raw count → axis fraction (before gap rescaling)

    # We work in normalized (gap-corrected) units
    # Need to scale individual flows by the same factor used in build_stacks
    def gap_scale(stk, counts, gap=0.012):
        """Return the scaling factor: axis_unit / fraction_of_total"""
        raw_total = sum((counts.get(g,0)/counts.sum()) for g in stk) + gap*(len(stk)-1)
        return 1 / raw_total

    gs_sp   = gap_scale(sp_stk,   sp_counts)
    gs_line = gap_scale(line_stk, line_counts)
    gs_wolb = gap_scale(wolb_stk, wolb_counts)
    gs_out  = gap_scale(out_stk,  out_counts)

    # Flow 1: Species → Line
    for _, row in (grp.groupby(["sp","line"])["n"].sum()
                      .reset_index(name="n").iterrows()):
        sp, line, n = row["sp"], row["line"], row["n"]
        h_sp   = (n / sp_counts.sum())   * gs_sp
        h_line = (n / line_counts.sum()) * gs_line
        y0b = sp_cur[sp];     sp_cur[sp]   += h_sp
        y1b = line_cur[line]; line_cur[line]+= h_line
        col = SP_COL.get(sp, "#888")
        bezier_ribbon(ax, xs[0], y0b, y0b+h_sp, xs[1], y1b, y1b+h_line,
                      color=col, alpha=0.40)

    # Reset cursors for Flow 2
    line_cur = {g: bot for g, (bot, top, col) in line_stk.items()}

    # Flow 2: Line → Wolbachia
    for _, row in (grp.groupby(["line","wolb"])["n"].sum()
                      .reset_index(name="n").iterrows()):
        line, wolb, n = row["line"], row["wolb"], row["n"]
        h_line = (n / line_counts.sum()) * gs_line
        h_wolb = (n / wolb_counts.sum()) * gs_wolb
        y0b = line_cur[line]; line_cur[line]+= h_line
        y1b = wolb_cur[wolb]; wolb_cur[wolb]+= h_wolb
        col = LINE_COL.get(line, "#888")
        bezier_ribbon(ax, xs[1], y0b, y0b+h_line, xs[2], y1b, y1b+h_wolb,
                      color=col, alpha=0.40)

    # Reset
    wolb_cur = {g: bot for g, (bot, top, col) in wolb_stk.items()}

    # Flow 3: Wolbachia → Outcome
    for _, row in (grp.groupby(["wolb","outcome"])["n"].sum()
                      .reset_index(name="n").iterrows()):
        wolb, out, n = row["wolb"], row["outcome"], row["n"]
        h_wolb = (n / wolb_counts.sum()) * gs_wolb
        h_out  = (n / out_counts.sum())  * gs_out
        y0b = wolb_cur[wolb]; wolb_cur[wolb]+= h_wolb
        y1b = out_cur[out];   out_cur[out]  += h_out
        col = OUT_COL.get(out, "#888")
        bezier_ribbon(ax, xs[2], y0b, y0b+h_wolb, xs[3], y1b, y1b+h_out,
                      color=col, alpha=0.55)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(save_path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print(f"Saved: {save_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# ALLUVIAL: Species → Fly Line → Wolbachia infection  (effort distribution)
# ═══════════════════════════════════════════════════════════════════════════════
def plot_alluvial(df: pd.DataFrame, save_path: str):
    df = df.dropna(subset=["sp", "line"])
    total = len(df)

    ALL_SP   = ["D. mel.", "D. sim.", "D. wil.", "F1 hybrid"]
    ALL_LINE = ["ub:Khc","mei-P26[1]","PlasMem:GFP","Jupiter:GFP",
                "Riv84","Mer23","w[501]","genome line","F1 hybrid"]
    ALL_WOLB = ["wMel_(DB)","wRi","wRi_(Riv84)","wRi_(Mer23)","wWil","Unknown/mixed"]

    SP_ORD   = [x for x in ALL_SP   if x in df["sp"].values]
    LINE_ORD = [x for x in ALL_LINE if x in df["line"].values]
    WOLB_ORD = [x for x in ALL_WOLB if x in df["wolb"].values]

    xs = [0.08, 0.45, 0.82]

    sp_counts   = df.groupby("sp").size()
    line_counts = df.groupby("line").size()
    wolb_counts = df.groupby("wolb").size()

    sp_stk   = build_stacks(SP_ORD,   sp_counts,   SP_COL,   gap=0.016)
    line_stk = build_stacks(LINE_ORD, line_counts, LINE_COL, gap=0.016)
    wolb_stk = build_stacks(WOLB_ORD, wolb_counts, WOLB_COL, gap=0.016)

    fig, ax = plt.subplots(figsize=(14, 9))
    ax.set_xlim(-0.02, 1.15)
    ax.set_ylim(-0.04, 1.10)
    ax.axis("off")
    fig.patch.set_facecolor("#FAFAFA")

    fig.text(0.5, 1.01,
             "Alluvial Diagram: Experimental Effort by Species, Fly Line & Wolbachia Infection",
             ha="center", va="top", fontsize=13, fontweight="bold",
             transform=ax.transAxes)
    fig.text(0.5, 0.97,
             f"Width ∝ number of flasks  (total n = {total})",
             ha="center", va="top", fontsize=9.5, color="#555",
             transform=ax.transAxes)

    for x, lbl in zip(xs, ["Species", "Fly Line", "Wolbachia Infection"]):
        ax.text(x, 1.055, lbl, ha="center", va="bottom", fontsize=11,
                fontweight="bold", color="#333", transform=ax.transAxes)

    NODE_W = 0.022
    # Species nodes — label left
    for g, (bot, top, col) in sp_stk.items():
        draw_node(ax, xs[0], bot, top, col, NODE_W,
                  label=f"{g}  (n={int(sp_counts.get(g,0))})",
                  label_side="left", fontsize=9)
    # Line nodes — label left
    for g, (bot, top, col) in line_stk.items():
        draw_node(ax, xs[1], bot, top, col, NODE_W,
                  label=g, label_side="left", fontsize=8.5)
    # Wolbachia nodes — label right
    for g, (bot, top, col) in wolb_stk.items():
        draw_node(ax, xs[2], bot, top, col, NODE_W,
                  label=f"{g}  (n={int(wolb_counts.get(g,0))})",
                  label_side="left", fontsize=8.5)

    # gap scaling factors
    def gap_scale(stk, counts, gap=0.016):
        raw_total = sum((counts.get(g,0)/counts.sum()) for g in stk) + gap*(len(stk)-1)
        return 1 / raw_total

    gs_sp   = gap_scale(sp_stk,   sp_counts)
    gs_line = gap_scale(line_stk, line_counts)
    gs_wolb = gap_scale(wolb_stk, wolb_counts)

    # Flow 1: Species → Line
    sp_cur   = {g: bot for g, (bot, top, _) in sp_stk.items()}
    line_cur = {g: bot for g, (bot, top, _) in line_stk.items()}

    for (sp, line), n in df.groupby(["sp","line"]).size().items():
        h_sp   = (n / sp_counts.sum())   * gs_sp
        h_line = (n / line_counts.sum()) * gs_line
        y0b = sp_cur[sp];     sp_cur[sp]   += h_sp
        y1b = line_cur[line]; line_cur[line]+= h_line
        col = SP_COL.get(sp, "#888")
        bezier_ribbon(ax, xs[0], y0b, y0b+h_sp, xs[1], y1b, y1b+h_line,
                      color=col, alpha=0.38)

    # Flow 2: Line → Wolbachia
    line_cur = {g: bot for g, (bot, top, _) in line_stk.items()}
    wolb_cur = {g: bot for g, (bot, top, _) in wolb_stk.items()}

    for (line, wolb), n in df.groupby(["line","wolb"]).size().items():
        h_line = (n / line_counts.sum()) * gs_line
        h_wolb = (n / wolb_counts.sum()) * gs_wolb
        y0b = line_cur[line]; line_cur[line]+= h_line
        y1b = wolb_cur[wolb]; wolb_cur[wolb]+= h_wolb
        col = LINE_COL.get(line, "#888")
        bezier_ribbon(ax, xs[1], y0b, y0b+h_line, xs[2], y1b, y1b+h_wolb,
                      color=col, alpha=0.38)

    # Percentage annotations on Wolbachia nodes
    for g, (bot, top, col) in wolb_stk.items():
        pct = wolb_counts.get(g, 0) / total * 100
        mid = (bot + top) / 2
        ax.text(xs[2] + NODE_W/2 + 0.25, mid,
                f"{pct:.0f}%", va="center", ha="left",
                fontsize=7.5, color="#555")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(save_path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print(f"Saved: {save_path}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():

    if len(sys.argv) < 3:
        print("Error: Please provide a file path and output directory.")
        print("Usage: python sankey_alluvial.py <input_csv> <output_dir>")
        sys.exit(1)

    path    = sys.argv[1]
    out_dir = sys.argv[2].rstrip("/") + "/"

    try:
        df = load_data(path)
    except FileNotFoundError:
        print(f"Error: file '{path}' not found.")
        sys.exit(1)

    # Species subsets — NaN sp/line rows are dropped inside each plot function
    SUBSETS = {
        "all":    None,                              # all four species
        "noF1":   ["D. mel.", "D. sim.", "D. wil."], # no F1 hybrid
        "melSim": ["D. mel.", "D. sim."],            # mel and sim only
    }

    for tag, sp_keep in SUBSETS.items():
        if sp_keep is not None:
            df_sub = df[df["sp"].isin(sp_keep)].copy()
        else:
            df_sub = df.copy()
        plot_sankey( df_sub, f"{out_dir}05_sankey_{tag}.pdf")
        plot_alluvial(df_sub, f"{out_dir}06_alluvial_{tag}.pdf")

    print("Done.")

if __name__ == "__main__":
    main()
