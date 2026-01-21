import os
from glob import glob
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

def axes_frac_to_data_y(ax, y_frac: float) -> float:
    """Convert an axes-fraction y (0..1) to a data y for this axis."""
    ypix = ax.transAxes.transform((0, y_frac))[1]
    return ax.transData.inverted().transform((0, ypix))[1]

def label_color(fp_val, rn_val, zero_tol=0.0):
    """
    Color by concordance of signs:
      - discordant (signs opposite): red
      - concordant (signs same): blue
      - neutral/zero: grey
    """
    if fp_val is None or rn_val is None or np.isnan(fp_val) or np.isnan(rn_val):
        return "0.5"
       # neutral grey
    return "#FF0000" if rn_val > zero_tol  else "#0000FF"  # red if discordant, blue if concordant

def place_list_labels_colored(ax, df, gene_col, idxs, colors,
                              x_edge, side="left",
                              y_start=0.04, y_end=0.22, x_pad=0.06,
                              fs=8):
    """
    Place labels in a non-overlapping vertical list between y_start..y_end (axes fraction),
    draw elbow leader lines from the column edge (x_edge) to the label.
    """
    if not idxs:
        return
    y_fracs = np.linspace(y_start, y_end, len(idxs))
    y_labels = [axes_frac_to_data_y(ax, yf) for yf in y_fracs]

    if side == "left":
        elbow_x = x_edge - 0.03
        label_x = x_edge - x_pad
        ha = "right"
    else:
        elbow_x = x_edge + 0.1
        label_x = x_edge + x_pad
        ha = "left"

    for i, y_lab, col in zip(idxs, y_labels, colors):
        row_y = i + 0.5
        # elbow leader line (horizontal then diagonal)
        ax.plot([x_edge, elbow_x], [row_y, row_y], color=col, lw=0.7, clip_on=False)
        # ax.plot([elbow_x, label_x], [row_y, y_lab], color=col, lw=0.7, clip_on=False)
        ax.text(label_x, row_y, str(df.loc[i, gene_col]), ha=ha, va="center",
                fontsize=fs, fontweight="bold", color=col, clip_on=False)

FILES = [
    "KLF15.Gn.1.csv",
    "KLF15.InN.5.csv",
    "KLF15.ExN.1.csv"
]
GLOB_PATTERN = None

GENE_COL       = "gene"
FOOTPRINT_COL  = "Footprint_log2fc"   # Footprint log2FC
RNA_COL        = "logFC"                      # RNA log2FC
OUTDIR         = "heatmaps"
FIG_SUFFIX     = ".pdf"

FIG_WIDTH_IN   = 3.8
ROW_HEIGHT_IN  = 0.16
MIN_HEIGHT_IN  = 4.0
DPI            = 300
FIXED_HEIGHT_IN = 6.0   # <- pick what looks good in your report (e.g., 5.5–7.0 in)
FIG_WIDTH_IN    = 3.8   # keep your current width (columns thinner if you reduce this)

SEP            = 0.35   # gap between the two columns (in column-width units)
OUTER_LW       = 2    # bold outer border linewidth
N_LABEL        = 5      # first/last N gene labels on the right

# Bottom colorbar sizing
CB_HEIGHT_FRAC = 0.06   # relative height of the colorbar (bigger = fatter)
CB_PAD         = 0.25   # padding between axes and colorbar
CB_SHRINK      = 0.60   # 0-1, shorter bar length (centered)
CB_HEIGHT_SCALE= 1.4    # scale up the colorbar height to make it fatter
# ==============================

def load_table(path: str, sheet_name=None) -> pd.DataFrame:
    ext = Path(path).suffix.lower()
    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(path, sheet_name=sheet_name)  # first sheet if None
    try:
        return pd.read_csv(path, sep=None, engine="python")
    except Exception:
        try:    return pd.read_csv(path, sep="\t")
        except: return pd.read_csv(path)

def draw_bold_outline(ax, x0, x1, y0, y1, lw=1.6, color="black"):
    z = 10
    ax.plot([x0, x1], [y0, y0], lw=lw, color=color, zorder=z, clip_on=False)  # bottom
    ax.plot([x0, x1], [y1, y1], lw=lw, color=color, zorder=z, clip_on=False)  # top
    ax.plot([x0, x0], [y0, y1], lw=lw, color=color, zorder=z, clip_on=False)  # left
    ax.plot([x1, x1], [y0, y1], lw=lw, color=color, zorder=z, clip_on=False)  # right

def plot_two_col_heatmap(df: pd.DataFrame, outfile: str, rna_label="RNA"):
    # coerce numeric & sort by Footprint high→low, keep ALL rows
    df = df.copy()

    if GENE_COL in df.columns:
        df = df.groupby(GENE_COL, as_index=False).mean(numeric_only=True)    

    for c in (FOOTPRINT_COL, RNA_COL):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df.sort_values(by=FOOTPRINT_COL, ascending=False, inplace=True, na_position="last")
    df.reset_index(drop=True, inplace=True)

    # ==========================================
    # USER-SELECTED GENES (only keep these rows)
    # ==========================================
    MY_GENES = [
        "HPCAL1",
        "KRAS",
        "CAMKV",
        "RAB6B",
        "SLC25A6",  
        "HOMER1"
    # add as many as you want
    ]

    # Keep only genes that appear in this dataset
    # find row indices of selected genes
    selected_idxs = df.index[df[GENE_COL].isin(MY_GENES)].tolist()

    # colors for selected genes
    selected_colors = [
    label_color(df.at[i, FOOTPRINT_COL], df.at[i, RNA_COL])
    for i in selected_idxs
    ]

    if df.empty:
        print(f"No selected genes found in {outfile}, skipping...")
        return

    fp = df[FOOTPRINT_COL].to_numpy(float).reshape(-1, 1)
    rn = df[RNA_COL].to_numpy(float).reshape(-1, 1)
    n_rows = len(df)
    y_edges = np.arange(0, n_rows + 1)

    # shared diverging map (white at 0)
    cmap = LinearSegmentedColormap.from_list("bwr_white_center", ["#0000FF", "white", "#FF0000"])
    vmax = np.nanmax(np.abs(np.concatenate([fp, rn]))) if n_rows > 0 else 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    # height_in = max(MIN_HEIGHT_IN, ROW_HEIGHT_IN * n_rows)
    # fig, ax = plt.subplots(figsize=(FIG_WIDTH_IN, height_in), dpi=DPI)

    height_in = FIXED_HEIGHT_IN        # <- uniform height for all datasets
    fig, ax = plt.subplots(figsize=(FIG_WIDTH_IN, height_in), dpi=DPI)

    # two SEPARATE columns (no row outlines)
    ax.pcolormesh([0, 1], y_edges, fp, cmap=cmap, norm=norm, shading="flat", edgecolors="none")
    ax.pcolormesh([1 + SEP, 2 + SEP], y_edges, rn, cmap=cmap, norm=norm, shading="flat", edgecolors="none")
    ax.invert_yaxis()

    # labels on top
    ax.set_xticks([0.5, 1.5 + SEP])
    ax.set_xticklabels(["Footprint", rna_label], rotation=55, ha="left", fontweight="bold")
    ax.xaxis.tick_top()
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    # bold outer borders per column
    draw_bold_outline(ax, 0, 1, 0, n_rows, lw=OUTER_LW)
    draw_bold_outline(ax, 1 + SEP, 2 + SEP, 0, n_rows, lw=OUTER_LW)

    # single bottom colorbar (shorter & fatter)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size=f"{int(CB_HEIGHT_FRAC*100)}%", pad=CB_PAD)
    bb = cax.get_position()
    new_w = bb.width * CB_SHRINK
    new_h = bb.height * CB_HEIGHT_SCALE
    new_x0 = bb.x0 + (bb.width - new_w) / 2.0
    cax.set_position([new_x0, bb.y0, new_w, new_h])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label("log2FC", fontsize=10, fontweight="bold")

    # first/last N gene labels on right
    # ---- labels: use two non-overlapping vertical lists with leader lines ----
    idx_first = list(range(n_rows))
    idx_last  = []

    # Colors per row based on discordance between FOOTPRINT_COL and RNA_COL
    colors_first = [label_color(df.at[i, FOOTPRINT_COL], df.at[i, RNA_COL]) for i in idx_first]
    colors_last  = [label_color(df.at[i, FOOTPRINT_COL], df.at[i, RNA_COL]) for i in idx_last]


    # make space on the right for labels
    x_edge_right = 2 + SEP
    ax.set_xlim(0.0, x_edge_right + 0.9)

    # top group (place in upper strip of the axis: 4%..22% of height)
    place_list_labels_colored(
    ax, df, gene_col=GENE_COL,
    idxs=selected_idxs, colors=selected_colors,
    x_edge=x_edge_right, side="right",
    y_start=0.96, y_end=0.06, x_pad=0.12, fs=8)

# Bottom group on the left (near bottom of the figure)
    fig.tight_layout()
    Path(OUTDIR).mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"Saved: {outfile}")

def main():
    # build the list of inputs
    inputs = []
    if GLOB_PATTERN:
        inputs.extend(sorted(glob(GLOB_PATTERN)))
    inputs.extend(FILES)

    if not inputs:
        raise SystemExit("No input files found. Add paths to FILES or set GLOB_PATTERN.")

    for path in inputs:
        if "InN.5" in path:
            rna_label = "InN.5"
        elif "GC.1" in path or "Gn.1" in path:   # some files use Gn.1, some GC.1
            rna_label = "GC.1"
        elif "ExN.1" in path:
            rna_label = "ExN.1"
        else:
            rna_label = "RNA"   # fallback    
        # skip empty strings (if both options used)
        if not path: 
            continue
        df = load_table(path)
        # multi-sheet Excel? loop sheets here (uncomment to process all):
        # if isinstance(df, dict):
        #     for sheet, df_sheet in df.items():
        #         out = Path(OUTDIR) / f"{Path(path).stem}_{sheet}{FIG_SUFFIX}"
        #         plot_two_col_heatmap(df_sheet, str(out))
        #     continue

        # normal case (single DataFrame)
        out = Path(OUTDIR) / f"{Path(path).stem}{FIG_SUFFIX}"
        plot_two_col_heatmap(df, str(out),rna_label)


if __name__ == "__main__":
    main()