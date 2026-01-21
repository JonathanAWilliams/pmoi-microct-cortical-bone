"""
Generate publication-ready slice-wise profiles and summary bar charts for:
Ct.Ar, TMD, pMOA, and pMOI.

What this script does
---------------------
1) Reads an Excel workbook where the first two sheets represent two experimental groups
   (e.g. 'SCI' and 'SHAM'), each containing slice-wise data for multiple specimens.

2) For each metric:
   - Extracts the slice-wise matrix [slices x specimens] for Group A and Group B
   - Interpolates each specimen to a common bone-length percentage axis
     (default: 58–62% of bone length) so specimens with different slice counts
     can be compared on the same x-axis
   - Plots group mean ± SE across the % axis
   - Computes Welch's t-test at each x-position and shows the p-values as a heatmap strip

3) Produces summary bar charts using per-specimen mean values across slices, with
   Welch's t-test between groups.

Expected Excel layout
---------------------
- The first two sheets in the workbook are interpreted as Group A and Group B.
- Each sheet is assumed to contain repeated metric columns for multiple specimens.
  The script detects metric columns by searching the top rows for strings like:
  'Ct.Ar', 'TMD', 'pMOA', 'pMOI' (excluding 'p-<metric>' columns).

Outputs
-------
Written to the same folder as this script:
- multi_panel_mean_se_<GroupA>_vs_<GroupB>.tiff
- bars_mean_<GroupA>_vs_<GroupB>.tiff
- bars_mean_<GroupA>_vs_<GroupB>.csv
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from scipy.interpolate import interp1d
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.gridspec import GridSpec

# -------- User settings --------
# Excel workbook to analyse.
# If a relative path is provided, it is resolved relative to this script's folder.
EXCEL_NAME   = "SCI_study_10W_slices_for_paper.xlsx"   # Excel in same folder as this script

# Order of metrics to plot (left to right): Ct.Ar, TMD, pMOA, pMOI
METRICS      = ["Ct.Ar", "TMD", "pMOA", "pMOI"]

SEARCH_ROWS  = 10
# Figure export settings
FIG_DPI      = 600   # high DPI for print

# Bone-length normalisation range (percent of total bone length)
START_PERCENT = 58.0
END_PERCENT   = 62.0
# --------------------------------

# Heatmap bins + colors for p-values
BOUNDS = [0.0, 0.001, 0.01, 0.05, 1.00001]
CMAP   = ListedColormap(["red", "yellow", "green", "blue"])
NORM   = BoundaryNorm(BOUNDS, CMAP.N)

# Group colors (consistent across all metrics)
GROUP_A_COLOR = "#1f4e79"  # dark blue
GROUP_B_COLOR = "#1b5e20"  # dark green

# Warm accent palette for borders/titles (per metric)
BORDER_TITLE_COLOR = {
    "Ct.Ar": "#C2410C",  # burnt orange
    "pMOA":  "#B45309",  # golden brown
    "TMD":   "#7C2D12",  # burgundy red
    "pMOI":  "#9A3412",  # deep orange-brown
}

# Publication font/line settings
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 12,
    "axes.titlesize": 16,
    "axes.titleweight": "bold",
    "axes.labelsize": 14,
    "axes.labelweight": "bold",
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "legend.frameon": False,
    "lines.linewidth": 2,
})

def safe_name(s: str) -> str:
    """Make a string safe for filenames."""
    return re.sub(r"[^A-Za-z0-9_\-]+", "_", s.strip())

def find_metric_columns_from_cells(df, label, search_rows=SEARCH_ROWS):
    """
    Scan the top rows (no header) for cells containing <label>, excluding 'p-<label>'.
    This picks up columns under the repeated label rows.
    """
    max_row = min(search_rows, len(df))
    lab = str(label).lower()
    cols = []
    for col in df.columns:
        vals = df[col].iloc[:max_row].astype(str).str.lower()
        if vals.str.contains(lab, na=False).any() and not vals.str.contains(f"p-{lab}", na=False).any():
            cols.append(col)
    return cols

def find_metric_columns_from_headers(df, label):
    """
    Fallback: if sheet is read with headers, find columns whose header contains <label>
    but not 'p-<label>'.
    """
    lab = str(label).lower()
    return [c for c in df.columns
            if (lab in str(c).lower()) and (f"p-{lab}" not in str(c).lower())]

def extract_metric_matrix(df_nohdr, label, excel_path, sheet_name_or_idx):
    """
    Identify metric columns on a sheet and return the numeric matrix only.
    - df_nohdr: sheet read with header=None.
    - label: 'Ct.Ar', 'pMOA', 'TMD', or 'pMOI'.

    Returns:
    - numeric: 2D array [slices x samples]
      (rows with all NaN removed)
    """
    metric_cols = find_metric_columns_from_cells(df_nohdr, label, SEARCH_ROWS)
    if not metric_cols:
        # Fallback: read with header and search column names
        df_hdr = pd.read_excel(excel_path, sheet_name=sheet_name_or_idx, header=0)
        header_cols = find_metric_columns_from_headers(df_hdr, label)
        if not header_cols:
            raise ValueError(f"No '{label}' columns detected on sheet '{sheet_name_or_idx}'.")
        numeric = df_hdr[header_cols].apply(pd.to_numeric, errors="coerce")
    else:
        numeric = df_nohdr[metric_cols].apply(pd.to_numeric, errors="coerce")

    # Remove rows that are all NaN (e.g. trailing slices)
    numeric = numeric.dropna(how="all")

    return numeric.to_numpy(float)

def interpolate_to_common_percent_grid(mat, start_percent, end_percent, num_points=None):
    """
    For a given [slices x samples] matrix 'mat', interpolate each sample (column)
    to a common % bone length grid between start_percent and end_percent.

    For each sample j:
    - Take its finite values y_j (length n_j).
    - Assume they are sampled uniformly from start_percent to end_percent:
        x_j = linspace(start_percent, end_percent, n_j)
    - Interpolate onto x_grid using:
        * cubic if n_j >= 4
        * linear if n_j == 2 or 3
        * constant if n_j == 1

    Returns:
    - x_grid: 1D array of length num_points
    - grid_vals: 2D array [num_points x samples]
    """
    num_slices, num_samples = mat.shape

    if num_points is None:
        # Use the maximum number of finite slices observed in any sample
        counts = np.sum(np.isfinite(mat), axis=0)
        max_n = int(np.nanmax(counts)) if np.nanmax(counts) > 1 else num_slices
        num_points = max_n

    x_grid = np.linspace(start_percent, end_percent, num_points)
    grid_vals = np.full((num_points, num_samples), np.nan, dtype=float)

    for j in range(num_samples):
        y = mat[:, j]
        mask = np.isfinite(y)
        y_valid = y[mask]
        n = y_valid.size
        if n == 0:
            continue

        x_j = np.linspace(start_percent, end_percent, n)

        if n >= 4:
            kind = "cubic"
        elif n >= 2:
            kind = "linear"
        else:
            # Only one point: treat as constant function
            grid_vals[:, j] = y_valid[0]
            continue

        f = interp1d(x_j, y_valid, kind=kind, bounds_error=False,
                     fill_value="extrapolate")
        grid_vals[:, j] = f(x_grid)

    return x_grid, grid_vals

def per_grid_pvalues(grid_a, grid_b):
    """
    Welch t-test at each position along the common % grid.

    - grid_a, grid_b: [num_points x samples]
    Returns:
    - pvals: length num_points, NaN where t-test not possible.
    """
    num_points = grid_a.shape[0]
    pvals = np.full(num_points, np.nan, dtype=float)
    for k in range(num_points):
        a = grid_a[k, :]
        b = grid_b[k, :]
        a = a[np.isfinite(a)]
        b = b[np.isfinite(b)]
        if len(a) >= 2 and len(b) >= 2:
            _, p = ttest_ind(a, b, equal_var=False)
            pvals[k] = p
    return pvals

def metric_ylabel(metric: str) -> str:
    """Return nicely formatted y-axis label with correct units."""
    if metric == "Ct.Ar":
        return "Ct.Ar (mm\u00b2)"           # mm²
    elif metric == "pMOA":
        return "pMOA (mm\u2074)"           # mm⁴
    elif metric == "TMD":
        return "TMD (g·cm\u207b\u00b2)"    # g·cm⁻²
    elif metric == "pMOI":
        return "pMOI (g\u00b7cm\u207b\u00b3\u00b7mm\u2074)"         # g·mm²
    else:
        return metric

def main():
    # Resolve Excel path relative to this script
    here = os.path.dirname(os.path.abspath(__file__))
    excel_path = os.path.join(here, EXCEL_NAME)
    if not os.path.exists(excel_path):
        raise FileNotFoundError(f"Excel not found: {excel_path}")

    # Load group names and data (first two sheets)
    xls = pd.ExcelFile(excel_path)
    group_a, group_b = xls.sheet_names[:2]
    group_a_safe, group_b_safe = safe_name(group_a), safe_name(group_b)

    df_a_nohdr = pd.read_excel(excel_path, sheet_name=group_a, header=None)
    df_b_nohdr = pd.read_excel(excel_path, sheet_name=group_b, header=None)

    # Collect per-metric stats
    results = []
    for metric in METRICS:
        # Get raw matrices [slices x samples]
        a_mat = extract_metric_matrix(df_a_nohdr, metric, excel_path, group_a)
        b_mat = extract_metric_matrix(df_b_nohdr, metric, excel_path, group_b)

        # Choose a grid resolution based on both groups
        counts_a = np.sum(np.isfinite(a_mat), axis=0)
        counts_b = np.sum(np.isfinite(b_mat), axis=0)
        max_n = int(np.nanmax(np.concatenate([counts_a, counts_b])))
        max_n = max(max_n, 2)  # at least 2 points for a sensible line

        # Interpolate each specimen to a common % bone-length grid
        x_grid, a_grid = interpolate_to_common_percent_grid(
            a_mat, START_PERCENT, END_PERCENT, num_points=max_n
        )
        _,      b_grid = interpolate_to_common_percent_grid(
            b_mat, START_PERCENT, END_PERCENT, num_points=max_n
        )

        # Group-level mean ± SE on the % grid
        # (SE uses number of specimens with finite values at that point)
        a_mean = np.nanmean(a_grid, axis=1)
        b_mean = np.nanmean(b_grid, axis=1)

        def se_along_axis(grid):
            se = np.full(grid.shape[0], np.nan, dtype=float)
            for i in range(grid.shape[0]):
                vals = grid[i, :]
                vals = vals[np.isfinite(vals)]
                n = len(vals)
                if n > 1:
                    se[i] = np.nanstd(vals, ddof=1) / np.sqrt(n)
            return se

        a_se = se_along_axis(a_grid)
        b_se = se_along_axis(b_grid)

        # Per-position p-values
        pvals = per_grid_pvalues(a_grid, b_grid)

        results.append(dict(
            metric=metric,
            group_a=group_a,
            group_b=group_b,
            x_grid=x_grid,
            a_mean=a_mean, a_se=a_se, a_grid=a_grid, a_mat=a_mat,
            b_mean=b_mean, b_se=b_se, b_grid=b_grid, b_mat=b_mat,
            pvals=pvals,
        ))

    # === Slice-by-slice (now % bone length) figure ===
    fig = plt.figure(figsize=(16, 6), constrained_layout=True)
    # Increase wspace to reduce overlap between y-axis labels and adjacent plots
    gs = GridSpec(2, len(METRICS),
                  height_ratios=[6, 0.7],
                  hspace=0.08,
                  wspace=0.35)

    for j, res in enumerate(results):
        metric = res["metric"]
        accent = BORDER_TITLE_COLOR.get(metric, "#444444")

        x_grid = res["x_grid"]
        a_mean, a_se = res["a_mean"], res["a_se"]
        b_mean, b_se = res["b_mean"], res["b_se"]
        pvals        = res["pvals"]

        # Main line plot vs % bone length
        ax = fig.add_subplot(gs[0, j])

        ax.plot(x_grid, a_mean, label=res["group_a"], color=GROUP_A_COLOR)
        ax.fill_between(x_grid, a_mean - a_se, a_mean + a_se,
                        color=GROUP_A_COLOR, alpha=0.25)
        ax.plot(x_grid, b_mean, label=res["group_b"], color=GROUP_B_COLOR)
        ax.fill_between(x_grid, b_mean - b_se, b_mean + b_se,
                        color=GROUP_B_COLOR, alpha=0.25)

        # Titles/borders (warm accents)
        ax.set_title(metric, color=accent, fontsize=16, fontweight="bold")
        for sp in ax.spines.values():
            sp.set_color(accent)
            sp.set_linewidth(2.5)

        ax.set_ylabel(metric_ylabel(metric))
        ax.tick_params(labelbottom=False)
        ax.set_xlim(START_PERCENT, END_PERCENT)

        if j == len(METRICS) - 1:
            ax.legend(loc="upper right")

        # Heatmap strip for per-position p-values across 58–62%
        ax_hm = fig.add_subplot(gs[1, j], sharex=ax)
        pv = pvals.copy()
        pv[np.isnan(pv)] = 1.0  # treat missing as non-significant

        if len(x_grid) > 0:
            ax_hm.imshow(
                pv.reshape(1, -1),
                aspect="auto",
                cmap=CMAP,
                norm=NORM,
                extent=[START_PERCENT, END_PERCENT, 0, 1]
            )

        ax_hm.set_yticks([])
        ax_hm.set_xlabel("Bone Length (%)", fontsize=12, fontweight="bold")
        for sp in ax_hm.spines.values():
            sp.set_color(accent)
            sp.set_linewidth(2.5)

        # e.g. ticks at 58, 59, 60, 61, 62
        ax_hm.set_xticks(np.linspace(START_PERCENT, END_PERCENT, 5))

    # Save + show slice-by-slice figure
    fig_w, fig_h = fig.get_size_inches()
    fig_name = f"multi_panel_mean_se_{group_a_safe}_vs_{group_b_safe}.tiff"
    plt.savefig(os.path.join(here, fig_name), dpi=FIG_DPI, format="tiff")
    plt.show()

    # === Summary bar charts (display AND save; same height, half width) ===
    # Based on per-sample means across original slices
    summary_rows = []
    for res in results:
        metric = res["metric"]

        a_vals = np.nanmean(res["a_mat"], axis=0)
        b_vals = np.nanmean(res["b_mat"], axis=0)
        a_vals = a_vals[np.isfinite(a_vals)]
        b_vals = b_vals[np.isfinite(b_vals)]

        nA, nB = len(a_vals), len(b_vals)
        meanA, meanB = np.nanmean(a_vals), np.nanmean(b_vals)
        seA = np.nanstd(a_vals, ddof=1) / np.sqrt(nA) if nA > 1 else np.nan
        seB = np.nanstd(b_vals, ddof=1) / np.sqrt(nB) if nB > 1 else np.nan

        _, p_val = ttest_ind(a_vals, b_vals, equal_var=False) if nA > 1 and nB > 1 else (None, np.nan)

        summary_rows.append(dict(
            metric=metric,
            group_a_mean=meanA, group_a_se=seA,
            group_b_mean=meanB, group_b_se=seB,
            p_value=p_val
        ))

    df_summary = pd.DataFrame(summary_rows)
    bars_csv = f"bars_mean_{group_a_safe}_vs_{group_b_safe}.csv"
    df_summary.to_csv(os.path.join(here, bars_csv), index=False)

    bars_w, bars_h = fig_w / 2.0, fig_h
    fig2, axes = plt.subplots(1, len(METRICS),
                              figsize=(bars_w, bars_h),
                              constrained_layout=True)

    if len(METRICS) == 1:
        axes = [axes]

    for i, row in enumerate(summary_rows):
        metric = row["metric"]
        accent = BORDER_TITLE_COLOR.get(metric, "#444444")
        ax = axes[i]

        means  = [row["group_a_mean"], row["group_b_mean"]]
        ses    = [row["group_a_se"],   row["group_b_se"]]
        colors = [GROUP_A_COLOR, GROUP_B_COLOR]

        # Plot bars with error bars
        ax.bar([0, 1], means, yerr=ses, capsize=4, color=colors, width=0.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([group_a, group_b], fontsize=12, fontweight="bold")

        ax.set_ylabel(metric_ylabel(metric))
        ax.set_title(metric, fontsize=16, fontweight="bold", color=accent)

        for sp in ax.spines.values():
            sp.set_color(accent)
            sp.set_linewidth(2.5)

        ax.yaxis.grid(True, linestyle=":", alpha=0.6)
        ax.set_axisbelow(True)

        # Overall p-value annotation with refined placement
        p = row["p_value"]
        if np.isfinite(p):
            txt = f"p = {p:.3g}" if p >= 0.001 else "p < 0.001"

            # Compute top of each bar including its error (if finite)
            bar_tops = []
            for m, s in zip(means, ses):
                if np.isfinite(m):
                    top = m + (s if np.isfinite(s) else 0.0)
                    bar_tops.append(top)

            if bar_tops:
                max_top = max(bar_tops)
                ymin, ymax = ax.get_ylim()
                span = ymax - ymin if ymax > ymin else 1.0

                # Ensure some headroom above the tallest bar+error
                # We want at least 8% of the current span above the bar tops
                required_ymax = max_top + 0.08 * span
                if ymax < required_ymax:
                    ymax = required_ymax
                    ax.set_ylim(ymin, ymax)
                    ymin, ymax = ax.get_ylim()
                    span = ymax - ymin

                # Place text slightly above the bar top, but also a bit below the upper border
                y_text = max_top + 0.03 * span       # 3% above bar+error
                upper_margin = 0.04 * span           # keep 4% gap below the top border
                max_y_text = ymax - upper_margin
                if y_text > max_y_text:
                    y_text = max_y_text

                ax.text(0.5, y_text, txt,
                        ha="center", va="bottom",
                        fontsize=12, fontweight="bold",
                        clip_on=True)

    bars_fig = f"bars_mean_{group_a_safe}_vs_{group_b_safe}.tiff"
    plt.savefig(os.path.join(here, bars_fig), dpi=FIG_DPI, format="tiff")
    plt.show()

    print(f"Saved slice-by-slice figure: {fig_name}")
    print(f"Saved summary bars figure : {bars_fig}")
    print(f"Saved summary CSV         : {bars_csv}")

if __name__ == "__main__":
    main()
