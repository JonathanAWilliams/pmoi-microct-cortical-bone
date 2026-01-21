"""
pMOA and pMOI contribution visualisations for TWO microCT samples with shared colour scales.

What you need to do
-------------------
1) Set SAMPLE_A_DIR and SAMPLE_B_DIR to point to your two slice folders.
   - Each can either contain the slice images directly, OR
   - contain a single subfolder that contains the images.

2) Set OUTPUT_ROOT to where you want PNG stacks and colourbars saved.

3) Confirm:
   - IMAGE_EXT matches your slice file type (e.g. .bmp / .tif / .png)
   - PIXEL_WIDTH_MICRONS matches the microCT pixel size
   - LOWER_THRESHOLD / UPPER_THRESHOLD match your cortical mask thresholds
   - Calibration values match your grayscale → BMD/TMD phantom calibration

What the script does
--------------------
1) First pass: scans both samples and computes per-pixel pMOA and pMOI contribution values
   for each slice, collecting global min/max for each metric across BOTH samples.

2) Second pass: writes transparent PNG heatmaps for each slice using the SAME colour scales
   (vmin/vmax) for pMOA and pMOI across both samples.

3) Saves one horizontal colourbar PNG for pMOA and one for pMOI in OUTPUT_ROOT.

Outputs
-------
OUTPUT_ROOT/
  <SAMPLE_A_LABEL>_pMOA/slice_001.png ...
  <SAMPLE_A_LABEL>_pMOI/slice_001.png ...
  <SAMPLE_B_LABEL>_pMOA/slice_001.png ...
  <SAMPLE_B_LABEL>_pMOI/slice_001.png ...
  colourbar_pMOA.png
  colourbar_pMOI.png

Notes / assumptions
-------------------
- Background is fully transparent (alpha=0) so stacks load nicely in 3D Slicer.
- pMOA contributions are computed about the geometric centroid.
- pMOI contributions are computed about the density-weighted centroid.
"""

import os
import numpy as np
from skimage.io import imread
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable


# ------------------------- USER CONFIG ---------------------------------
# Edit these paths and parameters for your dataset

# Two sample folders (edit these)
SAMPLE_A_DIR = r"D:\summer work\PI3K Study VOIs\PI3K-control\Sample31\VOI\cleaned_CORT_erosion"   # e.g. PI3K control
SAMPLE_B_DIR = r"D:\summer work\PI3K Study VOIs\PI3K\Sample13\VOI\cleaned_CORT_erosion"   # e.g. PI3K

# Short labels used to name the output folders
SAMPLE_A_LABEL = "PI3K_control"
SAMPLE_B_LABEL = "PI3K"


# Output folder where PNG stacks + colourbars will be written
OUTPUT_ROOT = r"C:\Users\JAW19\Desktop\Student Projects\2024-25 Projects\Arlo\PI3K Study"

# Image extension
IMAGE_EXT = ".bmp"     # or ".tif", ".png", ...

# Pixel width in microns (same for both samples in this study)
PIXEL_WIDTH_MICRONS = 9.041832
PIXEL_WIDTH_MM = PIXEL_WIDTH_MICRONS / 1000.0

# Threshold range used to define the cortical mask from grayscale
LOWER_THRESHOLD = 70
UPPER_THRESHOLD = 255

# Phantom BMD calibration (g/cm³)
BMD_0_75 = 0.75
BMD_0_25 = 0.25
CAL_GRAY_0_75 = 92.4032819
CAL_GRAY_0_25 = 43.0463712

# ----------------------------------------------------------------------


def calibrate_bmd(gray: np.ndarray) -> np.ndarray:
    """Convert grayscale intensities to BMD (g/cm^3)."""
    slope = (BMD_0_75 - BMD_0_25) / (CAL_GRAY_0_75 - CAL_GRAY_0_25)
    intercept = BMD_0_75 - slope * CAL_GRAY_0_75
    return slope * gray + intercept


def find_image_folder(sample_dir: str) -> str | None:
    """
    Given a sample directory, either:
    - return sample_dir if it already contains images, OR
    - search its subdirectories and return the first one containing images.

    Returns None if no images found.
    """
    # Images directly in sample_dir?
    direct_imgs = [f for f in os.listdir(sample_dir)
                   if f.lower().endswith(IMAGE_EXT)]
    if direct_imgs:
        return sample_dir

    # Otherwise, search subdirectories
    for d in os.listdir(sample_dir):
        sub = os.path.join(sample_dir, d)
        if os.path.isdir(sub):
            imgs = [f for f in os.listdir(sub)
                    if f.lower().endswith(IMAGE_EXT)]
            if imgs:
                return sub

    return None


def natural_sort_key(s: str):
    import re
    return [int(t) if t.isdigit() else t.lower()
            for t in re.split(r"(\d+)", s)]


def compute_contribution_maps(img: np.ndarray,
                              pixel_width_mm: float,
                              lower_thr: int,
                              upper_thr: int):
    """
    For a single slice:
        - segment bone
        - compute pMOA contribution per pixel (mm^4)
        - compute pMOI contribution per pixel (g/cm^3 * mm^4)

    Returns:
        bone_mask : (H, W) bool
        pmoa_contrib : (H, W) float
        pmoi_contrib : (H, W) float
    """

    # Ensure 2D grayscale
    if img.ndim == 3:
        img = img.mean(axis=2)
    img = img.astype(np.float32)

    bone_mask = (img >= lower_thr) & (img <= upper_thr)
    if not np.any(bone_mask):
        return bone_mask, np.zeros_like(img), np.zeros_like(img)

    # Pixel geometry
    a = pixel_width_mm
    A_pix = a**2  # mm^2

    # BMD in g/cm^3
    bmd_map = calibrate_bmd(img)
    rho_bone = bmd_map[bone_mask].astype(float)

    # Coordinates in mm
    coords = np.column_stack(np.nonzero(bone_mask))
    rows = coords[:, 0].astype(float)
    cols = coords[:, 1].astype(float)
    x_mm = cols * pixel_width_mm
    y_mm = rows * pixel_width_mm

    # Geometric centroid (for pMOA)
    cx_geom = x_mm.mean()
    cy_geom = y_mm.mean()
    dx_geom = x_mm - cx_geom
    dy_geom = y_mm - cy_geom
    r2_geom = dx_geom**2 + dy_geom**2

    # pMOA contribution: r^2 * A_pix
    pmoa_contrib = np.zeros_like(img, dtype=float)
    pmoa_contrib[coords[:, 0], coords[:, 1]] = r2_geom * A_pix  # mm^4

    # TMD-weighted centroid (for pMOI)
    rho_sum = np.sum(rho_bone)
    cx_rho = np.sum(x_mm * rho_bone) / rho_sum
    cy_rho = np.sum(y_mm * rho_bone) / rho_sum
    dx_rho = x_mm - cx_rho
    dy_rho = y_mm - cy_rho
    r2_rho = dx_rho**2 + dy_rho**2

    # pMOI contribution: rho * r^2 * A_pix   (g/cm^3 * mm^4)
    pmoi_contrib = np.zeros_like(img, dtype=float)
    pmoi_contrib[coords[:, 0], coords[:, 1]] = rho_bone * r2_rho * A_pix

    return bone_mask, pmoa_contrib, pmoi_contrib


def scan_sample_for_global_ranges(sample_dir: str):
    """
    First pass: for all slices in a sample, compute pMOA/pMOI contributions
    and return:
        (pmoa_min, pmoa_max, pmoi_min, pmoi_max)
    over bone pixels only.
    """
    img_dir = find_image_folder(sample_dir)
    if img_dir is None:
        print(f"[{sample_dir}] No images found.")
        return None

    files = sorted(
        [f for f in os.listdir(img_dir) if f.lower().endswith(IMAGE_EXT)],
        key=natural_sort_key
    )

    pmoa_vals = []
    pmoi_vals = []

    for fname in files:
        img_path = os.path.join(img_dir, fname)
        img = imread(img_path)

        bone_mask, pmoa_contrib, pmoi_contrib = compute_contribution_maps(
            img,
            pixel_width_mm=PIXEL_WIDTH_MM,
            lower_thr=LOWER_THRESHOLD,
            upper_thr=UPPER_THRESHOLD,
        )

        if np.any(bone_mask):
            pmoa_vals.append(pmoa_contrib[bone_mask])
            pmoi_vals.append(pmoi_contrib[bone_mask])

    if not pmoa_vals or not pmoi_vals:
        print(f"[{sample_dir}] No bone pixels found in any slice.")
        return None

    pmoa_all = np.concatenate(pmoa_vals)
    pmoi_all = np.concatenate(pmoi_vals)

    return (
        float(pmoa_all.min()),
        float(pmoa_all.max()),
        float(pmoi_all.min()),
        float(pmoi_all.max()),
    )


def visualize_sample(sample_dir: str,
                     sample_label: str,
                     out_root: str,
                     vmin_pmoa: float,
                     vmax_pmoa: float,
                     vmin_pmoi: float,
                     vmax_pmoi: float):
    """
    Second pass: for a single sample, recompute contribution maps and
    save transparent PNGs for pMOA and pMOI using the provided global
    colour ranges.
    """

    img_dir = find_image_folder(sample_dir)
    if img_dir is None:
        print(f"[{sample_dir}] No images found, skipping visualisation.")
        return

    files = sorted(
        [f for f in os.listdir(img_dir) if f.lower().endswith(IMAGE_EXT)],
        key=natural_sort_key
    )

    # Output folders
    out_pmoa = os.path.join(out_root, f"{sample_label}_pMOA")
    out_pmoi = os.path.join(out_root, f"{sample_label}_pMOI")
    os.makedirs(out_pmoa, exist_ok=True)
    os.makedirs(out_pmoi, exist_ok=True)

    for i, fname in enumerate(files, start=1):
        img_path = os.path.join(img_dir, fname)
        img = imread(img_path)

        bone_mask, pmoa_contrib, pmoi_contrib = compute_contribution_maps(
            img,
            pixel_width_mm=PIXEL_WIDTH_MM,
            lower_thr=LOWER_THRESHOLD,
            upper_thr=UPPER_THRESHOLD,
        )

        # pMOA heatmap (transparent background)
        if np.any(bone_mask):
            norm_pmoa = plt.Normalize(vmin=vmin_pmoa, vmax=vmax_pmoa)
            cmap = plt.get_cmap("hot")

            rgba_pmoa = cmap(norm_pmoa(pmoa_contrib))
            rgba_pmoa[..., 3] = bone_mask.astype(float)  # alpha mask

            fig, ax = plt.subplots(figsize=(8, 6))
            ax.imshow(rgba_pmoa)
            ax.axis("off")
            ax.set_aspect("equal")
            out_path_pmoa = os.path.join(out_pmoa, f"slice_{i:03d}.png")
            plt.savefig(out_path_pmoa, dpi=300, bbox_inches="tight",
                        pad_inches=0, transparent=True)
            plt.close(fig)

            # pMOI heatmap
            norm_pmoi = plt.Normalize(vmin=vmin_pmoi, vmax=vmax_pmoi)
            rgba_pmoi = cmap(norm_pmoi(pmoi_contrib))
            rgba_pmoi[..., 3] = bone_mask.astype(float)

            fig, ax = plt.subplots(figsize=(8, 6))
            ax.imshow(rgba_pmoi)
            ax.axis("off")
            ax.set_aspect("equal")
            out_path_pmoi = os.path.join(out_pmoi, f"slice_{i:03d}.png")
            plt.savefig(out_path_pmoi, dpi=300, bbox_inches="tight",
                        pad_inches=0, transparent=True)
            plt.close(fig)

    print(f"[{sample_label}] pMOA images saved in: {out_pmoa}")
    print(f"[{sample_label}] pMOI images saved in: {out_pmoi}")


def save_colourbar(vmin: float, vmax: float,
                   out_path: str,
                   label: str,
                   units: str = ""):
    """
    Save a horizontal colourbar for the given vmin/vmax using the 'hot' colormap,
    labelled with actual numeric range.
    """
    fig, ax = plt.subplots(figsize=(4, 0.6))
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(norm=norm, cmap="hot")
    sm.set_array([])

    cbar = fig.colorbar(sm, cax=ax, orientation="horizontal")
    if units:
        full_label = f"{label} ({units})"
    else:
        full_label = label
    cbar.set_label(full_label, fontsize=9)

    # Ticks at min, mid, max
    mid = 0.5 * (vmin + vmax)
    cbar.set_ticks([vmin, mid, vmax])
    cbar.set_ticklabels([f"{vmin:.2e}", f"{mid:.2e}", f"{vmax:.2e}"])

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0, transparent=True)
    plt.close(fig)


def main():
    os.makedirs(OUTPUT_ROOT, exist_ok=True)

    print("=== First pass: scanning ranges for Sample A ===")
    ranges_a = scan_sample_for_global_ranges(SAMPLE_A_DIR)
    if ranges_a is None:
        print("Sample A had no valid bone pixels; aborting.")
        return

    print("=== First pass: scanning ranges for Sample B ===")
    ranges_b = scan_sample_for_global_ranges(SAMPLE_B_DIR)
    if ranges_b is None:
        print("Sample B had no valid bone pixels; aborting.")
        return

    # Unpack
    pmoa_min_a, pmoa_max_a, pmoi_min_a, pmoi_max_a = ranges_a
    pmoa_min_b, pmoa_max_b, pmoi_min_b, pmoi_max_b = ranges_b

    # Global ranges across both samples
    vmin_pmoa = min(pmoa_min_a, pmoa_min_b)
    vmax_pmoa = max(pmoa_max_a, pmoa_max_b)
    vmin_pmoi = min(pmoi_min_a, pmoi_min_b)
    vmax_pmoi = max(pmoi_max_a, pmoi_max_b)

    print("\nGlobal pMOA range:", vmin_pmoa, "to", vmax_pmoa)
    print("Global pMOI range:", vmin_pmoi, "to", vmax_pmoi)

    print("\n=== Second pass: visualising Sample A ===")
    visualize_sample(
        SAMPLE_A_DIR,
        SAMPLE_A_LABEL,
        OUTPUT_ROOT,
        vmin_pmoa,
        vmax_pmoa,
        vmin_pmoi,
        vmax_pmoi,
    )

    print("\n=== Second pass: visualising Sample B ===")
    visualize_sample(
        SAMPLE_B_DIR,
        SAMPLE_B_LABEL,
        OUTPUT_ROOT,
        vmin_pmoa,
        vmax_pmoa,
        vmin_pmoi,
        vmax_pmoi,
    )

    # Save colourbars in OUTPUT_ROOT
    cb_pmoa_path = os.path.join(OUTPUT_ROOT, "colourbar_pMOA.png")
    cb_pmoi_path = os.path.join(OUTPUT_ROOT, "colourbar_pMOI.png")

    save_colourbar(
        vmin=vmin_pmoa,
        vmax=vmax_pmoa,
        out_path=cb_pmoa_path,
        label="pMOA contribution",
        units="mm\u2074",
    )
    save_colourbar(
        vmin=vmin_pmoi,
        vmax=vmax_pmoi,
        out_path=cb_pmoi_path,
        label="pMOI contribution",
        units="g·cm\u207b\u00b3·mm\u2074",
    )

    print("\nColourbars saved to:")
    print("  ", cb_pmoa_path)
    print("  ", cb_pmoi_path)
    print("\nDone.")


if __name__ == "__main__":
    main()
