#!/usr/bin/env python3
"""compute_slice_metrics_batch.py

Batch slice-by-slice cortical area (Ct.Ar), mean tissue mineral density (TMD),
geometric polar second moment of area (pMOA), and density-weighted polar second
moment (pMOI-like) from microCT images.

This is the reference Python implementation used in the accompanying study.

Overview
--------
- Expects a ROOT directory containing one subfolder per sample (e.g. Sample1, Sample2, ...).
- Within each sample folder, images may be nested arbitrarily (e.g., Sample35/VOI/cleaned_CORT_erosion).
- For each sample, the script:
    1) Recursively finds an image folder (prefers a folder named 'cleaned_CORT_erosion')
    2) Thresholds each slice to create a cortical bone mask
    3) Computes per-slice metrics and writes:
        - A per-sample CSV (full per-slice outputs)
        - A combined Excel workbook (one sheet per sample + a sheet of per-sample means)

Important conventions
---------------------
Coordinate convention:
- NumPy image indices: (row=y, col=x)
- Physical coordinates are computed in *mm* using the supplied pixel width.

Pixel self-term:
- The centroidal second moments include the square-pixel self-term (a^4/12) per pixel.

Units:
- Geometric moments (Ixx_geom, Iyy_geom, J_geom) are in mm^4.
- Density-weighted moments multiply by rho (TMD) in g/cm^3, so the units are:
    (g/cm^3) * mm^4
  These are *not* directly comparable in units to pMOA; they are a density-scaled second moment.

Outputs
-------
- Per-sample CSV written to ROOT:
    slice_metrics_full_<SampleName>.csv
- Combined Excel written to ROOT:
    <excel-name>.xlsx

Examples
--------
Dry run (show chosen image folders without computing):
    python compute_slice_metrics_batch.py --root ./data --dry-run

Process one sample only:
    python compute_slice_metrics_batch.py --root ./data --sample Sample35

"""

from __future__ import annotations

import argparse
import csv
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from skimage.io import imread


# ------------------------- Calibration -------------------------

@dataclass(frozen=True)
class TwoPointCalibration:
    """Linear calibration mapping grayscale -> BMD/TMD."""
    cal_gray_low: float
    cal_gray_high: float
    bmd_low: float
    bmd_high: float

    @property
    def slope(self) -> float:
        return (self.bmd_high - self.bmd_low) / (self.cal_gray_high - self.cal_gray_low)

    @property
    def intercept(self) -> float:
        return self.bmd_high - self.slope * self.cal_gray_high

    def apply(self, gray: np.ndarray) -> np.ndarray:
        return self.slope * gray + self.intercept


def principal_moments_and_angle(Ixx: float, Iyy: float, Ixy: float) -> Tuple[float, float, float]:
    """Return (Imax, Imin, theta) for a 2x2 centroidal second-moment tensor."""
    I_avg = (Ixx + Iyy) / 2.0
    I_diff = (Ixx - Iyy) / 2.0
    R = float(np.sqrt(I_diff * I_diff + Ixy * Ixy))
    Imax = I_avg + R
    Imin = I_avg - R
    theta = 0.5 * float(np.arctan2(2.0 * Ixy, Ixx - Iyy))
    return Imax, Imin, theta


def compute_slice_metrics(
    img: np.ndarray,
    pixel_width_mm: float,
    lower_thr: float,
    upper_thr: float,
    calib: TwoPointCalibration,
) -> Optional[Dict[str, float]]:
    """Compute Ct.Ar, TMD, geometric and density-weighted second moments for one 2D slice."""

    # Ensure 2D grayscale
    if img.ndim == 3:
        img = img.mean(axis=2)
    img = img.astype(np.float32)

    # Bone mask
    bone_mask = (img >= lower_thr) & (img <= upper_thr)
    if not np.any(bone_mask):
        return None

    # Pixel geometry
    a = float(pixel_width_mm)
    A_pix = a * a  # mm^2

    # Cortical area (Ct.Ar)
    N = int(np.count_nonzero(bone_mask))
    ct_ar_mm2 = N * A_pix

    # TMD over bone (g/cm^3)
    bmd_map = calib.apply(img)
    bmd_bone = bmd_map[bone_mask]
    tmd = float(np.nanmean(bmd_bone))

    # Coordinates (mm)
    coords = np.column_stack(np.nonzero(bone_mask))
    row = coords[:, 0].astype(float)
    col = coords[:, 1].astype(float)
    x_mm = col * a
    y_mm = row * a

    # ------------------ GEOMETRIC second moments ------------------
    cx_geom = float(x_mm.mean())
    cy_geom = float(y_mm.mean())
    dx = x_mm - cx_geom
    dy = y_mm - cy_geom

    # Include square-pixel self-term a^4/12 per pixel
    pix_self = (a ** 4) / 12.0
    Ixx_geom = float(np.sum(dy * dy * A_pix) + N * pix_self)
    Iyy_geom = float(np.sum(dx * dx * A_pix) + N * pix_self)
    Ixy_geom = float(-np.sum(dx * dy * A_pix))

    J_geom = Ixx_geom + Iyy_geom
    Imax_geom, Imin_geom, theta_geom = principal_moments_and_angle(Ixx_geom, Iyy_geom, Ixy_geom)

    # ------------------ DENSITY-WEIGHTED second moments ------------------
    rho = bmd_bone.astype(float)  # g/cm^3
    weight_sum = float(np.sum(rho))

    cx_rho = float(np.sum(x_mm * rho) / weight_sum)
    cy_rho = float(np.sum(y_mm * rho) / weight_sum)

    dx_rho = x_mm - cx_rho
    dy_rho = y_mm - cy_rho

    # Units: (g/cm^3) * mm^4
    Ixx_rho = float(np.sum(rho * (dy_rho ** 2) * A_pix) + weight_sum * pix_self)
    Iyy_rho = float(np.sum(rho * (dx_rho ** 2) * A_pix) + weight_sum * pix_self)
    Ixy_rho = float(-np.sum(rho * dx_rho * dy_rho * A_pix))

    J_rho = Ixx_rho + Iyy_rho
    Imax_rho, Imin_rho, theta_rho = principal_moments_and_angle(Ixx_rho, Iyy_rho, Ixy_rho)

    return {
        "Ct_Ar_mm2": float(ct_ar_mm2),
        "TMD_g_per_cm3": float(tmd),

        "cx_geom_mm": float(cx_geom),
        "cy_geom_mm": float(cy_geom),
        "Ixx_geom_mm4": float(Ixx_geom),
        "Iyy_geom_mm4": float(Iyy_geom),
        "Ixy_geom_mm4": float(Ixy_geom),
        "J_geom_mm4": float(J_geom),
        "Imax_geom_mm4": float(Imax_geom),
        "Imin_geom_mm4": float(Imin_geom),
        "theta_geom_rad": float(theta_geom),

        "cx_rho_mm": float(cx_rho),
        "cy_rho_mm": float(cy_rho),
        "Ixx_rho_g_per_cm3_mm4": float(Ixx_rho),
        "Iyy_rho_g_per_cm3_mm4": float(Iyy_rho),
        "Ixy_rho_g_per_cm3_mm4": float(Ixy_rho),
        "J_rho_g_per_cm3_mm4": float(J_rho),
        "Imax_rho_g_per_cm3_mm4": float(Imax_rho),
        "Imin_rho_g_per_cm3_mm4": float(Imin_rho),
        "theta_rho_rad": float(theta_rho),
    }


# ------------------------- IO helpers -------------------------

def natural_sort_key(s: str):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def find_image_folder(sample_dir: str, image_ext: str) -> Optional[str]:
    """Find the image folder within a sample directory.

    Preference:
      1) Any folder named 'cleaned_CORT_erosion' that contains images.
      2) Otherwise, the deepest folder that contains images.
    """
    candidates: List[str] = []

    for root, _dirs, files in os.walk(sample_dir):
        img_files = [f for f in files if f.lower().endswith(image_ext.lower())]
        if img_files:
            candidates.append(root)

    if not candidates:
        return None

    for path in candidates:
        if os.path.basename(path).lower() == "cleaned_cort_erosion":
            return path

    return max(candidates, key=lambda p: p.count(os.sep))


def iter_sample_dirs(root_dir: str, sample_filter: Optional[str]) -> List[str]:
    """Return list of sample directories under root_dir.

    If sample_filter is provided, it may be:
    - an exact folder name (e.g., 'Sample35')
    - a regex pattern (e.g., '^Sample(3|4)$')
    """
    subfolders = [
        os.path.join(root_dir, d)
        for d in os.listdir(root_dir)
        if os.path.isdir(os.path.join(root_dir, d))
    ]

    if sample_filter:
        # Try exact match first
        exact = [p for p in subfolders if os.path.basename(p) == sample_filter]
        if exact:
            return sorted(exact, key=os.path.basename)

        # Otherwise treat as regex
        pat = re.compile(sample_filter)
        subfolders = [p for p in subfolders if pat.search(os.path.basename(p))]

    return sorted(subfolders, key=os.path.basename)


# ------------------------- Processing -------------------------

def process_sample_folder(
    sample_dir: str,
    root_dir: str,
    image_ext: str,
    pixel_width_mm: float,
    lower_thr: float,
    upper_thr: float,
    calib: TwoPointCalibration,
    write_csv: bool = True,
) -> Tuple[Optional[pd.DataFrame], Optional[Dict[str, float]], Optional[str]]:
    """Process a sample folder. Returns (df_slice, sample_means, image_dir_used)."""

    sample_name = os.path.basename(os.path.normpath(sample_dir))

    image_dir = find_image_folder(sample_dir, image_ext)
    if image_dir is None:
        print(f"[{sample_name}] No images found — skipping.")
        return None, None, None

    files = [f for f in os.listdir(image_dir) if f.lower().endswith(image_ext.lower())]
    files.sort(key=natural_sort_key)

    slice_numbers: List[int] = []
    file_names: List[str] = []
    all_results: List[Dict[str, float]] = []

    for i, fname in enumerate(files, start=1):
        img_path = os.path.join(image_dir, fname)
        img = imread(img_path)

        metrics = compute_slice_metrics(
            img,
            pixel_width_mm=pixel_width_mm,
            lower_thr=lower_thr,
            upper_thr=upper_thr,
            calib=calib,
        )
        if metrics is None:
            continue

        slice_numbers.append(i)
        file_names.append(fname)
        all_results.append(metrics)

    if not all_results:
        print(f"[{sample_name}] No valid bone slices — skipping.")
        return None, None, image_dir

    # Full per-slice output table
    df_full = pd.DataFrame(all_results)
    df_full.insert(0, "File", file_names)
    df_full.insert(0, "Slice", slice_numbers)

    # Optional per-sample CSV (full)
    if write_csv:
        csv_name = f"slice_metrics_full_{sample_name}.csv"
        csv_path = os.path.join(root_dir, csv_name)
        df_full.to_csv(csv_path, index=False)
        print(f"[{sample_name}] CSV written: {csv_path}")

    # Minimal slice-by-slice table for Excel (paper-friendly)
    df_min = pd.DataFrame({
        "Slice": slice_numbers,
        "File": file_names,
        "Ct_Ar_mm2": df_full["Ct_Ar_mm2"].values,
        "TMD_g_per_cm3": df_full["TMD_g_per_cm3"].values,
        "pMOA_J_geom_mm4": df_full["J_geom_mm4"].values,
        "pMOI_J_rho_g_per_cm3_mm4": df_full["J_rho_g_per_cm3_mm4"].values,
    })

    sample_means = {
        "Sample": sample_name,
        "n_slices": int(len(df_min)),
        "Ct_Ar_mm2_mean": float(np.nanmean(df_min["Ct_Ar_mm2"])),
        "TMD_g_per_cm3_mean": float(np.nanmean(df_min["TMD_g_per_cm3"])),
        "pMOA_J_geom_mm4_mean": float(np.nanmean(df_min["pMOA_J_geom_mm4"])),
        "pMOI_J_rho_g_per_cm3_mm4_mean": float(np.nanmean(df_min["pMOI_J_rho_g_per_cm3_mm4"])),
    }

    return df_min, sample_means, image_dir


# ------------------------- CLI -------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Batch compute Ct.Ar, TMD, pMOA and density-weighted second moments from microCT slices.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument("--root", required=True, help="Root directory containing one subfolder per sample")
    p.add_argument("--ext", default=".bmp", help="Image file extension to process")
    p.add_argument("--pixel-width-microns", type=float, required=True, help="Pixel width (microns)")
    p.add_argument("--thr-low", type=float, default=70, help="Lower intensity threshold for bone mask")
    p.add_argument("--thr-high", type=float, default=255, help="Upper intensity threshold for bone mask")

    # Two-point calibration
    p.add_argument("--cal-gray-low", type=float, required=True, help="Grayscale value for low calibration point")
    p.add_argument("--cal-gray-high", type=float, required=True, help="Grayscale value for high calibration point")
    p.add_argument("--bmd-low", type=float, required=True, help="BMD/TMD value at low calibration point (g/cm^3)")
    p.add_argument("--bmd-high", type=float, required=True, help="BMD/TMD value at high calibration point (g/cm^3)")

    p.add_argument("--excel-name", default="TMD_pMOA_pMOI_summary.xlsx", help="Combined Excel workbook filename")
    p.add_argument("--sample", default=None, help="Process only one sample folder (exact name) or regex pattern")
    p.add_argument("--dry-run", action="store_true", help="List discovered samples and image folders without computing")
    p.add_argument("--no-csv", action="store_true", help="Do not write per-sample CSVs")

    return p


def main() -> int:
    args = build_arg_parser().parse_args()

    root_dir = os.path.abspath(args.root)
    if not os.path.isdir(root_dir):
        raise SystemExit(f"--root does not exist or is not a directory: {root_dir}")

    pixel_width_mm = float(args.pixel_width_microns) / 1000.0
    calib = TwoPointCalibration(
        cal_gray_low=float(args.cal_gray_low),
        cal_gray_high=float(args.cal_gray_high),
        bmd_low=float(args.bmd_low),
        bmd_high=float(args.bmd_high),
    )

    sample_dirs = iter_sample_dirs(root_dir, args.sample)
    if not sample_dirs:
        print(f"No sample subfolders found under: {root_dir}")
        return 1

    print(f"Found {len(sample_dirs)} sample folder(s) under: {root_dir}")

    if args.dry_run:
        for sd in sample_dirs:
            sn = os.path.basename(sd)
            img_dir = find_image_folder(sd, args.ext)
            if img_dir is None:
                print(f"[{sn}] no images with ext {args.ext}")
            else:
                rel = os.path.relpath(img_dir, sd)
                print(f"[{sn}] images: {rel}")
        return 0

    excel_sheets: Dict[str, pd.DataFrame] = {}
    sample_summaries: List[Dict[str, float]] = []

    for sd in sample_dirs:
        sn = os.path.basename(sd)
        print(f"\nProcessing: {sn}")

        df_min, means, img_dir = process_sample_folder(
            sd,
            root_dir=root_dir,
            image_ext=args.ext,
            pixel_width_mm=pixel_width_mm,
            lower_thr=float(args.thr_low),
            upper_thr=float(args.thr_high),
            calib=calib,
            write_csv=(not args.no_csv),
        )

        if df_min is None or means is None:
            continue

        # Excel worksheet names max 31 chars
        excel_sheets[sn[:31]] = df_min
        sample_summaries.append(means)

    if not excel_sheets:
        print("No valid samples processed; nothing to write.")
        return 2

    excel_path = os.path.join(root_dir, args.excel_name)
    with pd.ExcelWriter(excel_path) as writer:
        for sheet_name, df in excel_sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
        pd.DataFrame(sample_summaries).to_excel(writer, sheet_name="sample_means", index=False)

    print(f"\nCombined Excel workbook written: {excel_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
