#!/usr/bin/env python3
"""
Batch slice-by-slice TMD, cortical area (Ct.Ar), geometric inertia tensor,
and TMD-weighted inertia tensor from microCT images.

ROOT_DIR contains one subfolder per sample (e.g. Sample1, Sample2, ...).
Each sample folder contains a subfolder that holds the microCT slice images.

For each sample, this script:
  * Finds the image subfolder
  * Computes per-slice:
        - Ct.Ar (cortical bone area, mm^2)
        - TMD (g/cm^3)
        - geometric inertia (pMOA) and TMD-weighted inertia (pMOI)
  * Writes a per-sample CSV in ROOT_DIR:
        slice_metrics_TMD_full_inertia_<SampleName>.csv

Additionally, it creates a single Excel workbook in ROOT_DIR with one sheet
per sample, each containing only:
    Slice, File, Ct_Ar_mm2, TMD_g_per_cm3, pMOA_J_geom_mm4, pMOI_J_rho_cm
"""

import os
import csv
import numpy as np
import pandas as pd
from skimage.io import imread

# ------------------------- CONFIG ---------------------------------

# Root folder containing sample subfolders"
ROOT_DIR = r"C:\Users\JAW19\Desktop\Student Projects\2024-25 Projects\Arlo\SCI_SHAM_STUDY\SHAM data"   # <-- update this as needed

# File extension to process
IMAGE_EXT = ".bmp"   # change to ".tif", ".png", etc. if needed

# Pixel width in microns (change per dataset)
PIXEL_WIDTH_MICRONS = 6.89
PIXEL_WIDTH_MM = PIXEL_WIDTH_MICRONS / 1000.0  # pixel side length in mm

# Intensity thresholds for bone segmentation (change per dataset)
LOWER_THRESHOLD = 90
UPPER_THRESHOLD = 255

# Phantom BMD values (g/cm^3)
BMD_0_75 = 0.75
BMD_0_25 = 0.25

# Corresponding grayscale values from your calibration
CAL_GRAY_0_75 = 98.1016224
CAL_GRAY_0_25 = 47.6460718

# Name of combined Excel file (will be placed in ROOT_DIR)
EXCEL_NAME = "TMD_pMOA_pMOI_summary_SCI.xlsx"

# ------------------------------------------------------------------


def calibrate_bmd(gray: np.ndarray) -> np.ndarray:
    """Convert grayscale intensities to BMD (g/cm^3)."""
    slope = (BMD_0_75 - BMD_0_25) / (CAL_GRAY_0_75 - CAL_GRAY_0_25)
    intercept = BMD_0_75 - slope * CAL_GRAY_0_75
    return slope * gray + intercept


def principal_moments_and_angle(Ixx, Iyy, Ixy):
    """Compute principal moments and orientation angle from an inertia tensor."""
    I_avg = (Ixx + Iyy) / 2.0
    I_diff = (Ixx - Iyy) / 2.0
    R = np.sqrt(I_diff**2 + Ixy**2)
    Imax = I_avg + R
    Imin = I_avg - R
    theta = 0.5 * np.arctan2(2.0 * Ixy, Ixx - Iyy)
    return Imax, Imin, theta


def compute_slice_metrics(img, pixel_width_mm, lower_thr, upper_thr):
    """
    Compute Ct.Ar, TMD, geometric inertia tensor, and TMD-weighted inertia tensor
    for one slice.
    """

    # Ensure 2D grayscale
    if img.ndim == 3:
        img = img.mean(axis=2)
    img = img.astype(np.float32)

    # Bone mask
    bone_mask = (img >= lower_thr) & (img <= upper_thr)
    if not np.any(bone_mask):
        return None

    # Pixel geometry
    a = pixel_width_mm
    A_pix = a**2  # mm^2

    # Cortical area (Ct.Ar): number of bone pixels * pixel area (mm^2)
    N = int(np.count_nonzero(bone_mask))
    ct_ar_mm2 = N * A_pix

    # TMD over bone (g/cm^3)
    bmd_map = calibrate_bmd(img)
    bmd_bone = bmd_map[bone_mask]
    tmd = float(np.nanmean(bmd_bone))

    # Coordinates (mm)
    coords = np.column_stack(np.nonzero(bone_mask))
    row = coords[:, 0].astype(float)
    col = coords[:, 1].astype(float)
    x_mm = col * pixel_width_mm
    y_mm = row * pixel_width_mm

    # ------------------ GEOMETRIC inertia ------------------
    cx_geom = x_mm.mean()
    cy_geom = y_mm.mean()
    dx = x_mm - cx_geom
    dy = y_mm - cy_geom

    Ixx_geom = np.sum(dy**2 * A_pix) + N * (a**4 / 12.0)
    Iyy_geom = np.sum(dx**2 * A_pix) + N * (a**4 / 12.0)
    Ixy_geom = -np.sum(dx * dy * A_pix)

    J_geom = Ixx_geom + Iyy_geom
    Imax_geom, Imin_geom, theta_geom = principal_moments_and_angle(
        Ixx_geom, Iyy_geom, Ixy_geom
    )

    # ------------------ TMD-weighted (pMOI, 2D) ------------------
    rho = bmd_bone.astype(float)  # g/cm³
    weight_sum = np.sum(rho)

    cx_rho = np.sum(x_mm * rho) / weight_sum
    cy_rho = np.sum(y_mm * rho) / weight_sum

    dx_rho = x_mm - cx_rho
    dy_rho = y_mm - cy_rho

    Ixx_rho = np.sum(rho * dy_rho**2 * A_pix) + weight_sum * (a**4 / 12.0)
    Iyy_rho = np.sum(rho * dx_rho**2 * A_pix) + weight_sum * (a**4 / 12.0)
    Ixy_rho = -np.sum(rho * dx_rho * dy_rho * A_pix)

    J_rho = Ixx_rho + Iyy_rho
    Imax_rho, Imin_rho, theta_rho = principal_moments_and_angle(
        Ixx_rho, Iyy_rho, Ixy_rho
    )

    return dict(
        Ct_Ar_mm2=float(ct_ar_mm2),
        TMD_g_per_cm3=tmd,

        cx_geom_mm=float(cx_geom),
        cy_geom_mm=float(cy_geom),
        Ixx_geom_mm4=float(Ixx_geom),
        Iyy_geom_mm4=float(Iyy_geom),
        Ixy_geom_mm4=float(Ixy_geom),
        J_geom_mm4=float(J_geom),
        Imax_geom_mm4=float(Imax_geom),
        Imin_geom_mm4=float(Imin_geom),
        theta_geom_rad=float(theta_geom),

        cx_rho_mm=float(cx_rho),
        cy_rho_mm=float(cy_rho),
        Ixx_rho_cm=float(Ixx_rho),
        Iyy_rho_cm=float(Iyy_rho),
        Ixy_rho_cm=float(Ixy_rho),
        J_rho_cm=float(J_rho),
        Imax_rho_cm=float(Imax_rho),
        Imin_rho_cm=float(Imin_rho),
        theta_rho_rad=float(theta_rho),
    )


def natural_sort_key(s):
    import re
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def process_sample_folder(sample_dir, root_dir, excel_sheets):
    """Process images for each SampleX folder, supporting two layouts:
    
    1. Images directly in SampleX/
    2. Images inside a single subfolder (e.g. img_cort_*)
    """

    sample_name = os.path.basename(os.path.normpath(sample_dir))
    csv_name = f"slice_metrics_TMD_full_inertia_{sample_name}.csv"
    csv_path = os.path.join(root_dir, csv_name)

    print(f"\nSearching for microCT images inside: {sample_name}")

    # ---------------------------------------------------------
    # FIRST: check if images already exist directly in sample_dir
    # ---------------------------------------------------------
    direct_images = [
        f for f in os.listdir(sample_dir)
        if f.lower().endswith(IMAGE_EXT)
    ]

    if direct_images:
        image_dir = sample_dir
        print(f"[{sample_name}] Found images directly inside sample folder.")
    else:
        # ---------------------------------------------------------
        # SECOND: search subfolders (BCL3-style structure)
        # ---------------------------------------------------------
        subdirs = [
            os.path.join(sample_dir, d)
            for d in os.listdir(sample_dir)
            if os.path.isdir(os.path.join(sample_dir, d))
        ]

        image_dir = None
        for subdir in subdirs:
            img_files = [
                f for f in os.listdir(subdir)
                if f.lower().endswith(IMAGE_EXT)
            ]
            if img_files:
                image_dir = subdir
                print(f"[{sample_name}] Using image subfolder: {os.path.basename(subdir)}")
                break

        if image_dir is None:
            print(f"[{sample_name}] No images found — skipping.")
            return

    # ---------------------------------------------------------
    # Process images in image_dir
    # ---------------------------------------------------------
    files = [
        f for f in os.listdir(image_dir)
        if f.lower().endswith(IMAGE_EXT)
    ]
    files.sort(key=natural_sort_key)

    slice_numbers = []
    file_names = []
    all_results = []

    for i, fname in enumerate(files, start=1):
        img_path = os.path.join(image_dir, fname)
        img = imread(img_path)

        metrics = compute_slice_metrics(
            img,
            pixel_width_mm=PIXEL_WIDTH_MM,
            lower_thr=LOWER_THRESHOLD,
            upper_thr=UPPER_THRESHOLD,
        )

        if metrics is None:
            continue

        slice_numbers.append(i)
        file_names.append(fname)
        all_results.append(metrics)

    if not all_results:
        print(f"[{sample_name}] No valid bone slices — skipping.")
        return

    # ---------------------------------------------------------
    # Write per-sample CSV
    # ---------------------------------------------------------
    header = [
        "Slice", "File",
        "Ct_Ar_mm2", "TMD_g_per_cm3",
        "cx_geom_mm", "cy_geom_mm",
        "Ixx_geom_mm4", "Iyy_geom_mm4", "Ixy_geom_mm4",
        "J_geom_mm4", "Imax_geom_mm4", "Imin_geom_mm4", "theta_geom_rad",
        "cx_rho_mm", "cy_rho_mm",
        "Ixx_rho_cm", "Iyy_rho_cm", "Ixy_rho_cm",
        "J_rho_cm", "Imax_rho_cm", "Imin_rho_cm", "theta_rho_rad",
    ]

    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for i, fname, m in zip(slice_numbers, file_names, all_results):
            row = [i, fname] + [m[h] for h in header[2:]]
            w.writerow(row)

    print(f"[{sample_name}] CSV written: {csv_path}")

    # ---------------------------------------------------------
    # Add to Excel summary
    # ---------------------------------------------------------
    df = pd.DataFrame({
        "Slice": slice_numbers,
        "File": file_names,
        "Ct_Ar_mm2": [m["Ct_Ar_mm2"] for m in all_results],
        "TMD_g_per_cm3": [m["TMD_g_per_cm3"] for m in all_results],
        "pMOA_J_geom_mm4": [m["J_geom_mm4"] for m in all_results],
        "pMOI_J_rho_cm": [m["J_rho_cm"] for m in all_results],
    })

    excel_sheets[sample_name[:31]] = df



def main():
    # Find all immediate subfolders in ROOT_DIR
    subfolders = [
        os.path.join(ROOT_DIR, d)
        for d in os.listdir(ROOT_DIR)
        if os.path.isdir(os.path.join(ROOT_DIR, d))
    ]

    if not subfolders:
        print(f"No subfolders found in ROOT_DIR: {ROOT_DIR}")
        return

    print(f"Found {len(subfolders)} sample folders in {ROOT_DIR}")

    excel_sheets = {}

    for sample_dir in sorted(subfolders, key=os.path.basename):
        process_sample_folder(sample_dir, ROOT_DIR, excel_sheets)

    # -------- write combined Excel workbook --------
    if not excel_sheets:
        print("No Excel data collected (no valid samples).")
        return

    excel_path = os.path.join(ROOT_DIR, EXCEL_NAME)
    # Let pandas choose an available engine (e.g. openpyxl)
    with pd.ExcelWriter(excel_path) as writer:
        for sheet_name, df in excel_sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"\nCombined Excel workbook written to: {excel_path}")


if __name__ == "__main__":
    main()
