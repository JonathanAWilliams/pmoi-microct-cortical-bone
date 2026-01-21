#!/usr/bin/env python3
"""
SliceGeometry validation script (centroids + J only): compare Python vs BoneJ on two test images.

What we validate:
  - Geometric centroid (cx, cy)
  - Density-weighted centroid (cx, cy) using calibration m=1, c=0
  - Polar moment invariant J

What we do NOT validate:
  - Imin/Imax (not reported in manuscript; sensitive to convention details)

Important note on centroid conventions:
  BoneJ and this Python script may differ by a constant 0.5 px offset in reported centroids
  due to pixel-coordinate convention (pixel-center vs pixel-index reporting).
  This script reports both:
    - raw Python centroid
    - Python centroid shifted by (-0.5, -0.5) for direct comparison to BoneJ tables

J should agree closely regardless (rotational invariant; robust check).

Folder:
  C:\\Users\\JAW19\\Desktop\\Student Projects\\2024-25 Projects\\Arlo\\GitHub - BoneJ Implementation\\slicegeometry_validation_cases

Cases:
  - CaseA_2x2_block_4x4
  - CaseB_two_pixel_row_3x3_100_200
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

# Prefer imageio; fall back to skimage if needed
try:
    import imageio.v3 as iio
    def imread(path: str) -> np.ndarray:
        return iio.imread(path)
except Exception:
    from skimage.io import imread  # type: ignore


# ----------------------------
# CONFIG
# ----------------------------
ROOT_DIR = Path(__file__).resolve().parent / "slicegeometry_validation_cases"

CASE_BASENAMES = [
    "CaseA_2x2_block_4x4",
    "CaseB_two_pixel_row_3x3_100_200",
]

EXTS = [".tif", ".tiff", ".png", ".bmp", ".jpg", ".jpeg"]

# Foreground definition: pixels > 0 are included
FOREGROUND_THRESHOLD = 0.0

# Pixel size: keep 1.0 to match BoneJ px^n outputs
PIXEL_SIZE = 1.0

# Density calibration used in BoneJ for these checks
DENSITY_M = 1.0
DENSITY_C = 0.0

# Centroid convention adjustment:
# Observed: Python centroids are +0.5 px vs BoneJ, likely due to coordinate convention.
# We'll compare BoneJ to (Python - 0.5) to get PASS, while still reporting raw centroid.
CENTROID_SHIFT_TO_MATCH_BONEJ = (-0.5, -0.5)

# PASS/FAIL tolerances
ABS_TOL_CENTROID = 1e-3   # centroid comparisons after convention shift
ABS_TOL_J = 1e-3          # BoneJ table shown to 3 dp


# ----------------------------
# BoneJ reference values (from your screenshots)
# ----------------------------
BONEJ_REFERENCE: Dict[Tuple[str, str], Dict[str, float]] = {
    # -------- Case A --------
    ("CaseA_2x2_block_4x4", "geometric"): {
        "cx": 1.500,
        "cy": 1.500,
        "J": 2.667,
    },
    ("CaseA_2x2_block_4x4", "density"): {
        "cx": 1.500,
        "cy": 1.500,
        "J": 680.000,
    },

    # -------- Case B --------
    ("CaseB_two_pixel_row_3x3_100_200", "geometric"): {
        "cx": 1.500,
        "cy": 1.000,
        "J": 0.833,
    },
    ("CaseB_two_pixel_row_3x3_100_200", "density"): {
        "cx": 1.667,
        "cy": 1.000,
        "J": 116.667,
    },
}


# ----------------------------
# Helpers
# ----------------------------
def find_case_file(root: Path, basename: str) -> Path:
    """Find image file by stem (handles hidden extensions in Explorer)."""
    for ext in EXTS:
        p = root / f"{basename}{ext}"
        if p.exists():
            return p

    candidates: List[Path] = []
    for p in root.iterdir():
        if p.is_file() and p.suffix.lower() in EXTS and p.stem.lower() == basename.lower():
            candidates.append(p)

    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        return sorted(candidates)[0]

    raise FileNotFoundError(f"Could not find image for '{basename}' in {root}")


def to_2d_grayscale(img: np.ndarray) -> np.ndarray:
    """Ensure 2D grayscale float array."""
    if img.ndim == 2:
        arr = img
    elif img.ndim == 3:
        arr = img[..., 0]
    else:
        raise ValueError(f"Unsupported image shape: {img.shape}")
    return arr.astype(np.float64)


def compute_centroid_and_J(
    img: np.ndarray,
    mode: str,
    foreground_threshold: float,
    pixel_size: float,
    density_m: float,
    density_c: float,
) -> Dict[str, float]:
    """
    Compute centroid and polar moment J for a 2D image.

    mode:
      - "geometric": weight = 1 for foreground pixels
      - "density":   weight = (m*intensity + c) for foreground pixels; here m=1,c=0 -> intensity

    Includes pixel self-terms for square pixels of side length a:
      Ixx_self = Iyy_self = (a^4)/12 * weight
      Ixy_self = 0

    J is computed as Ixx + Iyy (equivalently Imin+Imax), and should match BoneJ.
    """
    arr = to_2d_grayscale(img)
    mask = arr > foreground_threshold

    if mode == "geometric":
        w = mask.astype(np.float64)
    elif mode == "density":
        dens = density_m * arr + density_c
        w = np.where(mask, dens, 0.0)
    else:
        raise ValueError("mode must be 'geometric' or 'density'")

    a = float(pixel_size)
    pixel_area = a * a

    W = float(np.sum(w) * pixel_area)
    if W == 0.0:
        return {"W": 0.0, "cx": np.nan, "cy": np.nan, "J": np.nan}

    h, ww = w.shape
    ys, xs = np.mgrid[0:h, 0:ww]

    # Pixel-center coordinate convention (matches your earlier Python implementation)
    x = (xs + 0.5) * a
    y = (ys + 0.5) * a

    cx = float(np.sum(w * x) * pixel_area / W)
    cy = float(np.sum(w * y) * pixel_area / W)

    xr = x - cx
    yr = y - cy

    # Second moments about centroid (point terms)
    Ixx_point = float(np.sum(w * (yr ** 2)) * pixel_area)
    Iyy_point = float(np.sum(w * (xr ** 2)) * pixel_area)

    # Self-terms for square pixels
    self_term = (a ** 4) / 12.0
    sumw = float(np.sum(w))
    Ixx = Ixx_point + sumw * self_term
    Iyy = Iyy_point + sumw * self_term

    J = float(Ixx + Iyy)

    return {"W": W, "cx": cx, "cy": cy, "J": J}


def check_against_bonej(
    case_name: str,
    mode: str,
    result: Dict[str, float],
    centroid_shift: Tuple[float, float],
    atol_centroid: float,
    atol_j: float,
) -> bool:
    """Print PASS/FAIL checks against BoneJ references for centroid + J only."""
    key = (case_name, mode)
    if key not in BONEJ_REFERENCE:
        print("  (No BoneJ reference available)")
        return False

    ref = BONEJ_REFERENCE[key]

    # Apply convention shift for centroid comparison
    sx, sy = centroid_shift
    cx_adj = result["cx"] + sx
    cy_adj = result["cy"] + sy

    print("  Validation vs BoneJ (centroid + J):")
    all_pass = True

    # Centroid
    for k, py_val_raw, py_val_adj, bj_val in [
        ("cx", result["cx"], cx_adj, ref["cx"]),
        ("cy", result["cy"], cy_adj, ref["cy"]),
    ]:
        diff = abs(py_val_adj - bj_val)
        passed = diff <= atol_centroid
        all_pass &= passed
        status = "PASS" if passed else "FAIL"
        print(
            f"    {k}: Python(raw)={py_val_raw:.6f}, Python(adj)={py_val_adj:.6f}, "
            f"BoneJ={bj_val:.6f}, Δ(adj)={diff:.3e} → {status}"
        )

    # J
    diff_j = abs(float(result["J"]) - float(ref["J"]))
    passed_j = diff_j <= atol_j
    all_pass &= passed_j
    status_j = "PASS" if passed_j else "FAIL"
    print(f"    J : Python={result['J']:.6f}, BoneJ={ref['J']:.6f}, Δ={diff_j:.3e} → {status_j}")

    print("  ✅ ALL CHECKS PASSED" if all_pass else "  ❌ VALIDATION FAILED")
    return all_pass


def main() -> None:
    print(f"ROOT_DIR: {ROOT_DIR}")
    if not ROOT_DIR.exists():
        raise FileNotFoundError(f"Folder not found: {ROOT_DIR}")

    sx, sy = CENTROID_SHIFT_TO_MATCH_BONEJ
    print(f"Using density calibration: m={DENSITY_M}, c={DENSITY_C} (BoneJ)")
    print(f"Centroid comparison uses shift: ({sx}, {sy}) to align conventions")
    print(f"Tolerances: centroid={ABS_TOL_CENTROID}, J={ABS_TOL_J}\n")

    overall_pass = True

    for base in CASE_BASENAMES:
        path = find_case_file(ROOT_DIR, base)
        img = imread(str(path))

        print("=" * 90)
        print(f"CASE: {base}")
        print(f"File: {path.name}")
        print(f"Shape: {img.shape}, dtype: {img.dtype}")

        for mode in ("geometric", "density"):
            res = compute_centroid_and_J(
                img=img,
                mode=mode,
                foreground_threshold=FOREGROUND_THRESHOLD,
                pixel_size=PIXEL_SIZE,
                density_m=DENSITY_M,
                density_c=DENSITY_C,
            )

            print(f"\n[{mode.upper()}]")
            print(f"  W   = {res['W']:.6g}")
            print(f"  c(raw) = ({res['cx']:.6g}, {res['cy']:.6g})")
            print(f"  c(adj) = ({res['cx']+sx:.6g}, {res['cy']+sy:.6g})   (adj = raw + shift)")
            print(f"  J   = {res['J']:.6g}")

            passed = check_against_bonej(
                case_name=base,
                mode=mode,
                result=res,
                centroid_shift=CENTROID_SHIFT_TO_MATCH_BONEJ,
                atol_centroid=ABS_TOL_CENTROID,
                atol_j=ABS_TOL_J,
            )
            overall_pass &= passed

        print()

    print("=" * 90)
    if overall_pass:
        print("✅ OVERALL: PASS (centroids + J match BoneJ within tolerance after convention alignment)")
        print("Note: centroid alignment uses a constant (-0.5, -0.5) shift due to coordinate convention.")
    else:
        print("❌ OVERALL: FAIL (one or more centroid/J checks exceeded tolerance)")
        print("Note: if fails occur only in centroid, re-check centroid shift assumption.")
    print("Done.")


if __name__ == "__main__":
    main()
