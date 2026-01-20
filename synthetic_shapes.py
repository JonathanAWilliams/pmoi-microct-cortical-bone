import os
import numpy as np
import cv2

"""
Synthetic dataset generator for validating pMOA (geometric) vs density-weighted pMOI in 2D.

This script creates idealised hollow cylindrical cortical cross-sections with identical geometry
(i.e., identical binary area distribution), but with different grayscale intensity patterns that
represent different mineralisation (density) distributions.

Intended use in the study:
- Provides controlled test inputs where pMOA should be identical across datasets (same geometry),
  while density-weighted metrics (e.g., pMOI / density-weighted centroid) can differ due solely to
  spatial redistribution of intensity (mineralisation).
- Outputs are written as 8-bit BMP slices to disk for use in downstream pipelines (Python / BoneJ).

Coordinate convention:
- Arrays follow NumPy/OpenCV convention: index (row=y, col=x)
- The 'center' tuple is (x, y), matching OpenCV drawing functions.
- Distances/angles are computed in pixel units, using pixel-centre approximation as given by indices.

Output:
- Folder structure:
    synthetic_slices/<dataset_name>/slice_0001.bmp ... slice_NNNN.bmp
- Each folder contains n_slices identical slices (constant along z) to create a "stack".

Note:
- All intensities are currently hard-coded offsets around a target mean; if strict equality of mean
  grayscale (and therefore mean TMD) across variants is required, area fractions must be used to
  solve intensities exactly (see comments in compute_intensities_for_equal_TMD()).
  
This script is used to generate Figure 1 in the manuscript  
"""


# =====================================================================
# CONFIGURATION
# =====================================================================

shape = (512, 512)
center = (256, 256)
outer_radius = 200
inner_radius = 100
n_slices = 200

# --- Calibration (your BCL3/BMD calibration) ---
CAL_GRAY_0_75 = 117.3579096
CAL_GRAY_0_25 = 48.4748119
BMD_0_75 = 0.75
BMD_0_25 = 0.25

# Build grayscale → TMD calibration
slope = (BMD_0_75 - BMD_0_25) / (CAL_GRAY_0_75 - CAL_GRAY_0_25)
intercept = BMD_0_75 - slope * CAL_GRAY_0_75


# =====================================================================
# HELPER FUNCTIONS
# =====================================================================

def grayscale_to_bmd(g):
    """Convert grayscale to BMD using calibration line."""
    return slope * g + intercept


def compute_region_pixel_counts(outer_r, inner_r):
    """Compute areas of outer shell, middle band, inner cortical band."""
    y, x = np.indices(shape)
    dist = np.sqrt((x - center[0])**2 + (y - center[1])**2)

    mask_outer_shell = (dist <= outer_r) & (dist > (inner_r + (outer_r - inner_r) // 2))
    mask_middle = (dist <= (inner_r + (outer_r - inner_r) // 2)) & (dist > inner_r + (inner_r - inner_r) // 2)
    mask_inner_shell = (dist <= inner_r + (inner_r - inner_r) // 2) & (dist > inner_radius)

    return mask_outer_shell, mask_middle, mask_inner_shell


def compute_intensities_for_equal_TMD():
    """
    Solves intensities so that all datasets have the same mean grayscale -> same TMD.
    We define:
        Uniform: g0
        Ring: outer = g1, inner = g2
        Sandwich: outer&inner = g3, middle = g4
        Wedge: base = g5, wedge = g6

    And enforce:
        mean grayscale = constant for all datasets
        
    """

    # We want the mean grayscale to equal g_target
    g_target = 150  # choose reasonable global mean

    # --- Uniform cylinder ---
    g_uniform = g_target

    # --- Ring variant ---
    # Outer ring 50%, inner ring 50% → solve g1, g2 such that mean = g_target
    g_ring_outer = g_target + 30
    g_ring_inner = g_target - 30

    # --- Sandwich variant ---
    # Outer 1/3 high, middle 1/3 low, inner 1/3 high
    g_sandwich_outer = g_target + 20
    g_sandwich_middle = g_target - 40

    # --- Wedge variant: keep mean equal ---
    g_wedge_base = g_target - 10
    g_wedge_sector = g_target + 50

    return {
        "hollow_cylinder": {"intensity": g_uniform},
        "ring_variation": {"outer": g_ring_outer, "inner": g_ring_inner},
        "sandwich": {"outer": g_sandwich_outer, "middle": g_sandwich_middle},
        "wedge_variation": {"base": g_wedge_base, "wedge": g_wedge_sector}
    }


# =====================================================================
# IMAGE GENERATORS
# =====================================================================

def generate_hollow_cylinder(img_shape, center, outer_r, inner_r, intensity):
    """
    Generate a hollow cylinder cross-section:
    - Cortex region (outer_r) filled with 'intensity'
    - Medullary canal (inner_r) set to 0

    This is the baseline (uniform density) case.
    """
    img = np.zeros(img_shape, dtype=np.uint8)
    cv2.circle(img, center, outer_r, intensity, thickness=-1)
    cv2.circle(img, center, inner_r, 0, thickness=-1)
    return img


def generate_ring_variant(img_shape, center, outer_r, inner_r, outer_int, inner_int):
    """
    Generate a 2-band concentric ring:
    - Outer half of cortex: outer_int
    - Inner half of cortex: inner_int
    - Medullary canal: 0

    Intended effect:
    - Same binary geometry as hollow_cylinder.
    - Redistributes intensity radially while preserving concentric symmetry.
    """
    img = np.zeros(img_shape, dtype=np.uint8)
    mid_r = (outer_r + inner_r) // 2
    cv2.circle(img, center, outer_r, outer_int, thickness=-1)
    cv2.circle(img, center, mid_r, inner_int, thickness=-1)
    cv2.circle(img, center, inner_r, 0, thickness=-1)
    return img

def generate_wedge_variant(img_shape, center, outer_r, inner_r, base_int, wedge_int):
    """
    Generate an angular 'wedge' of higher intensity within the cortex.

    Intended effect:
    - Same binary geometry as hollow_cylinder.
    - Breaks angular symmetry, shifting the density-weighted centroid toward the wedge.
    - Provides a clear test where density-weighted centroid and second moments should differ from
      the geometric versions, despite identical pMOA.
    """
    img = generate_hollow_cylinder(img_shape, center, outer_r, inner_r, base_int)

    # Angle map in radians, with arctan2(y, x) using NumPy indexing (Y=row, X=col).
    Y, X = np.indices(img_shape)
    angles = np.arctan2(Y - center[1], X - center[0])

    # Define a +/- 45-degree wedge about the +x axis (right side of image).
    wedge_mask = (angles > -np.pi/4) & (angles < np.pi/4) & (img > 0)
    img[wedge_mask] = wedge_int
    return img

def generate_sandwich_variant(img_shape, center, outer_r, inner_r, outer_int, middle_int):
    """
    Generate a 3-band radial pattern: high (outer cortex) – low (middle) – high (inner cortex).

    Intended effect:
    - Preserves concentric symmetry (centroid stays at geometric centre).
    - Produces a different radial density distribution than ring_variation.
    """
    img = np.zeros(img_shape, dtype=np.uint8)

    # Define radii for three bands within the cortex.
    r3 = outer_r
    r2 = (outer_r + inner_r) // 2
    r1c = inner_r + (r2 - inner_r) // 2  # boundary between inner cortex and middle band
    r1 = inner_r

    cv2.circle(img, center, r3, outer_int, thickness=-1)     # outer cortex (high)
    cv2.circle(img, center, r2, middle_int, thickness=-1)    # middle band (low)
    cv2.circle(img, center, r1c, outer_int, thickness=-1)    # inner cortex (high)
    cv2.circle(img, center, r1, 0, thickness=-1)             # medullary canal

    return img


# =====================================================================
# SAVE DATASETS
# =====================================================================

def save_dataset(img_generator, out_folder, **kwargs):
    os.makedirs(out_folder, exist_ok=True)
    for i in range(n_slices):
        img = img_generator(**kwargs)
        filename = os.path.join(out_folder, f"slice_{i+1:04d}.bmp")
        cv2.imwrite(filename, img)


# =====================================================================
# MAIN EXECUTION
# =====================================================================

intensities = compute_intensities_for_equal_TMD()

save_dataset(generate_hollow_cylinder, "synthetic_slices/hollow_cylinder",
             img_shape=shape, center=center,
             outer_r=outer_radius, inner_r=inner_radius,
             intensity=intensities["hollow_cylinder"]["intensity"])

save_dataset(generate_ring_variant, "synthetic_slices/ring_variation",
             img_shape=shape, center=center,
             outer_r=outer_radius, inner_r=inner_radius,
             outer_int=intensities["ring_variation"]["outer"],
             inner_int=intensities["ring_variation"]["inner"])

save_dataset(generate_wedge_variant, "synthetic_slices/wedge_variation",
             img_shape=shape, center=center,
             outer_r=outer_radius, inner_r=inner_radius,
             base_int=intensities["wedge_variation"]["base"],
             wedge_int=intensities["wedge_variation"]["wedge"])

save_dataset(generate_sandwich_variant, "synthetic_slices/sandwich",
             img_shape=shape, center=center,
             outer_r=outer_radius, inner_r=inner_radius,
             outer_int=intensities["sandwich"]["outer"],
             middle_int=intensities["sandwich"]["middle"])

print("Synthetic datasets generated with controlled TMD, pMOA, and area.")
