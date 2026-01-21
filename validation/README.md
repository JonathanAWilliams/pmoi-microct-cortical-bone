# Validation: Python vs BoneJ Slice Geometry

This folder contains minimal, reproducible validation scripts used in the manuscript to verify agreement between:

- a reference Python implementation of slice-wise geometric and density-weighted moments, and
- the Slice Geometry plugin implementation in BoneJ (ImageJ/Fiji).

The validation focuses on quantities explicitly reported and discussed in the paper.

---

## What is validated

For a set of synthetic 2D test images, we validate agreement between Python and BoneJ for:

- geometric centroid (cx, cy)  
- density-weighted centroid (cx, cy), using calibration m = 1, c = 0  
- polar moment invariant  

  $$
J = I_{xx} + I_{yy}
$$


These quantities are sufficient to validate the correctness of the density-weighted formulation and its implementation in Slice Geometry.

---

## What is not validated

- principal moments (Imin, Imax)  
- principal axes orientation  

These are intentionally excluded, as they are sensitive to convention choices and are not required for the manuscript’s conclusions.

---

## Centroid coordinate conventions

BoneJ and the Python reference implementation use different pixel coordinate conventions:

- Python uses pixel-centre coordinates  
  (x + 0.5, y + 0.5)

- BoneJ reports centroids using pixel-index coordinates

As a result, centroids differ by a constant offset of (−0.5, −0.5) pixels.

The validation script:
- reports raw Python centroids, and
- compares BoneJ values against Python centroids shifted by (−0.5, −0.5).

After this alignment, centroids agree within tolerance.

The polar moment invariant J is unaffected by this convention difference.

---

## Files in this folder

### validate_slicegeometry_vs_python.py

Main validation script.

- loads synthetic test images
- computes centroids and J using the Python reference implementation
- compares results against hard-coded BoneJ reference values
- reports PASS / FAIL per case and per metric
- uses strict absolute tolerances consistent with BoneJ table output

This script is intended as a methodological validation, not a general analysis tool.

---

### generate_slicegeometry_validation_cases.py

Utility script to generate the synthetic test images used for validation.

It reproduces the exact ImageJ macro test cases used during Slice Geometry development.

Running this script creates:


