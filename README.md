# Density-weighted polar moment analysis (pMOI vs pMOA)

This repository contains the reference Python implementation used in the manuscript:
A Tissue Mineral Density-Weighted Polar Moment for microCT Cortical Bone Analysis

The code implements slice-wise geometric (pMOA) and density-weighted (pMOI) polar moment
calculations from microCT cortical cross-sections, together with synthetic validation
cases demonstrating when and why these measures diverge.

---

## Repository structure

```text
.
├─ compute_slice_metrics_batch.py
├─ plot_slice_profiles_and_summary.py
├─ synthetic_shapes.py
├─ visualise_pmoa_i_contributions_two_samples.py
├─ validation/
│  ├─ README.md
│  ├─ validate_slicegeometry_vs_python.py
│  ├─ generate_slicegeometry_validation_cases.py
└─ README.md
```

Core computation:
- `compute_slice_metrics_batch.py`
- `slice_metrics_full_<Sample>.csv`
  Batch slice-wise computation of Ct.Ar, TMD, geometric pMOA, and density-weighted pMOI from 2D cortical microCT slices.

Synthetic validation data:
- `synthetic_shapes.py`  
  Generates the synthetic test cases used in the manuscript.

Plotting / figure reproduction:
- `plot_slice_profiles_and_summary.py`  
  Generates slice-wise (normalised bone length) profiles and summary plots/statistics from the exported results.

Contribution visualisations:
- `visualise_pmoa_i_contributions_two_samples.py`  
  Produces per-pixel pMOA and pMOI contribution heatmaps for two samples with shared colour scales (transparent background for 3D Slicer overlays).

Validation against BoneJ Slice Geometry:
- `validation/`  
  Scripts + minimal test images used to compare Python outputs to BoneJ Slice Geometry.


## Core computation: `compute_slice_metrics_batch.py`

### What it does

`compute_slice_metrics_batch.py` computes slice-wise geometric (pMOA) and density-weighted (pMOI)
polar moment metrics from segmented cortical microCT cross-sections.

For each slice it reports:
- Cortical area (Ct.Ar)
- Geometric centroid and geometric second moments (Ixx, Iyy, Ixy) and polar moment (J = pMOA)
- Density-weighted centroid (centre-of-mass analogue)
- Density-weighted second moments and polar moment (pMOI)

Outputs are saved as:
- One CSV per sample (slice-wise metrics)
- One combined Excel workbook containing all samples

---

### Script location

```text
scripts/compute_slice_metrics_batch.py
```
---

### Expected input folder structure

The script expects a root directory containing one or more sample folders, each
containing an image stack of segmented cortical cross-sections.

```text
PI3K/
├── Sample13/
│   └── VOI/
│       └── cleaned_CORT_erosion/
│           ├── slice_0001.bmp
│           ├── slice_0002.bmp
│           └── ...
├── Sample15/
└── Sample20/
```

Notes:

- --root should point to the directory containing the sample folders (e.g. PI3K/).
- Image slices must be 2D grayscale images (.bmp, .png, .tif, etc.).
- If a folder named cleaned_CORT_erosion is present, it is selected automatically.
- Otherwise, the script selects the deepest subfolder containing image files.
  
---

### Running the analysis
### 1. Dry run (recommended first)
A dry run verifies that:
- sample folders are detected correctly,
- the correct image subfolder is selected.
No computation is performed.

```bash
python scripts/compute_slice_metrics_batch.py ^
  --root "E:\summer work\PI3K Study VOIs\PI3K" ^
  --ext .bmp ^
  --pixel-width-microns 9.041832 ^
  --thr-low 70 ^
  --thr-high 255 ^
  --cal-gray-low 43.0463712 ^
  --cal-gray-high 92.4032819 ^
  --bmd-low 0.25 ^
  --bmd-high 0.75 ^
  --dry-run
```
Expected output (example)

```bash
Found 8 sample folder(s) under: E:\summer work\PI3K Study VOIs\PI3K
[Sample13] images: VOI\cleaned_CORT_erosion
...
Dry run complete. No computation performed.
```

### 2. Process a single sample

```bash
python scripts/compute_slice_metrics_batch.py ^
  --root "E:\summer work\PI3K Study VOIs\PI3K" ^
  --sample Sample13 ^
  --ext .bmp ^
  --pixel-width-microns 9.041832 ^
  --thr-low 70 ^
  --thr-high 255 ^
  --cal-gray-low 43.0463712 ^
  --cal-gray-high 92.4032819 ^
  --bmd-low 0.25 ^
  --bmd-high 0.75
```

### 3. Batch processing (all samples)

```bash
python scripts/compute_slice_metrics_batch.py ^
  --root "E:\summer work\PI3K Study VOIs\PI3K" ^
  --ext .bmp ^
  --pixel-width-microns 9.041832 ^
  --thr-low 70 ^
  --thr-high 255 ^
  --cal-gray-low 43.0463712 ^
  --cal-gray-high 92.4032819 ^
  --bmd-low 0.25 ^
  --bmd-high 0.75 ^
  --excel-name TMD_pMOA_pMOI_summary_PI3K.xlsx
```
---

### Outputs
All outputs are written to the directory specified by `--root`.

Per-sample output:

- `slice_metrics_full_<Sample>.csv` (slice-wise geometric and density-weighted metrics)

Combined output:
- One Excel workbook containing all samples (default name: `TMD_pMOA_pMOI_summary.xlsx`)
---

### Units and conventions
- Pixel width is specified in microns and internally converted to mm
- Geometric second moments and pMOA have units of mm⁴
- Density-weighted second moments and pMOI have units of (g/cm³)·mm⁴
- Density weighting uses a two-point grayscale → BMD/TMD calibration
- Slice-wise centroids are computed independently for each slice
---

### Relationship to BoneJ

This Python implementation is the reference implementation used in the manuscript.
An equivalent implementation has been integrated into the Slice Geometry plugin in
BoneJ and cross-validated using identical synthetic test cases.

---

### Citation

If you use this code, please cite:
Author(s). Title of paper. Journal, Year.




