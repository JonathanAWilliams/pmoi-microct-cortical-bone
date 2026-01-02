**Repository Overview**

- **Purpose:** Single-script batch processor that computes cortical area (Ct.Ar), TMD, geometric inertia (pMOA) and TMD-weighted inertia (pMOI) from microCT image slices.
- **Main entry:** `pmoi_code.py` — standalone script expecting a `ROOT_DIR` with one subfolder per sample.

**How To Run (quick)**

- Ensure Python 3.8+ and these packages are installed: `numpy`, `pandas`, `scikit-image`, `openpyxl` (pandas Excel engine). Example:

```powershell
python -m pip install numpy pandas scikit-image openpyxl
python pmoi_code.py
```

- The script is configured via constants at the top of `pmoi_code.py` — edit `ROOT_DIR`, `IMAGE_EXT`, `PIXEL_WIDTH_MICRONS`, thresholds and calibration values there. There are no CLI args.

**Big-picture architecture & data flow**

- Single-process, file-based pipeline: `pmoi_code.py` iterates immediate subfolders of `ROOT_DIR` (each is a sample), locates image files (either directly in the sample folder or inside one subfolder), computes per-slice metrics, writes per-sample CSV and a combined Excel workbook in `ROOT_DIR`.
- Per-slice steps (see `compute_slice_metrics`): convert grayscale -> threshold to bone mask -> calibrate grayscale to BMD -> compute centroids, geometric inertia tensor, TMD-weighted inertia tensor -> return numeric metrics.
- Output naming patterns:
  - Per-sample CSV: `slice_metrics_TMD_full_inertia_<SampleName>.csv`
  - Combined Excel workbook: filename is the `EXCEL_NAME` constant (default `TMD_pMOA_pMOI_summary_SCI.xlsx`).

**Project-specific conventions / gotchas**

- Units: pixel width specified in microns (`PIXEL_WIDTH_MICRONS`), converted to mm internally. Area in mm^2, TMD in g/cm^3, inertias in mm^4 or g·mm^4-like units as coded.
- Image discovery: supports two layouts — images directly in `SampleX/` OR inside a single child folder (common BCL3-style). The code searches both.
- Natural sorting: filenames are sorted with `natural_sort_key` to maintain slice order (handles numeric suffixes like slice1, slice2, ...).
- Calibration: grayscale→BMD mapping uses `CAL_GRAY_0_75` and `CAL_GRAY_0_25`; change these values only when you have calibration phantom results.
- Empty-mask behavior: if a slice produces no bone pixels (mask empty) it is skipped; if many/ all slices are skipped, no CSV/excel sheet is produced for that sample.

**Files of interest**

- `pmoi_code.py` — main script and configuration (see constants at top).

**Debugging & iterative work**

- To test quickly, duplicate a single sample folder containing a few images and point `ROOT_DIR` to it.
- If outputs are empty, try lowering `LOWER_THRESHOLD` or verify `IMAGE_EXT` matches your images (`.bmp`, `.tif`, `.png`).
- Add `print()` calls inside `process_sample_folder` and `compute_slice_metrics` to inspect intermediate arrays if needed.

**Integration points / external deps**

- Reads filesystem only (no network). External packages: `numpy`, `pandas`, `scikit-image` (for `imread`), and an Excel writer such as `openpyxl` for `pandas`.

**When editing or extending**

- Prefer modifying constants at the top of `pmoi_code.py` for dataset-specific changes. If adding CLI options, keep defaults matching current constants for backward compatibility.
- If you add engines that change Excel writing behavior, update the top-of-file comments and `EXCEL_NAME` usage.

**Feedback / Unknowns**

- If you want CLI flags, multi-threading, or unit tests added, specify desired behaviors and I can implement them.
