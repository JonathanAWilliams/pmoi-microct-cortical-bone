"""
Generate tiny validation images used for Slice Geometry vs Python checks.

These replicate the ImageJ macro test cases used in the manuscript validation.

Outputs (PNG):
  validation/slicegeometry_validation_cases/
    CaseA_2x2_block_4x4.png
    CaseB_two_pixel_row_4x3_100_200.png
"""


from pathlib import Path
import numpy as np
import imageio.v2 as imageio

def save_png(img: np.ndarray, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    assert img.dtype == np.uint8
    imageio.imwrite(out_path.as_posix(), img)

def make_case_a_2x2_block() -> np.ndarray:
    """
    ImageJ macro equivalent:
      newImage("CaseA_2x2_block", "8-bit black", 4, 4, 1);
      for (y=1; y<=2; y++) for (x=1; x<=2; x++) setPixel(x, y, 255);

    Returns: 4x4 uint8 image.
    """
    img = np.zeros((4, 4), dtype=np.uint8)
    # ImageJ uses (x,y) with origin at top-left; numpy uses [y, x]
    img[1:3, 1:3] = 255  # y=1..2, x=1..2
    return img

def make_case_b_two_pixel_row(g_left: int = 100, g_right: int = 200) -> np.ndarray:
    """
    Case B used in our validation:
      3x3 image, two adjacent pixels in a row:
        (x,y)=(1,1)->g_left
        (x,y)=(2,1)->g_right
    Returns: 4x3 uint8 image.
    """
    img = np.zeros((4, 3), dtype=np.uint8)
    img[1, 1] = np.uint8(g_left)   # (x=1,y=1)
    img[1, 2] = np.uint8(g_right)  # (x=2,y=1)
    return img

if __name__ == "__main__":
    out_dir = Path("slicegeometry_validation_cases")

    case_a = make_case_a_2x2_block()
    save_png(case_a, out_dir / "CaseA_2x2_block_4x4.png")

    case_b = make_case_b_two_pixel_row(100, 200)
    save_png(case_b, out_dir / "CaseB_two_pixel_row_3x3_100_200.png")

    print("Wrote:")
    print(" -", (out_dir / "CaseA_2x2_block_4x4.png").resolve())
    print(" -", (out_dir / "CaseB_two_pixel_row_3x3_100_200.png").resolve())
