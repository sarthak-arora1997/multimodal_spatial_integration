"""
Load landmark points, compute affine transformation, and evaluate registration.

Reads QuPath-exported landmark TSV files for both Xenium and DESI, matches
point pairs, converts pixel coordinates to physical micron coordinates,
computes a 2D affine transformation using SimpleITK, and evaluates
per-landmark registration errors.
"""

import warnings
warnings.filterwarnings("ignore")

from pathlib import Path

import numpy as np
import pandas as pd
import SimpleITK as sitk
import spatialdata as sd
import tifffile
from spatialdata.models import get_channel_names

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"

INPUT_ZARR = DATA_DIR / "01_create_spatial_data_object" / "sdata_object.zarr"
INPUT_DESI_FLIPPED = DATA_DIR / "02_channel_extraction" / "desi_flipped.ome.tif"

XENIUM_LANDMARKS_TSV = DATA_DIR / "04_landmark_selection" / "morphology_focus_0002.ome.tif-points.tsv"
DESI_LANDMARKS_TSV = DATA_DIR / "04_landmark_selection" / "3D_DESI_F5_Pos mode_40um_Flipped.ome.tif - Image0-points.tsv"

OUTPUT_LANDMARKS_CSV = DATA_DIR / "05_landmark_registration" / "landmarks.csv"

# Physical pixel sizes (microns per pixel)
XENIUM_PIXEL_SIZE = 0.2125
DESI_PIXEL_SIZE = 40.0

# DESI channel index for vessel reference (Heme B)
DESI_HEME_B_INDEX = 17


def main():
    # 1. Load landmark TSV files
    print("Loading landmark files...")
    xenium_df = pd.read_csv(XENIUM_LANDMARKS_TSV, sep="\t")
    desi_df = pd.read_csv(DESI_LANDMARKS_TSV, sep="\t")
    print(f"  Xenium landmarks: {len(xenium_df)}")
    print(f"  DESI landmarks:   {len(desi_df)}")

    # 2. Extract point numbers and compute physical coordinates
    xenium_df["point_num"] = xenium_df["class"].str.extract(r"_(\d+)$").astype(int)
    desi_df["point_num"] = desi_df["class"].str.extract(r"_(\d+)$").astype(int)

    xenium_df["x_microns"] = xenium_df["x"] * XENIUM_PIXEL_SIZE
    xenium_df["y_microns"] = xenium_df["y"] * XENIUM_PIXEL_SIZE
    desi_df["x_microns"] = desi_df["x"] * DESI_PIXEL_SIZE
    desi_df["y_microns"] = desi_df["y"] * DESI_PIXEL_SIZE

    # 3. Merge on point number
    landmarks_df = (
        xenium_df[["point_num", "x", "y", "x_microns", "y_microns"]]
        .merge(
            desi_df[["point_num", "x", "y", "x_microns", "y_microns"]],
            on="point_num",
            suffixes=("_xenium", "_desi"),
        )
        .sort_values("point_num")
        .reset_index(drop=True)
    )
    print(f"\n  Matched landmark pairs: {len(landmarks_df)}")
    print(landmarks_df.to_string(index=False))

    # 4. Save merged landmarks
    OUTPUT_LANDMARKS_CSV.parent.mkdir(parents=True, exist_ok=True)
    landmarks_df.to_csv(OUTPUT_LANDMARKS_CSV, index=False)
    print(f"\n  Saved landmarks to:\n    {OUTPUT_LANDMARKS_CSV}")

    # 5. Load reference channels for SimpleITK
    print("\nLoading reference channels...")
    sdata = sd.read_zarr(INPUT_ZARR)
    morph_array = sdata.images["morphology_focus"]["scale0"].to_dataset()["image"]
    alphasma_vim = morph_array.sel(c="AlphaSMA/Vimentin").values

    desi_flipped = tifffile.imread(INPUT_DESI_FLIPPED)
    heme_b = desi_flipped[DESI_HEME_B_INDEX]

    # 6. Convert to SimpleITK images
    xenium_sitk = sitk.GetImageFromArray(alphasma_vim)
    xenium_sitk.SetSpacing([XENIUM_PIXEL_SIZE, XENIUM_PIXEL_SIZE])
    xenium_sitk.SetOrigin([0.0, 0.0])

    desi_sitk = sitk.GetImageFromArray(heme_b)
    desi_sitk.SetSpacing([DESI_PIXEL_SIZE, DESI_PIXEL_SIZE])
    desi_sitk.SetOrigin([0.0, 0.0])

    print(f"  Xenium SimpleITK: size={xenium_sitk.GetSize()}, spacing={xenium_sitk.GetSpacing()}")
    print(f"  DESI SimpleITK:   size={desi_sitk.GetSize()}, spacing={desi_sitk.GetSpacing()}")

    # 7. Prepare landmark coordinates (flat list for SimpleITK)
    fixed_landmarks = [[float(r.x_microns_xenium), float(r.y_microns_xenium)] for _, r in landmarks_df.iterrows()]
    moving_landmarks = [[float(r.x_microns_desi), float(r.y_microns_desi)] for _, r in landmarks_df.iterrows()]

    fixed_flat = [c for pt in fixed_landmarks for c in pt]
    moving_flat = [c for pt in moving_landmarks for c in pt]

    # 8. Compute affine transformation
    print("\nComputing affine transformation...")
    transform = sitk.AffineTransform(2)
    transform = sitk.LandmarkBasedTransformInitializer(transform, fixed_flat, moving_flat)

    matrix = transform.GetMatrix()
    translation = transform.GetTranslation()
    rotation_deg = np.degrees(np.arctan2(matrix[2], matrix[0]))
    scale_x = np.sqrt(matrix[0] ** 2 + matrix[2] ** 2)
    scale_y = np.sqrt(matrix[1] ** 2 + matrix[3] ** 2)

    print(f"\n  Transformation Matrix:")
    print(f"    [{matrix[0]:8.5f}  {matrix[1]:8.5f}]")
    print(f"    [{matrix[2]:8.5f}  {matrix[3]:8.5f}]")
    print(f"  Translation: [{translation[0]:.2f}, {translation[1]:.2f}] um")
    print(f"  Rotation:    {rotation_deg:.2f} degrees")
    print(f"  Scale X:     {scale_x:.4f}")
    print(f"  Scale Y:     {scale_y:.4f}")

    # 9. Evaluate per-landmark errors
    transformed = np.array([transform.TransformPoint(pt) for pt in moving_landmarks])
    fixed_arr = np.array(fixed_landmarks)
    errors = np.linalg.norm(transformed - fixed_arr, axis=1)

    print(f"\n  Landmark Registration Errors (um):")
    for i, row in landmarks_df.iterrows():
        print(f"    Point {row.point_num:2f}: {errors[i]:.2f} um")

    print(f"\n  Mean error:   {errors.mean():.2f} um")
    print(f"  Median error: {np.median(errors):.2f} um")
    print(f"  Max error:    {errors.max():.2f} um")
    print(f"  Min error:    {errors.min():.2f} um")

    if errors.mean() < 50:
        print("\n  EXCELLENT registration (mean < 50 um)")
    elif errors.mean() < 100:
        print("\n  GOOD registration (mean < 100 um)")
    elif errors.mean() < 200:
        print("\n  FAIR registration (mean < 200 um)")
    else:
        print("\n  WARNING: Poor registration (mean > 200 um) -- consider re-placing landmarks")


if __name__ == "__main__":
    main()
