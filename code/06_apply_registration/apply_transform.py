"""
Apply the affine registration to resample DESI into Xenium coordinate space.

Loads the matched landmarks from step 05, recomputes the affine transformation,
resamples the DESI image into the Xenium coordinate system using SimpleITK,
and saves the registered DESI image as a TIFF file.
"""

import warnings
warnings.filterwarnings("ignore")

from pathlib import Path

import numpy as np
import pandas as pd
import SimpleITK as sitk
import spatialdata as sd
import tifffile

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"

INPUT_ZARR = DATA_DIR / "01_create_spatial_data_object" / "sdata_object.zarr"
INPUT_DESI_FLIPPED = DATA_DIR / "02_channel_extraction" / "desi_flipped.ome.tif"
INPUT_LANDMARKS_CSV = DATA_DIR / "05_landmark_registration" / "landmarks.csv"

OUTPUT_REGISTERED = DATA_DIR / "06_apply_registration" / "desi_registered.tif"

# Physical pixel sizes (microns per pixel)
XENIUM_PIXEL_SIZE = 0.2125
DESI_PIXEL_SIZE = 40.0

# DESI channel index for vessel reference (Heme B)
DESI_HEME_B_INDEX = 17


def main():
    # 1. Load landmarks
    print(f"Loading landmarks from:\n  {INPUT_LANDMARKS_CSV}")
    landmarks_df = pd.read_csv(INPUT_LANDMARKS_CSV)
    print(f"  {len(landmarks_df)} landmark pairs")

    # 2. Load reference channels
    print("\nLoading reference channels...")
    sdata = sd.read_zarr(INPUT_ZARR)
    morph_array = sdata.images["morphology_focus"]["scale0"].to_dataset()["image"]
    alphasma_vim = morph_array.sel(c="AlphaSMA/Vimentin").values

    desi_flipped = tifffile.imread(INPUT_DESI_FLIPPED)
    heme_b = desi_flipped[DESI_HEME_B_INDEX]

    # 3. Create SimpleITK images
    xenium_sitk = sitk.GetImageFromArray(alphasma_vim)
    xenium_sitk.SetSpacing([XENIUM_PIXEL_SIZE, XENIUM_PIXEL_SIZE])
    xenium_sitk.SetOrigin([0.0, 0.0])

    desi_sitk = sitk.GetImageFromArray(heme_b)
    desi_sitk.SetSpacing([DESI_PIXEL_SIZE, DESI_PIXEL_SIZE])
    desi_sitk.SetOrigin([0.0, 0.0])

    # 4. Compute affine transform from landmarks
    fixed_landmarks = [[float(r.x_microns_xenium), float(r.y_microns_xenium)] for _, r in landmarks_df.iterrows()]
    moving_landmarks = [[float(r.x_microns_desi), float(r.y_microns_desi)] for _, r in landmarks_df.iterrows()]

    fixed_flat = [c for pt in fixed_landmarks for c in pt]
    moving_flat = [c for pt in moving_landmarks for c in pt]

    transform = sitk.AffineTransform(2)
    transform = sitk.LandmarkBasedTransformInitializer(transform, fixed_flat, moving_flat)
    print(f"  Affine transform computed from {len(landmarks_df)} landmarks")

    # 5. Resample DESI into Xenium coordinate space
    print("\nResampling DESI into Xenium coordinate space...")
    desi_registered_sitk = sitk.Resample(
        desi_sitk,
        xenium_sitk,
        transform,
        sitk.sitkLinear,
        0.0,
        desi_sitk.GetPixelID(),
    )

    desi_registered = sitk.GetArrayFromImage(desi_registered_sitk)
    print(f"  Registered image shape: {desi_registered.shape}")
    print(f"  Pixel spacing: {desi_registered_sitk.GetSpacing()} um")

    # 6. Save registered image
    OUTPUT_REGISTERED.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(OUTPUT_REGISTERED, desi_registered)
    print(f"\n  Saved registered DESI to:\n    {OUTPUT_REGISTERED}")
    print("Done.")


if __name__ == "__main__":
    main()
