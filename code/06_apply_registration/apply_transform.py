"""
Apply the affine registration to resample DESI into Xenium coordinate space.

Loads the matched landmarks from step 05, recomputes the affine transformation,
resamples the 4 selected DESI channels into the Xenium coordinate system using
SimpleITK, and saves the registered DESI image as an OME-TIFF.
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
INPUT_DESI_COMPOSITE = DATA_DIR / "02_channel_extraction" / "desi_composite.ome.tif"
INPUT_LANDMARKS_CSV = DATA_DIR / "05_landmark_registration" / "landmarks.csv"

OUTPUT_REGISTERED = DATA_DIR / "06_apply_registration" / "registered_desi.ome.tif"

# Physical pixel sizes (microns per pixel)
XENIUM_PIXEL_SIZE = 0.2125
DESI_PIXEL_SIZE = 40.0


def main():
    # 1. Load landmarks
    print(f"Loading landmarks from:\n  {INPUT_LANDMARKS_CSV}")
    landmarks_df = pd.read_csv(INPUT_LANDMARKS_CSV)
    print(f"  {len(landmarks_df)} landmark pairs")

    # 2. Load Xenium morphology channels
    print("\nLoading Xenium morphology channels...")
    sdata = sd.read_zarr(INPUT_ZARR)
    morph_array = sdata.images["morphology_focus"]["scale0"].to_dataset()["image"]
    xenium_channels = [
        ch for ch in morph_array.coords["c"].values if ch != "dummy"
    ]
    xenium_arrays = []
    for ch_name in xenium_channels:
        arr = morph_array.sel(c=ch_name).values
        print(f"  {ch_name}: shape={arr.shape}, dtype={arr.dtype}")
        xenium_arrays.append(arr)

    reference = morph_array.sel(c="AlphaSMA/Vimentin").values
    assert reference.sum() != 0 and reference.mean() != 0, \
        "AlphaSMA/Vimentin reference channel is all zeros"

    # 3. Load DESI composite (4 channels)
    print("\nLoading DESI composite...")
    desi_composite = tifffile.imread(INPUT_DESI_COMPOSITE)  # (4, Y, X)
    print(f"  DESI composite: shape={desi_composite.shape}, dtype={desi_composite.dtype}")

    # 4. Create SimpleITK reference image (defines output coordinate space)
    xenium_sitk = sitk.GetImageFromArray(reference)
    xenium_sitk.SetSpacing([XENIUM_PIXEL_SIZE, XENIUM_PIXEL_SIZE])
    xenium_sitk.SetOrigin([0.0, 0.0])

    # 5. Compute affine transform from landmarks
    print("\nComputing affine transformation...")
    fixed_landmarks = [[float(r.x_microns_xenium), float(r.y_microns_xenium)] for _, r in landmarks_df.iterrows()]
    moving_landmarks = [[float(r.x_microns_desi), float(r.y_microns_desi)] for _, r in landmarks_df.iterrows()]

    fixed_flat = [c for pt in fixed_landmarks for c in pt]
    moving_flat = [c for pt in moving_landmarks for c in pt]

    transform = sitk.AffineTransform(2)
    transform = sitk.LandmarkBasedTransformInitializer(transform, fixed_flat, moving_flat)
    print(f"  Affine transform computed from {len(landmarks_df)} landmarks")

    # 6. Resample each DESI channel into Xenium coordinate space
    print("\nResampling DESI channels into Xenium coordinate space...")
    registered_channels = []
    for i in range(desi_composite.shape[0]):
        desi_ch_sitk = sitk.GetImageFromArray(desi_composite[i])
        desi_ch_sitk.SetSpacing([DESI_PIXEL_SIZE, DESI_PIXEL_SIZE])
        desi_ch_sitk.SetOrigin([0.0, 0.0])

        registered_sitk = sitk.Resample(
            desi_ch_sitk, xenium_sitk, transform,
            sitk.sitkLinear, 0.0, desi_ch_sitk.GetPixelID(),
        )
        registered_channels.append(sitk.GetArrayFromImage(registered_sitk))
        print(f"  Channel {i}: {registered_channels[-1].shape}")

    # 7. Stack registered DESI + all Xenium channels
    stacked = np.stack(registered_channels + xenium_arrays, axis=0)
    print(f"\n  Final registered image: shape={stacked.shape}, dtype={stacked.dtype}")

    OUTPUT_REGISTERED.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(OUTPUT_REGISTERED, stacked)
    print(f"Saved to:\n  {OUTPUT_REGISTERED}")
    print("Done.")


if __name__ == "__main__":
    main()
