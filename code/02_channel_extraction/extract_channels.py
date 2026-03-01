"""
Extract and flip DESI channels, and save a 4-channel DESI composite.

Reads the combined SpatialData object from step 01, flips the DESI image
along the x-axis to match Xenium orientation, saves the flipped DESI image
(all 88 channels), and saves a 4-channel composite of selected DESI channels.
"""

import warnings
warnings.filterwarnings("ignore")

from pathlib import Path

import numpy as np
import spatialdata as sd
import tifffile

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"

INPUT_ZARR = DATA_DIR / "01_create_spatial_data_object" / "sdata_object.zarr"
OUTPUT_DESI_FLIPPED = DATA_DIR / "02_channel_extraction" / "desi_flipped.ome.tif"
OUTPUT_DESI_COMPOSITE = DATA_DIR / "02_channel_extraction" / "desi_composite.ome.tif"


def main():
    # 1. Load SpatialData object
    print(f"Loading SpatialData from:\n  {INPUT_ZARR}")
    sdata = sd.read_zarr(INPUT_ZARR)

    # 2. Extract and flip DESI image
    print("\nExtracting and flipping DESI image...")
    desi = sdata.images["desi"]
    desi_full = desi.values  # (c, y, x)
    desi_flipped = np.flip(desi_full, axis=2)
    print(f"  DESI flipped shape: {desi_flipped.shape} (all {desi_flipped.shape[0]} channels)")

    # 3. Save flipped DESI image (all channels)
    OUTPUT_DESI_FLIPPED.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(OUTPUT_DESI_FLIPPED, desi_flipped)
    print(f"  Saved flipped DESI to:\n    {OUTPUT_DESI_FLIPPED}")

    # 4. Extract selected DESI channels and save composite
    desi_channel_names = list(desi.coords["c"].values)
    selected_indices = [3, 17, 23, 25]
    selected_channels = {desi_channel_names[i]: i for i in selected_indices}
    print(f"\n  Selected DESI channels: {list(selected_channels.keys())}")

    desi_composite = np.stack([desi_flipped[i] for i in selected_indices], axis=0)  # (4, Y, X)
    tifffile.imwrite(OUTPUT_DESI_COMPOSITE, desi_composite)
    print(f"  DESI composite: {desi_composite.shape}, dtype={desi_composite.dtype}")
    print(f"  Saved to:\n    {OUTPUT_DESI_COMPOSITE}")

    print("\nDone.")


if __name__ == "__main__":
    main()
