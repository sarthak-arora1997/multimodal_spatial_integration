"""
Extract morphology channels from Xenium and DESI, flip DESI, and create composites.

Reads the combined SpatialData object from step 01, extracts biologically
relevant channels from both modalities, flips the DESI image along the x-axis
to match Xenium orientation, and saves the flipped DESI image (all 88 channels).
"""

import warnings
warnings.filterwarnings("ignore")

from pathlib import Path

import numpy as np
import spatialdata as sd
import tifffile
from spatialdata.models import get_channel_names

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"

INPUT_ZARR = DATA_DIR / "01_create_spatial_data_object" / "sdata_object.zarr"
OUTPUT_DESI_FLIPPED = DATA_DIR / "02_channel_extraction" / "desi_flipped.ome.tif"


def normalize(img):
    """Normalize image intensities to [0, 1]."""
    img = img.astype(float)
    return (img - img.min()) / (img.max() - img.min() + 1e-10)


def main():
    # 1. Load SpatialData object
    print(f"Loading SpatialData from:\n  {INPUT_ZARR}")
    sdata = sd.read_zarr(INPUT_ZARR)

    # 2. Extract Xenium morphology channels at full resolution (scale0)
    print("\nExtracting Xenium morphology channels...")
    morph = sdata.images["morphology_focus"]
    channel_names = get_channel_names(morph)
    print(f"  Channels: {channel_names}")

    morph_full = morph["scale0"].to_dataset()
    morph_array = morph_full["image"]

    alphasma_vim = morph_array.sel(c="AlphaSMA/Vimentin").values
    atp1a1 = morph_array.sel(c="ATP1A1/CD45/E-Cadherin").values
    dapi = morph_array.sel(c="DAPI").values
    print(f"  AlphaSMA/Vimentin: {alphasma_vim.shape}")
    print(f"  ATP1A1/CD45/E-Cadherin: {atp1a1.shape}")
    print(f"  DAPI: {dapi.shape}")

    # 3. Extract and flip DESI image
    print("\nExtracting and flipping DESI image...")
    desi = sdata.images["desi"]
    desi_full = desi.values  # (c, y, x)
    desi_flipped = np.flip(desi_full, axis=2)
    print(f"  DESI flipped shape: {desi_flipped.shape} (all {desi_flipped.shape[0]} channels)")

    # 4. Save flipped DESI image (all channels)
    OUTPUT_DESI_FLIPPED.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(OUTPUT_DESI_FLIPPED, desi_flipped)
    print(f"  Saved flipped DESI to:\n    {OUTPUT_DESI_FLIPPED}")

    # 5. Extract specific DESI channels for composites
    desi_channel_names = list(desi.coords["c"].values)
    selected_channels = {
        "mz_204.13": 3,
        "heme_B_616.25": 17,
        "PE_PS_731.65": 23,
        "PE_PS_734.60": 25,
    }
    print(f"\n  Selected DESI channels: {list(selected_channels.keys())}")

    mz_204 = desi_flipped[selected_channels["mz_204.13"]]
    heme_b = desi_flipped[selected_channels["heme_B_616.25"]]
    pe_ps_731 = desi_flipped[selected_channels["PE_PS_731.65"]]
    pe_ps_734 = desi_flipped[selected_channels["PE_PS_734.60"]]

    # 6. Create RGB composite images
    print("\nCreating RGB composite images...")
    xenium_composite = np.stack([
        normalize(alphasma_vim),
        normalize(atp1a1),
        np.zeros_like(atp1a1, dtype=float),
    ], axis=-1)

    desi_composite = np.stack([
        normalize(heme_b),
        normalize(pe_ps_731),
        np.zeros_like(heme_b, dtype=float),
    ], axis=-1)

    print("  Red channel:   Vessels (AlphaSMA <-> Heme B)")
    print("  Green channel: Membranes (ATP1A1 <-> Lipid 731.65)")
    print(f"  Xenium composite: {xenium_composite.shape}")
    print(f"  DESI composite:   {desi_composite.shape}")
    print("\nDone.")


if __name__ == "__main__":
    main()
