"""
Visualize Xenium and DESI channels in napari.

Opens an interactive napari viewer with Xenium morphology channels and DESI
mass spectrometry channels overlaid at their respective physical scales.
Also displays RGB composite images for visual comparison of matched
biological features across modalities.
"""

import warnings
warnings.filterwarnings("ignore")

from pathlib import Path

import napari
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
INPUT_DESI_FLIPPED = DATA_DIR / "02_channel_extraction" / "desi_flipped.ome.tif"

# DESI channel indices for the 4 selected channels
DESI_CHANNEL_INDICES = {"mz_204.13": 3, "heme_B_616.25": 17, "PE_PS_731.65": 23, "PE_PS_734.60": 25}


def normalize(img):
    """Normalize image intensities to [0, 1]."""
    img = img.astype(float)
    return (img - img.min()) / (img.max() - img.min() + 1e-10)


def main():
    # 1. Load data
    print(f"Loading SpatialData from:\n  {INPUT_ZARR}")
    sdata = sd.read_zarr(INPUT_ZARR)

    print(f"Loading flipped DESI from:\n  {INPUT_DESI_FLIPPED}")
    desi_flipped = tifffile.imread(INPUT_DESI_FLIPPED)

    # 2. Extract Xenium channels
    morph = sdata.images["morphology_focus"]
    morph_array = morph["scale0"].to_dataset()["image"]
    alphasma_vim = morph_array.sel(c="AlphaSMA/Vimentin").values
    atp1a1 = morph_array.sel(c="ATP1A1/CD45/E-Cadherin").values
    dapi = morph_array.sel(c="DAPI").values

    # 3. Extract DESI channels
    mz_204 = desi_flipped[DESI_CHANNEL_INDICES["mz_204.13"]]
    heme_b = desi_flipped[DESI_CHANNEL_INDICES["heme_B_616.25"]]
    pe_ps_731 = desi_flipped[DESI_CHANNEL_INDICES["PE_PS_731.65"]]
    pe_ps_734 = desi_flipped[DESI_CHANNEL_INDICES["PE_PS_734.60"]]

    # 4. Open napari viewer with individual channels
    viewer = napari.Viewer()

    # Xenium channels (0.2125 um/pixel)
    viewer.add_image(atp1a1, name="ATP1A1 (membranes)", colormap="green",
                     scale=[0.2125, 0.2125], blending="additive", visible=True)
    viewer.add_image(alphasma_vim, name="AlphaSMA/Vimentin (vessels)", colormap="magenta",
                     scale=[0.2125, 0.2125], blending="additive", visible=True)
    viewer.add_image(dapi, name="DAPI (nuclei)", colormap="blue",
                     scale=[0.2125, 0.2125], blending="additive", visible=True)

    # DESI channels (40 um/pixel)
    viewer.add_image(mz_204, name="DESI mz_204.13", colormap="cyan",
                     scale=[40.0, 40.0], blending="additive", visible=False)
    viewer.add_image(heme_b, name="DESI Heme B (616.25)", colormap="yellow",
                     scale=[40.0, 40.0], blending="additive", visible=False)
    viewer.add_image(pe_ps_731, name="DESI PE/PS (731.65)", colormap="red",
                     scale=[40.0, 40.0], blending="additive", visible=False)
    viewer.add_image(pe_ps_734, name="DESI PE/PS (734.60)", colormap="bop orange",
                     scale=[40.0, 40.0], blending="additive", visible=False)

    # 5. RGB composites
    xenium_composite = np.stack([
        normalize(alphasma_vim), normalize(atp1a1),
        np.zeros_like(atp1a1, dtype=float),
    ], axis=-1)

    desi_composite = np.stack([
        normalize(heme_b), normalize(pe_ps_731),
        np.zeros_like(heme_b, dtype=float),
    ], axis=-1)

    viewer.add_image(xenium_composite, name="Xenium Composite (FIXED)", rgb=True,
                     scale=[0.2125, 0.2125], visible=False)
    viewer.add_image(desi_composite, name="DESI Composite (MOVING)", rgb=True,
                     scale=[40.0, 40.0], opacity=0.7, blending="additive", visible=False)

    print("\nnapari viewer opened. Close the window to exit.")
    napari.run()


if __name__ == "__main__":
    main()
