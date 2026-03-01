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
import spatialdata as sd
import tifffile

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"

INPUT_ZARR = DATA_DIR / "01_create_spatial_data_object" / "sdata_object.zarr"
INPUT_DESI_COMPOSITE = DATA_DIR / "02_channel_extraction" / "desi_composite.ome.tif"


def main():
    # 1. Load data
    print(f"Loading SpatialData from:\n  {INPUT_ZARR}")
    sdata = sd.read_zarr(INPUT_ZARR)

    print(f"Loading DESI composite from:\n  {INPUT_DESI_COMPOSITE}")
    desi_composite = tifffile.imread(INPUT_DESI_COMPOSITE)  # (4, Y, X)
    print(f"  DESI composite shape: {desi_composite.shape}, dtype={desi_composite.dtype}")

    # 2. Extract Xenium channels
    morph = sdata.images["morphology_focus"]
    morph_array = morph["scale0"].to_dataset()["image"]
    alphasma_vim = morph_array.sel(c="AlphaSMA/Vimentin").values
    atp1a1 = morph_array.sel(c="ATP1A1/CD45/E-Cadherin").values
    dapi = morph_array.sel(c="DAPI").values

    # 3. Open napari viewer with individual channels
    viewer = napari.Viewer()

    # Xenium channels (0.2125 um/pixel)
    viewer.add_image(atp1a1, name="ATP1A1 (membranes)", colormap="green",
                     scale=[0.2125, 0.2125], blending="additive", visible=True)
    viewer.add_image(alphasma_vim, name="AlphaSMA/Vimentin (vessels)", colormap="magenta",
                     scale=[0.2125, 0.2125], blending="additive", visible=True)
    viewer.add_image(dapi, name="DAPI (nuclei)", colormap="blue",
                     scale=[0.2125, 0.2125], blending="additive", visible=True)

    # DESI channels (40 um/pixel) from composite
    desi_channel_names = list(sdata.images["desi"].coords["c"].values)
    selected_indices = [3, 17, 23, 25]
    desi_colormaps = ["cyan", "yellow", "red", "bop orange"]
    for i, idx in enumerate(selected_indices):
        viewer.add_image(desi_composite[i], name=f"DESI {desi_channel_names[idx]}",
                         colormap=desi_colormaps[i], scale=[40.0, 40.0],
                         blending="additive", visible=False)

    print("\nnapari viewer opened. Close the window to exit.")
    napari.run()


if __name__ == "__main__":
    main()
