"""
Create a combined SpatialData object from Xenium and DESI raw data.

Loads 10x Genomics Xenium in situ transcriptomics output and a DESI mass
spectrometry imaging OME-TIFF (all 88 m/z channels), combines them into a
single SpatialData object, and writes it to a Zarr store.
"""

import warnings
warnings.filterwarnings("ignore")

import xml.etree.ElementTree as ET
from pathlib import Path

import tifffile
from spatialdata.models import Image2DModel
from spatialdata.transformations import Scale
from spatialdata_io import xenium

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"

XENIUM_RAW = DATA_DIR / "raw" / "xenium" / "output-XETG00169__0055588__55588_region_4__20250418__182706"
DESI_RAW = DATA_DIR / "raw" / "msi" / "3D_DESI_F5_Pos mode_40um_F5_5pos_40um_Features110425.ome.tif"

OUTPUT_ZARR = DATA_DIR / "01_create_spatial_data_object" / "sdata_object.zarr"


def main():
    # 1. Load Xenium data
    print(f"Loading Xenium data from:\n  {XENIUM_RAW}")
    sdata = xenium(XENIUM_RAW)
    print(f"  Xenium SpatialData loaded: {list(sdata.images.keys())}")

    # 2. Load DESI OME-TIFF (all 88 channels)
    print(f"\nLoading DESI image from:\n  {DESI_RAW}")
    desi_img = tifffile.imread(DESI_RAW)
    print(f"  DESI image shape: {desi_img.shape}  (channels, y, x)")

    # 3. Parse channel names from OME-XML metadata
    with tifffile.TiffFile(DESI_RAW) as tif:
        ome_xml = tif.ome_metadata
    root = ET.fromstring(ome_xml)
    ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
    channels = root.findall(".//ome:Channel", ns)
    channel_names = [ch.get("Name") for ch in channels]
    print(f"  Parsed {len(channel_names)} channel names from OME-XML")

    # 4. Create SpatialData image with 40 um/pixel scale
    desi_spatial = Image2DModel.parse(
        desi_img,
        dims=("c", "y", "x"),
        transformations={
            "global": Scale([40.0, 40.0], axes=("y", "x")),
        },
        c_coords=channel_names,
    )
    sdata.images["desi"] = desi_spatial
    print(f"\n  DESI image added to SpatialData: {sdata.images['desi'].shape}")
    print(f"  Available images: {list(sdata.images.keys())}")

    # 5. Write combined SpatialData object
    print(f"\nWriting SpatialData to:\n  {OUTPUT_ZARR}")
    OUTPUT_ZARR.parent.mkdir(parents=True, exist_ok=True)
    sdata.write(OUTPUT_ZARR, overwrite=True)
    print("  Done.")


if __name__ == "__main__":
    main()
