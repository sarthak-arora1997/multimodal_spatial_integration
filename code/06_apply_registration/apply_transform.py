"""
Apply the affine registration to resample DESI into Xenium coordinate space.

Loads the matched landmarks from step 05, recomputes the affine transformation,
resamples the 4 selected DESI channels into the Xenium coordinate system using
SimpleITK, stacks them with the Xenium morphology channels, and saves the
combined multi-channel image as an OME-TIFF.
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

OUTPUT_REGISTERED = DATA_DIR / "06_apply_registration" / "registered_desi_xenium.ome.tif"

# Physical pixel sizes (microns per pixel)
XENIUM_PIXEL_SIZE = 0.2125
DESI_PIXEL_SIZE = 40.0

# Selected DESI channels (index into the 88-channel flipped image)
DESI_CHANNELS = {
    "DESI_mz_204.13": 3,
    "DESI_heme_B_616.25": 17,
    "DESI_PE_PS_731.65": 23,
    "DESI_PE_PS_734.60": 25,
}

# Xenium morphology channels to include
XENIUM_CHANNELS = ["DAPI", "ATP1A1/CD45/E-Cadherin", "18S", "AlphaSMA/Vimentin"]


def main():
    # 1. Load landmarks
    print(f"Loading landmarks from:\n  {INPUT_LANDMARKS_CSV}")
    landmarks_df = pd.read_csv(INPUT_LANDMARKS_CSV)
    print(f"  {len(landmarks_df)} landmark pairs")

    # 2. Load Xenium morphology channels
    print("\nLoading Xenium morphology channels...")
    sdata = sd.read_zarr(INPUT_ZARR)
    morph_array = sdata.images["morphology_focus"]["scale0"].to_dataset()["image"]

    xenium_arrays = {}
    for ch_name in XENIUM_CHANNELS:
        xenium_arrays[ch_name] = morph_array.sel(c=ch_name).values
        print(f"  {ch_name}: {xenium_arrays[ch_name].shape}")

    # Use AlphaSMA/Vimentin as the reference image for SimpleITK
    reference = xenium_arrays["AlphaSMA/Vimentin"]

    # 3. Load DESI flipped channels
    print("\nLoading DESI flipped channels...")
    desi_flipped = tifffile.imread(INPUT_DESI_FLIPPED)

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
    registered_desi = {}
    for ch_name, ch_idx in DESI_CHANNELS.items():
        desi_ch_sitk = sitk.GetImageFromArray(desi_flipped[ch_idx])
        desi_ch_sitk.SetSpacing([DESI_PIXEL_SIZE, DESI_PIXEL_SIZE])
        desi_ch_sitk.SetOrigin([0.0, 0.0])

        registered_sitk = sitk.Resample(
            desi_ch_sitk, xenium_sitk, transform,
            sitk.sitkLinear, 0.0, desi_ch_sitk.GetPixelID(),
        )
        registered_desi[ch_name] = sitk.GetArrayFromImage(registered_sitk)
        print(f"  {ch_name}: {registered_desi[ch_name].shape}")

    # 7. Stack all channels: Xenium morphology + registered DESI
    all_channel_names = list(XENIUM_CHANNELS) + list(DESI_CHANNELS.keys())
    all_channels = [xenium_arrays[ch] for ch in XENIUM_CHANNELS] + \
                   [registered_desi[ch] for ch in DESI_CHANNELS]

    stacked = np.stack(all_channels, axis=0)  # (C, Y, X)
    print(f"\nCombined image: {stacked.shape}  ({len(all_channel_names)} channels)")
    for i, name in enumerate(all_channel_names):
        print(f"  Channel {i}: {name}")

    # 8. Save as OME-TIFF with channel metadata
    #    Build OME-XML manually so Bio-Formats/QuPath sees SizeC=N, SizeZ=1, SizeT=1
    OUTPUT_REGISTERED.parent.mkdir(parents=True, exist_ok=True)

    import uuid as _uuid
    import xml.etree.ElementTree as ET

    num_channels = len(all_channel_names)
    h, w = stacked.shape[1], stacked.shape[2]

    # Map numpy dtype → OME pixel type string (all ASCII, no Unicode)
    _ome_type = {
        "uint8": "uint8", "uint16": "uint16", "uint32": "uint32",
        "int8": "int8", "int16": "int16", "int32": "int32",
        "float32": "float", "float64": "double",
    }.get(str(stacked.dtype), "uint16")

    ome_ns = "http://www.openmicroscopy.org/Schemas/OME/2016-06"
    ome = ET.Element("OME", {
        "xmlns": ome_ns,
        "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
        "xsi:schemaLocation": f"{ome_ns} {ome_ns}/ome.xsd",
        "UUID": f"urn:uuid:{_uuid.uuid4()}",
    })
    image = ET.SubElement(ome, "Image", {"ID": "Image:0", "Name": "registered_desi_xenium"})
    pixels = ET.SubElement(image, "Pixels", {
        "ID": "Pixels:0",
        "DimensionOrder": "XYZCT",
        "Type": _ome_type,
        "SizeX": str(w),
        "SizeY": str(h),
        "SizeZ": "1",
        "SizeC": str(num_channels),
        "SizeT": "1",
        "PhysicalSizeX": str(XENIUM_PIXEL_SIZE),
        "PhysicalSizeY": str(XENIUM_PIXEL_SIZE),
    })
    for i, ch_name in enumerate(all_channel_names):
        ET.SubElement(pixels, "Channel", {
            "ID": f"Channel:0:{i}",
            "Name": ch_name,
            "SamplesPerPixel": "1",
        })
    for i in range(num_channels):
        ET.SubElement(pixels, "TiffData", {
            "IFD": str(i),
            "FirstC": str(i),
            "FirstZ": "0",
            "FirstT": "0",
            "PlaneCount": "1",
        })

    ome_xml = '<?xml version="1.0" encoding="UTF-8"?>' + ET.tostring(ome, encoding="unicode")

    print(f"\nWriting OME-TIFF ({num_channels} channels)...")
    with tifffile.TiffWriter(OUTPUT_REGISTERED, bigtiff=True) as tif:
        tif.write(
            stacked,           # (C, Y, X) — all channels in one call → one contiguous series
            photometric="minisblack",
            description=ome_xml,
            metadata=None,     # suppress tifffile's own OME/ImageJ metadata
        )

    print(f"\nSaved OME-TIFF to:\n  {OUTPUT_REGISTERED}")
    print("Done.")


if __name__ == "__main__":
    main()
