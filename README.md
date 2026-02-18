# Multimodal Spatial Integration

Cross-Modal Spatial Registration of Xenium, METASPACE, and PRIDE-derived Features for Joint Analysis of Tissue Microenvironments

## Project Overview

This project implements a **landmark-based spatial registration pipeline** to align low-resolution DESI (Desorption Electrospray Ionization) mass spectrometry images with high-resolution 10x Genomics Xenium in situ transcriptomics data. The goal is to enable joint analysis of metabolomic and transcriptomic features within the same tissue microenvironment.

### Datasets

| Modality | Sample | Resolution | Channels |
|----------|--------|-----------|----------|
| **Xenium** | FFPE Human Tissue (55588 region 4) | 0.2125 um/pixel | DAPI, ATP1A1/CD45/E-Cadherin, 18S, AlphaSMA/Vimentin |
| **DESI-MSI** | Same tissue section | 40 um/pixel | 88 m/z channels (4 selected: mz 204.13, Heme B 616.25, PE/PS 731.65, PE/PS 734.60) |

### Notebooks

- **`notebooks/Xenium_5k_data_analysis_journey_python.ipynb`** -- Single-cell spatial transcriptomics analysis pipeline (QC, clustering, spatial statistics) for lung and breast cancer Xenium datasets
- **`notebooks/DESI_Xenium_Registration.ipynb`** -- Cross-modal registration of DESI and Xenium images (described below)

## DESI-Xenium Registration Workflow

The registration notebook (`DESI_Xenium_Registration.ipynb`) performs the following steps:

### 1. Data Loading

- Load Xenium output into a `SpatialData` object via `spatialdata_io.xenium()`
- Load the 3D DESI OME-TIFF (88 m/z channels, 107 x 102 pixels at 40 um/pixel)
- Select biologically relevant DESI channels and parse them into a `spatialdata` `Image2DModel` with appropriate scale transforms

### 2. Channel Extraction and Visualization

Extract matched channel pairs from both modalities for visual comparison and registration:

- **Xenium**: AlphaSMA/Vimentin (vessels), ATP1A1/CD45/E-Cadherin (membranes), DAPI (nuclei)
- **DESI**: Heme B 616.25 (vessels), PE/PS 731.65 (membranes), mz 204.13, PE/PS 734.60

The DESI image is flipped along the x-axis to match Xenium orientation.

![DESI and Xenium channels visualized in napari](Channels_DESI_Xenium_Test.png)
*Individual DESI mass spectrometry imaging channels (top row from the left: Heme B 616.25, Xenium Morphology (AlphaSMA/Vimentin and ATP1A1/CD45/E-Cadherin), PE/PS 734.60. Bottom row from the left: m/z 204.13, PE/PS 731.65) visualized together in napari. Each channel is rendered with a distinct colormap to highlight tissue morphology across modalities.*

### 3. Composite Image Generation

RGB composite images are created for both modalities using matched biological features:
- **Red channel**: Vessel markers (AlphaSMA <-> Heme B)
- **Green channel**: Membrane markers (ATP1A1 <-> Lipid 731.65)

### 4. Landmark Selection

Corresponding landmark points are placed on both images using QuPath and exported as TSV files. Landmarks are matched by point number, and pixel coordinates are converted to physical micron coordinates (Xenium x 0.2125, DESI x 40.0).

### 5. SimpleITK Registration Setup

Both images are converted to SimpleITK format with correct physical spacing, and paired landmark coordinates are prepared for landmark-based affine registration.

## Project Structure

```
multimodal_spatial_integration/
├── notebooks/
│   ├── DESI_Xenium_Registration.ipynb
│   └── Xenium_5k_data_analysis_journey_python.ipynb
├── data/
│   ├── xenium/
│   │   ├── raw/          # Xenium instrument output
│   │   └── processed/    # Zarr files, clustering results
│   ├── msi/
│   │   ├── raw/          # DESI OME-TIFF files
│   │   └── processed/    # Flipped/preprocessed DESI images
│   ├── proteomics/       # (planned)
│   └── *.tsv             # Landmark point exports from QuPath
├── Channels_DESI_Xenium_Test.png
├── requirements.txt
├── xenium_scverse_py310.lock.txt
└── README.md
```

## Environment Setup (macOS, Python 3.10)

This project uses the **scverse / spatialdata ecosystem** for analysis of **10x Genomics Xenium In Situ** data. Because `spatialdata`, `spatialdata-io`, `anndata`, and `dask` are tightly coupled, the environment must be **carefully pinned** to ensure stability and reproducibility. The configuration below has been **fully tested and validated** on macOS (Apple Silicon and Intel).

**Python version:** 3.10.x
Python >=3.11 is intentionally not used because `anndata>=0.12` requires Python >=3.11 and is incompatible with this stack, while `scanpy 1.11.x` and `spatialdata 0.4.0` are most stable on Python 3.10.

### Environment creation and activation
```bash
brew install python@3.10
$(brew --prefix python@3.10)/bin/python3.10 -m venv ~/envs/xenium_scverse_py310
source ~/envs/xenium_scverse_py310/bin/activate
pip install --upgrade pip
pip install uv
```

### Dependency specification (requirements.txt)

```
setuptools<81

spatialdata[extra]==0.4.0
spatialdata-io==0.1.7

scanpy==1.11.2
squidpy==1.6.5
anndata==0.10.9

scikit-misc==0.5.1
harmonypy==0.0.10
geosketch==1.3
igraph==0.11.8

dask[dataframe]>=2024.1.0,<2024.8.0
distributed>=2024.1.0,<2024.8.0
dask-expr<1.1.0

```

### Install dependencies
```bash
uv pip install --force-reinstall -r requirements.txt
```

### Dask configuration (required)

Enable the new Dask DataFrame query-planning backend to avoid legacy behavior and future breakage. Add the following line to ~/.zshrc and reload the shell.

```bash
export DASK_DATAFRAME__QUERY_PLANNING=True
source ~/.zshrc
```

### Environment validation (smoke test)

```bash
python -c "
import spatialdata
import spatialdata_io
import scanpy
import squidpy
print('Environment OK')
"
```

### Environment summary

| Package | Version |
|---------|---------|
| Python | 3.10.x |
| spatialdata | 0.4.0 |
| spatialdata-io | 0.1.7 |
| scanpy | 1.11.2 |
| squidpy | 1.6.5 |
| anndata | 0.10.9 |

### Version compatibility notes

- `spatialdata-io` must match the installed spatialdata API
- `anndata>=0.12` requires Python >=3.11 and is not compatible with this setup
- Dask versions >=2024.8 may break spatialdata imports
- All versions listed above were validated together on macOS

### Reproducibility

```bash
uv pip freeze > xenium_scverse_py310.lock.txt
```

## Known Issues and Fixes

### Xenium Parquet Compatibility (Polars -> PyArrow)

Newer Xenium instrument runs (April 2025+) output parquet files written by **Polars**, which use parquet format features that **pyarrow cannot decode**. This causes an `OSError: Unexpected end of stream` when loading data with `spatialdata_io.xenium()`. The files are not corrupted -- they are simply incompatible with pyarrow's reader.

**Affected files:** `cells.parquet`, `transcripts.parquet`

#### Fix

Install Polars in your environment, then re-write the affected files as pyarrow-compatible parquet:

```bash
source ~/envs/xenium_scverse_py310/bin/activate
pip install polars
```

```python
import polars as pl
import shutil
from pathlib import Path

raw_data_path = Path("data/xenium/raw/<your_xenium_output_folder>")

for fname in ["cells.parquet", "transcripts.parquet"]:
    f = raw_data_path / fname
    df = pl.read_parquet(f)
    df.to_pandas().to_parquet(f.with_suffix(".parquet.bak"))
    shutil.move(f, f.with_suffix(".parquet.original"))  # preserve original
    shutil.move(f.with_suffix(".parquet.bak"), f)        # replace with compatible version
```

The original Polars-written files are preserved with a `.original` extension. After conversion, `xenium(raw_data_path)` will load without errors.
