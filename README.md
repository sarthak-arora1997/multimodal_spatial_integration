# multimodal_spatial_integration
Cross-Modal Spatial Registration of Xenium, METASPACE, and PRIDE-derived Features for Joint Analysis of Tissue Microenvironments

## 🧬 Xenium / scverse Python Environment Setup (macOS, Python 3.10)

This project uses the **scverse / spatialdata ecosystem** for analysis of **10x Genomics Xenium In Situ** data. Because `spatialdata`, `spatialdata-io`, `anndata`, and `dask` are tightly coupled, the environment must be **carefully pinned** to ensure stability and reproducibility. The configuration below has been **fully tested and validated** on macOS (Apple Silicon and Intel).

**Python version:** 3.10.x  
Python ≥3.11 is intentionally not used because `anndata>=0.12` requires Python ≥3.11 and is incompatible with this stack, while `scanpy 1.11.x` and `spatialdata 0.4.0` are most stable on Python 3.10.

### Environment creation and activation
```bash
brew install python@3.10
$(brew --prefix python@3.10)/bin/python3.10 -m venv ~/envs/xenium_scverse_py310
source ~/envs/xenium_scverse_py310/bin/activate
pip install --upgrade pip
pip install uv
```

## Dependency specification (requirements.txt)

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

## Install dependencies
```
uv pip install --force-reinstall -r requirements.txt
```

## Dask configuration (required)

Enable the new Dask DataFrame query-planning backend to avoid legacy behavior and future breakage. Add the following line to ~/.zshrc and reload the shell.

```
export DASK_DATAFRAME__QUERY_PLANNING=True
```

run the below to apply the changes

```
source ~/.zshrc
```

## Environment validation (smoke test)

python -c "
import spatialdata
import spatialdata_io
import scanpy
import squidpy
print('Environment OK')
"

## Version compatibility notes

spatialdata-io must match the installed spatialdata API

anndata>=0.12 requires Python ≥3.11 and is not compatible with this setup

Dask versions ≥2024.8 may break spatialdata imports

All versions listed above were validated together on macOS

## Reproducibility

uv pip freeze > xenium_scverse_py310.lock.txt

## Environment summary

Python: 3.10.x

spatialdata: 0.4.0

spatialdata-io: 0.1.7

scanpy: 1.11.2

squidpy: 1.6.5

anndata: 0.10.9

This environment supports Xenium data ingestion via spatialdata-io, unified multimodal spatial containers (SpatialData), single-cell analysis (scanpy), spatial statistics (squidpy), and interoperability with 10x Genomics Xenium Explorer.

## Fixing Xenium Parquet Compatibility (Polars → PyArrow)

Newer Xenium instrument runs (April 2025+) output parquet files written by **Polars**, which use parquet format features that **pyarrow cannot decode**. This causes an `OSError: Unexpected end of stream` when loading data with `spatialdata_io.xenium()`. The files are not corrupted — they are simply incompatible with pyarrow's reader.

**Affected files:** `cells.parquet`, `transcripts.parquet`

### Fix

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
