# rds_h5ad_convert_tool

Simple converters between Seurat `.rds` objects and AnnData `.h5ad` files.

## Environment

```bash
mamba create -n rds_h5ad_env -c conda-forge \
  python=3.11 \
  r-base=4.5 \
  r-seurat \
  r-matrix \
  rpy2 \
  anndata \
  numpy \
  scipy \
  pandas \
  h5py \
  pip \
  -y
```

## Smoke Test

```bash
conda activate rds_h5ad_env

python -c "import sys; print(sys.executable)"
python -c "import rpy2, anndata, numpy, scipy, pandas, h5py; print('python side ok')"
R -q -e "library(Seurat); sessionInfo()"
```

## Optional: Register Jupyter Kernel

```bash
python -m pip install ipykernel
python -m ipykernel install --user --name rds_h5ad_env --display-name "Python (rds_h5ad_env)"
```

## Usage

RDS to H5AD:

```bash
python rds_to_h5ad.py /path/to/input.rds --assay-name RNA
```

H5AD to RDS:

```bash
python h5ad_to_rds.py /path/to/input.h5ad --assay-name RNA
```
