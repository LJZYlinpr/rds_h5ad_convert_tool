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

# 冒烟测试
conda activate rds_h5ad_env

python -c "import sys; print(sys.executable)"
python -c "import rpy2, anndata, numpy, scipy, pandas, h5py; print('python side ok')"
R -q -e "library(Seurat); sessionInfo()"

# 注册进入内核

python -m pip install ipykernel
python -m ipykernel install --user --name rds_h5ad_env --display-name "Python (rds_h5ad_env)"

# 如何使用
python h5ad_to_rds.py input.h5ad --assay-name RNA
python rds_to_h5ad.py input.rds --assay-name RNA
