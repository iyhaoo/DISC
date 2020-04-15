
disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/JURKAT \
--min-expressed-cell=10 \
--library-size-factor=median


disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/293T \
--min-expressed-cell=10 \
--library-size-factor=median


CUDA_VISIBLE_DEVICES=5 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/JURKAT_293T \
--min-expressed-cell=10 \
--library-size-factor=median

