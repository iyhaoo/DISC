disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC_TISSUE/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/PBMC_TISSUE_1024 \
--min-expressed-cell=60 \
--library-size-factor=median \
--dimension-number=1024

disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC_TISSUE/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/PBMC_TISSUE_2048 \
--min-expressed-cell=60 \
--library-size-factor=median \
--dimension-number=2048


CUDA_VISIBLE_DEVICES=0 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/MELANOMA_1.1_test0 \
--min-expressed-cell=10 \
--library-size-factor=median \
--model-config-file=/home/yuanhao/github_repositories/DISC/disc/model_config_test.json

CUDA_VISIBLE_DEVICES=1 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/MELANOMA_1.1_test1 \
--min-expressed-cell=10 \
--library-size-factor=median \
--model-config-file=/home/yuanhao/github_repositories/DISC/disc/model_config.json



CUDA_VISIBLE_DEVICES=2 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/SSCORTEX/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/SSCORTEX_1.1 \
--min-expressed-cell=10 \
--library-size-factor=median \
--model-config-file=/home/yuanhao/github_repositories/DISC/disc/model_config.json


CUDA_VISIBLE_DEVICES=3 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1 \
--min-expressed-cell=10 \
--library-size-factor=median \
--model-config-file=/home/yuanhao/github_repositories/DISC/disc/model_config.json






