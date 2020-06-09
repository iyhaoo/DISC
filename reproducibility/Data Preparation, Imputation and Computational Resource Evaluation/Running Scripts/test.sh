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





CUDA_VISIBLE_DEVICES=0 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_0 \
--min-expressed-cell=49 \
--library-size-factor=median



CUDA_VISIBLE_DEVICES=1 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1 \
--min-expressed-cell=49 \
--library-size-factor=median




CUDA_VISIBLE_DEVICES=2 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_2 \
--min-expressed-cell=49 \
--library-size-factor=median





CUDA_VISIBLE_DEVICES=3 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1_0 \
--min-expressed-cell=49 \
--library-size-factor=median








R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_0/result/imputation.loom \
/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_0/result/feature.loom \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"


R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1/result/imputation.loom \
/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1/result/feature.loom \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"


R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_2/result/imputation.loom \
/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_2/result/feature.loom \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"






CUDA_VISIBLE_DEVICES=0 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1_1 \
--min-expressed-cell=49 \
--library-size-factor=median



CUDA_VISIBLE_DEVICES=1 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1_2 \
--min-expressed-cell=49 \
--library-size-factor=median



CUDA_VISIBLE_DEVICES=2 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1_3 \
--min-expressed-cell=49 \
--library-size-factor=median



R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_0/result/imputation.loom \
"" \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"


R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_1/result/imputation.loom \
"" \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"


R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_2/result/imputation.loom \
"" \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"


R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_3/result/imputation.loom \
"" \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"





CUDA_VISIBLE_DEVICES=0 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_0 \
--min-expressed-cell=49 \
--library-size-factor=median



CUDA_VISIBLE_DEVICES=2 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_1 \
--min-expressed-cell=49 \
--library-size-factor=median


CUDA_VISIBLE_DEVICES=3 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_2 \
--min-expressed-cell=49 \
--library-size-factor=median



R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_0/result/imputation.loom \
/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_0/result/feature.loom \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"


R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_1/result/imputation.loom \
/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_1/result/feature.loom \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"


R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_2/result/imputation.loom \
/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_2/result/feature.loom \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"



CUDA_VISIBLE_DEVICES=2 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1.0_0 \
--min-expressed-cell=49 \
--library-size-factor=median

CUDA_VISIBLE_DEVICES=3 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1.0_1 \
--min-expressed-cell=49 \
--library-size-factor=median

disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1.0_2 \
--min-expressed-cell=49 \
--library-size-factor=median






R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_0/result/imputation.loom \
/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_0/result/feature.loom \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"


R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_1/result/imputation.loom \
/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_1/result/feature.loom \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"


R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_2/result/imputation.loom \
/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_1_2/result/feature.loom \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"



R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1.0_0/result/imputation.loom \
"" \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"

R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1.0_1/result/imputation.loom \
"" \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"

R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1.0_2/result/imputation.loom \
"" \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"



disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1.0_3 \
--min-expressed-cell=49 \
--library-size-factor=median

R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1.0_3/result/imputation.loom \
"" \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"







CUDA_VISIBLE_DEVICES=2 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_2_0 \
--min-expressed-cell=49 \
--library-size-factor=median

R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_2_0/result/imputation.loom \
/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_2_0/result/feature.loom \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"




CUDA_VISIBLE_DEVICES=3 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_2_1 \
--min-expressed-cell=49 \
--library-size-factor=median

R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_2_1/result/imputation.loom \
/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_2_1/result/feature.loom \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"




CUDA_VISIBLE_DEVICES=3 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_2_2 \
--min-expressed-cell=49 \
--library-size-factor=median

R --slave < "/home/yuanhao/github_repositories/DISC/reproducibility/Down-stream Analysis Improvement/raw_scripts/cluster_fn_0.3.r" \
--args /home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_2_2/result/imputation.loom \
/home/yuanhao/DISC_imputation_result/RETINA_1.1_compressor_2_2/result/feature.loom \
/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/cell_type.rds \
30 1.4 "" "retina"



