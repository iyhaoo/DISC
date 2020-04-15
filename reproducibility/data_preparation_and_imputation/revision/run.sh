
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


R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom gene 10 1



R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom 16 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom 16 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom 16 10 1



python3 reproducibility/source/other_methods/MAGIC.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom \
--min-expressed-cell=10

python3 reproducibility/source/other_methods/MAGIC.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom \
--min-expressed-cell=10

python3 reproducibility/source/other_methods/MAGIC.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom \
--min-expressed-cell=10




python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom \
--min-expressed-cell=10



python3 reproducibility/source/other_methods/DCA.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom \
--min-expressed-cell=10

python3 reproducibility/source/other_methods/DCA.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom \
--min-expressed-cell=10

python3 reproducibility/source/other_methods/DCA.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom \
--min-expressed-cell=10



python3 reproducibility/source/other_methods/scVI.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom \
--min-expressed-cell=10

python3 reproducibility/source/other_methods/scVI.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom \
--min-expressed-cell=10

python3 reproducibility/source/other_methods/scVI.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom \
--min-expressed-cell=10




python3 reproducibility/source/other_methods/scScope.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom \
--min-expressed-cell=10

python3 reproducibility/source/other_methods/scScope.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom \
--min-expressed-cell=10

python3 reproducibility/source/other_methods/scScope.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom \
--min-expressed-cell=10

