
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


R --slave < /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/VIPER.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom gene 10

R --slave < /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/VIPER.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom gene 10

R --slave < /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/VIPER.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom gene 10



R --slave < /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scImpute.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom 16 10

R --slave < /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scImpute.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom 16 10

R --slave < /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scImpute.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom 16 10



python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/MAGIC.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom \
--min-expressed-cell=10

python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/MAGIC.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom \
--min-expressed-cell=10

python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/MAGIC.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom \
--min-expressed-cell=10




python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/deepImpute.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom \
--min-expressed-cell=10

python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/deepImpute.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom \
--min-expressed-cell=10

python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/deepImpute.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom \
--min-expressed-cell=10



python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/DCA.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom \
--min-expressed-cell=10

python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/DCA.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom \
--min-expressed-cell=10

python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/DCA.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom \
--min-expressed-cell=10



python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scVI.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom \
--min-expressed-cell=10

python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scVI.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom \
--min-expressed-cell=10

python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scVI.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom \
--min-expressed-cell=10




python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scScope.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom \
--min-expressed-cell=10

python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scScope.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom \
--min-expressed-cell=10

python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scScope.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom \
--min-expressed-cell=10


###   0.51
disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/JURKAT_0.51 \
--min-expressed-cell=10 \
--library-size-factor=median


disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/293T_0.51 \
--min-expressed-cell=10 \
--library-size-factor=median


CUDA_VISIBLE_DEVICES=5 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/JURKAT_293T_0.51 \
--min-expressed-cell=10 \
--library-size-factor=median
##########


CUDA_VISIBLE_DEVICES=5 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/10X_5CL/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/10X_5CL \
--min-expressed-cell=10 \
--library-size-factor=median


R --slave < /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/VIPER.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/10X_5CL/raw.loom gene 10


R --slave < /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scImpute.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/10X_5CL/raw.loom 16 10


python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/MAGIC.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/10X_5CL/raw.loom \
--min-expressed-cell=10


python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/deepImpute.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/10X_5CL/raw.loom \
--min-expressed-cell=10


python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/DCA.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/10X_5CL/raw.loom \
--min-expressed-cell=10


python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scVI.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/10X_5CL/raw.loom \
--min-expressed-cell=10


python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scScope.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/10X_5CL/raw.loom \
--min-expressed-cell=10





CUDA_VISIBLE_DEVICES=5 HDF5_USE_FILE_LOCKING=FALSE disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC_TISSUE/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/PBMC_TISSUE \
--min-expressed-cell=60 \
--library-size-factor=median


R --slave < /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/VIPER.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC_TISSUE/raw.loom gene 60


#R --slave < /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scImpute.r \
#--args /home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC_TISSUE/raw.loom 16 60


python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/MAGIC.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC_TISSUE/raw.loom \
--min-expressed-cell=60


python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/deepImpute.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC_TISSUE/raw.loom \
--min-expressed-cell=60


python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/DCA.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC_TISSUE/raw.loom \
--min-expressed-cell=60


python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scVI.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC_TISSUE/raw.loom \
--min-expressed-cell=60


python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/other_methods/scScope.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC_TISSUE/raw.loom \
--min-expressed-cell=60


