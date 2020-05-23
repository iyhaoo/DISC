## Run Imputation
We use a Linux CentOS 7 machine which has 2 Intel® Xeon® E5-2650 v4 CPUs, 128GB RAM and 1 NVIDIA® Tesla® V100 GPU.

For some datasets which have duplicated gene names, we run a python script in terminal like this 
    
    python3 reproducibility/source/loom_rename_duplicated.py \
    --input=matrix.loom

The renamed loom-formatted file will be generated in the same folder as `matrix.loom`.
#### SAVER
    R --slave < "./reproducibility/source/Running Scripts for Other Methods/SAVER.r" \
    --args /home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom 16 10
#### scVI
    python3 "./reproducibility/source/Running Scripts for Other Methods/scVI.py" \
    --loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
    --min-expressed-cell=10
#### MAGIC
    python3 "./reproducibility/source/Running Scripts for Other Methods/MAGIC.py" \
    --loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
    --min-expressed-cell=10
#### DCA
    python3 "./reproducibility/source/Running Scripts for Other Methods/DCA.py" \
    --loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
    --min-expressed-cell=10
#### scScope
    python3 "./reproducibility/source/Running Scripts for Other Methods/scScope.py" \
    --loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
    --min-expressed-cell=10
#### DeepImpute
    python3 "./reproducibility/source/Running Scripts for Other Methods/DeepImpute.py" \
    --loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
    --min-expressed-cell=10
#### VIPER
    R --slave < "./reproducibility/source/Running Scripts for Other Methods/VIPER.r" \
    --args /home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom gene 10
#### scImpute
    R --slave < "./reproducibility/source/Running Scripts for Other Methods/scImpute.r" \
    --args /home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom 16 10
    
Here, all results will be saved in `/home/yuanhao/data/fn/melanoma/imputation` with genes filtered.
We can easily run

    python3 reproducibility/source/resume_dim.py \
    --raw-loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
    --impute-h5=output_file
and get `output_file_resume_dim.loom` with complete dimension (only genes as we don't resume cells here) as `/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom`.

#### DISC
For `DISC`, the imputation result (DISC outputs a matrix of the same dimension with `input.loom`, in which genes selected by "gene selection" are updated using DISC imputaion result) is saved in `out_dir/result/imputation.loom`, the low dimensional (default 50) of  DISC imputaion result is saved in `out_dir/result/feature.loom`. 

    disc \
    --dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
    --out-dir=/home/yuanhao/DISC_imputation_result/melanoma \
    --min-expressed-cell=10 \
    --library-size-factor=median