## Run imputation
We use a Linux CentOS 7 machine which has 2 Intel® Xeon® E5-2650 v4 CPUs, 128GB RAM and 1 NVIDIA® Tesla® V100 GPU.

For Some datasets which have duplicated gene names, we can run a python script in terminal as 
    
    python3 reproducibility/source/loom_rename_duplicated.py \
    --input=matrix.loom

The renamed loom-formatted file will be generated in the same directory of `matrix.loom`.
#### SAVER
Here, 16 means using 16 cores and 10 means removing genes that expressed in less than 10 cells.

    R --slave < reproducibility/source/other_methods/SAVER.r \
    --args /home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom 16 10
#### MAGIC
    python3 reproducibility/source/other_methods/MAGIC.py \
    --loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
    --min-expressed-cell=10

#### DCA
    python3 reproducibility/source/other_methods/DCA.py \
    --loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
    --min-expressed-cell=10
#### scScope
    python3 reproducibility/source/other_methods/scScope.py \
    --loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
    --min-expressed-cell=10
#### scVI
    python3 reproducibility/source/other_methods/scVI.py \
    --loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
    --min-expressed-cell=10
Here, all results will be saved in `/home/yuanhao/data/fn/melanoma/imputation` with genes filtered.
We can easily run

    python3 reproducibility/source/resume_dim.py \
    --raw-loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
    --impute-h5=output_file    
and get `output_file_resume_dim.loom` with complete dimension (only genes as we don't resume cells here) as `/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom`.

#### DISC
For `DISC`, the loom-formatted result with complete genes is saved in `out_dir/result/imputation.loom`, we also provide a low dimensional recovery result which is saved in `out_dir/result/feature.loom` and the reduced dimension can be difined as `50` here. 

    disc \
    --datase=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
    --out-dir=/home/yuanhao/DISC_imputation_result/melanoma \
    --min-expressed-cell=10 \
    --compress-dimensions=50 \
    --library-size-factor=median