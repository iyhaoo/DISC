
disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/MELANOMA \
--min-expressed-cell=10 \
--library-size-factor=median

disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/SSCORTEX/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/SSCORTEX \
--min-expressed-cell=10 \
--library-size-factor=median

disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/RETINA \
--min-expressed-cell=49 \
--library-size-factor=median

disc \
--dataset=/home/yuanhao/github_repositories/DISC/reproducibility/data/BRAIN_SPLiT/raw.loom \
--out-dir=/home/yuanhao/DISC_imputation_result/BRAIN_SPLiT \
--min-expressed-cell=150 \
--library-size-factor=median







R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/SSCORTEX/raw.loom gene 10 1



python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/github_repositories/DISC/reproducibility/data/SSCORTEX/raw.loom \
--min-expressed-cell=10







R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_1/dropseq_filt_ls_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_2/dropseq_filt_ls_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_3/dropseq_filt_ls_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_4/dropseq_filt_ls_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_5/dropseq_filt_ls_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_1/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_2/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_3/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_4/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_5/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom gene 10 1



python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_1/dropseq_filt_ls_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_2/dropseq_filt_ls_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_3/dropseq_filt_ls_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_4/dropseq_filt_ls_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_5/dropseq_filt_ls_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_1/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_2/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_3/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_4/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_5/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom \
--min-expressed-cell=10


