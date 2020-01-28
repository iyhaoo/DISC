
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


R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/raw.loom 16 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/github_repositories/DISC/reproducibility/data/SSCORTEX/raw.loom 16 10 1







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


R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_1/GSE100866_CBMC_8K_filtered_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_2/GSE100866_CBMC_8K_filtered_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_3/GSE100866_CBMC_8K_filtered_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_4/GSE100866_CBMC_8K_filtered_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_5/GSE100866_CBMC_8K_filtered_ds_0.5.loom gene 10 1


R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_1/pbmc3k_filtered_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_2/pbmc3k_filtered_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_3/pbmc3k_filtered_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_4/pbmc3k_filtered_ds_0.5.loom gene 10 1

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_5/pbmc3k_filtered_ds_0.5.loom gene 10 1



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


python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_1/GSE100866_CBMC_8K_filtered_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_2/GSE100866_CBMC_8K_filtered_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_3/GSE100866_CBMC_8K_filtered_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_4/GSE100866_CBMC_8K_filtered_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_5/GSE100866_CBMC_8K_filtered_ds_0.5.loom \
--min-expressed-cell=10


python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_1/pbmc3k_filtered_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_2/pbmc3k_filtered_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_3/pbmc3k_filtered_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_4/pbmc3k_filtered_ds_0.5.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_5/pbmc3k_filtered_ds_0.5.loom \
--min-expressed-cell=10



python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_1/pbmc3k_filtered_ds_0.3.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_2/pbmc3k_filtered_ds_0.3.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_3/pbmc3k_filtered_ds_0.3.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_4/pbmc3k_filtered_ds_0.3.loom \
--min-expressed-cell=10

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_5/pbmc3k_filtered_ds_0.3.loom \
--min-expressed-cell=10





python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_0/merged_set_unique_rename_ds_0.5.loom \
--min-expressed-cell=49

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_1/merged_set_unique_rename_ds_0.5.loom \
--min-expressed-cell=49

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_2/merged_set_unique_rename_ds_0.5.loom \
--min-expressed-cell=49

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_3/merged_set_unique_rename_ds_0.5.loom \
--min-expressed-cell=49

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_4/merged_set_unique_rename_ds_0.5.loom \
--min-expressed-cell=49

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_5/merged_set_unique_rename_ds_0.5.loom \
--min-expressed-cell=49

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_6/merged_set_ds_0.5.loom \
--min-expressed-cell=49





python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_0/merged_set_unique_rename_ds_0.3.loom \
--min-expressed-cell=49

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_1/merged_set_unique_rename_ds_0.3.loom \
--min-expressed-cell=49

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_2/merged_set_unique_rename_ds_0.3.loom \
--min-expressed-cell=49

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_3/merged_set_unique_rename_ds_0.3.loom \
--min-expressed-cell=49

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_4/merged_set_unique_rename_ds_0.3.loom \
--min-expressed-cell=49

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_5/merged_set_unique_rename_ds_0.3.loom \
--min-expressed-cell=49

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_6/merged_set_ds_0.3.loom \
--min-expressed-cell=49






python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_1/GSM3017261_150000_CNS_nuclei_ds_0.5.loom \
--min-expressed-cell=150

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_2/GSM3017261_150000_CNS_nuclei_ds_0.5.loom \
--min-expressed-cell=150

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_3/GSM3017261_150000_CNS_nuclei_ds_0.5.loom \
--min-expressed-cell=150

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_4/GSM3017261_150000_CNS_nuclei_ds_0.5.loom \
--min-expressed-cell=150

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_5/GSM3017261_150000_CNS_nuclei_ds_0.5.loom \
--min-expressed-cell=150



python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_1/GSM3017261_150000_CNS_nuclei_ds_0.3.loom \
--min-expressed-cell=150

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_2/GSM3017261_150000_CNS_nuclei_ds_0.3.loom \
--min-expressed-cell=150

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_3/GSM3017261_150000_CNS_nuclei_ds_0.3.loom \
--min-expressed-cell=150

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_4/GSM3017261_150000_CNS_nuclei_ds_0.3.loom \
--min-expressed-cell=150

python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_5/GSM3017261_150000_CNS_nuclei_ds_0.3.loom \
--min-expressed-cell=150





disc \
--dataset=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_1/GSM3017261_150000_CNS_nuclei_ds_0.3.loom \
--out-dir=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_1/DISC_ds_0.3 \
--min-expressed-cell=150 \
--library-size-factor=median

disc \
--dataset=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_2/GSM3017261_150000_CNS_nuclei_ds_0.3.loom \
--out-dir=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_2/DISC_ds_0.3 \
--min-expressed-cell=150 \
--library-size-factor=median

disc \
--dataset=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_3/GSM3017261_150000_CNS_nuclei_ds_0.3.loom \
--out-dir=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_3/DISC_ds_0.3 \
--min-expressed-cell=150 \
--library-size-factor=median

disc \
--dataset=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_4/GSM3017261_150000_CNS_nuclei_ds_0.3.loom \
--out-dir=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_4/DISC_ds_0.3 \
--min-expressed-cell=150 \
--library-size-factor=median

disc \
--dataset=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_5/GSM3017261_150000_CNS_nuclei_ds_0.3.loom \
--out-dir=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_5/DISC_ds_0.3 \
--min-expressed-cell=150 \
--library-size-factor=median




disc \
--dataset=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_1/GSM3017261_150000_CNS_nuclei_ds_0.5.loom \
--out-dir=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_1/DISC_ds_0.5 \
--min-expressed-cell=150 \
--library-size-factor=median

disc \
--dataset=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_2/GSM3017261_150000_CNS_nuclei_ds_0.5.loom \
--out-dir=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_2/DISC_ds_0.5 \
--min-expressed-cell=150 \
--library-size-factor=median

disc \
--dataset=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_3/GSM3017261_150000_CNS_nuclei_ds_0.5.loom \
--out-dir=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_3/DISC_ds_0.5 \
--min-expressed-cell=150 \
--library-size-factor=median

disc \
--dataset=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_4/GSM3017261_150000_CNS_nuclei_ds_0.5.loom \
--out-dir=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_4/DISC_ds_0.5 \
--min-expressed-cell=150 \
--library-size-factor=median

disc \
--dataset=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_5/GSM3017261_150000_CNS_nuclei_ds_0.5.loom \
--out-dir=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_5/DISC_ds_0.5 \
--min-expressed-cell=150 \
--library-size-factor=median








R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_1/dropseq_filt_ls_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_2/dropseq_filt_ls_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_3/dropseq_filt_ls_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_4/dropseq_filt_ls_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_5/dropseq_filt_ls_ds_0.5.loom 16


R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_1/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_2/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_3/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_4/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_5/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_0.5.loom 16


R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_1/GSE100866_CBMC_8K_filtered_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_2/GSE100866_CBMC_8K_filtered_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_3/GSE100866_CBMC_8K_filtered_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_4/GSE100866_CBMC_8K_filtered_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_5/GSE100866_CBMC_8K_filtered_ds_0.5.loom 16


R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_1/pbmc3k_filtered_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_2/pbmc3k_filtered_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_3/pbmc3k_filtered_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_4/pbmc3k_filtered_ds_0.5.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_5/pbmc3k_filtered_ds_0.5.loom 16




python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/retina/merged_set.loom \
--min-expressed-cell=49



python3 /home/yuanhao/github_repositories/DISC/reproducibility/source/resume_dim.py \
--raw-loom=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_5/GSM3017261_150000_CNS_nuclei_ds_0.5.loom \
--impute-h5=/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_5/imputation/GSM3017261_150000_CNS_nuclei_ds_0.5_deepImpute_mc_150_mce_1.hdf5


R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/cluster_evaluation_retina1_with_outliers.r \
--args /home/yuanhao/DISC_imputation_result/0.45/RETINA2/result/imputation.loom \
"" \
"" \
/home/yuanhao/data/fn/retina/retina_cell_type_main.tsv



R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/cluster_fn_0.2.r \
--args /home/yuanhao/DISC_imputation_result/0.45/PBMC_transfer/result/imputation.loom \
"" /home/yuanhao/data/fn/pbmc3k/cluster_evaluation/pbmc3k_filtered/pca/cell_type.rds \
10 0.5 "" "pbmc"






R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_0/merged_set_unique_rename_ds_0.5.loom gene 49

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_1/merged_set_unique_rename_ds_0.5.loom gene 49

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_2/merged_set_unique_rename_ds_0.5.loom gene 49

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_3/merged_set_unique_rename_ds_0.5.loom gene 49

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_4/merged_set_unique_rename_ds_0.5.loom gene 49

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_5/merged_set_unique_rename_ds_0.5.loom gene 49

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_6/merged_set_ds_0.5.loom gene 49




R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_0/merged_set_unique_rename_ds_0.5.loom 16 49


R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_1/pbmc3k_filtered_ds_0.3.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_2/pbmc3k_filtered_ds_0.3.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_3/pbmc3k_filtered_ds_0.3.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_4/pbmc3k_filtered_ds_0.3.loom 16

R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_5/pbmc3k_filtered_ds_0.3.loom 16





