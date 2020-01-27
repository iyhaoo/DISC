#!/bin/bash
#PBS -N CBMC_8K_ds0.5_r5
#PBS -l nodes=cu16
#PBS -l walltime=7200:00:00


R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_5/GSE100866_CBMC_8K_filtered_ds_0.5.loom gene 10 1



