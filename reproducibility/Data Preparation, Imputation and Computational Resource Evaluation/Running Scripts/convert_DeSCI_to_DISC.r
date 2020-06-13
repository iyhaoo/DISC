setwd("E:/DISC/reproducibility")
utilities_path = "./source/utilities.r"
source(utilities_path)
DeSCI_imputation = readh5_imputation("./data/MELANOMA/DISC.hdf5", with_outliers = T)
save_h5("./data/MELANOMA/DISC.loom", t(DeSCI_imputation))
DeSCI_imputation = readh5_imputation("./data/MELANOMA/ds_0.5/DISC.hdf5", with_outliers = T)
save_h5("./data/MELANOMA/ds_0.5/DISC.loom", t(DeSCI_imputation))
DeSCI_imputation = readh5_imputation("./data/PBMC/ds_0.5/DISC.hdf5", with_outliers = T)
save_h5("./data/PBMC/ds_0.5/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_2/DeSCI_2.7.4.33/ds_0.5/results/epoch_1170/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_2/DeSCI_2.7.4.33/ds_0.5/results/epoch_1170/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_1/DeSCI_2.7.4.33/ds_0.5/results/epoch_1170/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_1/DeSCI_2.7.4.33/ds_0.5/results/epoch_1170/DISC.loom", t(DeSCI_imputation))

get_optimal_point33("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_3/DeSCI_2.7.4.33/ds_0.5/log.txt")


DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_3/DeSCI_2.7.4.33/ds_0.5/results/epoch_1170/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_3/DeSCI_2.7.4.33/ds_0.5/results/epoch_1170/DISC.loom", t(DeSCI_imputation))


get_optimal_point33("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_4/DeSCI_2.7.4.33/ds_0.5/log.txt")
DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_4/DeSCI_2.7.4.33/ds_0.5/results/epoch_1170/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_4/DeSCI_2.7.4.33/ds_0.5/results/epoch_1170/DISC.loom", t(DeSCI_imputation))



get_optimal_point33("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_5/DeSCI_2.7.4.33/ds_0.5/log.txt")
DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_5/DeSCI_2.7.4.33/ds_0.5/results/epoch_1170/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_5/DeSCI_2.7.4.33/ds_0.5/results/epoch_1170/DISC.loom", t(DeSCI_imputation))



get_optimal_point33("/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_1/DeSCI_2.7.4.33/ds_0.5/log.txt")

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_1/DeSCI_2.7.4.33/ds_0.5/results/epoch_1248/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/CBMC/ds_0.5/r1/DISC.loom", t(DeSCI_imputation))




get_optimal_point33("/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_2/DeSCI_2.7.4.33/ds_0.5/log.txt")

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_2/DeSCI_2.7.4.33/ds_0.5/results/epoch_1248/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/CBMC/ds_0.5/r2/DISC.loom", t(DeSCI_imputation))


get_optimal_point33("/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_3/DeSCI_2.7.4.33/ds_0.5/log.txt")

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_3/DeSCI_2.7.4.33/ds_0.5/results/epoch_1248/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/CBMC/ds_0.5/r3/DISC.loom", t(DeSCI_imputation))


get_optimal_point33("/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_4/DeSCI_2.7.4.33/ds_0.5/log.txt")

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_4/DeSCI_2.7.4.33/ds_0.5/results/epoch_1248/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/CBMC/ds_0.5/r4/DISC.loom", t(DeSCI_imputation))


get_optimal_point33("/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_5/DeSCI_2.7.4.33/ds_0.5/log.txt")

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/CITE-seq/downsampling_first_repeat_5/DeSCI_2.7.4.33/ds_0.5/results/epoch_1248/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/CBMC/ds_0.5/r5/DISC.loom", t(DeSCI_imputation))




DeSCI_imputation = readh5_imputation(get_optimal_point33("/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_1/DeSCI_2.7.4.33/ds_0.5/log.txt"), with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/SSCORTEX/ds_0.5/r1/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation(get_optimal_point33("/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_2/DeSCI_2.7.4.33/ds_0.5/log.txt"), with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/SSCORTEX/ds_0.5/r2/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation(get_optimal_point33("/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_3/DeSCI_2.7.4.33/ds_0.5/log.txt"), with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/SSCORTEX/ds_0.5/r3/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation(get_optimal_point33("/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_4/DeSCI_2.7.4.33/ds_0.5/log.txt"), with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/SSCORTEX/ds_0.5/r4/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation(get_optimal_point33("/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds/downsampling_first_repeat_5/DeSCI_2.7.4.33/ds_0.5/log.txt"), with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/SSCORTEX/ds_0.5/r5/DISC.loom", t(DeSCI_imputation))





DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_1/DeSCI_2.7.4.17/ds_0.5/results/epoch_8/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/BRAIN_SPLiT/ds_0.5/r1/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_2/DeSCI_2.7.4.17/ds_0.5/results/epoch_8/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/BRAIN_SPLiT/ds_0.5/r2/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_3/DeSCI_2.7.4.17/ds_0.5/results/epoch_8/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/BRAIN_SPLiT/ds_0.5/r3/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_4/DeSCI_2.7.4.17/ds_0.5/results/epoch_8/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/BRAIN_SPLiT/ds_0.5/r4/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_5/DeSCI_2.7.4.17/ds_0.5/results/epoch_8/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/BRAIN_SPLiT/ds_0.5/r5/DISC.loom", t(DeSCI_imputation))




DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_1/DeSCI_2.7.4.17/ds_0.3/results/epoch_8/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/BRAIN_SPLiT/ds_0.3/r1/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_2/DeSCI_2.7.4.17/ds_0.3/results/epoch_8/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/BRAIN_SPLiT/ds_0.3/r2/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_3/DeSCI_2.7.4.17/ds_0.3/results/epoch_8/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/BRAIN_SPLiT/ds_0.3/r3/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_4/DeSCI_2.7.4.17/ds_0.3/results/epoch_8/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/BRAIN_SPLiT/ds_0.3/r4/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_5/DeSCI_2.7.4.17/ds_0.3/results/epoch_8/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/BRAIN_SPLiT/ds_0.3/r5/DISC.loom", t(DeSCI_imputation))




DeSCI_imputation = readh5_imputation(get_optimal_point33("/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_1/DeSCI_2.7.4.33/ds_0.5/log.txt"), with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/ds_0.5/r1/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation(get_optimal_point33("/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_0/DeSCI_2.7.4.33/ds_0.5/log.txt"), with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/ds_0.5/r2/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_loom("/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_3/ds_0.5_new/result/imputation.loom")
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/ds_0.5/r3/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_6/DeSCI_2.7.4.17/ds_0.5/results/epoch_102/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/ds_0.5/r4/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_5/DeSCI_2.7.4.17/ds_0.5/results/epoch_60/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/ds_0.5/r5/DISC.loom", t(DeSCI_imputation))




DeSCI_imputation = readh5_imputation(get_optimal_point33("/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_1/DeSCI_2.7.4.33/ds_0.3/log.txt"), with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/ds_0.3/r1/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation(get_optimal_point33("/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_0/DeSCI_2.7.4.33/ds_0.3/log.txt"), with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/ds_0.3/r2/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_loom("/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_3/ds_0.3_new/result/imputation.loom")
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/ds_0.3/r3/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_6/DeSCI_2.7.4.17/ds_0.3/results/epoch_108/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/ds_0.3/r4/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/retina/ds/downsampling_first_repeat_5/DeSCI_2.7.4.17/ds_0.3/results/epoch_120/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/RETINA/ds_0.3/r5/DISC.loom", t(DeSCI_imputation))





source("/home/yuanhao/github_repositories/DISC/reproducibility/source/utilities.r")



DeSCI_imputation = readh5_imputation(get_optimal_point33("/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_1/DeSCI_2.7.4.33/ds_0.5/log.txt"), with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC/ds_0.5/r1/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation(get_optimal_point33("/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_2_showed_0.3/DeSCI_2.7.4.33/ds_0.5/log.txt"), with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC/ds_0.5/r2/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation(get_optimal_point33("/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_3_showed_0.5/DeSCI_2.7.4.33/ds_0.5/log.txt"), with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC/ds_0.5/r3/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation(get_optimal_point33("/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_4/DeSCI_2.7.4.33/ds_0.5/log.txt"), with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC/ds_0.5/r4/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation(get_optimal_point33("/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_5/DeSCI_2.7.4.33/ds_0.5/log.txt"), with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC/ds_0.5/r5/DISC.loom", t(DeSCI_imputation))



DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_1/DeSCI_2.7.4.33/ds_0.3/results/epoch_3515/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC/ds_0.3/r1/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_2_showed_0.3/DeSCI_2.7.4.33/ds_0.3/results/epoch_3515/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC/ds_0.3/r2/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_3_showed_0.5/DeSCI_2.7.4.33/ds_0.3/results/epoch_3515/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC/ds_0.3/r3/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_4/DeSCI_2.7.4.33/ds_0.3/results/epoch_3515/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC/ds_0.3/r4/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_5/DeSCI_2.7.4.33/ds_0.3/results/epoch_3515/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/PBMC/ds_0.3/r5/DISC.loom", t(DeSCI_imputation))





DISC_imputation = readh5_loom("/home/yuanhao/DISC_imputation_result/BONE_MARROW/result/imputation.loom")
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/BONE_MARROW/DISC.loom", t(DISC_imputation))



DISC_imputation = readh5_loom("/home/yuanhao/DISC_imputation_result/SSCORTEX/result/imputation.loom")
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/SSCORTEX/DISC.loom", t(DISC_imputation))



DeSCI_imputation = readh5_imputation("/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/DISC.hdf5", with_outliers = T)
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/DISC.loom", t(DeSCI_imputation))



gene_bc_mat = readh5_loom("E:/DISC/reproducibility/data/JURKAT_293T/raw.loom")
DISC_imputation = readh5_loom("E:/DISC/reproducibility/data/JURKAT_293T/DISC_JURKAT_293T.loom")
print(sum(rownames(DISC_imputation) != rownames(gene_bc_mat)))
rownames(DISC_imputation) = rownames(gene_bc_mat)
save_h5("E:/DISC/reproducibility/data/JURKAT_293T/DISC_JURKAT_293T.loom", t(DISC_imputation))



gene_bc_filt = gene_bc_mat[gene_selection(gene_bc_mat, 10), ]
for(ii in c("scVI", "MAGIC", "DCA", "scScope", "DeepImpute", "VIPER", "scImpute")){
  imputation_result = readh5_imputation(paste0("E:/DISC/reproducibility/data/JURKAT_293T/", ii, "_JURKAT_293T.hdf5"))
  print(sum(rownames(imputation_result) != rownames(gene_bc_filt)))
  rownames(imputation_result) = rownames(gene_bc_filt)
  save_h5(paste0("E:/DISC/reproducibility/data/JURKAT_293T/", ii, "_JURKAT_293T.loom"), t(imputation_result))
}

imputation_result = readh5_loom(paste0("E:/DISC/reproducibility/data/JURKAT_293T/gene_selection_JURKAT_293T.loom"))
print(sum(rownames(imputation_result) != rownames(gene_bc_filt)))
rownames(imputation_result) = rownames(gene_bc_filt)
save_h5(paste0("E:/DISC/reproducibility/data/JURKAT_293T/gene_selection_JURKAT_293T.loom"), t(imputation_result))




DISC_imputation = readh5_loom("/home/yuanhao/DISC_imputation_result/10X_5CL/result/imputation.loom")
save_h5("/home/yuanhao/github_repositories/DISC/reproducibility/data/10X_5CL/DISC.loom", t(DISC_imputation))


gene_bc_mat = readh5_loom("E:/DISC/reproducibility/data/10X_5CL/raw.loom")
for(ii in c("DISC", "gene_selection")){
  imputation_result = readh5_loom(paste0("E:/DISC/reproducibility/data/10X_5CL/", ii, ".loom"))
  colnames(imputation_result) = colnames(gene_bc_mat)
  save_h5(paste0("E:/DISC/reproducibility/data/10X_5CL/", ii, ".loom"), t(imputation_result))
  print(ii)
}

for(ii in c("scVI", "MAGIC", "DCA", "scScope", "DeepImpute", "VIPER", "scImpute")){
  imputation_result = readh5_imputation(paste0("E:/DISC/reproducibility/data/10X_5CL/", ii, ".hdf5"))
  gene_bc_mat = readh5_loom("E:/DISC/reproducibility/data/10X_5CL/raw.loom")
  gene_bc_mat[rownames(imputation_result), ] = imputation_result
  save_h5(paste0("E:/DISC/reproducibility/data/10X_5CL/", ii, ".loom"), t(gene_bc_mat))
}

for(ii in c("scVI", "MAGIC", "DCA", "scScope", "DeepImpute", "VIPER", "scImpute")){
  imputation_result = readh5_loom(paste0("E:/DISC/reproducibility/data/JURKAT_293T/", ii, "_JURKAT_293T.loom"))
  gene_bc_mat = readh5_loom("E:/DISC/reproducibility/data/JURKAT_293T/raw.loom")
  gene_bc_mat[rownames(imputation_result), ] = imputation_result
  save_h5(paste0("E:/DISC/reproducibility/data/JURKAT_293T/", ii, "_JURKAT_293T.loom"), t(gene_bc_mat))
}


