setwd("/home/yuanhao/github_repositories/DISC/reproducibility")
utilities_path = "./source/utilities.r"
source(utilities_path)
#  Use the same color as Figure 4
method_name = c("Observed", "DISC", "scImpute", "VIPER", "MAGIC", "DCA", "DeepImpute", "scScope", "scVI")
method_color = c("#A5A5A5", "#E83828", "#278BC4", "#EADE36", "#198B41", "#920783", "#F8B62D", "#8E5E32", "#1E2B68")
names(method_color) = method_name
text_color = rep("black", length(method_name))
### Load raw data and downsampling
raw_data = readh5_loom("./data/PBMC/raw.loom")
compared_genes = rownames(raw_data)
cell_number = ncol(raw_data)
data_list = list(Raw = raw_data)
rm(raw_data)
ds_dir = "./data/PBMC/ds_0.5"
dir.create(ds_dir, showWarnings = F, recursive = T)
output_dir = "./results/PBMC/ds_0.5"
dir.create(output_dir, showWarnings = F, recursive = T)
repeats = c("r1", "r2", "r3", "r4", "r5")
data_list[["Observed"]] = list()
for(ii in repeats){
  observed_path = paste(ds_dir, ii, "gene_selection.loom", sep = "/")
  data_list[["Observed"]][[ii]] = readh5_loom(observed_path)
  compared_genes = intersect(compared_genes, rownames(data_list[["Observed"]][[ii]]))
}
gene_number = length(compared_genes)
data_list[["Raw"]] = data_list[["Raw"]][compared_genes, ]
### load vst file
vst_file = paste0(output_dir, "/vst_gene.tsv")
if(file.exists(vst_file)){
  hvg_info = read.table(vst_file)
  print("load vst_file")
}else{
  hvg_info = FindVariableFeatures_vst_by_genes(data_list[["Raw"]])
  hvg_info = hvg_info[order(hvg_info$variance.standardized, decreasing = T), ]
  write.table(hvg_info, paste0(output_dir, "/vst_gene.tsv"), sep = "\t", quote = F, row.names = T, col.names = T)
}
hvg_genes = rownames(hvg_info)[rownames(hvg_info) %in% compared_genes]
cat("Use ", gene_number, " genes for comparison.\n")
#  Read imputation results
for(ii in repeats){
  data_list[["Observed"]][[ii]] = data_list[["Observed"]][[ii]][compared_genes, ]
  for(jj in setdiff(method_name, "Observed")){
    if(jj == "DISC"){
      file_format = "loom"
    }else{
      file_format = "hdf5"
    }
    data_list[[jj]][[ii]] = get_gene_bc_mat(paste0(ds_dir, "/", ii, "/", jj, ".", file_format))[compared_genes, ]
    print(dim(data_list[[jj]][[ii]]))
  }
}
result_list = list()
### MAE
#  Here, we only evaluate the expressed genes in RAW (before down-sampling) dataset.
MAE_O_mat = matrix(nrow = cell_number, ncol = length(method_name), dimnames = list(c(), method_name)) # Overall
MAE_DID_mat = matrix(nrow = cell_number, ncol = length(method_name), dimnames = list(c(), method_name)) # Dropouts Introduced by Down-sampling
MAE_RGE_mat = matrix(nrow = cell_number, ncol = length(method_name), dimnames = list(c(), method_name)) # Retained Gene Expressions
top_1000_genes = hvg_genes[1:1000]
top_1000_index = which(rownames(data_list[["Raw"]]) %in% top_1000_genes)
top_MAE_O_mat = matrix(nrow = cell_number, ncol = length(method_name), dimnames = list(c(), method_name)) # Overall
top_MAE_DID_mat = matrix(nrow = cell_number, ncol = length(method_name), dimnames = list(c(), method_name)) # Dropouts Introduced by Down-sampling
top_MAE_RGE_mat = matrix(nrow = cell_number, ncol = length(method_name), dimnames = list(c(), method_name)) # Retained Gene Expressions
ls_raw = colSums(data_list[["Raw"]])
for(ii in method_name){
  for(jj in repeats){
    ls_this = colSums(data_list[[ii]][[jj]])
    scale_factor = ls_raw / ls_this
    MAE_cell = sapply(seq(cell_number), function(x){
      expressed_index = which(data_list[["Raw"]][, x] > 0)
      DID_mask_in_expressed_index = expressed_index %in% which(data_list[["Observed"]][[jj]][, x] == 0)
      RGE_mask_in_expressed_index = expressed_index %in% which(data_list[["Observed"]][[jj]][, x] > 0)
      top_expressed_mask = expressed_index %in% top_1000_index
      top_DID_mask = DID_mask_in_expressed_index & top_expressed_mask
      top_RGE_mask = RGE_mask_in_expressed_index & top_expressed_mask
      #  Results - All genes
      error_O = data_list[["Raw"]][expressed_index, x] - (data_list[[ii]][[jj]][expressed_index, x] * scale_factor[x])
      MAE_O = sum(abs(error_O)) / length(expressed_index)
      MAE_DID = sum(abs(error_O[DID_mask_in_expressed_index])) / sum(DID_mask_in_expressed_index)
      MAE_RGE = sum(abs(error_O[RGE_mask_in_expressed_index])) / sum(RGE_mask_in_expressed_index)
      #  Results - Top genes
      top_MAE_O = sum(abs(error_O[top_expressed_mask])) / sum(top_expressed_mask)
      top_MAE_DID = sum(abs(error_O[top_DID_mask])) / sum(top_DID_mask)
      top_MAE_RGE = sum(abs(error_O[top_RGE_mask])) / sum(top_RGE_mask)
      return(c(MAE_O, MAE_DID, MAE_RGE, top_MAE_O, top_MAE_DID, top_MAE_RGE))
    })
    if(jj == repeats[1]){
      MAE_O = MAE_cell[1, ]
      MAE_DID = MAE_cell[2, ]
      MAE_RGE = MAE_cell[3, ]
      top_MAE_O = MAE_cell[4, ]
      top_MAE_DID = MAE_cell[5, ]
      top_MAE_RGE = MAE_cell[6, ]
    }else{
      MAE_O = cbind(MAE_O, MAE_cell[1, ])
      MAE_DID = cbind(MAE_DID, MAE_cell[2, ])
      MAE_RGE = cbind(MAE_RGE, MAE_cell[3, ])
      top_MAE_O = cbind(top_MAE_O, MAE_cell[4, ])
      top_MAE_DID = cbind(top_MAE_DID, MAE_cell[5, ])
      top_MAE_RGE = cbind(top_MAE_RGE, MAE_cell[6, ])
    }
    cat("Finish: ", ii, " - ", jj, "\n")
  }
  MAE_O_mat[, ii] = rowMeans(MAE_O)
  MAE_DID_mat[, ii] = rowMeans(MAE_DID)
  MAE_RGE_mat[, ii] = rowMeans(MAE_RGE)
  top_MAE_O_mat[, ii] = rowMeans(top_MAE_O)
  top_MAE_DID_mat[, ii] = rowMeans(top_MAE_DID)
  top_MAE_RGE_mat[, ii] = rowMeans(top_MAE_RGE)
  print(ii)
}
MAE_O_mat = MAE_O_mat[rowSums(is.na(MAE_O_mat)) < 1, ]
MAE_DID_mat = MAE_DID_mat[rowSums(is.na(MAE_DID_mat)) < 1, ]
MAE_RGE_mat = MAE_RGE_mat[rowSums(is.na(MAE_RGE_mat)) < 1, ]
top_MAE_O_mat = top_MAE_O_mat[rowSums(is.na(top_MAE_O_mat)) < 1, ]
top_MAE_DID_mat = top_MAE_DID_mat[rowSums(is.na(top_MAE_DID_mat)) < 1, ]
top_MAE_RGE_mat = top_MAE_RGE_mat[rowSums(is.na(top_MAE_RGE_mat)) < 1, ]
result_list[["MAE_O"]] = MAE_O_mat
result_list[["MAE_DID"]] = MAE_DID_mat
result_list[["MAE_RGE"]] = MAE_RGE_mat
result_list[["top_MAE_O"]] = top_MAE_O_mat
result_list[["top_MAE_DID"]] = top_MAE_DID_mat
result_list[["top_MAE_RGE"]] = top_MAE_RGE_mat
MAE_O_df = melt(t(MAE_O_mat))
MAE_DID_df = melt(t(MAE_DID_mat))
MAE_RGE_df = melt(t(MAE_RGE_mat))
top_MAE_O_df = melt(t(top_MAE_O_mat))
top_MAE_DID_df = melt(t(top_MAE_DID_mat))
top_MAE_RGE_df = melt(t(top_MAE_RGE_mat))
MAE_O_levels = colnames(MAE_O_mat)
MAE_DID_levels = colnames(MAE_DID_mat)
MAE_RGE_levels = colnames(MAE_RGE_mat)
top_MAE_O_levels = colnames(top_MAE_O_mat)
top_MAE_DID_levels = colnames(top_MAE_DID_mat)
top_MAE_RGE_levels = colnames(top_MAE_RGE_mat)
### MAE - Overall
p = ggplot(MAE_O_df, aes(x = factor(Var1, levels = MAE_O_levels), y = value, fill = factor(Var1, levels = MAE_O_levels))) +
  geom_boxplot(outlier.shape = NA) + stat_boxplot(geom = "errorbar", width = 0.3) +
  ylim(min(apply(MAE_O_mat, 2, quantile, 0.1)), max(apply(MAE_O_mat, 2, quantile, 0.9))) + theme_classic() +
  labs(x="", y="MAE") + scale_fill_manual(values=method_color) +
  theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position = "none")
ggsave(paste0(output_dir, "/MAE_Overall.pdf"), p, height = 6, width = 5)
### MAE - Dropouts Introduced by Down-sampling
p = ggplot(MAE_DID_df, aes(x = factor(Var1, levels = MAE_DID_levels), y = value, fill = factor(Var1, levels = MAE_DID_levels))) +
  geom_boxplot(outlier.shape = NA) + stat_boxplot(geom = "errorbar", width = 0.3) +
  ylim(min(apply(MAE_DID_mat, 2, quantile, 0.1)), max(apply(MAE_DID_mat, 2, quantile, 0.9))) + theme_classic() +
  labs(x="", y="MAE") + scale_fill_manual(values=method_color) +
  theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position = "none")
ggsave(paste0(output_dir, "/MAE_Sampled_Zero.pdf"), p, height = 6, width = 5)
### MAE - Retained Gene Expressions
p = ggplot(MAE_RGE_df, aes(x = factor(Var1, levels = MAE_RGE_levels), y = value, fill = factor(Var1, levels = MAE_RGE_levels))) +
  geom_boxplot(outlier.shape = NA) + stat_boxplot(geom = "errorbar", width = 0.3) +
  ylim(min(apply(MAE_RGE_mat, 2, quantile, 0.1)), max(apply(MAE_RGE_mat, 2, quantile, 0.9))) + theme_classic() +
  labs(x="", y="MAE") + scale_fill_manual(values=method_color) +
  theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position = "none")
ggsave(paste0(output_dir, "/MAE_Retained_Expressions.pdf"), p, height = 6, width = 5)
### top MAE - Overall
p = ggplot(top_MAE_O_df, aes(x = factor(Var1, levels = top_MAE_O_levels), y = value, fill = factor(Var1, levels = top_MAE_O_levels))) +
  geom_boxplot(outlier.shape = NA) + stat_boxplot(geom = "errorbar", width = 0.3) +
  ylim(min(apply(top_MAE_O_mat, 2, quantile, 0.1)), max(apply(top_MAE_O_mat, 2, quantile, 0.9))) + theme_classic() +
  labs(x="", y="MAE") + scale_fill_manual(values=method_color) +
  theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position = "none")
ggsave(paste0(output_dir, "/top_MAE_Overall.pdf"), p, height = 6, width = 5)
### top MAE - Dropouts Introduced by Down-sampling
p = ggplot(top_MAE_DID_df, aes(x = factor(Var1, levels = top_MAE_DID_levels), y = value, fill = factor(Var1, levels = top_MAE_DID_levels))) +
  geom_boxplot(outlier.shape = NA) + stat_boxplot(geom = "errorbar", width = 0.3) +
  ylim(min(apply(top_MAE_DID_mat, 2, quantile, 0.1)), max(apply(top_MAE_DID_mat, 2, quantile, 0.9))) + theme_classic() +
  labs(x="", y="MAE") + scale_fill_manual(values=method_color) +
  theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position = "none")
ggsave(paste0(output_dir, "/top_MAE_Sampled_Zero.pdf"), p, height = 6, width = 5)
### top MAE - Retained Gene Expressions
p = ggplot(top_MAE_RGE_df, aes(x = factor(Var1, levels = top_MAE_RGE_levels), y = value, fill = factor(Var1, levels = top_MAE_RGE_levels))) +
  geom_boxplot(outlier.shape = NA) + stat_boxplot(geom = "errorbar", width = 0.3) +
  ylim(min(apply(top_MAE_RGE_mat, 2, quantile, 0.1)), max(apply(top_MAE_RGE_mat, 2, quantile, 0.9))) + theme_classic() +
  labs(x="", y="MAE") + scale_fill_manual(values=method_color) +
  theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position = "none")
ggsave(paste0(output_dir, "/top_MAE_Retained_Expressions.pdf"), p, height = 6, width = 5)
### CMD
top_300_genes = hvg_genes[1:300]
cor_all = list()
for(ii in names(data_list)){
  if(ii == "Raw"){
    cor_all[[ii]] = calc_cor_mat(data_list[[ii]][top_300_genes, ])
  }else{
    cor_all[[ii]] = list()
    for(jj in repeats){
      cor_all[[ii]][[jj]] = calc_cor_mat(delete_lt0.5(data_list[[ii]][[jj]])[top_300_genes, ])
    }
  }
  print(ii)
}
saveRDS(cor_all, paste0(output_dir, "/cor_all.rds"))
cmd_mat = matrix(nrow = length(method_name), ncol = length(repeats), dimnames = list(method_name, repeats))
for(ii in method_name){
  for(jj in repeats){
    cmd_mat[ii, jj] = calc_cmd(cor_all[["Raw"]], cor_all[[ii]][[jj]])
  }
}
pdf(paste0(output_dir, "/CMD.pdf"), height = 6, width = 5)
barplot_usage(rowMeans(cmd_mat), standard_error = apply(cmd_mat, 1, ste), main = "", cex.main = 1.5, bar_color = method_color, text_color = text_color, ylab = "CMD", cex.lab = 1.5, font.main = 1, ylim = c(-0.1, 1), use_border = F)
dev.off()
result_list[["CMD"]] = cmd_mat
### Gene correlation
gene_corr_mat = matrix(nrow = gene_number, ncol = length(method_name), dimnames = list(c(), method_name))
for(ii in method_name){
  for(jj in repeats){
    if(jj == repeats[1]){
      cor_mat = calc_corr(data_list[["Raw"]], data_list[[ii]][[jj]], "gene")
    }else{
      cor_mat = cbind(cor_mat, calc_corr(data_list[["Raw"]], data_list[[ii]][[jj]], "gene"))
    }
  }
  gene_corr_mat[, ii] = rowMeans(cor_mat)
  print(ii)
}
gene_corr_mat = gene_corr_mat[rowSums(is.na(gene_corr_mat)) < 1, ]
gene_corr_df = melt(t(gene_corr_mat))
gene_corr_levels = colnames(gene_corr_mat)
p = ggplot(gene_corr_df, aes(x = factor(Var1, levels = gene_corr_levels), y = value, fill = factor(Var1, levels = gene_corr_levels))) +
  geom_boxplot(outlier.shape = NA) + stat_boxplot(geom = "errorbar", width = 0.3) +
  ylim(0, 1) + theme_classic() +
  labs(x="", y="Gene correlation with reference") + scale_fill_manual(values=method_color) +
  theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position = "none")
ggsave(paste0(output_dir, "/CORR_GENE.pdf"), p, height = 6, width = 5)
result_list[["CORR_GENE"]] = gene_corr_mat
### Gene correlation top 1000 genes
top_1000_genes = hvg_genes[1:1000]
gene_corr_mat = matrix(nrow = length(top_1000_genes), ncol = length(method_name), dimnames = list(c(), method_name))
for(ii in method_name){
  for(jj in repeats){
    if(jj == repeats[1]){
      cor_mat = calc_corr(data_list[["Raw"]][top_1000_genes, ], data_list[[ii]][[jj]][top_1000_genes, ], "gene")
    }else{
      cor_mat = cbind(cor_mat, calc_corr(data_list[["Raw"]][top_1000_genes, ], data_list[[ii]][[jj]][top_1000_genes, ], "gene"))
    }
  }
  gene_corr_mat[, ii] = rowMeans(cor_mat)
  print(ii)
}
gene_corr_mat = gene_corr_mat[rowSums(is.na(gene_corr_mat)) < 1, ]
gene_corr_df = melt(t(gene_corr_mat))
gene_corr_levels = colnames(gene_corr_mat)
p = ggplot(gene_corr_df, aes(x = factor(Var1, levels = gene_corr_levels), y = value, fill = factor(Var1, levels = gene_corr_levels))) +
  geom_boxplot(outlier.shape = NA) + stat_boxplot(geom = "errorbar", width = 0.3) +
  ylim(0, 1) + theme_classic() +
  labs(x="", y="Gene correlation with reference") + scale_fill_manual(values=method_color) +
  theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position = "none")
ggsave(paste0(output_dir, "/CORR_GENE_top_1000_genes.pdf"), p, height = 6, width = 5)
result_list[["CORR_GENE_top_1000_genes"]] = gene_corr_mat
### Cell correlation
cell_corr_mat = matrix(nrow = cell_number, ncol = length(method_name), dimnames = list(c(), method_name))
for(ii in method_name){
  for(jj in repeats){
    if(jj == repeats[1]){
      cor_mat = calc_corr(data_list[["Raw"]], data_list[[ii]][[jj]], "cell")
    }else{
      cor_mat = cbind(cor_mat, calc_corr(data_list[["Raw"]], data_list[[ii]][[jj]], "cell"))
    }
  }
  cell_corr_mat[, ii] = rowMeans(cor_mat)
  print(ii)
}
cell_corr_mat = cell_corr_mat[rowSums(is.na(cell_corr_mat)) < 1, ]
cell_corr_df = melt(t(cell_corr_mat))
cell_corr_levels = colnames(cell_corr_mat)
p = ggplot(cell_corr_df, aes(x = factor(Var1, levels = cell_corr_levels), y = value, fill = factor(Var1, levels = cell_corr_levels))) +
  geom_boxplot(outlier.shape = NA) + stat_boxplot(geom = "errorbar", width = 0.3) +
  ylim(0, 1) + theme_classic() +
  labs(x="", y="Cell correlation with reference") + scale_fill_manual(values=method_color) +
  theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position = "none")
ggsave(paste0(output_dir, "/CORR_CELL.pdf"), p, height = 6, width = 5)
result_list[["CORR_CELL"]] = cell_corr_mat
### Cell correlation top 1000 genes
cell_corr_mat = matrix(nrow = cell_number, ncol = length(method_name), dimnames = list(c(), method_name))
for(ii in method_name){
  for(jj in repeats){
    if(jj == repeats[1]){
      cor_mat = calc_corr(data_list[["Raw"]][top_1000_genes, ], data_list[[ii]][[jj]][top_1000_genes, ], "cell")
    }else{
      cor_mat = cbind(cor_mat, calc_corr(data_list[["Raw"]][top_1000_genes, ], data_list[[ii]][[jj]][top_1000_genes, ], "cell"))
    }
  }
  cell_corr_mat[, ii] = rowMeans(cor_mat)
  print(ii)
}
cell_corr_mat = cell_corr_mat[rowSums(is.na(cell_corr_mat)) < 1, ]
cell_corr_df = melt(t(cell_corr_mat))
cell_corr_levels = colnames(cell_corr_mat)
p = ggplot(cell_corr_df, aes(x = factor(Var1, levels = cell_corr_levels), y = value, fill = factor(Var1, levels = cell_corr_levels))) +
  geom_boxplot(outlier.shape = NA) + stat_boxplot(geom = "errorbar", width = 0.3) +
  ylim(0, 1) + theme_classic() +
  labs(x="", y="Cell correlation with reference") + scale_fill_manual(values=method_color) +
  theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position = "none")
ggsave(paste0(output_dir, "/CORR_CELL_top_1000_genes.pdf"), p, height = 6, width = 5)
result_list[["CORR_CELL_top_1000_genes"]] = cell_corr_mat
saveRDS(result_list, paste0(output_dir, "/result.rds"))








