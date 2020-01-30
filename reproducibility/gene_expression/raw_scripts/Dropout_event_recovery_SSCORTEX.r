utilities_path = "/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/utilities.r"
source(utilities_path)
min_expressed_cell = 10
min_expressed_cell_average_expression = 1
ds_mode = 0.5
method_names = c("DISC", "SAVER", "scImpute", "VIPER", "MAGIC", "DCA", "deepImpute", "scScope", "scVI")
### Load raw data and downsampling
raw_data = readh5_loom("/home/yuanhao/github_repositories/DISC/reproducibility/data/SSCORTEX/raw.loom")
compare_gene = rownames(raw_data)
cell_number = ncol(raw_data)
data_list = list(Raw = raw_data)
rm(raw_data)
ds_dir = "/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/ds"
dir.create(ds_dir, showWarnings = F, recursive = T)
output_dir = paste0(ds_dir, "/ds_", ds_mode, "_results")
dir.create(output_dir, showWarnings = F, recursive = T)
repeats = paste0("downsampling_first_repeat_", seq(5))
observed_file = paste0("L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_ds_", ds_mode, ".loom")
observed_name = paste(delete_last_element(unlist(strsplit(observed_file, ".", fixed = T))), collapse = ".")
data_list[["Observed"]] = list()
for(ii in repeats){
  observed_path = paste(ds_dir, ii, observed_file, sep = "/")
  if(!file.exists(observed_path)){
    cat("Generate ", observed_path, " ...")
    data_list[["Observed"]][[ii]] = downsampling_cell(ds_mode, raw_data)
    save_h5(observed_path, t(data_list[["Observed"]][[ii]]))
  }else{
    data_list[["Observed"]][[ii]] = readh5_loom(observed_path)
  }
  expressed_cell = rowSums(data_list[["Observed"]][[ii]] > 0)
  gene_expression = rowSums(data_list[["Observed"]][[ii]])
  gene_filter = expressed_cell >= min_expressed_cell & gene_expression > expressed_cell * min_expressed_cell_average_expression
  this_used_genes = rownames(data_list[["Observed"]][[ii]])[gene_filter]
  compare_gene = intersect(compare_gene, this_used_genes)
}
gene_number = length(compare_gene)
cat("Use ", gene_number, " genes for comparison.\n")
data_list[["Raw"]] = data_list[["Raw"]][compare_gene, ]
for(ii in method_names){
  data_list[[ii]] = list()
}
for(ii in repeats){
  data_list[["Observed"]][[ii]] = data_list[["Observed"]][[ii]][compare_gene, ]
  ### Load downsampling data and imputation results
  data_list[["DISC"]][[ii]] = readh5_imputation(get_optimal_point33(paste(ds_dir, ii, paste0("DeSCI_2.7.4.33/ds_", ds_mode, "/log.txt"), sep = "/")), with_outliers = T)[compare_gene, ]
  print(dim(data_list[["DISC"]][[ii]]))
  for(jj in setdiff(method_names, "DISC")){
    if(jj == "VIPER"){
      tmp_name = "VIPER_gene"
    }else{
      tmp_name = jj
    }
    data_list[[jj]][[ii]] = readh5_imputation(paste0(paste(ds_dir, ii, "imputation", paste(observed_name, tmp_name, "mc", min_expressed_cell, "mce", min_expressed_cell_average_expression, sep = "_"), sep = "/"), ".hdf5"))[compare_gene, ]
    print(dim(data_list[[jj]][[ii]]))
  }
}
### Settings
method_names = setdiff(names(data_list), "Raw")
method_color = c("gray80", "#FF0000", "#000080", "#BFBF00", "#408000", "#804000", "#00FF00", "#FF8000", "#FF00FF", "#00FFFF")
names(method_color) = method_names
text_color = rep("black", length(method_names))
names(text_color) = method_names
text_color["DISC"] = "red"
bar_color = rep("gray50", length(method_names))
names(bar_color) = method_names
bar_color["Observed"] = "gray80"
bar_color["DISC"] = "red"
### MAE
mae_eq0 = matrix(nrow = length(method_names), ncol = length(repeats), dimnames = list(method_names, repeats))
scale_factor = 1 / ds_mode
for(ii in method_names){
  for(jj in repeats){
    eq0 = sum(data_list[["Raw"]] > 0 & data_list[["Observed"]][[jj]] == 0)
    sae_eq0 = sum(sapply(compare_gene, function(x){
      expressed_mask = data_list[["Raw"]][x, ] > 0 & data_list[["Observed"]][[jj]][x, ] == 0
      return(sum(abs(data_list[["Raw"]][x, expressed_mask] - (data_list[[ii]][[jj]][x, expressed_mask] * scale_factor))))
    }))
    mae_eq0[ii, jj] = sae_eq0 / eq0
  }
  print(ii)
}
pdf(paste0(output_dir, "/MAE.pdf"), height = 6, width = 5)
barplot_usage(rowMeans(mae_eq0), standard_error = apply(mae_eq0, 1, ste), main = "Zero entries", cex.main = 1.5, bar_color = bar_color, text_color = text_color, use_data_order = T, ylab = "Log (MAE + 1)", use_log1p = T, cex.lab = 1.5, font.main = 1)
dev.off()
### CMD
CMD_output_dir = paste0(output_dir, "/CMD")
dir.create(CMD_output_dir, showWarnings = F, recursive = T)
vst_file = paste0(CMD_output_dir, "/vst_gene.tsv")
if(file.exists(vst_file)){
  hvg_info = read.table(vst_file)
  print("load vst_file")
}else{
  hvg_info = FindVariableFeatures_vst_by_genes(data_list[["Raw"]])
  hvg_info = hvg_info[order(hvg_info$variance.standardized, decreasing = T), ]
  write.table(hvg_info, paste0(CMD_output_dir, "/vst_gene.tsv"), sep = "\t", quote = F, row.names = T, col.names = T)
}
used_feature_genes = rownames(hvg_info)[1:300]
cor_all = list()
for(ii in names(data_list)){
  if(ii == "Raw"){
    cor_all[[ii]] = calc_cor_mat(data_list[[ii]][used_feature_genes, ])
  }else{
    cor_all[[ii]] = list()
    for(jj in repeats){
      cor_all[[ii]][[jj]] = calc_cor_mat(delete_lt0.5(data_list[[ii]][[jj]])[used_feature_genes, ])
    }
  }
  print(ii)
}
saveRDS(cor_all, paste0(CMD_output_dir, "/cor_all.rds"))
cmd_mat = matrix(nrow = length(method_names), ncol = length(repeats), dimnames = list(method_names, repeats))
for(ii in method_names){
  for(jj in repeats){
    cmd_mat[ii, jj] = calc_cmd(cor_all[["Raw"]], cor_all[[ii]][[jj]])
  }
}
pdf(paste0(output_dir, "/CMD.pdf"), height = 6, width = 5)
barplot_usage(rowMeans(cmd_mat), standard_error = apply(cmd_mat, 1, ste), main = "", cex.main = 1.5, bar_color = bar_color, text_color = text_color, use_data_order = T, ylab = "CMD", cex.lab = 1.5, font.main = 1, ylim = c(-0.1, 1))
dev.off()
for(ii in method_names){
  for(jj in repeats){
    dimnames(data_list[[ii]][[jj]]) = dimnames(data_list[["Raw"]])
  }
  print(ii)
}
### Gene correlation
gene_corr_mat = matrix(nrow = length(method_names), ncol = length(repeats), dimnames = list(method_names, repeats))
for(ii in method_names){
  for(jj in repeats){
    gene_corr_mat[ii, jj] = mean(calc_corr(data_list[["Raw"]], data_list[[ii]][[jj]], "gene"), na.rm = T)
  }
  print(ii)
}
pdf(paste0(output_dir, "/CORR_GENE.pdf"), height = 6, width = 5)
barplot_usage(rowMeans(gene_corr_mat), standard_error = apply(gene_corr_mat, 1, ste), main = "", cex.main = 1.5, bar_color = bar_color, text_color = text_color, use_data_order = T, decreasing = T, ylab = "Gene correlation with reference", cex.lab = 1.5, font.main = 1, ylim = c(-0.1, 1))
dev.off()
### Cell correlation
cell_corr_mat = matrix(nrow = length(method_names), ncol = length(repeats), dimnames = list(method_names, repeats))
for(ii in method_names){
  for(jj in repeats){
    cell_corr_mat[ii, jj] = mean(calc_corr(data_list[["Raw"]], data_list[[ii]][[jj]], "cell"), na.rm = T)
  }
  print(ii)
}
pdf(paste0(output_dir, "/CORR_CELL.pdf"), height = 6, width = 5)
barplot_usage(rowMeans(cell_corr_mat), standard_error = apply(cell_corr_mat, 1, ste), main = "", cex.main = 1.5, bar_color = bar_color, text_color = text_color, use_data_order = T, decreasing = T, ylab = "Gene correlation with reference", cex.lab = 1.5, font.main = 1, ylim = c(-0.1, 1))
dev.off()








