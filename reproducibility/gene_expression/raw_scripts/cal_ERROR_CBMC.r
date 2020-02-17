utilities_path = "/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/utilities.r"
source(utilities_path)
min_expressed_cell = 10
min_expressed_cell_average_expression = 1
ds_mode = 0.5
method_names = c("DISC", "scImpute", "VIPER", "MAGIC", "DCA", "deepImpute", "scScope", "scVI")
### Load raw data and downsampling
raw_data = readh5_loom("/home/yuanhao/github_repositories/DISC/reproducibility/data/CBMC/raw.loom")
compare_gene = rownames(raw_data)
cell_number = ncol(raw_data)
data_list = list(Raw = raw_data)
rm(raw_data)
ds_dir = "/home/yuanhao/data/fn/CITE-seq"
dir.create(ds_dir, showWarnings = F, recursive = T)
output_dir = paste0(ds_dir, "/ds_", ds_mode, "_results_new_strategy3")
dir.create(output_dir, showWarnings = F, recursive = T)
repeats = paste0("downsampling_first_repeat_", seq(5))
observed_file = paste0("GSE100866_CBMC_8K_filtered_ds_", ds_mode, ".loom")
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

MAE_mat = matrix(nrow = cell_number, ncol = length(method_names), dimnames = list(c(), method_names))

ls_raw = colSums(data_list[["Raw"]])
for(ii in method_names){
  for(jj in repeats){
    ls_this = colSums(data_list[[ii]][[jj]])
    scale_factor = ls_raw / ls_this
    ERROR_MAE = sapply(seq(cell_number), function(x){
      expressed_mask = data_list[["Raw"]][, x] > 0
      expressed_number = sum(expressed_mask)
      error = data_list[["Raw"]][expressed_mask, x] - (data_list[[ii]][[jj]][expressed_mask, x] * scale_factor[x])
      return(sum(abs(error)) / expressed_number)
    })
    if(jj == repeats[1]){
      MAE_cell = ERROR_MAE
    }else{
      MAE_cell = cbind(MAE_cell, ERROR_MAE)
    }
  }
  MAE_mat[, ii] = rowMeans(MAE_cell)
  print(ii)
}
MAE_mat = MAE_mat[rowSums(is.na(MAE_mat)) < 1, ]

boxplot_usage = function(data_matrix, main, bar_color, text_color=NULL, use_data_order=F, decreasing=F, cex.main=2, axis_by=0.25, ...){
  data_means = colMeans(data_matrix, na.rm = T)
  data_order = c(1, order(data_means[-1], decreasing = decreasing) + 1)
  if(use_data_order){
    data_matrix = data_matrix[, data_order]
    bar_color = bar_color[colnames(data_matrix)][data_order]
    if(!is.null(text_color)){
      text_color = text_color[colnames(data_matrix)][data_order]
    }
  }
  bp = boxplot(data_matrix, main = main, at = seq(ncol(data_matrix)), names = rep("", ncol(data_matrix)), las = 2, names.arg="", col = bar_color, cex.axis = 1.2, cex.main = cex.main, outline = F, ...)
  if(is.null(text_color)){
    text_color = rep("black", ncol(data_matrix))
  }
  for(ii in seq(ncol(data_matrix))){
    mtext(colnames(data_matrix)[ii], side = 1, line = -0.25, at = ii, las = 2, font = 1, col = text_color[ii], cex=1)
  }
}

MAE_df = melt(t(MAE_mat))
MAE_levels = colnames(MAE_mat)[c(1, order(colMeans(MAE_mat)[-1], decreasing = F) + 1)]

p = ggplot(MAE_df, aes(x = factor(Var1, levels = MAE_levels), y = value, fill = factor(Var1, levels = MAE_levels))) +
  geom_boxplot(outlier.shape = NA) + stat_boxplot(geom = "errorbar", width = 0.3) +
  ylim(min(apply(MAE_mat, 2, quantile, 0.1)), max(apply(MAE_mat, 2, quantile, 0.9))) + theme_classic() + 
  theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"),
        legend.position = "none")
ggsave(paste0(output_dir, "/ERROR_cell.pdf"), p, height = 6, width = 5)
saveRDS(MAE_mat, paste0(output_dir, "/MAE_mat.rds"))



