utilities_path = "/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/utilities.r"
source(utilities_path)
result_list = list()
method_names = c("Observed", "DISC", "SAVER", "scImpute", "VIPER", "MAGIC", "DCA", "deepImpute", "scScope", "scVI")
text_color = rep("black", length(method_names))
names(text_color) = method_names
text_color["DISC"] = "red"
bar_color = rep("gray50", length(method_names))
names(bar_color) = method_names
bar_color["Observed"] = "gray80"
bar_color["DISC"] = "red"
for(ii in method_names){
  result_list[[ii]] = list()
  for(jj in c(1, 2, 3, 4, 5)){
    if(ii %in% c("Observed", "DISC", "SAVER", "MAGIC", "DCA", "scScope", "scVI")){
      if(ii %in% c("SAVER", "MAGIC", "DCA", "scScope", "scVI")){
        use_name = ii
      }else if(ii == "Observed"){
        use_name = "RAW"
      }else if(ii == "DISC"){
        use_name = "DeSCI"
      }else{
        stop("Invaild use_name!")
      }
      cell_type_result = as.matrix(read.table(paste0("/home/wucheng/imputation/pbmc/ds0.3/repeat", jj, "/", use_name, "/Recall1.txt"), sep = "\t"))
      overall_result = as.matrix(read.table(paste0("/home/wucheng/imputation/pbmc/ds0.3/repeat", jj, "/", use_name, "/Accuracy.txt"), sep = "\t", header = T))
      this_result = list(Jaccard = cell_type_result["Jaac", ], ACC = overall_result[, "Accuracy"])
      result_list[[ii]][[jj]] = this_result
    }else{
      if(ii == "VIPER"){
        use_name = "VIPER_gene"
      }else{
        use_name = ii
      }
      all_result = readRDS(paste0("/home/yuanhao/data/fn/pbmc3k/ds/downsampling_first_repeat_", jj, "/imputation/cluster_evaluation_0.3/3_", use_name, "_mc_10_mce_1_resume_dim/pca/identification_result.rds"))
      this_result = list(Jaccard = all_result[["cell_type_result"]][, "Jaccard"], ACC = all_result[["summary"]]["ACC"])
      result_list[[ii]][[jj]] = this_result
    }
  }
}
output_dir = paste0("/home/yuanhao/data/fn/pbmc3k/ds/ds_0.3_results/clustering")
dir.create(output_dir, showWarnings = F, recursive = T)


acc_result = t(sapply(result_list, function(x){
  return(as.numeric(sapply(x, function(y){
    return(y[["ACC"]])
  })))
}))
jaccard_mat = t(sapply(result_list, function(x){
  return(rowMeans(sapply(x, function(y){
    return(y[["Jaccard"]])
  })))
}))

pdf(paste0(output_dir, "/ACC.pdf"), height = 6, width = 5)
barplot_usage(rowMeans(acc_result), standard_error = apply(acc_result, 1, ste), main = "", cex.main = 1.5, bar_color = bar_color, text_color = text_color, use_data_order = T, ylab = "ACC", cex.lab = 1.5, font.main = 1, ylim = c(-0.1, 1), decreasing = TRUE)
dev.off()


pdf(paste0(output_dir, "/Jaccard.pdf"), height = 5, width = 10)
cell_type_heatmap(jaccard_mat, "Jaccard")
dev.off()











