setwd("/home/yuanhao/github_repositories/DISC/reproducibility")
utilities_path = "./source/utilities.r"
source(utilities_path)
used_datasets = c("MELANOMA", "SSCORTEX", "CBMC", "PBMC", "RETINA", "BRAIN_SPLiT")
method_names = c("DISC", "scScope")
output_dir = "./results/Identification_of_true_zeros"
dir.create(output_dir, showWarnings = F, recursive = T)
data_list = list()
for(ii in used_datasets){
  data_list[[ii]] = list()
  raw_data = readh5_loom(paste0("./data/", ii, "/raw.loom"))
  vst_file = paste0(output_dir, "/", ii, "_vst_gene.tsv")
  if(file.exists(vst_file)){
    hvg_info = read.table(vst_file)
    print("load vst_file")
  }else{
    hvg_info = FindVariableFeatures_vst_by_genes(raw_data)
    hvg_info = hvg_info[order(hvg_info$variance.standardized, decreasing = T), ]
    write.table(hvg_info, vst_file, sep = "\t", quote = F, row.names = T, col.names = T)
  }
  ds_dir = paste0("./data/", ii, "/ds_0.5/r1")
  observed_data = readh5_loom(paste0(ds_dir, "/gene_selection.loom"))
  compared_genes = rownames(observed_data)
  top_1000_genes = rownames(hvg_info)[rownames(hvg_info) %in% compared_genes][1:1000]
  raw_data = raw_data[compared_genes, ]
  for(jj in method_names){
    if(jj == "DISC"){
      imputation = readh5_loom(paste0(ds_dir, "/", jj, ".loom"))[compared_genes, ]
    }else{
      imputation = readh5_imputation(paste0(ds_dir, "/", jj, ".hdf5"))[compared_genes, ]
    }
    imputation_cp10k = sweep(imputation, 2, 10000 / colSums(imputation), "*")
    raw_zero = raw_data == 0 & observed_data == 0
    induced_zero = raw_data != 0 & observed_data == 0
    kept_values = raw_data != 0 & observed_data != 0
    data_list[[ii]][[jj]] = list(raw_zero_imputation = imputation_cp10k[raw_zero],
                                 induced_zero_imputation = imputation_cp10k[induced_zero],
                                 top_raw_zero_imputation = imputation_cp10k[top_1000_genes, ][raw_zero[top_1000_genes, ]],
                                 top_induced_zero_imputation = imputation_cp10k[top_1000_genes, ][induced_zero[top_1000_genes, ]])
    cat("Finish: ", ii, " - ", jj, "\n")
  }
  cat("Finish: ", ii, "\n")
}
#  All genes
pdf(file = paste0(output_dir, "/all_genes.pdf"), height = 5, width = 2 * length(used_datasets))
par(mfrow = c(1, length(used_datasets)), mar = c(5, 0, 0.5, 0), oma = c(0, 5, 6, 2), mgp = c(3.5, 1, 0), cex.axis = 1.1, cex.lab = 1.5, font.lab = 2, cex.main = 1.5)
mtext_position = seq(from = 0, to = 1, by = 1 / length(used_datasets) / 2)[seq(length(used_datasets)) * 2]
names(mtext_position) = used_datasets
for(ii in used_datasets){
  max_value = 0
  for(jj in method_names){
    max_value = max(c(max_value, data_list[[ii]][[jj]][["raw_zero_imputation"]], data_list[[ii]][[jj]][["induced_zero_imputation"]]))
  }
  n = 2 ^ (round(log(max_value, base = 2)) + 7)
  bw = (max_value * (n + 10)) / (n ^ 2)
  to = max_value + (bw * 5)
  from = -(bw * 5)
  x_max = 0
  for(jj in method_names){
    data_list[[ii]][[jj]][["raw_zero_density"]] = density(data_list[[ii]][[jj]][["raw_zero_imputation"]], from = from, to = to, n = n, bw = bw)
    data_list[[ii]][[jj]][["induced_zero_density"]] = density(data_list[[ii]][[jj]][["induced_zero_imputation"]], from = from, to = to, n = n, bw = bw)
    data_list[[ii]][[jj]][["y_max"]] = max(c(data_list[[ii]][[jj]][["raw_zero_density"]]$y, data_list[[ii]][[jj]][["induced_zero_density"]]$y))
    cutoff = data_list[[ii]][[jj]][["y_max"]] * 0.005
    test_vector = !(data_list[[ii]][[jj]][["raw_zero_density"]]$y < cutoff & data_list[[ii]][[jj]][["induced_zero_density"]]$y < cutoff)
    for(kk in seq(from = length(test_vector), to = 1)){
      if(test_vector[kk]){
        x_max = max(c(x_max, data_list[[ii]][[jj]][["raw_zero_density"]]$x[kk + 1]))
        break()
      }
    }
  }
  data_list[[ii]][["x_max"]] = x_max
  cat("Finish: ", ii, "\n")
}
for(ii in method_names){
  for(jj in used_datasets){
    plot(0, type = "n", xlim = c(-min(c(1, data_list[[jj]][["x_max"]] * 0.05)), data_list[[jj]][["x_max"]]), ylim = c(0, data_list[[jj]][[ii]][["y_max"]]), axes = FALSE, ann = FALSE, frame.plot = TRUE)
    lines(data_list[[jj]][[ii]][["raw_zero_density"]], col = "black", lwd = 2)
    lines(data_list[[jj]][[ii]][["induced_zero_density"]], col = "red", lwd = 2)
    axis(1)
    mtext(jj, outer = TRUE, cex = 1.2, font = 2, line = 0.5, at = mtext_position[jj])
    cat("Finish: ", ii, " - ", jj, "\n")
  }
  mtext(paste0(ii, " - All genes"), outer = TRUE, cex = 1.4, font = 2, line = 3)
  legend("topright", c("Raw zero", "Sampled zero"), lty = 1, col = c("black", "red"), lwd = 2)
}
dev.off()

#  Top 1000 genes
pdf(file = paste0(output_dir, "/top_1000_genes.pdf"), height = 5, width = 2 * length(used_datasets))
par(mfrow = c(1, length(used_datasets)), mar = c(5, 0, 0.5, 0), oma = c(0, 5, 6, 2), mgp = c(3.5, 1, 0), cex.axis = 1.1, cex.lab = 1.5, font.lab = 2, cex.main = 1.5)
mtext_position = seq(from = 0, to = 1, by = 1 / length(used_datasets) / 2)[seq(length(used_datasets)) * 2]
names(mtext_position) = used_datasets
for(ii in used_datasets){
  max_value = 0
  for(jj in method_names){
    max_value = max(c(max_value, data_list[[ii]][[jj]][["top_raw_zero_imputation"]], data_list[[ii]][[jj]][["top_induced_zero_imputation"]]))
  }
  n = 2 ^ (round(log(max_value, base = 2)) + 7)
  bw = (max_value * (n + 10)) / (n ^ 2)
  to = max_value + (bw * 5)
  from = -(bw * 5)
  x_max = 0
  for(jj in method_names){
    data_list[[ii]][[jj]][["top_raw_zero_density"]] = density(data_list[[ii]][[jj]][["top_raw_zero_imputation"]], from = from, to = to, n = n, bw = bw)
    data_list[[ii]][[jj]][["top_induced_zero_density"]] = density(data_list[[ii]][[jj]][["top_induced_zero_imputation"]], from = from, to = to, n = n, bw = bw)
    data_list[[ii]][[jj]][["top_y_max"]] = max(c(data_list[[ii]][[jj]][["top_raw_zero_density"]]$y, data_list[[ii]][[jj]][["top_induced_zero_density"]]$y))
    cutoff = data_list[[ii]][[jj]][["top_y_max"]] * 0.005
    test_vector = !(data_list[[ii]][[jj]][["top_raw_zero_density"]]$y < cutoff & data_list[[ii]][[jj]][["top_induced_zero_density"]]$y < cutoff)
    for(kk in seq(from = length(test_vector), to = 1)){
      if(test_vector[kk]){
        x_max = max(c(x_max, data_list[[ii]][[jj]][["top_raw_zero_density"]]$x[kk + 1]))
        break()
      }
    }
  }
  data_list[[ii]][["top_x_max"]] = x_max
  cat("Finish: ", ii, "\n")
}
for(ii in method_names){
  for(jj in used_datasets){
    plot(0, type = "n", xlim = c(-min(c(1, data_list[[jj]][["top_x_max"]] * 0.05)), data_list[[jj]][["top_x_max"]]), ylim = c(0, data_list[[jj]][[ii]][["top_y_max"]]), axes = FALSE, ann = FALSE, frame.plot = TRUE)
    lines(data_list[[jj]][[ii]][["top_raw_zero_density"]], col = "black", lwd = 2)
    lines(data_list[[jj]][[ii]][["top_induced_zero_density"]], col = "red", lwd = 2)
    axis(1)
    mtext(jj, outer = TRUE, cex = 1.2, font = 2, line = 0.5, at = mtext_position[jj])
    cat("Finish: ", ii, " - ", jj, "\n")
  }
  mtext(paste0(ii, " - Top 1000 genes"), outer = TRUE, cex = 1.4, font = 2, line = 3)
  legend("topright", c("Raw zero", "Sampled zero"), lty = 1, col = c("black", "red"), lwd = 2)
}
dev.off()





