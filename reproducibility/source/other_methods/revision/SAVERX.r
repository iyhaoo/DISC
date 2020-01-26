library(devtools)
install_github("jingshuw/SAVERX")



library(SAVERX)
file <- saverx("D:/study/single_cell/saverx/shekhar_downsampled.csv")
denoised.data <- readRDS(file)


file <- saverx("./testdata/shekhar_downsampled.csv", data.species = "Mouse", use.pretrain = T, pretrained.weights.file = "./mouse_retina.hdf5", model.species = "Mouse")
denoised.data <- readRDS(file)











generic_functions_path = "/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/generic_functions.r"
source(generic_functions_path)
####################main###############
dataset_list = list()
dataset_list[["FISH"]] = readh5_loom("/home/yuanhao/data/fn/sscortex/osmFISH_SScortex_mouse_sscortex_cells.loom")
raw_input_data = "/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/imputation/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_mc_10_mce_1.loom"
use_genes = intersect(rownames(dataset_list[["FISH"]]), get_loom_gene(raw_input_data))
dataset_list[["Raw"]] = readh5_loom(raw_input_data, use_genes)
used_cells = colnames(dataset_list[["Raw"]])
### our imputation
our_result = get_optimal_point1("/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/DeSCI_2.7.4.17/log.txt")
dataset_list[["DeSCI"]] = readh5_imputation(our_result, use_genes, with_outliers=T)
### theirs
dataset_list[["SAVER"]] = readRDS("/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/imputation/SAVER_tmp/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_SAVER_mc_10_mce_1.rds")
dataset_list[["SAVER_gamma"]] = gamma_result(dataset_list[["SAVER"]], num_of_obs=1)[use_genes, used_cells]
dataset_list[["MAGIC"]] = readh5_imputation("/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/imputation/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_MAGIC_mc_10_mce_1.hdf5", use_genes, used_cells)
dataset_list[["DCA"]] = readh5_imputation("/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/imputation/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_DCA_mc_10_mce_1.hdf5", use_genes, used_cells)
dataset_list[["scScope"]] = readh5_imputation("/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/imputation/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_scScope_mc_10_mce_1.hdf5", use_genes, used_cells)
dataset_list[["scVI"]] = readh5_imputation("/home/yuanhao/data/fn/sscortex/filt_gene_500_5000/merge/imputation/L1_Cortex2_filt_ls_merged_s1_s2_unique_rename_scVI_mc_10_mce_1.hdf5", use_genes, used_cells)
### loaded
method_names = setdiff(names(dataset_list), c("FISH", "SAVER_gamma"))
method_color = c("gray80", "red", "blue4", "yellow4", "green", "purple", "cyan")
names(method_color) = method_names
bar_color = rep("gray80", length(method_names))
names(bar_color) = method_names
bar_color["Raw"] = "black"
bar_color["DeSCI"] = "red"
text_color = rep("black", length(method_names))
names(text_color) = method_names
text_color["DeSCI"] = "red"
path_fragments = unlist(strsplit(our_result, "/", fixed = T))
length_path_fragments = length(path_fragments)
### make output dir
outdir = paste(c(path_fragments[1:(length_path_fragments - 1)], "quantitative_measures_round", path_fragments[length_path_fragments]), collapse = "/")
dir.create(outdir, showWarnings = F, recursive = T)
### calculate correlation of FISH
cor_mat = matrix(ncol = 3, nrow = 0, dimnames = list(c(), c("Gene x", "Gene y", "FISH Correlation")))
cor_all = list()
for(ii in c("FISH", method_names)){
  cor_all[[ii]] = matrix(nrow = length(use_genes), ncol = length(use_genes), dimnames = list(use_genes, use_genes))
}
for(ii in use_genes){
  for(jj in use_genes){
    cat(ii, "\t", jj, "\n")
    for(kk in c("FISH", method_names)){
      if(kk == "SAVER"){
        cor_all[["SAVER"]][ii, jj] = cor(dataset_list[["SAVER_gamma"]][ii, ], dataset_list[["SAVER_gamma"]][jj, ], use = "pairwise.complete.obs")
      }else{
        cor_all[[kk]][ii, jj] = cor(dataset_list[[kk]][ii, ], dataset_list[[kk]][jj, ], use = "pairwise.complete.obs")
      }
    }
  }
}
###correlation_map
fish_mask_mat = !is.na(cor_all[["FISH"]])
for(ii in 1:length(use_genes)){
  fish_mask_mat[ii, seq(ii)] = FALSE
}
fish_mask = as.vector(fish_mask_mat)
fish_gene_mat = apply(fish_mask_mat, 1, names)
outfile = paste0(outdir, "/correlation_map.pdf")
pdf(outfile, width = 5.75, height = 6)
for(ii in setdiff(method_names, "Raw")){
  fish_corr = as.vector(cor_all[["FISH"]])[fish_mask]
  raw_corr = as.vector(cor_all[["Raw"]])[fish_mask]
  this_corr = as.vector(cor_all[[ii]])[fish_mask]
  plot(fish_corr, raw_corr, xlim = c(-0.2, 1), ylim = c(-0.2, 1), col = "gray80", pch = 16, xlab = "FISH", ylab = "scRNA-seq", main = paste0("Gene - Gene Correlation\n", sum(fish_mask), " pairs"))
  points(fish_corr, this_corr, col = "red")
  abline(c(0, 1))
  legend("bottomright", c(paste0("RMSE(", ii, ") - ", round(rmse(fish_corr, this_corr), 6)),
                      paste0("RMSE(Raw) - ", round(rmse(fish_corr, raw_corr), 6))),
         col = c("red", "gray80"), pch=c(1, 16))
}
dev.off()
#  heatmap & CMD
outfile = paste0(outdir, "/correlation_heatmap.pdf")
use_order = order(colMeans(matrix(apply(cor_all[["FISH"]], 2, function(x){
  x[is.na(x)] = 0
  return(x)
})[t(t(fish_gene_mat) != colnames(fish_gene_mat))], nrow = nrow(cor_all[["FISH"]]) - 1, ncol = ncol(cor_all[["FISH"]]), dimnames = list(c(), colnames(cor_all[["FISH"]])))), decreasing = T)
pdf(outfile, height = 7, width = 12.5)
cmd_vector = layout_correlogram_plot(cor_all, use_order=use_order)
dev.off()

names(cmd_vector) = names(cor_all)
cmd_vector = cmd_vector[setdiff(names(cmd_vector), "FISH")]
barplot_usage(cmd_vector, main = "CMD", bar_color = bar_color, text_color = text_color, use_data_order = T)
dev.off()
###
###two dimensional distributions compare
rescale_mean_list = list()
for(ii in method_names){
  if(ii != "SAVER"){
    rescale_mean_list[[ii]] = mean_norm_fun(dataset_list[[ii]], dataset_list[["FISH"]])
  }else{
    rescale_mean_list[[ii]] = mean_norm_fun(dataset_list[["SAVER_gamma"]], dataset_list[["FISH"]])
  }
}
max_points = ncol(dataset_list[["FISH"]])
for(ii in rescale_mean_list){
  max_points = max(c(max_points, ncol(ii)))
}
if("GAPDH" %in% rownames(dataset_list[["FISH"]])){
  saver_style_filt_norm_list = list()
  for(ii in method_names){
    if(ii != "SAVER"){
      saver_style_filt_norm_list[[ii]] = saver_density_norm_method(dataset_list[[ii]], dataset_list[["FISH"]])
    }else{
      saver_style_filt_norm_list[[ii]] = saver_norm_fun(saver_filter_fun(binom_result(dataset_list[["SAVER"]], dataset_list[["FISH"]])))
    }
  }
  saver_style_filt_norm_list[["FISH"]] = saver_norm_fun(saver_filter_fun(dataset_list[["FISH"]]))
}else{
  saver_style_filt_norm_list = rescale_mean_list
  saver_style_filt_norm_list[["SAVER"]] = binom_result_mean(dataset_list[["SAVER"]], dataset_list[["FISH"]])
  saver_style_filt_norm_list[["FISH"]] = dataset_list[["FISH"]]
}
#  Gini
gini_result_list = list()
for(ii in names(saver_style_filt_norm_list)){
  if(ii != "SAVER"){
    gini_result_list[[ii]] = apply(dataset_list[[ii]], 1, function(x){gini(x[complete.cases(x)])})
  }else{
    gini_result_list[[ii]] = apply(dataset_list[["SAVER_gamma"]], 1, function(x){gini(x[complete.cases(x)])})
  }
}
outfile = paste0(outdir, "/Gini.pdf")
mfrow = c(2, 4)
pdf(outfile, height = mfrow[1] * 3.75, width = mfrow[2] * 3.5)
gini_rmse = layout_gini(gini_result_list, use_genes)
dev.off()
###density
plot_genes = setdiff(use_genes, "GAPDH")
ks_matrix = matrix(nrow = length(plot_genes), ncol = length(method_names), dimnames = list(plot_genes, method_names))
mfrow = c(3, 6)
outfile = paste0(outdir, "/density.pdf")
pdf(outfile, height = mfrow[1] * 3, width = mfrow[2] * 2.75)
par(mfrow = mfrow)
for(ii in plot_genes){
  this_fish = saver_style_filt_norm_list[["FISH"]][ii, ]
  this_fish = this_fish[!is.na(this_fish)]
  fish_density = density(this_fish)
  zero_proportion = round(100 * (1 - sum(dataset_list[["Raw"]][ii, ] > 0) / length(used_cells)), digits = 4)
  xlim_max = as.numeric(quantile(this_fish, 0.90)) * 2
  ylim_max = max(fish_density$y) * 2
  plot(fish_density, lwd = 4, col = "gray80",
       xlim = c(min(fish_density$x, 0), xlim_max + 5),
       ylim = c(0, ylim_max), yaxt = "n", bty="n",
       main = paste0(firstup(ii), " (", zero_proportion, "%)"),
       sub = "", ylab = "Density", xlab = "mRNA Counts")
  axis(2, labels = FALSE, lwd.ticks = 0)
  par(las = 0)
  dens.bw <- fish_density$bw
  for(method_index in 1:length(method_names)){
    this_method_expression = saver_style_filt_norm_list[[method_names[method_index]]][ii, ]
    if(length(unique(this_method_expression)) == 1){
      this_method_expression[1] = this_method_expression[1] * 1.0001
    }
    if(method_index <= 2){
      lines(density(this_method_expression, bw = dens.bw), lwd = 3, col=method_color[method_index])
    }
    ks_matrix[ii, method_index] = ks.test(delect_lt0.5(this_method_expression), delect_lt0.5(this_fish))$statistic
  }
  if(ii %in% use_genes[mfrow[2] + (mfrow[1] * mfrow[2] * seq(0, floor(length(use_genes) / mfrow[1] * mfrow[2])))]){
    legend("topright", c("RNA FISH", method_names[1:2]), lty = rep(1, length(method_names[1:2])),
           lwd = rep(3, length(method_names[1:2])), col = c("gray80", method_color[1:length(method_names[1:2])]), box.lty = 0, xjust = 1, yjust = 1)
  }
}
dev.off()
outfile = paste0(outdir, "/density_summary.pdf")
pdf(outfile, height = 6, width = ncol(ks_matrix) * 0.6)
ks_mean = Matrix::colMeans(ks_matrix)
standard_error_ks = apply(ks_matrix, 2, function(x) sqrt(var(x)/length(x)))
par(mfrow = c(1, 1))
barplot_usage(ks_mean, main = "K-S Statistic", bar_color = bar_color, text_color = text_color, use_data_order = T, standard_error = standard_error_ks)
dev.off()
#####2d_distribution
# Make the plots
dist_outdir = paste0(outdir, "/2d_distribution")
dir.create(dist_outdir, showWarnings = F)
mfrow = c(2, 5)
pairs_2d_distribution = cor_mat[order(abs(cor_mat[, 3]), decreasing = TRUE), ]

library(parallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterExport(cl, varlist = c("dist_outdir", "text_color", "bar_color", "max_points", "rescale_mean_list", "method_names", "used_cells", "mfrow", "fish_gene_mat", "fish_mask_mat", "dataset_list", "generic_functions_path"))
return_list = parLapply(cl, 1:sum(fish_mask), function(ii){
  source(generic_functions_path)
  gene_x = fish_gene_mat[t(fish_mask_mat)][ii]
  gene_y = t(fish_gene_mat)[t(fish_mask_mat)][ii]
  x_dropout_rate = round(100 * (1 - sum(dataset_list[["Raw"]][gene_x, ] > 0) / length(used_cells)), digits = 4)
  y_dropout_rate = round(100 * (1 - sum(dataset_list[["Raw"]][gene_y, ] > 0) / length(used_cells)), digits = 4)
  x_fish_raw = dataset_list[["FISH"]][gene_x, ]
  y_fish_raw = dataset_list[["FISH"]][gene_y, ]
  select_cell = !is.na(x_fish_raw) & !is.na(y_fish_raw)
  x_fish = x_fish_raw[select_cell]
  y_fish = y_fish_raw[select_cell]
  fish_pair_mat = matrix(c(x_fish, y_fish), ncol = 2)
  ks_stat = c()
  corr_score = cor(x_fish, y_fish)
  x_i = list(FISH = x_fish)
  y_i = list(FISH = y_fish)
  for(jj in 1:length(method_names)){
    method_name = method_names[jj]
    x_i[[method_name]] = rescale_mean_list[[method_name]][gene_x, ]
    y_i[[method_name]] = rescale_mean_list[[method_name]][gene_y, ]
    ks_result = unlist(ks2d2s(np_array(round(x_fish)),
                              np_array(round(y_fish)), 
                              np_array(round(x_i[[method_name]])),
                              np_array(round(y_i[[method_name]])), extra=T))
    ks_stat = c(ks_stat, ks_result[2])
    corr_score = c(corr_score, cor(x_i[[method_name]], y_i[[method_name]]))
  }
  names(corr_score) = c("FISH", method_names)
  names(ks_stat) = method_names
  this_name = c(paste(gene_x, gene_y, sep = "_"))
  pdf(paste0(dist_outdir, "/", this_name, ".pdf"),
      height = mfrow[1] * 4,
      width = mfrow[2] * 3.75)
  par(mfrow = mfrow)
  nbin = 128
  x_fish_95 = quantile(x_fish, 0.95) + 1### R is from 1 to max + 1
  y_fish_95 = quantile(y_fish, 0.95) + 1
  bandwidth = c(max(x_fish) / nbin, max(y_fish) / nbin)
  for(jj in c("Raw", "FISH", "DeSCI", "SAVER")){
    if(jj == "DeSCI"){
      col.main = "red"
    }else{
      col.main = "black"
    }
    smoothScatter(x = x_i[[jj]], y = y_i[[jj]],
                  ylab = paste0(gene_y, " (", y_dropout_rate, "%)"),
                  xlab = paste0(gene_x, " (", x_dropout_rate, "%)"), cex = 1.5,
                  xlim = c(0, x_fish_95), ylim = c(0, max(y_fish_95)),
                  lwd = 2, main = jj, nrpoints = 0, col.main = col.main, nbin = nbin,
                  bandwidth = bandwidth)
  }
  barplot_usage(ks_stat, main = "Fasano and Franceschini's Test", bar_color = bar_color, text_color = text_color, use_data_order = T)
  for(jj in c("MAGIC", "scVI", "DCA", "scScope")){
    smoothScatter(x = x_i[[jj]], y = y_i[[jj]],
                  ylab = paste0(gene_y, " (", y_dropout_rate, "%)"),
                  xlab = paste0(gene_x, " (", x_dropout_rate, "%)"), cex = 1.5,
                  xlim = c(0, x_fish_95), ylim = c(0, max(y_fish_95)),
                  lwd = 2, main = jj, nrpoints = 0, col.main = col.main, nbin = nbin,
                  bandwidth = bandwidth)
  }
  dev.off()
  return(list("ks_stat" = matrix(ks_stat, nrow = 1, dimnames = list(paste(gene_x, gene_y, sep = " - "), c())),
              "corr_score" = matrix(corr_score, nrow = 1, dimnames = list(paste(gene_x, gene_y, sep = " - "), c()))))
})
stopCluster(cl)
ks_stat_mat = matrix(nrow = 0, ncol = length(method_names), dimnames = list(c(), method_names))
corr_mat = matrix(nrow = 0, ncol = length(method_names) + 1, dimnames = list(c(), c("FISH", method_names)))

for(ii in return_list){
  ks_stat_mat = rbind(ks_stat_mat, ii$ks_stat)
  corr_mat = rbind(corr_mat, ii$corr_score)
}
###all_compare
mean_ks_stat = Matrix::colMeans(ks_stat_mat)
standard_error_ks_stat = apply(ks_stat_mat, 2, function(x) sqrt(var(x)/length(x)))
outfile = paste(outdir, "/score_compare.pdf", sep = "")
pdf(outfile, height = 5, width = 9)
par(mfrow = c(1, 2))
barplot_usage(mean_ks_stat, main = "Fasano and\nFranceschini's Test", cex.main = 1.5,bar_color = bar_color, text_color = text_color, use_data_order = T, standard_error = standard_error_ks_stat)
corr_rmse = sapply(method_names, function(x) rmse(corr_mat[, "FISH"], corr_mat[, x]))
barplot_usage(corr_rmse, main = "FISH - Impute\nCorrelation RMSE", cex.main = 1.5, bar_color = bar_color, text_color = text_color, use_data_order = T)
dev.off()

outfile = paste(outdir, "/Bar_plot.pdf", sep = "")
pdf(outfile, height = 3, width = 8)
plot_height = 4
plot_width = 3.5
plot_region = matrix(seq(3), nrow = 1)
this_height = rep(plot_height, nrow(plot_region))
this_width = rep(plot_width, ncol(plot_region))
this_index = max(plot_region) + 1
layout_mat = plot_region
xlab_region = matrix(rep(this_index, ncol(plot_region)), nrow = 1)
layout_mat = rbind(layout_mat, xlab_region)
this_height = c(this_height, 0.5)
layout(mat = layout_mat, heights = this_height, widths = this_width)
par(mar = c(1, 4.1, 4.1, 2.1))
barplot_usage_new(cmd_vector, main = "CMD", method_color = method_color, use_data_order = T)
barplot_usage_new(gini_rmse, main = "Gini RMSE", method_color = method_color, use_data_order = T)
barplot_usage_new(mean_ks_stat, main = "Fasano and\nFranceschini's Test", cex.main = 1.5, method_color = method_color, use_data_order = T, standard_error = standard_error_ks_stat)
par(mar = rep(0, 4))
plot(1, type = "n", axes = FALSE, xlab="", ylab="")
legend(x = "top",inset = 0, legend = names(method_color), fill = method_color, horiz = TRUE, border = NA, bty = "n")
dev.off()