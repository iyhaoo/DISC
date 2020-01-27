utilities_path = "/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/utilities.r"
source(utilities_path)
####################main###############
dataset_list = list()
dataset_list[["FISH"]] = readh5_loom("/home/yuanhao/data/fn/melanoma/fishSubset.loom")
raw_input_data = "/home/yuanhao/data/fn/melanoma/imputation/dropseq_filt_ls_mc_10_mce_1.loom"
use_genes = intersect(rownames(dataset_list[["FISH"]]), get_loom_gene(raw_input_data))
print(length(use_genes))
print(use_genes)
dataset_list[["Raw"]] = readh5_loom(raw_input_data, use_genes)
used_cells = colnames(dataset_list[["Raw"]])
### DISC
our_result = "/home/yuanhao/DISC_imputation_result/MELANOMA/result/imputation.loom"
dataset_list[["DISC"]] = readh5_loom(our_result, use_genes)
### Other methods
dataset_list[["SAVER"]] = readRDS("/home/yuanhao/data/fn/melanoma/imputation/SAVER_tmp/dropseq_filt_ls_SAVER_mc_10_mce_1.rds")
set.seed(42)
dataset_list[["SAVER_gamma"]] = gamma_result(dataset_list[["SAVER"]], num_of_obs=1)[use_genes, used_cells]
dataset_list[["MAGIC"]] = readh5_imputation("/home/yuanhao/data/fn/melanoma/imputation/dropseq_filt_ls_MAGIC_mc_10_mce_1.hdf5", use_genes, used_cells)
dataset_list[["DCA"]] = readh5_imputation("/home/yuanhao/data/fn/melanoma/imputation/dropseq_filt_ls_DCA_mc_10_mce_1.hdf5", use_genes, used_cells)
dataset_list[["deepImpute"]] = readh5_imputation("/home/yuanhao/github_repositories/DISC/reproducibility/data/MELANOMA/imputation/raw_deepImpute_mc_10_mce_1.hdf5", use_genes, used_cells)
dataset_list[["scScope"]] = readh5_imputation("/home/yuanhao/data/fn/melanoma/imputation/dropseq_filt_ls_scScope_mc_10_mce_1.hdf5", use_genes, used_cells)
dataset_list[["scVI"]] = readh5_imputation("/home/yuanhao/data/fn/melanoma/imputation/dropseq_filt_ls_scVI_mc_10_mce_1.hdf5", use_genes, used_cells)
### Output settings
method_names = c("Raw", "DISC", "SAVER", "MAGIC", "DCA", "deepImpute", "scScope", "scVI")
method_color = c("gray80", "red", "blue4", "yellow4", "green", "orange", "purple", "cyan")
names(method_color) = method_names
bar_color = rep("gray50", length(method_names))
names(bar_color) = method_names
bar_color["Raw"] = "gray80"
bar_color["DISC"] = "red"
text_color = rep("black", length(method_names))
names(text_color) = method_names
text_color["DISC"] = "red"
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
#### correlation map
fish_mask_mat = !is.na(cor_all[["FISH"]])
for(ii in 1:length(use_genes)){
  fish_mask_mat[ii, seq(ii)] = FALSE
}
fish_mask = as.vector(fish_mask_mat)
fish_gene_mat = apply(fish_mask_mat, 1, names)
outfile = paste0(outdir, "/correlation_map.pdf")
mfrow = c(2, 4)
pdf(outfile, height = mfrow[1] * 3.75, width = mfrow[2] * 3.25)
layout_scatter(cor_all, method_names, fish_mask, this_xlab = "FISH", this_ylab = "scRNA-seq", xlim = c(-0.2, 1), ylim = c(-0.2, 1))
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
saveRDS(cor_all, paste0(outdir, "/cor_all.rds"))
#### Normalization
example_gene = c("WNT5A", "SOX10")
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
saver_style_filt_norm_list = rescale_mean_list
saver_style_filt_norm_list[["SAVER"]] = binom_result_mean(dataset_list[["SAVER"]], dataset_list[["FISH"]])
saver_style_filt_norm_list[["FISH"]] = dataset_list[["FISH"]]
norm_for_density = saver_style_filt_norm_list
for(ii in names(norm_for_density)){
  norm_for_density[[ii]] = norm_for_density[[ii]][example_gene,]
}
saveRDS(norm_for_density, paste0(outdir, "/norm_for_density.rds"))

gene1 = example_gene[1]
gene2 = example_gene[2]
gene_x = list()
gene_y = list()
gene_x_ws = list()
gene_y_ws = list()
for(ii in method_names){
  gene_x[[ii]] = rescale_mean_list[[ii]][gene1, ]
  gene_y[[ii]] = rescale_mean_list[[ii]][gene2, ]
  if(ii != "SAVER"){
    gene_x_ws[[ii]] = dataset_list[[ii]][gene1, ]
    gene_y_ws[[ii]] = dataset_list[[ii]][gene2, ]
  }else{
    gene_x_ws[[ii]] = dataset_list[["SAVER_gamma"]][gene1, ]
    gene_y_ws[[ii]] = dataset_list[["SAVER_gamma"]][gene2, ]
  }
}
x_fish_raw = dataset_list[["FISH"]][gene1, ]
y_fish_raw = dataset_list[["FISH"]][gene2, ]
select_cell = !is.na(x_fish_raw) & !is.na(y_fish_raw)
x_fish = x_fish_raw[select_cell]
y_fish = y_fish_raw[select_cell]
gene_x[["FISH"]] = x_fish
gene_y[["FISH"]] = y_fish
gene_x_ws[["FISH"]] = x_fish
gene_y_ws[["FISH"]] = y_fish
saveRDS(gene_x, file = paste0(outdir, "/", gene1, ".rds"))
saveRDS(gene_y, file = paste0(outdir, "/", gene2, ".rds"))
saveRDS(gene_x_ws, file = paste0(outdir, "/", gene1, "_ws.rds"))
saveRDS(gene_y_ws, file = paste0(outdir, "/", gene2, "_ws.rds"))

#  Gini
gini_result_list = list()
for(ii in names(saver_style_filt_norm_list)){
  if(ii != "SAVER"){
    gini_result_list[[ii]] = apply(dataset_list[[ii]], 1, function(x){gini(x[complete.cases(x)])})
  }else{
    gini_result_list[[ii]] = apply(dataset_list[["SAVER_gamma"]], 1, function(x){gini(x[complete.cases(x)])})
  }
}
color_point = c("#31a354", "#a63603")
names(color_point) = example_gene
outfile = paste0(outdir, "/Gini.pdf")
mfrow = c(2, 4)
pdf(outfile, height = mfrow[1] * 3.75, width = mfrow[2] * 3.25)
gini_rmse = layout_scatter(gini_result_list, method_names, use_genes, color_point = color_point, this_xlab = "FISH Gini", this_ylab = "scRNA-seq Gini", xlim = c(0, 1), ylim = c(0, 1))
dev.off()
saveRDS(gini_result_list, paste0(outdir, "/gini_result_list.rds"))
###density
#plot_genes = setdiff(use_genes, "GAPDH")
plot_genes = use_genes
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
  dens.bw = fish_density$bw
  ylim_max = max(fish_density$y)
  use_density = list()
  for(method_index in 1:length(method_names)){
    this_method_expression = saver_style_filt_norm_list[[method_names[method_index]]][ii, ]
    if(length(unique(this_method_expression)) == 1){
      this_method_expression[1] = this_method_expression[1] * 1.0001
    }
    #if(method_index <= 2){
    if(T){
      this_density = density(this_method_expression, bw = dens.bw)
      use_density[[method_names[method_index]]] = this_density
      ylim_max = max(c(ylim_max, this_density$y))
    }
    ks_matrix[ii, method_index] = ks.test(delete_lt0.5(this_method_expression), delete_lt0.5(this_fish))$statistic
  }
  plot(fish_density, lwd = 2, col = "black", lty = 1,
       xlim = c(min(fish_density$x, 0), xlim_max + 5),
       ylim = c(0, ylim_max), yaxt = "s", bty="n",
       main = paste0(firstup(ii), " (", zero_proportion, "%)"),
       sub = "", ylab = "Density", xlab = "mRNA Counts")
  par(las = 0)
  for(this_density_name in names(use_density)){
    lines(use_density[[this_density_name]], lwd = 3, col=method_color[this_density_name])
  }
  if(ii %in% plot_genes[mfrow[2] + (mfrow[1] * mfrow[2] * seq(0, floor(length(plot_genes) / mfrow[1] * mfrow[2])))]){
    legend("topright", c("FISH", names(use_density)), lty = rep(1, 1 + length(names(use_density))),
           lwd = rep(3, 1 + length(names(use_density))), col = c("black", method_color[names(use_density)]), box.lty = 0, xjust = 1, yjust = 1)
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
mfrow = c(2, 4)
pairs_2d_distribution = cor_mat[order(abs(cor_mat[, 3]), decreasing = TRUE), ]
library(parallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterExport(cl, varlist = c("dist_outdir", "text_color", "bar_color", "max_points", "rescale_mean_list", "method_names", "used_cells", "mfrow", "fish_gene_mat", "fish_mask_mat", "dataset_list", "utilities_path"))
return_list = parLapply(cl, 1:sum(fish_mask), function(ii){
  source(utilities_path)
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
    ks_stat = c(ks_stat, ks2d2s(round(x_fish), round(y_fish), round(x_i[[method_name]]), round(y_i[[method_name]])))
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
  for(jj in method_names){
    if(jj == "DISC"){
      col.main = "red"
    }else{
      col.main = "black"
    }
    if(jj == "Raw"){
      x_use = dataset_list[[jj]][gene_x, ]
      y_use = dataset_list[[jj]][gene_y, ]
      xlim = c(0, max(x_use))
      ylim = c(0, max(y_use))
      bandwidth = c(xlim[2] / nbin, ylim[2] / nbin)
    }else{
      x_use = x_i[[jj]]
      y_use = y_i[[jj]]
      xlim = c(0, x_fish_95)
      ylim = c(0, y_fish_95)
      bandwidth = c(max(x_fish) / nbin, max(y_fish) / nbin)
    }
    smoothScatter1(x = x_use, y = y_use,
                   xlab = paste0(gene_x, " (", x_dropout_rate, "%)"),
                   ylab = paste0(gene_y, " (", y_dropout_rate, "%)"),
                   cex = 1.5, xlim = xlim, ylim = ylim,
                   lwd = 2, main = paste0(jj, " - FF = ", round(ks_stat[jj], 4)),
                   nrpoints = 0, col.main = col.main, nbin = nbin, bandwidth = bandwidth)
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
saveRDS(ks_stat_mat, paste(outdir, "/ks_stat_mat.rds", sep = ""))
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
output_list = list()
output_list[["cmd"]] = cmd_vector
output_list[["gini_rmse"]] = gini_rmse
output_list[["mean_ks_stat"]] = mean_ks_stat
output_list[["standard_error_ks_stat"]] = standard_error_ks_stat
saveRDS(output_list, paste(outdir, "/Bar_stat.rds", sep = ""))
