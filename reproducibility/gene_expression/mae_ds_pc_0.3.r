source("D:/study/single_cell/our_pipeline/qm_functions.r")
method_names = c("Observed", "DeSCI", "SAVER", "MAGIC", "DCA", "scScope", "scVI")
method_color = c("gray80", "red", "blue4", "yellow4", "green", "purple", "cyan")
names(method_color) = method_names
MAE_result_dir = "E:/DeSCI/fn/pbmc3k/MAE_ds_0.3"
mae_eq0 = as.matrix(read.table(paste0(MAE_result_dir, "/mae_eq0_with_outliers.tsv"), sep = "\t"))
mae_gt0 = as.matrix(read.table(paste0(MAE_result_dir, "/mae_gt0_with_outliers.tsv"), sep = "\t"))
outfile = paste0(MAE_result_dir, "/MAE_ds_0.3_with_outliers.pdf")
pdf(outfile, height = 3.5, width = 7)
this_width = c(4, 4)
this_height = c(4, 0.5)
layout_mat = t(matrix(c(1, 2, 3, 3), ncol = 2))
layout(mat = layout_mat, heights = this_height, widths = this_width)
par(mar = c(1, 5, 4.1, 1))
if(is.na(max(mae_eq0[, "ds_0.3_ste"]))){
  barplot_usage_new(mae_eq0[, "ds_0.3"], main = "Zero entries", method_color = method_color, use_data_order = T, use_log1p = T, ylab = "Log1p (MAE)", cex.lab = 1.5, font.lab = 2, cex.main = 1.5)
  print("Without standard error")
}else{
  barplot_usage_new(mae_eq0[, "ds_0.3"], main = "Zero entries", method_color = method_color, use_data_order = T, standard_error = mae_eq0[, "ds_0.3_ste"], use_log1p = T, ylab = "Log1p (MAE)", cex.lab = 1.5, font.lab = 2, cex.main = 1.5)
}
if(is.na(max(mae_gt0[, "ds_0.3_ste"]))){
  barplot_usage_new(mae_gt0[, "ds_0.3"], main = "Non-zero entries", method_color = method_color, use_data_order = T, use_log1p = T, ylab = "Log1p (MAE)", cex.lab = 1.5, font.lab = 2, cex.main = 1.5)
  print("Without standard error")
}else{
  barplot_usage_new(mae_gt0[, "ds_0.3"], main = "Non-zero entries", method_color = method_color, use_data_order = T, standard_error = mae_gt0[, "ds_0.3_ste"], use_log1p = T, ylab = "Log1p (MAE)", cex.lab = 1.5, font.lab = 2, cex.main = 1.5)
}
par(mar = rep(0, 4))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0, legend = names(method_color), fill = method_color, horiz = TRUE, border = NA, bty = "n", col = 2)
dev.off()






