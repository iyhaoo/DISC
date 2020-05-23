single_plot_line = function(data_mat, method_color, title = "1000 Genes", ylab = "Duration in hours", use_log=T, ymax=NULL){
  data_df = data.frame()
  for(ii in 1:ncol(data_mat)){
    data_df[seq(nrow(data_mat)) + nrow(data_mat) * (ii - 1), "method"] <- rep(colnames(data_mat)[ii], nrow(data_mat))
    data_df[seq(nrow(data_mat)) + nrow(data_mat) * (ii - 1), "method_color"] <- rep(method_color[ii], nrow(data_mat))
    data_df[seq(nrow(data_mat)) + nrow(data_mat) * (ii - 1), "cells"] <- factor(rownames(data_mat), levels = rownames(data_mat))
    data_df[seq(nrow(data_mat)) + nrow(data_mat) * (ii - 1), "value"] <- as.numeric(data_mat[, ii])
  }
  data_df$method = factor(data_df$method, levels = colnames(data_mat))
  this_plot = ggplot(data = data_df, aes(x = cells, y = value, color = method, group = method)) + geom_line(size = 1.5) + geom_point(size = 4) +
    scale_color_manual(values = method_color) + 
    scale_size_manual(values = 8) + 
    labs(title = title) +  xlab("Cell Count") + ylab(ylab) +
    scale_x_discrete() + 
    theme(panel.border = element_blank(), axis.line = element_line(), legend.background = element_blank(),
          panel.background = element_blank(), legend.spacing = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA), legend.title = element_blank(),
          axis.title = element_text(size = 15, face = "bold"), axis.text = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
  if(use_log){
    this_plot = this_plot + scale_y_continuous(trans = 'log10')
  }
  if(!is.null(ymax)){
    this_plot = this_plot + ylim(c(0, ymax * 1.08)) + geom_hline(yintercept = ymax, color = "gray", size = 1, linetype=2)
  }
  return(this_plot)
}

single_plot_bar = function(data_mat, method_color, title = "10000 Genes", ylab = "Duration in hours", use_log=T, ymax=NULL){
  data_df = data.frame()
  for(ii in 1:ncol(data_mat)){
    data_df[seq(nrow(data_mat)) + nrow(data_mat) * (ii - 1), "method"] <- rep(colnames(data_mat)[ii], nrow(data_mat))
    data_df[seq(nrow(data_mat)) + nrow(data_mat) * (ii - 1), "method_color"] <- rep(method_color[ii], nrow(data_mat))
    data_df[seq(nrow(data_mat)) + nrow(data_mat) * (ii - 1), "cells"] <- factor(rownames(data_mat), levels = rownames(data_mat))
    data_df[seq(nrow(data_mat)) + nrow(data_mat) * (ii - 1), "value"] <- as.numeric(data_mat[, ii])
  }
  data_df$method = factor(data_df$method, levels = colnames(data_mat))
  this_plot = ggplot(data = data_df, aes(x = method, y = value, fill = cells, group = cells)) + geom_bar(stat="identity", position=position_dodge()) +
    scale_color_manual(values = method_color) + 
    scale_size_manual(values = 8) + 
    labs(title = title) +  xlab("Cell Count") + ylab(ylab) +
    scale_x_discrete() + 
    theme(panel.border = element_blank(), axis.line = element_line(), legend.background = element_blank(),
          panel.background = element_blank(), legend.spacing = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA), legend.title = element_blank(),
          axis.title = element_text(size = 15, face = "bold"), axis.text = element_text(size = 12, face = "bold"), axis.text.x.bottom = element_text(vjust = 0.5, angle = 45),
          legend.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
  if(use_log){
    this_plot = this_plot + scale_y_continuous(trans = 'log10')
  }
  if(!is.null(ymax)){
    this_plot = this_plot + ylim(c(0, ymax * 1.08)) + geom_hline(yintercept = ymax, color = "gray", size = 1, linetype=2)
  }
  return(this_plot)
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

library(ggplot2)
library(gridExtra)
method_name = c("DISC", "SAVER", "scImpute", "VIPER", "MAGIC", "DCA", "deepImpute", "scScope", "scVI")
cell_number = c("10k", "50k", "100k", "500k", "1.3m", "2.6m")
method_color = c("#FF0000", "#000080", "#BFBF00", "#408000", "#804000", "#00FF00", "#FF8000", "#FF00FF", "#00FFFF")
method_lty = c(rep(1, 4), rep(4, 4), 1)


performance_mat = as.matrix(read.table("E:/DISC/reproducibility/data_preparation_and_imputation/raw_scripts/summary.txt", sep = "\t", header = F))[, -1]
runtime = performance_mat[2:14, ][-1, -1]
ram = performance_mat[17:29, ][-1, -1]


runtime1k = runtime[1:6, ]
colnames(runtime1k) = method_name
rownames(runtime1k) = cell_number

runtime10k = runtime[7:12, ]
colnames(runtime10k) = method_name
rownames(runtime10k) = cell_number

ram1k = ram[1:6, ]
colnames(ram1k) = method_name
rownames(ram1k) = cell_number

ram10k = ram[7:12, ]
colnames(ram10k) = method_name
rownames(ram10k) = cell_number

runtime1k_plot = single_plot_line(runtime1k, method_color, "1000 Genes", "Duration in hours", use_log = F, ymax = 24)
runtime10k_plot = single_plot_bar(runtime10k, method_color, "10000 Genes", "Duration in hours", use_log = F)

ram1k_plot = single_plot_line(ram1k, method_color, "1000 Genes", "Peak RAM (GB)", use_log = F, ymax = 128)
ram10k_plot = single_plot_bar(ram10k, method_color, "10000 Genes", "Peak RAM (GB)", use_log = F, ymax = 128)


pdf("E:/DISC/reproducibility/data_preparation_and_imputation/raw_scripts/performance_ram_1k.pdf", height = 4, width = 12)
grid.arrange(arrangeGrob(runtime1k_plot + theme(legend.position = "none", plot.title = element_blank(), axis.title.x = element_blank()),
                         ram1k_plot + theme(legend.position = "none", plot.title = element_blank(), axis.title.x = element_blank()),
                         top = gridExtra:::text_grob(label = "1000 Genes", fontsize = 15, fontface = 2, hjust = 0.2), nrow = 1),
             g_legend(runtime1k_plot), bottom = gridExtra:::text_grob(label = "Cell Count", fontsize = 15, fontface = 2, hjust = 1), ncol = 2, widths = c(4, 2))

dev.off()

pdf("E:/DISC/reproducibility/data_preparation_and_imputation/raw_scripts/performance_ram_10k.pdf", height = 4, width = 12)
grid.arrange(arrangeGrob(runtime10k_plot + theme(legend.position = "none", plot.title = element_blank(), axis.title.x = element_blank()),
                         ram10k_plot + theme(legend.position = "none", plot.title = element_blank(), axis.title.x = element_blank()),
                         top = gridExtra:::text_grob(label = "10000 Genes", fontsize = 15, fontface = 2, hjust = 0.2), nrow = 1),
             g_legend(runtime10k_plot), bottom = gridExtra:::text_grob(label = "Cell Count", fontsize = 15, fontface = 2, hjust = 1), ncol = 2, widths = c(8, 1))

dev.off()











