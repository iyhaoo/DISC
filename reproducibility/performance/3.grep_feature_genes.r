args = commandArgs(trailingOnly=TRUE)
if(length(args) < 3){
  stop("R --slave < this_code.r --args <gene-bc imputation h5/loom> <min_expressed_cell> <min_expressed_cell_average_expression>")
}
source("/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/generic_functions.r")
loom_path = args[1]
min_expressed_cell = as.integer(args[2])
min_expressed_cell_average_expression = as.numeric(args[3])
output_tsv = paste0(paste(delete_last_element(unlist(strsplit(loom_path, ".", fixed = T))), collapse = "/"), "_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression, "_hvg.tsv")
gene_bc_mat = t(h5read(loom_path, "matrix"))
colnames(gene_bc_mat) = h5read(loom_path, "col_attrs/CellID")
row.names(gene_bc_mat) = h5read(loom_path, "row_attrs/Gene")
expressed_cell = Matrix::rowSums(gene_bc_mat > 0)
gene_expression = Matrix::rowSums(gene_bc_mat)
gene_filter = expressed_cell >= min_expressed_cell & gene_expression > expressed_cell * min_expressed_cell_average_expression
hvf.info = FindVariableFeatures_vst_by_genes(gene_bc_mat[gene_filter, ])
hvg = rownames(hvf.info)[order(hvf.info$variance.standardized, decreasing = T)]
write.table(hvf.info[order(hvf.info$variance.standardized, decreasing = T), ], output_tsv, sep = "\t", quote = F, row.names = T, col.names = T)





