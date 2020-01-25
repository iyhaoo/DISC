args<-commandArgs(trailingOnly=TRUE)
if(length(args) != 1){
  stop("R --slave < this_code.r --args <filt-loom>")
}
#install.packages("devtools")
#library(devtools)
#install_github("ChenMengjie/VIPER")
library(VIPER)
library(loomR)
library(rhdf5)

delete_last_element <- function(x){
  return(x[1: (length(x) - 1)])
}
get_last_element <- function(x){
  return(x[length(x)])
}

starttime = Sys.time()
loom_path = args[1]
dir_path = paste(delete_last_element(unlist(strsplit(loom_path, "/", fixed = T))), collapse = "/")
output_dir = paste0(dir_path, "/imputation")
method_dir = paste0(dir_path, "/imputation/VIPER_tmp")
dir.create(output_dir, recursive = T, showWarnings = F)
dir.create(method_dir, showWarnings = F)
print(loom_path)
min_expressed_cell = 10
min_expressed_cell_average_expression = 1
output_h5 = paste(output_dir, sub(".loom", paste0("_VIPER_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression, ".hdf5"), get_last_element(unlist(strsplit(loom_path, "/", fixed = T))), fixed = TRUE), sep = "/")
print(output_h5)
if(file.exists(output_h5)){
  print("pass")
  next()
}
raw_data_loompy = connect(loom_path)
gene_bc_mat = t(as.matrix(raw_data_loompy$matrix[, ]))
gene_name = raw_data_loompy[["row_attrs/Gene"]][]
cell_id = raw_data_loompy[["col_attrs/CellID"]][]
raw_data_loompy$close_all()
expressed_cell = Matrix::rowSums(gene_bc_mat > 0)
gene_expression = Matrix::rowSums(gene_bc_mat)
gene_filter = expressed_cell >= min_expressed_cell & gene_expression > expressed_cell * min_expressed_cell_average_expression
filt_data = gene_bc_mat[gene_filter, ]
gene_name_filt = gene_name[gene_filter]
rownames(filt_data) = paste0("Gene_", 1:length(gene_name_filt))
colnames(filt_data) = paste0("Cell_", 1:length(cell_id))
out_dir = paste(method_dir, sub(".loom", paste0("_VIPER_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression), get_last_element(unlist(strsplit(loom_path, "/", fixed = T))), fixed = TRUE), sep = "/")
dir.create(out_dir, showWarnings = F, recursive = T)
print(filt_data[1:5, 1:5])
print(dim(filt_data))
res <- VIPER(gene.expression = filt_data, outdir = out_dir)
row.names(res$imputed_log) = gene_name_filt
row.names(res$imputed) = gene_name_filt
colnames(res$imputed_log) = cell_id
colnames(res$imputed) = cell_id
output_rds_name = sub(".loom", paste0("_VIPER_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression, ".rds"), get_last_element(unlist(strsplit(loom_path, "/", fixed = T))), fixed = TRUE)
saveRDS(res, paste(method_dir, output_rds_name, sep = "/"))
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(length(gene_name_filt), length(cell_id)),
                storage.mode = "double",
                chunk=c(length(gene_name_filt), 1))
h5write(cell_id, output_h5,"cell_id")
h5write(gene_name_filt, output_h5,"gene_name")
h5write(res$imputed, file=output_h5, name="imputation")
print(Sys.time() - starttime)

