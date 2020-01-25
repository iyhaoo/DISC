args<-commandArgs(trailingOnly=TRUE)
if(length(args) != 1){
  stop("R --slave < this_code.r --args <filt-loom>")
}
#library(devtools)
#install_github('gongx030/DrImpute')
library(DrImpute)
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
dir.create(output_dir, recursive = T, showWarnings = F)
print(loom_path)
min_expressed_cell = 10
min_expressed_cell_average_expression = 1
output_h5 = paste(output_dir, sub(".loom", paste0("_DrImpute_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression, ".hdf5"), get_last_element(unlist(strsplit(loom_path, "/", fixed = T))), fixed = TRUE), sep = "/")
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
colnames(gene_bc_mat) = cell_id
rownames(gene_bc_mat) = gene_name
expressed_cell = Matrix::rowSums(gene_bc_mat > 0)
gene_expression = Matrix::rowSums(gene_bc_mat)
gene_filter = expressed_cell >= min_expressed_cell & gene_expression > expressed_cell * min_expressed_cell_average_expression
filt_data = gene_bc_mat[gene_filter, ]
gene_name_filt = gene_name[gene_filter]
logdat_imp <- DrImpute(log1p(filt_data))
resume_impute = expm1(logdat_imp)
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(length(gene_name_filt), length(cell_id)),
                storage.mode = "double",
                chunk=c(length(gene_name_filt), 1))
h5write(cell_id, output_h5,"cell_id")
h5write(gene_name_filt, output_h5,"gene_name")
h5write(resume_impute, file=output_h5, name="imputation")
print(Sys.time() - starttime)