args<-commandArgs(trailingOnly=TRUE)
if(length(args) < 2){
  stop("R --slave < this_code.r --args <loom> <use_core> <min_expressed_cell> <min_expressed_cell_average_expression>")
}
library(doParallel)
library(Matrix)
library(SAVER)
library(scales)
library(rhdf5)
library(loomR)

delete_last_element <- function(x){
  return(x[1: (length(x) - 1)])
}
get_last_element <- function(x){
  return(x[length(x)])
}

gamma_result = function(saver_result, num_of_obs=5){
  saver_est = saver_result$estimate
  saver_se = saver_result$se
  saver_sf = saver_result$info$size.factor
  result_mat = t(sapply(1:nrow(saver_est), FUN = function(ii){
    rgamma(n = num_of_obs * ncol(saver_est),
           shape = saver_est[ii, ]^2/saver_se[ii, ]^2,
           rate = saver_est[ii, ]/saver_se[ii, ]^2/saver_sf)
  }))
  rownames(result_mat) = rownames(saver_est)
  colnames(result_mat) = colnames(saver_est)
  return(result_mat)
}

starttime = Sys.time()
loom_path = args[1]
use_core = as.integer(args[2])
print(use_core)
dir_path = paste(delete_last_element(unlist(strsplit(loom_path, "/", fixed = T))), collapse = "/")
output_dir = paste0(dir_path, "/imputation")
method_dir = paste0(dir_path, "/imputation/SAVER_tmp")
dir.create(output_dir, recursive = T, showWarnings = F)
dir.create(method_dir, showWarnings = F)
print(loom_path)
if(length(args) >= 3){
  min_expressed_cell = as.integer(args[3])
}else{
  min_expressed_cell = 10
}
if(length(args) >= 4){
  min_expressed_cell_average_expression = as.numeric(args[4])
}else{
  min_expressed_cell_average_expression = 1
}
output_h5 = paste(output_dir, sub(".loom", paste0("_SAVER_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression, ".hdf5"), get_last_element(unlist(strsplit(loom_path, "/", fixed = T))), fixed = TRUE), sep = "/")
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
print(dim(filt_data))
gene_name_filt = gene_name[gene_filter]
cl <- makeCluster(use_core, outfile = "")
registerDoParallel(cl)
saver_result <- saver(filt_data)
stopCluster(cl)
output_rds_name = sub(".loom", paste0("_SAVER_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression, ".rds"), get_last_element(unlist(strsplit(loom_path, "/", fixed = T))), fixed = TRUE)
saveRDS(saver_result, paste(method_dir, output_rds_name, sep = "/"))
gamma_res = gamma_result(saver_result, 1)
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(length(gene_name_filt), length(cell_id)),
                storage.mode = "double",
                chunk=c(length(gene_name_filt),1))
h5write(cell_id, output_h5,"cell_id")
h5write(gene_name_filt, output_h5,"gene_name")
h5write(gamma_res, file=output_h5, name="imputation")
print(Sys.time() - starttime)


