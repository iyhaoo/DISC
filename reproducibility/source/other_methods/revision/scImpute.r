args<-commandArgs(trailingOnly=TRUE)
if(length(args) < 2){
  stop("R --slave < this_code.r --args <loom> <use_core> <min_expressed_cell> <min_expressed_cell_average_expression>")
}
source("/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/utilities.r")
library(scImpute)
starttime = Sys.time()
loom_path = args[1]
use_core = as.integer(args[2])
print(use_core)
dir_path = paste(delete_last_element(unlist(strsplit(loom_path, "/", fixed = T))), collapse = "/")
output_dir = paste0(dir_path, "/imputation")
method_dir = paste0(dir_path, "/imputation/scImpute_tmp")
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
output_h5 = paste(output_dir, sub(".loom", paste0("_scImpute_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression, ".hdf5"), get_last_element(unlist(strsplit(loom_path, "/", fixed = T))), fixed = TRUE), sep = "/")
print(output_h5)
if(file.exists(output_h5)){
  print("pass")
  next()
}
gene_bc_mat = readh5_loom(loom_path)
gene_name = rownames(gene_bc_mat)
cell_id = colnames(gene_bc_mat)
expressed_cell = rowSums(gene_bc_mat > 0)
gene_expression = rowSums(gene_bc_mat)
gene_filter = expressed_cell >= min_expressed_cell & gene_expression > expressed_cell * min_expressed_cell_average_expression
filt_data = gene_bc_mat[gene_filter, ]
print(dim(filt_data))
print(filt_data[1:5, 1:5])
gene_name_filt = gene_name[gene_filter]
rds_file = paste(method_dir, sub(".loom", paste0("_scImpute_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression, ".rds"), get_last_element(unlist(strsplit(loom_path, "/", fixed = T))), fixed = TRUE), sep = "/")
saveRDS(filt_data, rds_file)
out_dir = sub(".rds", "/scImpute/", rds_file, fixed = T)
scimpute(# full path to raw count matrix
  count_path = rds_file, 
  infile = "rds",           # format of input file
  outfile = "rds",          # format of output file
  out_dir = out_dir,           # full path to output directory
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = 1,
  ncores = use_core)              # number of cores used in parallel computation   3g/cpu
saveRDS(filt_data, rds_file)
impute_path = paste(out_dir, "scimpute_count.rds", sep = "/")
imputed_mat = readRDS(impute_path)
row.names(imputed_mat) = gene_name_filt
colnames(imputed_mat) = cell_id
saveRDS(imputed_mat, impute_path)
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(length(gene_name_filt), length(cell_id)),
                storage.mode = "double",
                chunk=c(length(gene_name_filt), 1))
h5write(cell_id, output_h5,"cell_id")
h5write(gene_name_filt, output_h5,"gene_name")
h5write(imputed_mat, file=output_h5, name="imputation")
print(Sys.time() - starttime)

