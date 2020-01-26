args<-commandArgs(trailingOnly=TRUE)
if(length(args) != 4){
  stop("R --slave < this_code.r --args <loom> <cell/gene> <min_expressed_cell> <min_expressed_cell_average_expression>")
}
source("/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/utilities.r")
library(VIPER)
starttime = Sys.time()
loom_path = args[1]
dir_path = paste(delete_last_element(unlist(strsplit(loom_path, "/", fixed = T))), collapse = "/")
output_dir = paste0(dir_path, "/imputation")
method_dir = paste0(dir_path, "/imputation/VIPER_tmp")
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
if(args[2] == "cell"){
  output_h5 = paste(output_dir, sub(".loom", paste0("_VIPER_cell_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression, ".hdf5"), get_last_element(unlist(strsplit(loom_path, "/", fixed = T))), fixed = TRUE), sep = "/")
}else if(args[2] == "gene"){
  output_h5 = paste(output_dir, sub(".loom", paste0("_VIPER_gene_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression, ".hdf5"), get_last_element(unlist(strsplit(loom_path, "/", fixed = T))), fixed = TRUE), sep = "/")
}else{
  stop("args[2] must be cell or gene")
}
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
rownames(filt_data) = paste0("Gene_", 1:length(gene_name_filt))
colnames(filt_data) = paste0("Cell_", 1:length(cell_id))
if(args[2] == "cell"){
  res <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1, report = FALSE, outdir = out_dir, prefix = NULL)
  rownames(res$imputed_log) = gene_name_filt
  rownames(res$imputed) = gene_name_filt
  colnames(res$imputed_log) = cell_id
  colnames(res$imputed) = cell_id
  imputed_mat = res$imputed
  out_dir = paste(method_dir, sub(".loom", paste0("_VIPER_cell_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression), get_last_element(unlist(strsplit(loom_path, "/", fixed = T))), fixed = TRUE), sep = "/")
  output_rds_name = sub(".loom", paste0("_VIPER_cell_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression, ".rds"), get_last_element(unlist(strsplit(loom_path, "/", fixed = T))), fixed = TRUE)
}else{
  res = VIPER(t(filt_data), num = 1000, percentage.cutoff = 0.5, minbool = FALSE, alpha = 0.5, report = FALSE, outdir = out_dir, prefix = NULL)
  colnames(res$imputed_log) = gene_name_filt
  colnames(res$imputed) = gene_name_filt
  rownames(res$imputed_log) = cell_id
  rownames(res$imputed) = cell_id
  imputed_mat = t(res$imputed)
  out_dir = paste(method_dir, sub(".loom", paste0("_VIPER_gene_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression), get_last_element(unlist(strsplit(loom_path, "/", fixed = T))), fixed = TRUE), sep = "/")
  output_rds_name = sub(".loom", paste0("_VIPER_gene_mc_", min_expressed_cell, "_mce_", min_expressed_cell_average_expression, ".rds"), get_last_element(unlist(strsplit(loom_path, "/", fixed = T))), fixed = TRUE)
}
dir.create(out_dir, showWarnings = F, recursive = T)
saveRDS(res, paste(method_dir, output_rds_name, sep = "/"))
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

