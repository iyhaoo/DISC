setwd("/home/yuanhao/github_repositories/DISC/reproducibility")
source_path = "./Down-stream Analysis Improvement/raw_scripts/pseudotemporal_analysis_source.r"
source(source_path)


###
this_method = "raw"
gene_bc_mat = readh5_loom(paste0("./data/BONE_MARROW/", this_method, ".loom"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells])
result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = this_output_dir, type_level = type_level)
print(this_method)
print(result_list)
saveRDS(cds, paste0(this_output_dir, "/all_cds.rds"))
saveRDS(result_list, paste0(this_output_dir, "/all_result_list.rds"))


###
this_method = "DISC"
gene_bc_mat = readh5_loom(paste0("./data/BONE_MARROW/", this_method, ".loom"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells])
result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = this_output_dir, type_level = type_level)
print(this_method)
print(result_list)
saveRDS(cds, paste0(this_output_dir, "/all_cds.rds"))
saveRDS(result_list, paste0(this_output_dir, "/all_result_list.rds"))


###
this_method = "scVI"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells])
result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = this_output_dir, type_level = type_level)
print(this_method)
print(result_list)
saveRDS(cds, paste0(this_output_dir, "/all_cds.rds"))
saveRDS(result_list, paste0(this_output_dir, "/all_result_list.rds"))


###
this_method = "DCA"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells])
result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = this_output_dir, type_level = type_level)
print(this_method)
print(result_list)
saveRDS(cds, paste0(this_output_dir, "/all_cds.rds"))
saveRDS(result_list, paste0(this_output_dir, "/all_result_list.rds"))


###
this_method = "scScope"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells])
result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = this_output_dir, type_level = type_level)
print(this_method)
print(result_list)
saveRDS(cds, paste0(this_output_dir, "/all_cds.rds"))
saveRDS(result_list, paste0(this_output_dir, "/all_result_list.rds"))


###
this_method = "MAGIC"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells])
result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = this_output_dir, type_level = type_level)
print(this_method)
print(result_list)
saveRDS(cds, paste0(this_output_dir, "/all_cds.rds"))
saveRDS(result_list, paste0(this_output_dir, "/all_result_list.rds"))


###
this_method = "VIPER"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells])
result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = this_output_dir, type_level = type_level)
print(this_method)
print(result_list)
saveRDS(cds, paste0(this_output_dir, "/all_cds.rds"))
saveRDS(result_list, paste0(this_output_dir, "/all_result_list.rds"))


###
this_method = "scImpute"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells])
result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = this_output_dir, type_level = type_level)
print(this_method)
print(result_list)
saveRDS(cds, paste0(this_output_dir, "/all_cds.rds"))
saveRDS(result_list, paste0(this_output_dir, "/all_result_list.rds"))


###
this_method = "DeepImpute"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells])
result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = this_output_dir, type_level = type_level)
print(this_method)
print(result_list)
saveRDS(cds, paste0(this_output_dir, "/all_cds.rds"))
saveRDS(result_list, paste0(this_output_dir, "/all_result_list.rds"))





##########
method_names = c("raw", "DISC", "scImpute", "VIPER", "MAGIC", "DCA", "DeepImpute", "scScope", "scVI")
result_list = list()
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates")
for(ii in method_names){
  result_list[[ii]] = readRDS(paste0(this_output_dir, "/", ii, "/all_result_list.rds"))
  print(ii)
}
print(result_list)

sapply(result_list, function(x){
  return(x[which.max(x[, 2]), ])
})

###   revised   ###
setwd("/home/yuanhao/github_repositories/DISC/reproducibility")
source_path = "./Down-stream Analysis Improvement/raw_scripts/pseudotemporal_analysis_source.r"
source(source_path)


method_names = c("raw", "DISC", "scImpute", "VIPER", "MAGIC", "DCA", "DeepImpute", "scScope", "scVI")
output_dir_1 = paste0(output_dir, "/calculate_all_available_rootstates_revised")
result_list = list()
for(ii in method_names){
  cds = readRDS(paste0(output_dir, "/calculate_all_available_rootstates", "/", ii, "/all_cds.rds"))
  this_output_dir = paste0(output_dir_1, "/", ii)
  dir.create(this_output_dir, recursive = T, showWarnings = F)
  result_list[[ii]] = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order_list = wrong_order_list, output_dir = this_output_dir, type_level = type_level)
  saveRDS(result_list[[ii]], paste0(this_output_dir, "/all_result_list.rds"))
}



setwd("/home/yuanhao/github_repositories/DISC/reproducibility")
output_dir = "./results/BONE_MARROW"
method_names = c("raw", "DISC", "scImpute", "VIPER", "MAGIC", "DCA", "DeepImpute", "scScope", "scVI")
result_list = list()
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates_revised")
for(ii in method_names){
  result_list[[ii]] = readRDS(paste0(this_output_dir, "/", ii, "/all_result_list.rds"))
  print(ii)
}
print(result_list)

result_mat = sapply(result_list, function(x){
  return(x[which.max(x[, 1]), ])
})





