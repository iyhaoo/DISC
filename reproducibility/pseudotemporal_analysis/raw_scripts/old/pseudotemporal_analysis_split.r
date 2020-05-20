setwd("/home/yuanhao/github_repositories/DISC/reproducibility")
source_path = "./pseudotemporal_analysis/raw_scripts/pseudotemporal_analysis_source.r"
source(source_path)


###
this_method = "raw"
gene_bc_mat = readh5_loom(paste0("./data/BONE_MARROW/", this_method, ".loom"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
for(ii in path_name){
  used_cells = names(cell_type)[cell_type %in% order_list[["cell_type"]][[ii]]]
  cds = get_cds_monocle2(gene_bc_filt[, used_cells])
  result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = paste0(this_output_dir, "/", ii), type_level = type_level)
  cat(this_method, " - ", ii, "\n")
  print(result_list)
  saveRDS(cds, paste0(this_output_dir, "/", ii, "_cds.rds"))
  saveRDS(result_list, paste0(this_output_dir, "/", ii, "_result_list.rds"))
}


###
this_method = "DISC"
gene_bc_mat = readh5_loom(paste0("./data/BONE_MARROW/", this_method, ".loom"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
for(ii in path_name){
  used_cells = names(cell_type)[cell_type %in% order_list[["cell_type"]][[ii]]]
  cds = get_cds_monocle2(gene_bc_filt[, used_cells])
  result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = paste0(this_output_dir, "/", ii), type_level = type_level)
  cat(this_method, " - ", ii, "\n")
  print(result_list)
  saveRDS(cds, paste0(this_output_dir, "/", ii, "_cds.rds"))
  saveRDS(result_list, paste0(this_output_dir, "/", ii, "_result_list.rds"))
}


###
this_method = "scVI"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
for(ii in path_name){
  used_cells = names(cell_type)[cell_type %in% order_list[["cell_type"]][[ii]]]
  cds = get_cds_monocle2(gene_bc_filt[, used_cells])
  result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = paste0(this_output_dir, "/", ii), type_level = type_level)
  cat(this_method, " - ", ii, "\n")
  print(result_list)
  saveRDS(cds, paste0(this_output_dir, "/", ii, "_cds.rds"))
  saveRDS(result_list, paste0(this_output_dir, "/", ii, "_result_list.rds"))
}


###
this_method = "DCA"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
for(ii in path_name){
  used_cells = names(cell_type)[cell_type %in% order_list[["cell_type"]][[ii]]]
  cds = get_cds_monocle2(gene_bc_filt[, used_cells])
  result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = paste0(this_output_dir, "/", ii), type_level = type_level)
  cat(this_method, " - ", ii, "\n")
  print(result_list)
  saveRDS(cds, paste0(this_output_dir, "/", ii, "_cds.rds"))
  saveRDS(result_list, paste0(this_output_dir, "/", ii, "_result_list.rds"))
}


###
this_method = "scScope"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
for(ii in path_name){
  used_cells = names(cell_type)[cell_type %in% order_list[["cell_type"]][[ii]]]
  cds = get_cds_monocle2(gene_bc_filt[, used_cells])
  result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = paste0(this_output_dir, "/", ii), type_level = type_level)
  cat(this_method, " - ", ii, "\n")
  print(result_list)
  saveRDS(cds, paste0(this_output_dir, "/", ii, "_cds.rds"))
  saveRDS(result_list, paste0(this_output_dir, "/", ii, "_result_list.rds"))
}


###
this_method = "MAGIC"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
for(ii in path_name){
  used_cells = names(cell_type)[cell_type %in% order_list[["cell_type"]][[ii]]]
  cds = get_cds_monocle2(gene_bc_filt[, used_cells])
  result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = paste0(this_output_dir, "/", ii), type_level = type_level)
  cat(this_method, " - ", ii, "\n")
  print(result_list)
  saveRDS(cds, paste0(this_output_dir, "/", ii, "_cds.rds"))
  saveRDS(result_list, paste0(this_output_dir, "/", ii, "_result_list.rds"))
}


###
this_method = "VIPER"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
for(ii in path_name){
  used_cells = names(cell_type)[cell_type %in% order_list[["cell_type"]][[ii]]]
  cds = get_cds_monocle2(gene_bc_filt[, used_cells])
  result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = paste0(this_output_dir, "/", ii), type_level = type_level)
  cat(this_method, " - ", ii, "\n")
  print(result_list)
  saveRDS(cds, paste0(this_output_dir, "/", ii, "_cds.rds"))
  saveRDS(result_list, paste0(this_output_dir, "/", ii, "_result_list.rds"))
}


###
this_method = "scImpute"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
for(ii in path_name){
  used_cells = names(cell_type)[cell_type %in% order_list[["cell_type"]][[ii]]]
  cds = get_cds_monocle2(gene_bc_filt[, used_cells])
  result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = paste0(this_output_dir, "/", ii), type_level = type_level)
  cat(this_method, " - ", ii, "\n")
  print(result_list)
  saveRDS(cds, paste0(this_output_dir, "/", ii, "_cds.rds"))
  saveRDS(result_list, paste0(this_output_dir, "/", ii, "_result_list.rds"))
}


###
this_method = "DeepImpute"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
for(ii in path_name){
  used_cells = names(cell_type)[cell_type %in% order_list[["cell_type"]][[ii]]]
  cds = get_cds_monocle2(gene_bc_filt[, used_cells])
  result_list = get_score_monocle2(cds, cell_type, correct_order = correct_order_all, wrong_order = wrong_order_all, output_dir = paste0(this_output_dir, "/", ii), type_level = type_level)
  cat(this_method, " - ", ii, "\n")
  print(result_list)
  saveRDS(cds, paste0(this_output_dir, "/", ii, "_cds.rds"))
  saveRDS(result_list, paste0(this_output_dir, "/", ii, "_result_list.rds"))
}





##########
method_names = c("raw", "DISC", "scImpute", "VIPER", "MAGIC", "DCA", "DeepImpute", "scScope", "scVI")
result_list = list()
this_output_dir = paste0(output_dir, "/calculate_all_available_rootstates")
for(ii in path_name){
  result_list[[ii]] = list()
}
for(ii in method_names){
  for(jj in path_name){
    result_list[[ii]][[jj]] = readRDS(paste0(this_output_dir, "/", ii, "/", jj, "_result_list.rds"))
    cat(ii, " - ", jj, "\n")
  }
}
print(result_list)
















