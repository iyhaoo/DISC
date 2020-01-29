args = commandArgs(trailingOnly=TRUE)
if(length(args) < 7){
  stop("R --slave < this_code.r --args <gene-bc imputation h5/loom> <gene-bc feature h5/loom> <cell_type.rds> <pca_dim> <res> <gene_name_rds> <mapping function>")
}
source("/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/utilities.r")
if(args[4] == ""){
  pca_dim = 50  #  useless when feature is provided
}else{
  pca_dim = as.integer(args[4])
}
if(args[5] == ""){
  res = 1.4
}else{
  res = as.numeric(args[5])
}
if(args[2] != ""){
  feature_bc_mat = readh5_loom(args[2], is_feature = TRUE)
  output_dir = paste(c(delete_last_element(unlist(strsplit(args[1], "/", fixed = T))), "cluster_evaluation_0.3", get_last_element(delete_last_element(unlist(strsplit(args[1], "[/\\.]", perl = T)))), "feature"), collapse = "/")
}else{
  output_dir = paste(c(delete_last_element(unlist(strsplit(args[1], "/", fixed = T))), "cluster_evaluation_0.3", get_last_element(delete_last_element(unlist(strsplit(args[1], "[/\\.]", perl = T)))), "pca"), collapse = "/")
  feature_bc_mat = NULL
}
gene_bc_mat = get_gene_bc_mat(args[1])
if(args[6] != ""){
  gene_bc_mat = gene_bc_mat[gsub("_", "-", rownames(gene_bc_mat)) %in% readRDS(args[6]), ] 
}
if(args[7] != ""){
  switch(args[7],
         "pbmc"={
           mapping_function = cell_type_identification_pbmc
         },
         "retina"={
           mapping_function = cell_type_identification_retina
         },{
           stop("no such mapping")
         })
  seurat_classification(gene_bc_mat = gene_bc_mat, feature_bc_mat = feature_bc_mat, cell_type = readRDS(args[3]),
                        output_dir = output_dir, pca_dim = pca_dim, res = res, min_pct = 0.25, show_plots = F, cell_type_identification_fun = mapping_function)
}





