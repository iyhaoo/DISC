args = commandArgs(trailingOnly=TRUE)
if(length(args) < 6){
  stop("R --slave < this_code.r --args <gene-bc imputation h5/loom> <gene-bc feature h5/loom> <cell_type.rds> <pca_dim> <res> <gene_name_rds>")
}
#args = c("/home/yuanhao/data/fn/pbmc3k/pbmc3k_filtered.loom", "")
source("/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/generic_functions.r")
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
  cell.embeddings = t(get_gene_bc_mat(args[2], is_feature = TRUE))
  output_dir = paste(c(delete_last_element(unlist(strsplit(args[1], "/", fixed = T))), "cluster_evaluation", get_last_element(delete_last_element(unlist(strsplit(args[1], "[/\\.]", perl = T)))), "feature"), collapse = "/")
}else{
  output_dir = paste(c(delete_last_element(unlist(strsplit(args[1], "/", fixed = T))), "cluster_evaluation", get_last_element(delete_last_element(unlist(strsplit(args[1], "[/\\.]", perl = T)))), "pca"), collapse = "/")
}
gene_bc_imputation = get_gene_bc_mat(args[1])
if(args[6] != ""){
  gene_bc_imputation = gene_bc_imputation[gsub("_", "-", rownames(gene_bc_imputation)) %in% readRDS(args[6]), ] 
}
this_data = CreateSeuratObject(as.data.frame(gene_bc_imputation))
this_data <- NormalizeData(object = this_data)
dir.create(output_dir, showWarnings = F, recursive = T)
if(args[2] != ""){
  this_data@reductions = list()
  this_data@reductions[["pca"]] = Seurat:::DimReduc(cell.embeddings = cell.embeddings, assay.used = "RNA", key = "feature_")
  pca_dim = ncol(cell.embeddings)
}else{
  this_data <- FindVariableFeatures(this_data)
  this_data <- ScaleData(object = this_data)
  this_data <- RunPCA(object = this_data, npcs = pca_dim * 2)
  tmp_plot = ElbowPlot(this_data, ndims = pca_dim * 2)
  ggsave(plot = tmp_plot, paste0(output_dir, "/elbow.pdf"), height = 8, width = 11)
}
this_data <- RunTSNE(this_data, dims = 1:pca_dim, tsne.method = "Rtsne", check_duplicates = FALSE)
this_data = RunUMAP(this_data, dims = 1:pca_dim)
if(args[3] != ""){
  cell_type = readRDS(args[3])[colnames(gene_bc_imputation)]
  this_data@active.ident = factor(cell_type, levels = sort(unique(cell_type)))
  tmp_plot = TSNEPlot(this_data, cells = names(cell_type))
  ggsave(plot = tmp_plot, filename = paste0(output_dir, "/tsne_cell_type.pdf"), height = 8, width = 11)
  tmp_plot = UMAPPlot(this_data)
  ggsave(plot = tmp_plot, paste0(output_dir, "/umap_cell_type.pdf"), height = 8, width = 11)
}
this_data <- FindNeighbors(object = this_data, dims = 1:pca_dim, force.recalc = TRUE)
this_data <- FindClusters(object = this_data, resolution = res)
tmp_plot = TSNEPlot(this_data)
ggsave(plot = tmp_plot, paste0(output_dir, "/tsne_cluster.pdf"), height = 8, width = 11)
tmp_plot = UMAPPlot(this_data)
ggsave(plot = tmp_plot, paste0(output_dir, "/umap_cluster.pdf"), height = 8, width = 11)
this_metadata = this_data@meta.data
write.table(this_metadata, paste0(output_dir, "/metadata.txt"), row.names = T, col.names = T, quote = F)
this_markers <- FindAllMarkers(object = this_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(this_markers, paste0(output_dir, "/markers.txt"), row.names = T, col.names = T, quote = F)
write.table(this_data@active.ident, paste0(output_dir, "/cluster_cell_type.txt"), row.names = T, col.names = F, quote = F)
cluster_cell_type = this_data@active.ident
save(this_markers, this_metadata, cluster_cell_type, file = paste0(output_dir, "/tmp.rdata"))

