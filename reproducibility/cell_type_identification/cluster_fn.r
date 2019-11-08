args = commandArgs(trailingOnly=TRUE)
if(length(args) < 2){
  stop("R --slave < this_code.r --args <gene-bc imputation h5/loom> <gene-bc feature h5/loom>")
}
#args = c("/home/yuanhao/data/fn/pbmc3k/pbmc3k_filtered.loom", "")
source("/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/generic_functions.r")
pca_dim = 50 #  useless when feature is provided
res = 1.4
if(args[2] != ""){
  cell.embeddings = t(readh5_loom(args[2], is_feature = TRUE))
  output_dir = paste(c(delete_last_element(unlist(strsplit(args[1], "/", fixed = T))), "cluster_evaluation", get_last_element(delete_last_element(unlist(strsplit(args[1], "[/\\.]", perl = T)))), "feature"), collapse = "/")
}else{
  output_dir = paste(c(delete_last_element(unlist(strsplit(args[1], "/", fixed = T))), "cluster_evaluation", get_last_element(delete_last_element(unlist(strsplit(args[1], "[/\\.]", perl = T)))), "pca"), collapse = "/")
}
gene_bc_imputation = readh5_loom(args[1])
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
this_data <- FindNeighbors(object = this_data, dims = 1:pca_dim, force.recalc = TRUE)
this_data <- FindClusters(object = this_data, resolution = res)
tmp_plot = TSNEPlot(this_data)
ggsave(plot = tmp_plot, paste0(output_dir, "/tsne_cluster.pdf"), height = 8, width = 11)
this_metadata = this_data@meta.data
write.table(this_metadata, paste0(output_dir, "/metadata.txt"), row.names = T, col.names = T, quote = F)
this_markers <- FindAllMarkers(object = this_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(this_markers, paste0(output_dir, "/markers.txt"), row.names = T, col.names = T, quote = F)
write.table(this_data@active.ident, paste0(output_dir, "/cluster_cell_type.txt"), row.names = T, col.names = F, quote = F)
cluster_cell_type = this_data@active.ident
save(this_markers, this_metadata, cluster_cell_type, file = paste0(output_dir, "/tmp.rdata"))

