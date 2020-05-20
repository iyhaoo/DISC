setwd("/home/yuanhao/github_repositories/DISC/reproducibility")
utilities_path = "./source/utilities.r"
source(utilities_path)
#### STEP 1
#Here, we use BONE_MARROW dataset. The detail information of this dataset can be seen at https://raw.githack.com/iyhaoo/DISC/master/reproducibility/data_preparation_and_imputation/data_preprocessing_BONE_MARROW.nb.html.</br>
#  We used the raw data after gene selection for cell identification.
gene_bc_mat = readh5_loom("./data/BONE_MARROW/raw.loom")
gene_bc_filt = gene_bc_mat[gene_selection(gene_bc_mat, 10), ]
dim(gene_bc_filt) #  13813, 6939
used_genes = rownames(gene_bc_filt)
output_dir = "./results/BONE_MARROW"
dir.create(output_dir, showWarnings = F, recursive = T)
#### STEP 2
#Following this script (https://github.com/Winnie09/imputationBenchmark/blob/master/data/code/process/07_hca_assign_celltype.R), we use the bulk-sequence data (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74246) of 13 normal hematopoietic cell types and 3 acute myeloid leukemia cell types for cell identification, the file is downloaded from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE74246&format=file&file=GSE74246%5FRNAseq%5FAll%5FCounts%2Etxt%2Egz.
if(!file.exists("./data/BONE_MARROW/cell_type.rds")){
  library(scran)
  scran_normalization = function(gene_bc_mat){
    #  The source of this function is https://github.com/Winnie09/imputationBenchmark/blob/master/data/code/process/06_make_hca_MantonBM6.R
    dimnames_gene_bc_mat = dimnames(gene_bc_mat)
    dimnames(gene_bc_mat) = list()
    sce = SingleCellExperiment(list(counts = gene_bc_mat))
    no_cores = max(c(detectCores() - 1, 1))
    if(ncol(gene_bc_mat) < 21){
      sce = computeSumFactors(sce, BPPARAM = MulticoreParam(workers = no_cores), sizes = c(5, 10, 15, 20))
    } else {
      sce = computeSumFactors(sce, BPPARAM = MulticoreParam(workers = no_cores))  
    }
    sf = sizeFactors(sce)
    dimnames(gene_bc_mat) = dimnames_gene_bc_mat
    return(log2(sweep(gene_bc_mat, 2, sf, "/") + 1))
  }
  scalematrix = function(data){
    cm = rowMeans(data)
    csd = apply(data, 1, sd)
    (data - cm) / csd
  }
  corfunc = function(m1, m2){
    scalematrix(t(m1)) %*% t(scalematrix(t(m2))) / (nrow(m1) - 1)
  }
  gene_bulk_all = as.matrix(read.table("./data/BONE_MARROW/original_data/GSE74246_RNAseq_All_Counts.txt.gz", header = T, row.names = 1))
  gene_bulk_mat = gene_bulk_all[, grep("^X", colnames(gene_bulk_all))]
  #  Use annotation information
  gz_path = "./data/hg19/Homo_sapiens.GRCh37.87.gtf.gz"
  annotation_mat = get_map(gz_path)
  tgngl = tapply(annotation_mat[, "gene_length"] / 1000, annotation_mat[, "gene_name"], max)
  gngl = as.vector(tgngl)
  names(gngl) = names(tgngl)
  gene_bulk_filt = gene_bulk_mat[row.names(gene_bulk_mat) %in% names(gngl),]
  gene_bulk_norm = sweep(gene_bulk_filt / gngl[rownames(gene_bulk_filt)], 2, colSums(gene_bulk_filt) / 1e6, "/")
  bulk_data = log2(gene_bulk_norm[rowSums(gene_bulk_norm) > 0,] + 1)
  bulk_cell_type = sapply(colnames(bulk_data), function(x){
    strsplit(x,"\\.")[[1]][2]
  }, USE.NAMES = F)
  sc_data = scran_normalization(gene_bc_filt)
  rownames(sc_data) = sub(".*:", "", rownames(gene_bc_filt))
  used_genes = intersect(rownames(bulk_data), rownames(sc_data))
  bulk_filt = bulk_data[used_genes, ]
  sc_filt = sc_data[used_genes, ]
  #  The expression level for each cell type in bulk sequencing
  bulk_mean = sapply(unique(bulk_cell_type),function(x) {
    rowMeans(bulk_filt[ , bulk_cell_type == x])
  })
  #  Find 100 top postive differentially expressed genes for each celltype pair.
  DEG_list = list()
  top_number = 100
  unique_celltype_pairs = combn(ncol(bulk_mean), 2)
  for(ii in seq(ncol(unique_celltype_pairs))){
    celltype_1 = colnames(bulk_mean)[unique_celltype_pairs[1, ii]]
    celltype_2 = colnames(bulk_mean)[unique_celltype_pairs[2, ii]]
    sort_result = sort(bulk_mean[ , celltype_1] - bulk_mean[ , celltype_2], decreasing = FALSE)
    DEG_list[[paste(celltype_2, celltype_1, sep = "-")]] = names(sort_result)[seq(top_number)]
    DEG_list[[paste(celltype_1, celltype_2, sep = "-")]] = names(sort_result)[seq(from = length(sort_result), to = length(sort_result) - (top_number - 1))]
  }
  #  Calculate the mean expression of these top-gene combinations across cell types (bulk) or cells (single-cell).
  expression_mean_function = function(gene_bc_norm, DEG_list){
    return(t(sapply(DEG_list, function(x){
      colMeans(gene_bc_norm[x, ])
    })))
  }
  bulk_DEG_expression_mean = expression_mean_function(bulk_mean, DEG_list)
  sc_DEG_expression_mean = expression_mean_function(sc_filt, DEG_list)
  #  Calculate the expression variation of these top-gene combinations across cell types (bulk) or cells (single-cell).
  expression_variation_function = function(x){
    return((x - rowMeans(x)) / apply(x, 1, sd))
  }
  bulk_DEG_expression_variation = expression_variation_function(bulk_DEG_expression_mean)
  sc_DEG_expression_variation = expression_variation_function(sc_DEG_expression_mean)
  #  Each top-gene combination correspond a cell type.
  bulk_DEG_combination_rank = apply(bulk_DEG_expression_variation, 2, rank)
  sc_DEG_combination_rank = apply(sc_DEG_expression_variation, 2, rank)
  #  Cell type identification.
  maxcorcut = 0.6
  difcorcut = 0
  cormat = corfunc(sc_DEG_combination_rank, bulk_DEG_combination_rank)
  maxcor = apply(cormat, 1, max)
  max2cor = apply(cormat, 1, function(x){
    sort(x, decreasing = T)[2]
  })
  cell_type = colnames(cormat)[apply(cormat, 1, which.max)]
  cell_type[maxcor < maxcorcut] = NA
  cell_type[maxcor - max2cor < difcorcut] = NA
  names(cell_type) = colnames(sc_data)
  saveRDS(cell_type, "./data/BONE_MARROW/cell_type.rds")
}else{
  cell_type = readRDS("./data/BONE_MARROW/cell_type.rds")
}
print("Cell Type ... OK!")
### Trajectory evaluation
#After cell identification, we evaluate the trajectory performance using monocle following these scripts(https://github.com/Winnie09/imputationBenchmark/blob/93f27e890a86fdc732257a4036bf38a52faf9f33/trajectory/code/hca/monocle2/01_get_score.R, https://github.com/Winnie09/imputationBenchmark/blob/93f27e890a86fdc732257a4036bf38a52faf9f33/trajectory/code/hca/tscan/01_get_score.R).
###  monocle2
library(monocle)
get_cds_monocle2 = function(gene_bc_mat, cell_type, cell_level_df){
  #  Firstly make a new CDS and use DDRTree for dimension reduction.
  pd = new("AnnotatedDataFrame", data = data.frame(row.names = colnames(gene_bc_mat), cell = colnames(gene_bc_mat)))
  fd = new("AnnotatedDataFrame", data = data.frame(row.names = rownames(gene_bc_mat), gene_short_name = rownames(gene_bc_mat)))
  cds = newCellDataSet(gene_bc_mat, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  print("Reducing dimension...")
  cds = reduceDimension(cds)
  cds = orderCells(cds)
  print("Looking for the root state...")
  cell_states = as.numeric(as.character(pData(cds)$State))
  names(cell_states) = colnames(gene_bc_mat)
  unique_states = unique(cell_states)
  checkroot = sapply(unique_states, function(x){
    cds = orderCells(cds, root_state = x)
    return(length(cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell))
  })
  unique_states_filt = unique_states[checkroot > 0]
  root_state = unique_states_filt[which.min(sapply(unique_states_filt, function(x){
    return(mean(cell_level_df[cell_type[names(cell_states)[cell_states == x]], "level"], na.rm = T))
  }))]
  #  Reorder the cells.
  return(orderCells(cds, root_state = root_state))
}
get_order_pairs_monocle2 = function(cds, cell_type){
  df = data.frame(pData(cds),stringsAsFactors = F)
  df = df[df$cell %in% names(cell_type)[!is.na(cell_type)], ]
  df = df[order(df$Pseudotime), ]
  ordered_cell = as.character(df[, 1])
  ordered_celltype = cell_type[ordered_cell]
  index_pair = combn(length(ordered_cell), 2)
  if(min(index_pair[2,] - index_pair[1,]) < 0){
    stop("index_pair error")
  }
  cell_order = sprintf("%s_%s",ordered_celltype[index_pair[1, ]], ordered_celltype[index_pair[2, ]])
  cell_pair = sprintf("%s_%s",ordered_cell[index_pair[1, ]], ordered_cell[index_pair[2, ]])
  names(cell_order) = cell_pair
  return(cell_order)
}
acc_function = function(cell_order, correctorder, wrongorder, compare_pairs = NULL){
  if(is.null(compare_pairs)){
    score = c(sum(cell_order %in% correctorder), sum(cell_order %in% wrongorder))
  }else{
    compare_order = cell_order[names(cell_order) %in% compare_pairs]
    score = c(sum(compare_order %in% correctorder), sum(compare_order %in% wrongorder))
  }
  return(score[1] / sum(score))
}
get_compare_pairs = function(pair_list){
  flag = TRUE
  compare_pairs = NULL
  for(ii in names(pair_list)){
    this_pair = names(pair_list[[ii]])
    cat(ii, " - ", length(this_pair), " pairs\n")
    this_unique_pair = unique(this_pair)
    this_reverse_pair = apply(stri_list2matrix(stri_split_fixed(this_unique_pair, "_"))[c(2, 1), ], 2, stri_join, collapse = "_")
    this_candidate_pair = union(this_unique_pair, this_reverse_pair)
    if(flag){
      compare_pairs = this_candidate_pair
      flag = FALSE
    }else{
      compare_pairs = intersect(compare_pairs, this_candidate_pair)
    }
    
  }
  cat("Use ", length(compare_pairs), " pairs to compare.\n")
  return(compare_pairs)
}
###
cell_level_df = data.frame(level = c(1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5), 
                           immunepath = c(1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0), 
                           monopath = c(1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0), 
                           erypath = c(1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1), 
                           stringsAsFactors = F)
row.names(cell_level_df) = c("HSC", "MPP", "LMPP", "CMP", "CLP", "GMP", "MEP", "Bcell", "CD4Tcell", "CD8Tcell", "NKcell", "Mono", "Ery")
path_name = c("immunepath", "monopath", "erypath")
order_list = list(correct = list(), wrong = list(), cell_type = list())
for(ii in path_name){
  path_celltype = rownames(cell_level_df)[cell_level_df[, ii] == 1]
  order_list[["cell_type"]][[ii]] = path_celltype
  cell_type_pair = as.matrix(apply(expand.grid(path_celltype, path_celltype), 2, as.character))
  cell_type_pair = cell_type_pair[cell_type_pair[, 1] != cell_type_pair[, 2], ]
  correct_mat = cell_type_pair[cell_level_df[cell_type_pair[, 1], "level"] < cell_level_df[cell_type_pair[, 2], "level"], ]
  wrong_mat = cell_type_pair[cell_level_df[cell_type_pair[, 1], "level"] > cell_level_df[cell_type_pair[, 2], "level"], ]
  order_list[["correct"]][[ii]] = apply(correct_mat, 1, paste, collapse = "_")
  order_list[["wrong"]][[ii]] = apply(wrong_mat, 1, paste, collapse = "_")
}


###
this_method = "raw"
gene_bc_mat = readh5_loom(paste0("./data/BONE_MARROW/", this_method, ".loom"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/all_path/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells], cell_type, cell_level_df)
cell_order = get_order_pairs_monocle2(cds, cell_type)
acc = acc_function(cell_order, unique(unlist(order_list[["correct"]])), unique(unlist(order_list[["wrong"]])))
p = plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_pseudotime.pdf"), p) 
pData(cds)$CellType = cell_type[used_cells]
p = plot_cell_trajectory(cds, color_by = "CellType")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_celltype.pdf"), p) 
pData(cds)$level = as.character(cell_level_df[cell_type[used_cells], "level"])
p = plot_cell_trajectory(cds, color_by = "level")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_level.pdf"), p)
cat(ii, " - ", length(cell_order), " pairs\n")
cat(ii, " - ", acc, "\n")
result_list = list(cds = cds, cell_order = cell_order, acc = acc)
saveRDS(result_list, paste0(output_dir, "/", this_method, "_monocle2_all_result.rds"))


###
this_method = "DISC"
gene_bc_mat = readh5_loom(paste0("./data/BONE_MARROW/", this_method, ".loom"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/all_path/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells], cell_type, cell_level_df)
cell_order = get_order_pairs_monocle2(cds, cell_type)
acc = acc_function(cell_order, unique(unlist(order_list[["correct"]])), unique(unlist(order_list[["wrong"]])))
p = plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_pseudotime.pdf"), p) 
pData(cds)$CellType = cell_type[used_cells]
p = plot_cell_trajectory(cds, color_by = "CellType")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_celltype.pdf"), p) 
pData(cds)$level = as.character(cell_level_df[cell_type[used_cells], "level"])
p = plot_cell_trajectory(cds, color_by = "level")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_level.pdf"), p)
cat(ii, " - ", length(cell_order), " pairs\n")
cat(ii, " - ", acc, "\n")
result_list = list(cds = cds, cell_order = cell_order, acc = acc)
saveRDS(result_list, paste0(output_dir, "/", this_method, "_monocle2_all_result.rds"))



###
this_method = "scVI"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/all_path/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells], cell_type, cell_level_df)
cell_order = get_order_pairs_monocle2(cds, cell_type)
acc = acc_function(cell_order, unique(unlist(order_list[["correct"]])), unique(unlist(order_list[["wrong"]])))
p = plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_pseudotime.pdf"), p) 
pData(cds)$CellType = cell_type[used_cells]
p = plot_cell_trajectory(cds, color_by = "CellType")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_celltype.pdf"), p) 
pData(cds)$level = as.character(cell_level_df[cell_type[used_cells], "level"])
p = plot_cell_trajectory(cds, color_by = "level")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_level.pdf"), p)
cat(ii, " - ", length(cell_order), " pairs\n")
cat(ii, " - ", acc, "\n")
result_list = list(cds = cds, cell_order = cell_order, acc = acc)
saveRDS(result_list, paste0(output_dir, "/", this_method, "_monocle2_all_result.rds"))



###
this_method = "DCA"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/all_path/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells], cell_type, cell_level_df)
cell_order = get_order_pairs_monocle2(cds, cell_type)
acc = acc_function(cell_order, unique(unlist(order_list[["correct"]])), unique(unlist(order_list[["wrong"]])))
p = plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_pseudotime.pdf"), p) 
pData(cds)$CellType = cell_type[used_cells]
p = plot_cell_trajectory(cds, color_by = "CellType")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_celltype.pdf"), p) 
pData(cds)$level = as.character(cell_level_df[cell_type[used_cells], "level"])
p = plot_cell_trajectory(cds, color_by = "level")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_level.pdf"), p)
cat(ii, " - ", length(cell_order), " pairs\n")
cat(ii, " - ", acc, "\n")
result_list = list(cds = cds, cell_order = cell_order, acc = acc)
saveRDS(result_list, paste0(output_dir, "/", this_method, "_monocle2_all_result.rds"))


###
this_method = "scScope"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/all_path/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells], cell_type, cell_level_df)
cell_order = get_order_pairs_monocle2(cds, cell_type)
acc = acc_function(cell_order, unique(unlist(order_list[["correct"]])), unique(unlist(order_list[["wrong"]])))
p = plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_pseudotime.pdf"), p) 
pData(cds)$CellType = cell_type[used_cells]
p = plot_cell_trajectory(cds, color_by = "CellType")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_celltype.pdf"), p) 
pData(cds)$level = as.character(cell_level_df[cell_type[used_cells], "level"])
p = plot_cell_trajectory(cds, color_by = "level")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_level.pdf"), p)
cat(ii, " - ", length(cell_order), " pairs\n")
cat(ii, " - ", acc, "\n")
result_list = list(cds = cds, cell_order = cell_order, acc = acc)
saveRDS(result_list, paste0(output_dir, "/", this_method, "_monocle2_all_result.rds"))


###
this_method = "MAGIC"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/all_path/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells], cell_type, cell_level_df)
cell_order = get_order_pairs_monocle2(cds, cell_type)
acc = acc_function(cell_order, unique(unlist(order_list[["correct"]])), unique(unlist(order_list[["wrong"]])))
p = plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_pseudotime.pdf"), p) 
pData(cds)$CellType = cell_type[used_cells]
p = plot_cell_trajectory(cds, color_by = "CellType")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_celltype.pdf"), p) 
pData(cds)$level = as.character(cell_level_df[cell_type[used_cells], "level"])
p = plot_cell_trajectory(cds, color_by = "level")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_level.pdf"), p)
cat(ii, " - ", length(cell_order), " pairs\n")
cat(ii, " - ", acc, "\n")
result_list = list(cds = cds, cell_order = cell_order, acc = acc)
saveRDS(result_list, paste0(output_dir, "/", this_method, "_monocle2_all_result.rds"))



###
this_method = "VIPER"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/all_path/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells], cell_type, cell_level_df)
cell_order = get_order_pairs_monocle2(cds, cell_type)
acc = acc_function(cell_order, unique(unlist(order_list[["correct"]])), unique(unlist(order_list[["wrong"]])))
p = plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_pseudotime.pdf"), p) 
pData(cds)$CellType = cell_type[used_cells]
p = plot_cell_trajectory(cds, color_by = "CellType")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_celltype.pdf"), p) 
pData(cds)$level = as.character(cell_level_df[cell_type[used_cells], "level"])
p = plot_cell_trajectory(cds, color_by = "level")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_level.pdf"), p)
cat(ii, " - ", length(cell_order), " pairs\n")
cat(ii, " - ", acc, "\n")
result_list = list(cds = cds, cell_order = cell_order, acc = acc)
saveRDS(result_list, paste0(output_dir, "/", this_method, "_monocle2_all_result.rds"))




###
this_method = "scImpute"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/all_path/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells], cell_type, cell_level_df)
cell_order = get_order_pairs_monocle2(cds, cell_type)
acc = acc_function(cell_order, unique(unlist(order_list[["correct"]])), unique(unlist(order_list[["wrong"]])))
p = plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_pseudotime.pdf"), p) 
pData(cds)$CellType = cell_type[used_cells]
p = plot_cell_trajectory(cds, color_by = "CellType")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_celltype.pdf"), p) 
pData(cds)$level = as.character(cell_level_df[cell_type[used_cells], "level"])
p = plot_cell_trajectory(cds, color_by = "level")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_level.pdf"), p)
cat(ii, " - ", length(cell_order), " pairs\n")
cat(ii, " - ", acc, "\n")
result_list = list(cds = cds, cell_order = cell_order, acc = acc)
saveRDS(result_list, paste0(output_dir, "/", this_method, "_monocle2_all_result.rds"))



###
this_method = "DeepImpute"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/all_path/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells], cell_type, cell_level_df)
cell_order = get_order_pairs_monocle2(cds, cell_type)
acc = acc_function(cell_order, unique(unlist(order_list[["correct"]])), unique(unlist(order_list[["wrong"]])))
p = plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_pseudotime.pdf"), p) 
pData(cds)$CellType = cell_type[used_cells]
p = plot_cell_trajectory(cds, color_by = "CellType")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_celltype.pdf"), p) 
pData(cds)$level = as.character(cell_level_df[cell_type[used_cells], "level"])
p = plot_cell_trajectory(cds, color_by = "level")
ggsave(paste0(this_output_dir, "/", this_method, "_monocle2_", ii, "_level.pdf"), p)
cat(ii, " - ", length(cell_order), " pairs\n")
cat(ii, " - ", acc, "\n")
result_list = list(cds = cds, cell_order = cell_order, acc = acc)
saveRDS(result_list, paste0(output_dir, "/", this_method, "_monocle2_all_result.rds"))







###
method_names = c("raw", "DISC", "scImpute", "VIPER", "MAGIC", "DCA", "DeepImpute", "scScope", "scVI")
pair_list = list()
for(ii in method_names){
  result_list = readRDS(paste0(output_dir, "/", ii, "_monocle2_all_result.rds"))
  pair_list[[ii]] = result_list[["cell_order"]]
  cat(ii, " - ", result_list[["acc"]], "\n")
}
compare_pairs = get_compare_pairs(pair_list)
acc_all = c()
for(ii in names(pair_list)){
  acc_all = c(acc_all, acc_function(pair_list[[ii]], unique(unlist(order_list[["correct"]])), unique(unlist(order_list[["wrong"]])), compare_pairs))
}
names(acc_all) = names(pair_list)
print(acc_all)
pair_list1 = list(DISC = pair_list[["DISC"]], raw = pair_list[["raw"]])
compare_pairs1 = get_compare_pairs(pair_list1)
acc_all1 = c()
for(ii in names(pair_list1)){
  acc_all1 = c(acc_all1, acc_function(pair_list1[[ii]], unique(unlist(order_list[["correct"]])), unique(unlist(order_list[["wrong"]])), compare_pairs1))
}
names(acc_all1) = names(pair_list1)
print(acc_all1)










