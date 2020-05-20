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
ct <- data.frame(cell=names(cell_type),cluster=NA,celltype=cell_type,stringsAsFactors = F)
ctlevel <- data.frame(ct=c('HSC','MPP','LMPP','CMP','CLP','GMP','MEP',"Bcell","CD4Tcell","CD8Tcell",'NKcell','Mono','Ery'),level=c(1,2,3,3,4,4,4,5,5,5,5,5,5),immunepath=c(1,1,1,0,1,0,0,1,1,1,1,0,0),monopath=c(1,1,1,1,0,1,0,0,0,0,0,1,0),erypath=c(1,1,0,1,0,0,1,0,0,0,0,0,1),stringsAsFactors = F)
row.names(ctlevel) <- ctlevel[,1]

correctorder <- wrongorder <- NULL
for(pid in c('immunepath','monopath','erypath')) {
  evct <- ctlevel[ctlevel[,pid]==1,1]
  pair <- expand.grid(evct,evct)
  pair[,1] <- as.character(pair[,1])
  pair[,2] <- as.character(pair[,2])
  pair <- pair[pair[,1]!=pair[,2],]
  corid <- which(ctlevel[pair[,1],'level'] < ctlevel[pair[,2],'level'])
  wroid <- which(ctlevel[pair[,1],'level'] > ctlevel[pair[,2],'level'])
  correctorder <- c(correctorder,sapply(corid,function(si) paste0(pair[si,],collapse = '_')))
  wrongorder <- c(wrongorder,sapply(wroid,function(si) paste0(pair[si,],collapse = '_')))
}
correctorder <- unique(correctorder)
wrongorder <- unique(wrongorder)


scorefunc <- function(cds) {
  if (length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points) > 0) {
    sl <- NULL
    for(i in 1:length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points)){
      tryCatch({cds_reduced <- buildBranchCellDataSet(cds,branch_point=i)},error=function(e) {})
      df = data.frame(pData(cds_reduced),stringsAsFactors = F)
      df <- df[order(df$Pseudotime),]
      sl <- rbind(sl,rowSums(sapply(unique(df$Branch),function(ub) {
        so <- as.character(df[df[,'Branch']==ub,1])
        soct <- ct[match(so,ct[,1]),3]
        eid <- expand.grid(1:length(soct),1:length(soct))
        eid <- eid[eid[,1]<eid[,2],]
        eid <- sprintf('%s_%s',soct[eid[,1]],soct[eid[,2]])
        c(sum(eid %in% correctorder),sum(eid %in% wrongorder))
      })))
    }
    print(sl)
    outpur_sl = sl[which.max(sl[,1]/rowSums(sl)),]
    pair_number = sum(outpur_sl)
    acc = outpur_sl[1] / pair_number
    return(c(pair_number, acc))
  } else {
    return(NA)
  }
}

######
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
result_list = list(cds = cds, acc = acc)
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
result_list = list(cds = cds, acc = acc)
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
result_list = list(cds = cds, acc = acc)
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
result_list = list(cds = cds, acc = acc)
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
result_list = list(cds = cds, acc = acc)
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
result_list = list(cds = cds, acc = acc)
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
result_list = list(cds = cds, acc = acc)
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
result_list = list(cds = cds, acc = acc)
saveRDS(result_list, paste0(output_dir, "/", this_method, "_monocle2_all_result.rds"))



###
this_method = "DeepImpute"
gene_bc_mat = readh5_imputation(paste0("./data/BONE_MARROW/", this_method, ".hdf5"))
gene_bc_filt = gene_bc_mat[used_genes, ]
this_output_dir = paste0(output_dir, "/all_path/", this_method)
dir.create(this_output_dir, recursive = T, showWarnings = F)
used_cells = names(cell_type)[!is.na(cell_type)]
cds = get_cds_monocle2(gene_bc_filt[, used_cells], cell_type, cell_level_df)
acc = get_score_monocle2(cds, cell_type, unique(unlist(order_list[["correct"]])), unique(unlist(order_list[["wrong"]])))
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
result_list = list(cds = cds, acc = acc)
saveRDS(result_list, paste0(output_dir, "/", this_method, "_monocle2_all_result.rds"))







###
method_names = c("raw", "DISC", "scImpute", "VIPER", "MAGIC", "DCA", "DeepImpute", "scScope", "scVI")
pair_list = list()
for(ii in method_names){
  result_list = readRDS(paste0(output_dir, "/", ii, "_monocle2_all_result.rds"))
  pair_num_acc = scorefunc(result_list[["cds"]])
  cat(ii, " - ", pair_num_acc[1], " - ", pair_num_acc[2], "\n")
}














result_list = readRDS(paste0(output_dir, "/", "DISC", "_monocle2_all_result.rds"))
p = plot_cell_trajectory(result_list[["cds"]], color_by = "Pseudotime")
ggsave(paste0(output_dir, "/test_monocle2_pseudotime.pdf"), p) 
p = plot_cell_trajectory(result_list[["cds"]], color_by = "CellType")
ggsave(paste0(output_dir, "/test_monocle2_celltype.pdf"), p) 
p = plot_cell_trajectory(result_list[["cds"]], color_by = "level")
ggsave(paste0(output_dir, "/test_monocle2_level.pdf"), p)
all_branch_points = result_list[["cds"]]@auxOrderingData[[result_list[["cds"]]@dim_reduce_type]]$branch_points
for(i in 1:length(all_branch_points)){
  cds_reduced <- buildBranchCellDataSet(result_list[["cds"]], branch_point = i)
  print(all_branch_points[i])
  pData(result_list[["cds"]])$Branch = pData(cds_reduced)[pData(result_list[["cds"]])$cell, "Branch"]
  p = plot_cell_trajectory(result_list[["cds"]], color_by = "Branch")
  ggsave(paste0(output_dir, "/test_", all_branch_points[i], "_monocle2_level.pdf"), p)
}





method_names = c("raw", "DISC", "scImpute", "VIPER", "MAGIC", "DCA", "DeepImpute", "scScope", "scVI")
for(ii in method_names){
  result_list = readRDS(paste0(output_dir, "/", ii, "_monocle2_all_result.rds"))
  p = plot_cell_trajectory(result_list[["cds"]], color_by = "Pseudotime")
  ggsave(paste0(output_dir, "/", ii, "_monocle2_pseudotime.pdf"), p) 
  p = plot_cell_trajectory(result_list[["cds"]], color_by = "CellType")
  ggsave(paste0(output_dir, "/", ii, "_monocle2_celltype.pdf"), p) 
  p = plot_cell_trajectory(result_list[["cds"]], color_by = "level")
  ggsave(paste0(output_dir, "/", ii, "_monocle2_level.pdf"), p)
  p = plot_cell_trajectory(result_list[["cds"]], color_by = "State")
  ggsave(paste0(output_dir, "/", ii, "_monocle2_state.pdf"), p)
  all_branch_points = result_list[["cds"]]@auxOrderingData[[result_list[["cds"]]@dim_reduce_type]]$branch_points
  for(jj in 1:length(all_branch_points)){
    cds_tmp = result_list[["cds"]]
    cds_reduced = buildBranchCellDataSet(cds_tmp, branch_point = jj)
    df = data.frame(pData(cds_reduced),stringsAsFactors = F)
    cells = as.character(pData(cds_tmp)$cell)
    print(all_branch_points[jj])
    pData(cds_tmp)$Pseudotime = df[cells, "Pseudotime"]
    pData(cds_tmp)$Branch = df[cells, "Branch"]
    pData(cds_tmp)$State = df[cells, "State"]
    p = plot_cell_trajectory(cds_tmp, color_by = "Pseudotime")
    ggsave(paste0(output_dir, "/", ii, "_", all_branch_points[jj], "_monocle2_pseudotime.pdf"), p)
    p = plot_cell_trajectory(cds_tmp, color_by = "Branch")
    ggsave(paste0(output_dir, "/", ii, "_", all_branch_points[jj], "_monocle2_branch.pdf"), p)
    p = plot_cell_trajectory(cds_tmp, color_by = "State")
    ggsave(paste0(output_dir, "/", ii, "_", all_branch_points[jj], "_monocle2_state.pdf"), p)
    cat(ii, " - ", jj, "\n")
  }
}






ii = "raw"
result_list = readRDS(paste0(output_dir, "/", ii, "_monocle2_all_result.rds"))
p = plot_cell_trajectory(result_list[["cds"]], color_by = "Pseudotime")
ggsave(paste0(output_dir, "/", ii, "_monocle2_pseudotime.pdf"), p) 
p = plot_cell_trajectory(result_list[["cds"]], color_by = "CellType")
ggsave(paste0(output_dir, "/", ii, "_monocle2_celltype.pdf"), p) 
p = plot_cell_trajectory(result_list[["cds"]], color_by = "level")
ggsave(paste0(output_dir, "/", ii, "_monocle2_level.pdf"), p)
all_branch_points = result_list[["cds"]]@auxOrderingData[[result_list[["cds"]]@dim_reduce_type]]$branch_points
for(i in 1:length(all_branch_points)){
  cds_reduced <- buildBranchCellDataSet(result_list[["cds"]],branch_point = i)
  print(all_branch_points[i])
  pData(result_list[["cds"]])$Branch = pData(cds_reduced)[pData(result_list[["cds"]])$cell, "Branch"]
  p = plot_cell_trajectory(result_list[["cds"]], color_by = "Branch")
  ggsave(paste0(output_dir, "/", ii, "_", all_branch_points[i], "_monocle2_level.pdf"), p)
}




ii = "DCA"
result_list = readRDS(paste0(output_dir, "/", ii, "_monocle2_all_result.rds"))
p = plot_cell_trajectory(result_list[["cds"]], color_by = "Pseudotime")
ggsave(paste0(output_dir, "/", ii, "_monocle2_pseudotime.pdf"), p) 
p = plot_cell_trajectory(result_list[["cds"]], color_by = "CellType")
ggsave(paste0(output_dir, "/", ii, "_monocle2_celltype.pdf"), p) 
p = plot_cell_trajectory(result_list[["cds"]], color_by = "level")
ggsave(paste0(output_dir, "/", ii, "_monocle2_level.pdf"), p)
all_branch_points = result_list[["cds"]]@auxOrderingData[[result_list[["cds"]]@dim_reduce_type]]$branch_points
for(i in 1:length(all_branch_points)){
  cds_reduced <- buildBranchCellDataSet(result_list[["cds"]],branch_point = i)
  print(all_branch_points[i])
  pData(result_list[["cds"]])$Branch = pData(cds_reduced)[pData(result_list[["cds"]])$cell, "Branch"]
  p = plot_cell_trajectory(result_list[["cds"]], color_by = "Branch")
  ggsave(paste0(output_dir, "/", ii, "_", all_branch_points[i], "_monocle2_level.pdf"), p)
}










get_order_pairs_monocle2 = function(cds, cell_type){
  pair_list = list()
  if(length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points) > 0){
    for(ii in seq(length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points))){
      cds_reduced = buildBranchCellDataSet(cds, branch_point = ii)
      df = data.frame(pData(cds_reduced),stringsAsFactors = F)[names(cell_type), ]
      df = df[order(df$Pseudotime), ]
      compare_order = c()
      for(jj in unique(df$Branch)){
        branch_cell = as.character(df[df$Branch == jj, 1])
        branch_celltype = cell_type[branch_cell]
        index_pair = combn(length(branch_cell), 2)
        if(min(index_pair[2,] - index_pair[1,]) < 0){
          stop("index_pair error")
        }
        branch_order = sprintf('%s_%s',branch_celltype[index_pair[1, ]], branch_celltype[index_pair[2, ]])
        branch_pair = sprintf('%s_%s',branch_cell[index_pair[1, ]], branch_cell[index_pair[2, ]])
        names(branch_order) = branch_pair
        compare_order = c(compare_order, branch_order)
      }
      pair_list[[ii]] = compare_order
    }
  }else{
    print("Trajectory has only 1 branch.")
    df = data.frame(pData(cds),stringsAsFactors = F)[names(cell_type), ]
    df = df[order(df$Pseudotime), ]
    branch_cell = as.character(df[, 1])
    branch_celltype = cell_type[branch_cell]
    index_pair = combn(length(branch_cell), 2)
    if(min(index_pair[2,] - index_pair[1,]) < 0){
      stop("index_pair error")
    }
    branch_order = sprintf('%s_%s',branch_celltype[index_pair[1, ]], branch_celltype[index_pair[2, ]])
    branch_pair = sprintf('%s_%s',branch_cell[index_pair[1, ]], branch_cell[index_pair[2, ]])
    names(branch_order) = branch_pair
    pair_list[[1]] = branch_order
  }
  return(pair_list)
}
acc_function = function(pair_list, correctorder, wrongorder, compare_pairs = NULL){
  acc = 0
  for(ii in pair_list){
    if(is.null(compare_pairs)){
      score = c(sum(ii %in% correctorder), sum(ii %in% wrongorder))
    }else{
      compare_order = ii[intersect(names(ii), compare_pairs)]
      score = c(sum(compare_order %in% correctorder), sum(compare_order %in% wrongorder))
    }
    acc = max(c(acc, score[1] / sum(score)))
  }
  return(acc)
}

cell_type = readRDS("./data/BONE_MARROW/cell_type.rds")
ctlevel <- data.frame(cell_type=c('HSC','MPP','LMPP','CMP','CLP','GMP','MEP',"Bcell","CD4Tcell","CD8Tcell",'NKcell','Mono','Ery'),level=c(1,2,3,3,4,4,4,5,5,5,5,5,5),immunepath=c(1,1,1,0,1,0,0,1,1,1,1,0,0),monopath=c(1,1,1,1,0,1,0,0,0,0,0,1,0),erypath=c(1,1,0,1,0,0,1,0,0,0,0,0,1),stringsAsFactors = F)
row.names(ctlevel) <- ctlevel[,1]
correctorder <- wrongorder <- NULL
for(pid in c('immunepath','monopath','erypath')) {
  evct <- ctlevel[ctlevel[,pid]==1,1]
  pair <- expand.grid(evct,evct)
  pair[,1] <- as.character(pair[,1])
  pair[,2] <- as.character(pair[,2])
  pair <- pair[pair[,1]!=pair[,2],]
  corid <- which(ctlevel[pair[,1],'level'] < ctlevel[pair[,2],'level'])
  wroid <- which(ctlevel[pair[,1],'level'] > ctlevel[pair[,2],'level'])
  correctorder <- c(correctorder,sapply(corid,function(si) paste0(pair[si,],collapse = '_')))
  wrongorder <- c(wrongorder,sapply(wroid,function(si) paste0(pair[si,],collapse = '_')))
}
correctorder <- unique(correctorder)
wrongorder <- unique(wrongorder)
scVI_res = readRDS("/home/yuanhao/github_repositories/DISC/reproducibility/results/BONE_MARROW/scVI_monocle2_all_result.rds")
pair_list = get_order_pairs_monocle2(scVI_res[["cds"]], cell_type)
acc = acc_function(pair_list, correctorder, wrongorder)
print(acc)

DISC_res = readRDS("/home/yuanhao/github_repositories/DISC/reproducibility/results/BONE_MARROW/DISC_monocle2_all_result.rds")
pair_list = get_order_pairs_monocle2(DISC_res[["cds"]], cell_type)
acc = acc_function(pair_list, correctorder, wrongorder)
print(acc)






