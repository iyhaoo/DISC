library(reticulate)
source_python("/home/wucheng/imputation/readloom.py")
readloom = function(loom_path){
   print(loom_path)
   return_list = readloom_py(loom_path)
   gene_bc_mat = as.matrix(return_list[["gene_bc_mat"]])
   gene_bc_mat[gene_bc_mat == -1] = NA
   rownames(gene_bc_mat) = return_list[["gene_name"]]
   colnames(gene_bc_mat) = return_list[["cell_id"]]
   print(dim(gene_bc_mat))
   return(gene_bc_mat)
 }
library(Seurat) 
library(rhdf5)
readh5 = function(h5_path, use_genes=NULL, used_cells=NULL){
  print(h5_path)
  gene_name = h5read(h5_path, "gene_name")
  cell_id = h5read(h5_path, "cell_id")
  if(!is.null(used_cells)){
    cell_grasp_index = which(cell_id %in% used_cells)
  }else{
    cell_grasp_index = c(1: length(cell_id))
  }
  if(!is.null(use_genes)){
    gene_grasp_index = which(gene_name %in% use_genes)
  }else{
    gene_grasp_index = c(1: length(gene_name))
  }
  gene_bc_mat = h5read(h5_path, "imputation", index = list(gene_grasp_index, cell_grasp_index))
  gene_bc_mat[gene_bc_mat == -1] = NA
  rownames(gene_bc_mat) = gene_name[gene_grasp_index]
  colnames(gene_bc_mat) = cell_id[cell_grasp_index]
  print(dim(gene_bc_mat))
  return(gene_bc_mat)
}

#####

pbmc.data = readloom("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_4/imputation/GSM3017261_150000_CNS_nuclei_ds_0.3_SAVER_mc_150_mce_1_resume_dim.loom")
dim(pbmc.data)
rownames(pbmc.data) = as.character(rownames(pbmc.data))
colnames(pbmc.data) = as.character(colnames(pbmc.data))

#celltype <-read.table("/home/wucheng/imputation/split-seq/celltype.tsv",sep = "\t")
celltype1 <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
celltype2 <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")

pbmc.data1 <-pbmc.data[,rownames(celltype1)]
dim(pbmc.data1)
pbmc.data2 <-pbmc.data[,rownames(celltype2)]
dim(pbmc.data2)
output_h5= "/home/wucheng/imputation/split-seq/ds/repeat4/0.3/MAGIC/neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data1), ncol(pbmc.data1)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data1), 1))
h5write(colnames(pbmc.data1), output_h5,"cell_id")
h5write(rownames(pbmc.data1), output_h5,"gene_name")
h5write(pbmc.data1, output_h5, "imputation")

output_h5= "/home/wucheng/imputation/split-seq/ds/repeat4/0.3/MAGIC/non_neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data2), ncol(pbmc.data2)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data2), 1))
h5write(colnames(pbmc.data2), output_h5,"cell_id")
h5write(rownames(pbmc.data2), output_h5,"gene_name")
h5write(pbmc.data2, output_h5, "imputation")
#######
pbmc.data = readloom("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_4/imputation/GSM3017261_150000_CNS_nuclei_ds_0.5_SAVER_mc_150_mce_1_resume_dim.loom")
dim(pbmc.data)
rownames(pbmc.data) = as.character(rownames(pbmc.data))
colnames(pbmc.data) = as.character(colnames(pbmc.data))

#celltype <-read.table("/home/wucheng/imputation/split-seq/celltype.tsv",sep = "\t")
celltype1 <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
celltype2 <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")

pbmc.data1 <-pbmc.data[,rownames(celltype1)]
dim(pbmc.data1)
pbmc.data2 <-pbmc.data[,rownames(celltype2)]
dim(pbmc.data2)
output_h5= "/home/wucheng/imputation/split-seq/ds/repeat5/0.5/SAVER/neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data1), ncol(pbmc.data1)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data1), 1))
h5write(colnames(pbmc.data1), output_h5,"cell_id")
h5write(rownames(pbmc.data1), output_h5,"gene_name")
h5write(pbmc.data1, output_h5, "imputation")

output_h5= "/home/wucheng/imputation/split-seq/ds/repeat5/0.5/SAVER/non_neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data2), ncol(pbmc.data2)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data2), 1))
h5write(colnames(pbmc.data2), output_h5,"cell_id")
h5write(rownames(pbmc.data2), output_h5,"gene_name")
h5write(pbmc.data2, output_h5, "imputation")

###
pbmc.data = readloom("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_4/imputation/GSM3017261_150000_CNS_nuclei_ds_0.7_SAVER_mc_150_mce_1_resume_dim.loom")
dim(pbmc.data)
rownames(pbmc.data) = as.character(rownames(pbmc.data))
colnames(pbmc.data) = as.character(colnames(pbmc.data))

#celltype <-read.table("/home/wucheng/imputation/split-seq/celltype.tsv",sep = "\t")
celltype1 <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
celltype2 <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")

pbmc.data1 <-pbmc.data[,rownames(celltype1)]
dim(pbmc.data1)
pbmc.data2 <-pbmc.data[,rownames(celltype2)]
dim(pbmc.data2)
output_h5= "/home/wucheng/imputation/split-seq/ds/repeat5/0.7/SAVER/neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data1), ncol(pbmc.data1)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data1), 1))
h5write(colnames(pbmc.data1), output_h5,"cell_id")
h5write(rownames(pbmc.data1), output_h5,"gene_name")
h5write(pbmc.data1, output_h5, "imputation")

output_h5= "/home/wucheng/imputation/split-seq/ds/repeat5/0.7/SAVER/non_neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data2), ncol(pbmc.data2)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data2), 1))
h5write(colnames(pbmc.data2), output_h5,"cell_id")
h5write(rownames(pbmc.data2), output_h5,"gene_name")
h5write(pbmc.data2, output_h5, "imputation")