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
pbmc.data = readloom("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_1/GSM3017261_150000_CNS_nuclei_ds_0.3.loom")
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
output_h5= "/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data1), ncol(pbmc.data1)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data1), 1))
h5write(colnames(pbmc.data1), output_h5,"cell_id")
h5write(rownames(pbmc.data1), output_h5,"gene_name")
h5write(pbmc.data1, output_h5, "imputation")

output_h5= "/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data2), ncol(pbmc.data2)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data2), 1))
h5write(colnames(pbmc.data2), output_h5,"cell_id")
h5write(rownames(pbmc.data2), output_h5,"gene_name")
h5write(pbmc.data2, output_h5, "imputation")
############
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
library(Seurat) 

pbmc.data = readh5("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous.h5")
rownames(pbmc.data) = as.character(rownames(pbmc.data))
colnames(pbmc.data) = as.character(colnames(pbmc.data))
dim(pbmc.data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "neurous",min.cells = 10)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
##this_data@active.ident = cell_type_factor
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc),npcs = 100)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/PCA.pdf")
DimPlot(object = pbmc, reduction = "pca")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/PC_score.pdf",width=15,height=10)
ElbowPlot(object = pbmc, ndims = 100)
dev.off()
pbmc <- FindNeighbors(object = pbmc, dims = 1:50)
pbmc <- FindClusters(object = pbmc, resolution = 1.4)
pbmc <- RunUMAP(object = pbmc, dims = 1:50)
pbmc <- RunTSNE(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/UMAP.pdf")
DimPlot(object = pbmc, reduction = "umap")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/tsne.pdf")
DimPlot(object = pbmc, reduction = "tsne")
dev.off()
write.table(pbmc@meta.data,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data.txt",row.names=TRUE,col.names=TRUE)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pbmc.markers,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/markers.txt",row.names=TRUE,col.names=TRUE)
saveRDS(pbmc, file = "/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/data.rds")
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
pbmc@meta.data$res2 <-celltype[,2]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/tsne1.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res2")
dev.off()
pbmc@meta.data$res3 <-celltype[,3]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/tsne2.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res3")
dev.off()

data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/markers.txt",header=T,row.names=1)
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))
Others <-setdiff(A,hebing)
Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/Accuracy.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/Recall.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:10,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)
Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/Accuracy_10.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/Recall_10.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:20,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])






Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)

Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/Accuracy_20.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/Recall_20.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:50,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)

Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/Accuracy_50.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/neurous/Recall_50.txt",row.names=T,col.names=T,sep="\t")

#####
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
library(Seurat) 

pbmc.data = readh5("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous.h5")
rownames(pbmc.data) = as.character(rownames(pbmc.data))
colnames(pbmc.data) = as.character(colnames(pbmc.data))
dim(pbmc.data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "neurous",min.cells = 10)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
##this_data@active.ident = cell_type_factor
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc),npcs = 100)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/PCA.pdf")
DimPlot(object = pbmc, reduction = "pca")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/PC_score.pdf",width=15,height=10)
ElbowPlot(object = pbmc, ndims = 100)
dev.off()
pbmc <- FindNeighbors(object = pbmc, dims = 1:50)
pbmc <- FindClusters(object = pbmc, resolution = 1.4)
pbmc <- RunUMAP(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/UMAP.pdf")
DimPlot(object = pbmc, reduction = "umap")
dev.off()
pbmc <- RunTSNE(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/tsne.pdf")
DimPlot(object = pbmc, reduction = "tsne")
dev.off()
write.table(pbmc@meta.data,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.data.txt",row.names=TRUE,col.names=TRUE)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pbmc.markers,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/markers.txt",row.names=TRUE,col.names=TRUE)
saveRDS(pbmc, file = "/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/data.rds")
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
pbmc@meta.data$res2 <-celltype[,2]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/tsne1.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res2")
dev.off()
pbmc@meta.data$res3 <-celltype[,3]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/tsne2.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res3")
dev.off()

###
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/markers.txt",header=T,row.names=1)

ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/Accuracy.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/Recall.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:10,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/Accuracy_10.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/Recall_10.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:20,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/Accuracy_20.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/Recall_20.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:50,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/Accuracy_50.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat1/RAW/non_neurous/Recall_50.txt",row.names=T,col.names=T,sep="\t")


##########
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
pbmc.data = readloom("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_2/GSM3017261_150000_CNS_nuclei_ds_0.3.loom")
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
output_h5= "/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data1), ncol(pbmc.data1)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data1), 1))
h5write(colnames(pbmc.data1), output_h5,"cell_id")
h5write(rownames(pbmc.data1), output_h5,"gene_name")
h5write(pbmc.data1, output_h5, "imputation")

output_h5= "/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data2), ncol(pbmc.data2)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data2), 1))
h5write(colnames(pbmc.data2), output_h5,"cell_id")
h5write(rownames(pbmc.data2), output_h5,"gene_name")
h5write(pbmc.data2, output_h5, "imputation")
############
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
library(Seurat) 

pbmc.data = readh5("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous.h5")
rownames(pbmc.data) = as.character(rownames(pbmc.data))
colnames(pbmc.data) = as.character(colnames(pbmc.data))
dim(pbmc.data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "neurous",min.cells = 10)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
##this_data@active.ident = cell_type_factor
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc),npcs = 100)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/PCA.pdf")
DimPlot(object = pbmc, reduction = "pca")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/PC_score.pdf",width=15,height=10)
ElbowPlot(object = pbmc, ndims = 100)
dev.off()
pbmc <- FindNeighbors(object = pbmc, dims = 1:50)
pbmc <- FindClusters(object = pbmc, resolution = 1.4)
pbmc <- RunUMAP(object = pbmc, dims = 1:50)
pbmc <- RunTSNE(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/UMAP.pdf")
DimPlot(object = pbmc, reduction = "umap")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/tsne.pdf")
DimPlot(object = pbmc, reduction = "tsne")
dev.off()
write.table(pbmc@meta.data,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data.txt",row.names=TRUE,col.names=TRUE)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pbmc.markers,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/markers.txt",row.names=TRUE,col.names=TRUE)
saveRDS(pbmc, file = "/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/data.rds")
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
pbmc@meta.data$res2 <-celltype[,2]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/tsne1.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res2")
dev.off()
pbmc@meta.data$res3 <-celltype[,3]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/tsne2.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res3")
dev.off()

data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/markers.txt",header=T,row.names=1)
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))
Others <-setdiff(A,hebing)
Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/Accuracy.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/Recall.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:10,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)
Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/Accuracy_10.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/Recall_10.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:20,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])






Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)

Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/Accuracy_20.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/Recall_20.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:50,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)

Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/Accuracy_50.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/neurous/Recall_50.txt",row.names=T,col.names=T,sep="\t")

#####
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
library(Seurat) 

pbmc.data = readh5("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous.h5")
rownames(pbmc.data) = as.character(rownames(pbmc.data))
colnames(pbmc.data) = as.character(colnames(pbmc.data))
dim(pbmc.data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "neurous",min.cells = 10)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
##this_data@active.ident = cell_type_factor
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc),npcs = 100)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/PCA.pdf")
DimPlot(object = pbmc, reduction = "pca")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/PC_score.pdf",width=15,height=10)
ElbowPlot(object = pbmc, ndims = 100)
dev.off()
pbmc <- FindNeighbors(object = pbmc, dims = 1:50)
pbmc <- FindClusters(object = pbmc, resolution = 1.4)
pbmc <- RunUMAP(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/UMAP.pdf")
DimPlot(object = pbmc, reduction = "umap")
dev.off()
pbmc <- RunTSNE(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/tsne.pdf")
DimPlot(object = pbmc, reduction = "tsne")
dev.off()
write.table(pbmc@meta.data,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.data.txt",row.names=TRUE,col.names=TRUE)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pbmc.markers,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/markers.txt",row.names=TRUE,col.names=TRUE)
saveRDS(pbmc, file = "/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/data.rds")
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
pbmc@meta.data$res2 <-celltype[,2]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/tsne1.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res2")
dev.off()
pbmc@meta.data$res3 <-celltype[,3]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/tsne2.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res3")
dev.off()

###
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/markers.txt",header=T,row.names=1)

ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/Accuracy.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/Recall.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:10,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/Accuracy_10.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/Recall_10.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:20,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/Accuracy_20.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/Recall_20.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:50,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/Accuracy_50.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat2/RAW/non_neurous/Recall_50.txt",row.names=T,col.names=T,sep="\t")

############
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
pbmc.data = readloom("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_3/GSM3017261_150000_CNS_nuclei_ds_0.3.loom")
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
output_h5= "/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data1), ncol(pbmc.data1)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data1), 1))
h5write(colnames(pbmc.data1), output_h5,"cell_id")
h5write(rownames(pbmc.data1), output_h5,"gene_name")
h5write(pbmc.data1, output_h5, "imputation")

output_h5= "/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data2), ncol(pbmc.data2)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data2), 1))
h5write(colnames(pbmc.data2), output_h5,"cell_id")
h5write(rownames(pbmc.data2), output_h5,"gene_name")
h5write(pbmc.data2, output_h5, "imputation")
############
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
library(Seurat) 

pbmc.data = readh5("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous.h5")
rownames(pbmc.data) = as.character(rownames(pbmc.data))
colnames(pbmc.data) = as.character(colnames(pbmc.data))
dim(pbmc.data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "neurous",min.cells = 10)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
##this_data@active.ident = cell_type_factor
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc),npcs = 100)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/PCA.pdf")
DimPlot(object = pbmc, reduction = "pca")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/PC_score.pdf",width=15,height=10)
ElbowPlot(object = pbmc, ndims = 100)
dev.off()
pbmc <- FindNeighbors(object = pbmc, dims = 1:50)
pbmc <- FindClusters(object = pbmc, resolution = 1.4)
pbmc <- RunUMAP(object = pbmc, dims = 1:50)
pbmc <- RunTSNE(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/UMAP.pdf")
DimPlot(object = pbmc, reduction = "umap")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/tsne.pdf")
DimPlot(object = pbmc, reduction = "tsne")
dev.off()
write.table(pbmc@meta.data,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data.txt",row.names=TRUE,col.names=TRUE)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pbmc.markers,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/markers.txt",row.names=TRUE,col.names=TRUE)
saveRDS(pbmc, file = "/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/data.rds")
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
pbmc@meta.data$res2 <-celltype[,2]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/tsne1.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res2")
dev.off()
pbmc@meta.data$res3 <-celltype[,3]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/tsne2.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res3")
dev.off()

data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/markers.txt",header=T,row.names=1)
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))
Others <-setdiff(A,hebing)
Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/Accuracy.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/Recall.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:10,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)
Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/Accuracy_10.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/Recall_10.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:20,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])






Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)

Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/Accuracy_20.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/Recall_20.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:50,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)

Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/Accuracy_50.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/neurous/Recall_50.txt",row.names=T,col.names=T,sep="\t")

#####
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
library(Seurat) 

pbmc.data = readh5("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous.h5")
rownames(pbmc.data) = as.character(rownames(pbmc.data))
colnames(pbmc.data) = as.character(colnames(pbmc.data))
dim(pbmc.data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "neurous",min.cells = 10)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
##this_data@active.ident = cell_type_factor
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc),npcs = 100)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/PCA.pdf")
DimPlot(object = pbmc, reduction = "pca")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/PC_score.pdf",width=15,height=10)
ElbowPlot(object = pbmc, ndims = 100)
dev.off()
pbmc <- FindNeighbors(object = pbmc, dims = 1:50)
pbmc <- FindClusters(object = pbmc, resolution = 1.4)
pbmc <- RunUMAP(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/UMAP.pdf")
DimPlot(object = pbmc, reduction = "umap")
dev.off()
pbmc <- RunTSNE(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/tsne.pdf")
DimPlot(object = pbmc, reduction = "tsne")
dev.off()
write.table(pbmc@meta.data,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.data.txt",row.names=TRUE,col.names=TRUE)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pbmc.markers,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/markers.txt",row.names=TRUE,col.names=TRUE)
saveRDS(pbmc, file = "/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/data.rds")
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
pbmc@meta.data$res2 <-celltype[,2]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/tsne1.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res2")
dev.off()
pbmc@meta.data$res3 <-celltype[,3]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/tsne2.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res3")
dev.off()

###
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/markers.txt",header=T,row.names=1)

ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/Accuracy.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/Recall.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:10,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/Accuracy_10.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/Recall_10.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:20,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/Accuracy_20.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/Recall_20.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:50,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/Accuracy_50.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat3/RAW/non_neurous/Recall_50.txt",row.names=T,col.names=T,sep="\t")


#########
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
pbmc.data = readloom("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_4/GSM3017261_150000_CNS_nuclei_ds_0.3.loom")
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
output_h5= "/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data1), ncol(pbmc.data1)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data1), 1))
h5write(colnames(pbmc.data1), output_h5,"cell_id")
h5write(rownames(pbmc.data1), output_h5,"gene_name")
h5write(pbmc.data1, output_h5, "imputation")

output_h5= "/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data2), ncol(pbmc.data2)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data2), 1))
h5write(colnames(pbmc.data2), output_h5,"cell_id")
h5write(rownames(pbmc.data2), output_h5,"gene_name")
h5write(pbmc.data2, output_h5, "imputation")
############
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
library(Seurat) 

pbmc.data = readh5("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous.h5")
rownames(pbmc.data) = as.character(rownames(pbmc.data))
colnames(pbmc.data) = as.character(colnames(pbmc.data))
dim(pbmc.data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "neurous",min.cells = 10)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
##this_data@active.ident = cell_type_factor
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc),npcs = 100)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/PCA.pdf")
DimPlot(object = pbmc, reduction = "pca")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/PC_score.pdf",width=15,height=10)
ElbowPlot(object = pbmc, ndims = 100)
dev.off()
pbmc <- FindNeighbors(object = pbmc, dims = 1:50)
pbmc <- FindClusters(object = pbmc, resolution = 1.4)
pbmc <- RunUMAP(object = pbmc, dims = 1:50)
pbmc <- RunTSNE(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/UMAP.pdf")
DimPlot(object = pbmc, reduction = "umap")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/tsne.pdf")
DimPlot(object = pbmc, reduction = "tsne")
dev.off()
write.table(pbmc@meta.data,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data.txt",row.names=TRUE,col.names=TRUE)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pbmc.markers,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/markers.txt",row.names=TRUE,col.names=TRUE)
saveRDS(pbmc, file = "/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/data.rds")
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
pbmc@meta.data$res2 <-celltype[,2]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/tsne1.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res2")
dev.off()
pbmc@meta.data$res3 <-celltype[,3]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/tsne2.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res3")
dev.off()

data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/markers.txt",header=T,row.names=1)
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))
Others <-setdiff(A,hebing)
Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/Accuracy.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/Recall.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:10,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)
Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/Accuracy_10.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/Recall_10.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:20,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])






Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)

Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/Accuracy_20.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/Recall_20.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:50,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)

Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/Accuracy_50.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/neurous/Recall_50.txt",row.names=T,col.names=T,sep="\t")

#####
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
library(Seurat) 

pbmc.data = readh5("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous.h5")
rownames(pbmc.data) = as.character(rownames(pbmc.data))
colnames(pbmc.data) = as.character(colnames(pbmc.data))
dim(pbmc.data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "neurous",min.cells = 10)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
##this_data@active.ident = cell_type_factor
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc),npcs = 100)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/PCA.pdf")
DimPlot(object = pbmc, reduction = "pca")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/PC_score.pdf",width=15,height=10)
ElbowPlot(object = pbmc, ndims = 100)
dev.off()
pbmc <- FindNeighbors(object = pbmc, dims = 1:50)
pbmc <- FindClusters(object = pbmc, resolution = 1.4)
pbmc <- RunUMAP(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/UMAP.pdf")
DimPlot(object = pbmc, reduction = "umap")
dev.off()
pbmc <- RunTSNE(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/tsne.pdf")
DimPlot(object = pbmc, reduction = "tsne")
dev.off()
write.table(pbmc@meta.data,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.data.txt",row.names=TRUE,col.names=TRUE)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pbmc.markers,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/markers.txt",row.names=TRUE,col.names=TRUE)
saveRDS(pbmc, file = "/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/data.rds")
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
pbmc@meta.data$res2 <-celltype[,2]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/tsne1.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res2")
dev.off()
pbmc@meta.data$res3 <-celltype[,3]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/tsne2.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res3")
dev.off()

###
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/markers.txt",header=T,row.names=1)

ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/Accuracy.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/Recall.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:10,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/Accuracy_10.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/Recall_10.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:20,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/Accuracy_20.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/Recall_20.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:50,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/Accuracy_50.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat4/RAW/non_neurous/Recall_50.txt",row.names=T,col.names=T,sep="\t")

#######
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
pbmc.data = readloom("/home/yuanhao/data/fn/split_seq/downsampling_first_repeat_5/GSM3017261_150000_CNS_nuclei_ds_0.3.loom")
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
output_h5= "/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data1), ncol(pbmc.data1)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data1), 1))
h5write(colnames(pbmc.data1), output_h5,"cell_id")
h5write(rownames(pbmc.data1), output_h5,"gene_name")
h5write(pbmc.data1, output_h5, "imputation")

output_h5= "/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous.h5"
h5createFile(output_h5)
h5createDataset(file = output_h5,
                dataset = "imputation", 
                dims = c(nrow(pbmc.data2), ncol(pbmc.data2)),
                storage.mode = "double",
                chunk=c(nrow(pbmc.data2), 1))
h5write(colnames(pbmc.data2), output_h5,"cell_id")
h5write(rownames(pbmc.data2), output_h5,"gene_name")
h5write(pbmc.data2, output_h5, "imputation")
############
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
library(Seurat) 

pbmc.data = readh5("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous.h5")
rownames(pbmc.data) = as.character(rownames(pbmc.data))
colnames(pbmc.data) = as.character(colnames(pbmc.data))
dim(pbmc.data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "neurous",min.cells = 10)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
##this_data@active.ident = cell_type_factor
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc),npcs = 100)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/PCA.pdf")
DimPlot(object = pbmc, reduction = "pca")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/PC_score.pdf",width=15,height=10)
ElbowPlot(object = pbmc, ndims = 100)
dev.off()
pbmc <- FindNeighbors(object = pbmc, dims = 1:50)
pbmc <- FindClusters(object = pbmc, resolution = 1.4)
pbmc <- RunUMAP(object = pbmc, dims = 1:50)
pbmc <- RunTSNE(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/UMAP.pdf")
DimPlot(object = pbmc, reduction = "umap")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/tsne.pdf")
DimPlot(object = pbmc, reduction = "tsne")
dev.off()
write.table(pbmc@meta.data,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data.txt",row.names=TRUE,col.names=TRUE)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pbmc.markers,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/markers.txt",row.names=TRUE,col.names=TRUE)
saveRDS(pbmc, file = "/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/data.rds")
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
pbmc@meta.data$res2 <-celltype[,2]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/tsne1.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res2")
dev.off()
pbmc@meta.data$res3 <-celltype[,3]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/tsne2.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res3")
dev.off()

data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/markers.txt",header=T,row.names=1)
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))
Others <-setdiff(A,hebing)
Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/Accuracy.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/Recall.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:10,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)
Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/Accuracy_10.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/Recall_10.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:20,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])






Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)

Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/Accuracy_20.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/Recall_20.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:50,])
}
ind1 <-which(data[,7]=="Deptor")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Rarb")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Satb2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tfap2d")##736
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fign")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Arap1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pax3")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Ntn1")
h <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pax2")
i <-data[ind1,][,6]
ind1 <-which(data[,7]=="Slc6a3") #### 47
j <-data[ind1,][,6]
ind1 <-which(data[,7]=="Fn1")
k <-data[ind1,][,6]
ind1 <-which(data[,7]=="Tspan18")
l <-data[ind1,][,6]
ind1 <-which(data[,7]=="Pde11a")
m <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dlx6os1")
n <-data[ind1,][,6]

Olfa <-a
Stri <-b
Cort <-c
Rost <-d
Thal <-setdiff(e,union(f,m))
Cere <-union(union(f,g),setdiff(setdiff(h,b),j))
Medu <-setdiff(i,h)
Basal <-intersect(h,j)
Hipp <-union(k,setdiff(l,m))
Spin <-m
Mirg <-n

uni <-c(Olfa,Stri,Cort,Rost,Thal,Cere,Medu,Basal,Hipp,Spin,Mirg)
dup <- unique(uni[duplicated(uni)])

Olfa1 <-setdiff(Olfa,dup)
Stri1 <-setdiff(Stri,dup)
Cort1 <-setdiff(Cort,dup)
Rost1 <-setdiff(Rost,dup)
Thal1 <-setdiff(Thal,dup)
Cere1 <-setdiff(Cere,dup)
Medu1 <-setdiff(Medu,dup)
Basal1 <-setdiff(Basal,dup)
Hipp1 <-setdiff(Hipp,dup)
Spin1 <-setdiff(Spin,dup)
Mirg1 <-setdiff(Mirg,dup)
hebing <-unique(c(Olfa1,Stri1,Cort1,Rost1,Thal1,Cere1,Medu1,Basal1,Hipp1,Spin1,Mirg1))
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)

Olfa2 <-cbind(Olfa1,rep("Olfactory Bulb",length(Olfa1)))
Stri2 <-cbind(Stri1,rep("Striatum",length(Stri1)))
Cort2 <-cbind(Cort1,rep("Cortex",length(Cort1)))
Rost2 <-cbind(Rost1,rep("Rostral Midbrain",length(Rost1)))
Thal2 <-cbind(Thal1,rep("Thalamus",length(Thal1)))
Cere2 <-cbind(Cere1,rep("Cerebellum",length(Cere1)))
Medu2 <-cbind(Medu1,rep("Medulla",length(Medu1)))
Basal2 <-cbind(Basal1,rep("Basal Ganglia",length(Basal1)))
Hipp2 <-cbind(Hipp1,rep("Hippocampus",length(Hipp1)))
Spin2 <-cbind(Spin1,rep("Spinalcord",length(Spin1)))
Mirg2 <-cbind(Mirg1,rep("Mirgrating Interneurous",length(Mirg1)))
Others1 <-cbind(Others,rep("Unresolved",length(Others)))
mat  <-rbind(Olfa2,Stri2,Cort2,Rost2,Thal2,Cere2,Medu2,Basal2,Hipp2,Spin2,Mirg2,Others1)

x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/meta.data_10.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Neuronal/cell_Neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]

index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Olfa1),length(Stri1),length(Cort1),length(Rost1),length(Thal1),
length(Cere1),length(Medu1),length(Basal1),length(Hipp1),length(Spin1),length(Mirg1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/Accuracy_50.txt",col.names=F,sep="\t")

#########
type <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus",
"Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","Mirgrating Interneurous")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[9])
ind1 <-which(data[,6]==type[9])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
i <- length(index)/length(dat[,7])
i1 <-length(ind)
i2 <-length(ind1)
i3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[10])
ind1 <-which(data[,6]==type[10])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
j <- length(index)/length(dat[,7])
j1 <-length(ind)
j2 <-length(ind1)
j3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[11])
ind1 <-which(data[,6]==type[11])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
l <- length(index)/length(dat[,7])
l1 <-length(ind)
l2 <-length(ind1)
l3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h,i,j,l)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,l1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,l2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,l3)

res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/neurous/Recall_50.txt",row.names=T,col.names=T,sep="\t")

#####
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
library(Seurat) 

pbmc.data = readh5("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous.h5")
rownames(pbmc.data) = as.character(rownames(pbmc.data))
colnames(pbmc.data) = as.character(colnames(pbmc.data))
dim(pbmc.data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "neurous",min.cells = 10)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
##this_data@active.ident = cell_type_factor
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc),npcs = 100)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/PCA.pdf")
DimPlot(object = pbmc, reduction = "pca")
dev.off()
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/PC_score.pdf",width=15,height=10)
ElbowPlot(object = pbmc, ndims = 100)
dev.off()
pbmc <- FindNeighbors(object = pbmc, dims = 1:50)
pbmc <- FindClusters(object = pbmc, resolution = 1.4)
pbmc <- RunUMAP(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/UMAP.pdf")
DimPlot(object = pbmc, reduction = "umap")
dev.off()
pbmc <- RunTSNE(object = pbmc, dims = 1:50)
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/tsne.pdf")
DimPlot(object = pbmc, reduction = "tsne")
dev.off()
write.table(pbmc@meta.data,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.data.txt",row.names=TRUE,col.names=TRUE)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pbmc.markers,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/markers.txt",row.names=TRUE,col.names=TRUE)
saveRDS(pbmc, file = "/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/data.rds")
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
pbmc@meta.data$res2 <-celltype[,2]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/tsne1.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res2")
dev.off()
pbmc@meta.data$res3 <-celltype[,3]
pdf(file="/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/tsne2.pdf")
DimPlot(object = pbmc, reduction = "tsne", group.by = "res3")
dev.off()

###
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/markers.txt",header=T,row.names=1)

ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/Accuracy.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/Recall.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:10,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/Accuracy_10.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/Recall_10.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:20,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/Accuracy_20.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/Recall_20.txt",row.names=T,col.names=T,sep="\t")

marker <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/markers.txt",header=T,row.names=1)
dat <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))-1
data <-NULL
for(i in 0:clu){
ind <-which(marker[,6]==i)
data <-rbind(data,marker[ind,][1:50,])
}
ind1 <-which(data[,7]=="Mbp")
a <- data[ind1,][,6]
ind1 <-which(data[,7]=="Pdgfra")
b <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dock2")
c <-data[ind1,][,6]
ind1 <-which(data[,7]=="Rgs5")
d <-data[ind1,][,6]
ind1 <-which(data[,7]=="Col1a2")
e <-data[ind1,][,6]
ind1 <-which(data[,7]=="Aldh1l1")
f <-data[ind1,][,6]
ind1 <-which(data[,7]=="Dnah11")
g <-data[ind1,][,6]
ind1 <-which(data[,7]=="Mybpc1")
h <- data[ind1,][,6]

Oligo <-a
OPC <-setdiff(b,e)
Immune <-c
VLMC <-d
Vasc <-e
Astrocyte <-setdiff(f,g)
Epend <-g
OEC <-h
uni <-c(Oligo,OPC,Immune,VLMC,Vasc,Astrocyte,Epend,OEC)
dup <- unique(uni[duplicated(uni)])
Oligo1 <-setdiff(Oligo,dup)
OPC1 <-setdiff(OPC,dup)
Immune1 <-setdiff(Immune,dup)
VLMC1 <-setdiff(VLMC,dup)
Vasc1 <-setdiff(Vasc,dup)
Astrocyte1 <-setdiff(Astrocyte,dup)
Epend1 <-setdiff(Epend,dup)
OEC1<-setdiff(OEC,dup)
hebing <-unique(c(Oligo1,OPC1,Immune1,VLMC1,Vasc1,Astrocyte1,Epend1,OEC1))
meta <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))

Oligo2 <-cbind(Oligo1,rep("Oligo",length(Oligo1)))
OPC2 <-cbind(OPC1,rep("OPC",length(OPC1)))
Immune2 <-cbind(Immune1,rep("Immune",length(Immune1)))
VLMC2 <-cbind(VLMC1,rep("VLMC",length(VLMC1)))
Vasc2 <-cbind(Vasc1,rep("Vasc",length(Vasc1)))
Astrocyte2 <-cbind(Astrocyte1,rep("Astrocyte",length(Astrocyte1)))
Epend2 <-cbind(Epend1,rep("Epend",length(Epend1)))
OEC2 <-cbind(OEC1,rep("OEC",length(OEC1)))
Others1 <-cbind(Others,rep("Others",length(Others)))
mat  <-rbind(Oligo2,OPC2,Immune2,VLMC2,Vasc2,Astrocyte2,Epend2,OEC2,Others1)
x<-mat[order(as.numeric(mat[,1])),2]
data <-read.table("/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
head(data)
celltype <-read.table("/home/wucheng/imputation/split-seq/raw1/Non_neuronal/cell_non_neuronal_type.txt")
cell <-intersect(rownames(data),rownames(celltype))
data <-data[cell,]
celltype <-celltype[cell,]
data$V8 <-celltype[,3]
index <-which(data[,6]==data[,7])
a <-length(index)/length(data[,7])
library(mclust)
b <-adjustedRandIndex(data[,6], data[,7])
c <-length(y)
d <-c(length(Astrocyte1),length(Epend1),length(Immune1),length(OEC1),length(Oligo1),length(OPC1),length(Vasc1),length(VLMC1))
d <-sum(d!=0)
e <-length(Others)
aa <-c(a,b,c,d,e)
bb <-c("Accuracy","ARI","Clusters","Types","Unknow")
cc <-rbind(bb,aa)
write.table(cc,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/Accuracy_50.txt",col.names=F,sep="\t")
#####

type <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
ind <-which(data[,7]==type[1])
ind1 <-which(data[,6]==type[1])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
a <- length(index)/length(dat[,7])
a1 <-length(ind)
a2 <-length(ind1)
a3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[2])
ind1 <-which(data[,6]==type[2])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
b <- length(index)/length(dat[,7])
b1 <-length(ind)
b2 <-length(ind1)
b3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[3])
ind1 <-which(data[,6]==type[3])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
c <- length(index)/length(dat[,7])
c1 <-length(ind)
c2 <-length(ind1)
c3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[4])
ind1 <-which(data[,6]==type[4])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
d <- length(index)/length(dat[,7])
d1 <-length(ind)
d2 <-length(ind1)
d3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[5])
ind1 <-which(data[,6]==type[5])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
e <- length(index)/length(dat[,7])
e1 <-length(ind)
e2 <-length(ind1)
e3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[6])
ind1 <-which(data[,6]==type[6])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
f <- length(index)/length(dat[,7])
f1 <-length(ind)
f2 <-length(ind1)
f3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[7])
ind1 <-which(data[,6]==type[7])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
g <- length(index)/length(dat[,7])
g1 <-length(ind)
g2 <-length(ind1)
g3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

ind <-which(data[,7]==type[8])
ind1 <-which(data[,6]==type[8])
dat <-data[ind,]
index <-which(dat[,6]==dat[,7])
h <- length(index)/length(dat[,7])
h1 <-length(ind)
h2 <-length(ind1)
h3 <-length(intersect(ind,ind1))/length(union(ind,ind1))

Recall <-c(a,b,c,d,e,f,g,h)
real <-c(a1,b1,c1,d1,e1,f1,g1,h1)
predict <-c(a2,b2,c2,d2,e2,f2,g2,h2)
Jaac <-c(a3,b3,c3,d3,e3,f3,g3,h3)


res <-rbind(Recall,real,predict,Jaac)
rownames(res) <-c("Recall","real","predict","Jaac")
colnames(res) <-type
write.table(res,"/home/wucheng/imputation/split-seq/ds/repeat5/RAW/non_neurous/Recall_50.txt",row.names=T,col.names=T,sep="\t")




