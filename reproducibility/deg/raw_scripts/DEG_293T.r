data =read.table("/home/wucheng/imputation/bulk/GSE129240_rsem_expected_counts.tsv",header=T,row.names=1,sep="\t")
data <-data[1:58440,c(1,2,13,14)]
data <-data[rowSums(data) > 0,]
gen <-matrix(unlist(strsplit(rownames(data),"\\.")),2)[1,]
infer <-as.matrix(read.table("/home/wucheng/imputation/deg/reference.gtf"))
infer <-infer[,c(1,6)]
gene <-NULL
for(i in 1:length(gen)){
index <-which(infer[,1]==gen[i])
gene <-c(gene,infer[index,2])
}
rs <-rowSums(data)
kid <- sapply(unique(gene),function(sid) {
            tmp <- which(gene==sid)
            if (length(tmp)==1) {
                  tmp
            } else {
                  tmp[which.max(rs[tmp])]
            }
      })
	  
data <-data[kid,]
rownaems(data) <-gene[kid]
data <- data[!grepl('^MT-',row.names(data)),]
data = round(data)
saveRDS(data,"/home/wucheng/imputation/deg/293T_Jurkat/RSEM_gene.rds")

#######
cnt <- readRDS("/home/wucheng/imputation/deg/293T_Jurkat/RSEM_gene.rds")
colnames(cnt) <-c("HEK","HEK","Jurkat","Jurkat")
ct <- colnames(cnt)
for (i in unique(ct)){colnames(cnt)[ct == i] <- paste0(i,'_',1:sum(ct==i))
}
for (n1 in 1:(length(unique(ct))-1)){
  i = unique(ct)[n1]
  for (n2 in ((n1+1):length(unique(ct)))){
      j = unique(ct)[n2]
      print(paste0(i,'_',j))
      expr = cnt[, ct %in% c(i,j)]
      expr = expr[rowSums(expr>0)>0,]
      sct <- sub('_.*','',colnames(expr))
      library(limma)
      des <- cbind(1,ifelse(sct==i,1,0))
      fit <- eBayes(lmFit(voom(expr,des),design=des))
      res <- topTable(fit,coef=2,number=nrow(expr))
#     res <- res[res[,'adj.P.Val']<0.05,]
	  ind <-intersect(c(which(res[,1]>=2),which(res[,1]<=(-2))),which(res[,'adj.P.Val']<=0.05))
      res <- res[ind,]
      gs <- rownames(res)
      saveRDS(gs,paste0('/home/wucheng/imputation/deg/293T_Jurkat/bulk1/',i,'_',j,'_diffgene.rds'))
  }
}

########RAW
pbmc.data = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_mc_10_mce_1.loom")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

metadata <-as.data.frame(as.matrix(pbmc@meta.data))
pbmc@active.ident <-metadata[,1]
cluster1.markers <- FindMarkers(pbmc, ident.1 = "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "wilcox")
saveRDS(cluster1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/Raw/wilcox/HEK_Jurkat.rds")
clu1.markers <- FindMarkers(pbmc, ident.1 =  "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "MAST")
saveRDS(clu1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/Raw/MAST/HEK_Jurkat.rds")

####DISC
pbmc.data = readloom("/home/yuanhao/DISC_imputation_result/JURKAT_293T/result/imputation.loom")
raw = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_mc_10_mce_1.loom")
pbmc.data <-pbmc.data[rownames(raw),]
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
metadata <-as.data.frame(as.matrix(pbmc@meta.data))
pbmc@active.ident <-metadata[,1]
cluster1.markers <- FindMarkers(pbmc, ident.1 = "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "wilcox")
saveRDS(cluster1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/DISC/wilcox/HEK_Jurkat.rds")
clu1.markers <- FindMarkers(pbmc, ident.1 =  "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "MAST")
saveRDS(clu1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/DISC/MAST/HEK_Jurkat.rds")
####MAGIC
pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_MAGIC_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
metadata <-as.data.frame(as.matrix(pbmc@meta.data))
pbmc@active.ident <-metadata[,1]
cluster1.markers <- FindMarkers(pbmc, ident.1 = "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "wilcox")
saveRDS(cluster1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/MAGIC/wilcox/HEK_Jurkat.rds")
clu1.markers <- FindMarkers(pbmc, ident.1 =  "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "MAST")
saveRDS(clu1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/MAGIC/MAST/HEK_Jurkat.rds")
####DCA
pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_DCA_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
metadata <-as.data.frame(as.matrix(pbmc@meta.data))
pbmc@active.ident <-metadata[,1]
cluster1.markers <- FindMarkers(pbmc, ident.1 = "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "wilcox")
saveRDS(cluster1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/Method/DCA/wilcox/HEK_Jurkat.rds")
clu1.markers <- FindMarkers(pbmc, ident.1 =  "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "MAST")
saveRDS(clu1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/DCA/MAST/HEK_Jurkat.rds")
###scScope
pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_scScope_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
metadata <-as.data.frame(as.matrix(pbmc@meta.data))
pbmc@active.ident <-metadata[,1]
cluster1.markers <- FindMarkers(pbmc, ident.1 = "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "wilcox")
saveRDS(cluster1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/scScope/wilcox/HEK_Jurkat.rds")
clu1.markers <- FindMarkers(pbmc, ident.1 =  "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "MAST")
saveRDS(clu1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/scScope/MAST/HEK_Jurkat.rds")
##scVI
pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_scVI_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
metadata <-as.data.frame(as.matrix(pbmc@meta.data))
pbmc@active.ident <-metadata[,1]
cluster1.markers <- FindMarkers(pbmc, ident.1 = "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "wilcox")
saveRDS(cluster1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/scVI/wilcox/HEK_Jurkat.rds")
clu1.markers <- FindMarkers(pbmc, ident.1 =  "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "MAST")
saveRDS(clu1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/scVI/MAST/HEK_Jurkat.rds")
##scImpute
pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_scImpute_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
metadata <-as.data.frame(as.matrix(pbmc@meta.data))
pbmc@active.ident <-metadata[,1]
cluster1.markers <- FindMarkers(pbmc, ident.1 = "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "wilcox")
saveRDS(cluster1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/scImpute/wilcox/HEK_Jurkat.rds")
clu1.markers <- FindMarkers(pbmc, ident.1 =  "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "MAST")
saveRDS(clu1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/scImpute/MAST/HEK_Jurkat.rds")
###deepImpute
pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_deepImpute_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
metadata <-as.data.frame(as.matrix(pbmc@meta.data))
pbmc@active.ident <-metadata[,1]
cluster1.markers <- FindMarkers(pbmc, ident.1 = "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "wilcox")
saveRDS(cluster1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/deepImpute/wilcox/HEK_Jurkat.rds")
clu1.markers <- FindMarkers(pbmc, ident.1 =  "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "MAST")
saveRDS(clu1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/deepImpute/MAST/HEK_Jurkat.rds")
###VIPER
pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_VIPER_gene_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
metadata <-as.data.frame(as.matrix(pbmc@meta.data))
pbmc@active.ident <-metadata[,1]
cluster1.markers <- FindMarkers(pbmc, ident.1 = "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "wilcox")
saveRDS(cluster1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/VIPER/wilcox/HEK_Jurkat.rds")
clu1.markers <- FindMarkers(pbmc, ident.1 =  "293T", ident.2 = "JURKAT", min.pct = 0.1,logfc.threshold=0,test.use = "MAST")
saveRDS(clu1.markers,"/home/wucheng/imputation/deg/293T_Jurkat/Method/VIPER/MAST/HEK_Jurkat.rds")

#####overlap
allmtd = list.files('/home/wucheng/imputation/deg/293T_Jurkat/Method/')
mtd = allmtd[2]
allf = list.files(paste0('/home/wucheng/imputation/deg/293T_Jurkat/Method/',mtd,'/wilcox/'))
ove <- sapply(allmtd, function(mtd){
  sapply(allf,function(f) {
      print(f)
      if (file.exists(paste0('/home/wucheng/imputation/deg/293T_Jurkat/Method/',mtd,'/wilcox/',f))){
        res = readRDS(paste0('/home/wucheng/imputation/deg/293T_Jurkat/Method/',mtd,'/wilcox/',f))
        res = res[order(res[,'p_val_adj']),]
        gs = readRDS(paste0('/home/wucheng/imputation/deg/293T_Jurkat/bulk/', sub('.rds','',f),'_diffgene.rds'))
        tmp <- mean(sapply(c(1:100)*10,function(i) {
          mean(rownames(res)[1:i] %in% gs)  ## discuss
        }))  
      } else {
        return(NA)
      }
  })
})
ove = t(ove)
colnames(ove) = sub('.rds','', colnames(ove))
saveRDS(ove, '/home/wucheng/imputation/deg/293T_Jurkat/wilcox_bulk_sc_overlaps.rds')


allmtd = list.files('/home/wucheng/imputation/deg/293T_Jurkat/Method/')
mtd = allmtd[2]
allf = list.files(paste0('/home/wucheng/imputation/deg/293T_Jurkat/Method/',mtd,'/MAST/'))
ove <- sapply(allmtd, function(mtd){
  sapply(allf,function(f) {
      print(f)
      if (file.exists(paste0('/home/wucheng/imputation/deg/293T_Jurkat/Method/',mtd,'/MAST/',f))){
        res = readRDS(paste0('/home/wucheng/imputation/deg/293T_Jurkat/Method/',mtd,'/MAST/',f))
        res = res[order(res[,'p_val_adj']),]
        gs = readRDS(paste0('/home/wucheng/imputation/deg/293T_Jurkat/bulk/', sub('.rds','',f),'_diffgene.rds'))
        tmp <- mean(sapply(c(1:100)*10,function(i) {
          mean(rownames(res)[1:i] %in% gs)  ## discuss
        }))  
      } else {
        return(NA)
      }
  })
})
ove = t(ove)
colnames(ove) = sub('.rds','', colnames(ove))
saveRDS(ove, '/home/wucheng/imputation/deg/293T_Jurkat/MAST_bulk_sc_overlaps.rds')

#####bulk1
allmtd = list.files('/home/wucheng/imputation/deg/293T_Jurkat/Method/')
mtd = allmtd[2]
allf = list.files(paste0('/home/wucheng/imputation/deg/293T_Jurkat/Method/',mtd,'/wilcox/'))
ove <- sapply(allmtd, function(mtd){
  sapply(allf,function(f) {
      print(f)
      if (file.exists(paste0('/home/wucheng/imputation/deg/293T_Jurkat/Method/',mtd,'/wilcox/',f))){
        res = readRDS(paste0('/home/wucheng/imputation/deg/293T_Jurkat/Method/',mtd,'/wilcox/',f))
        res = res[order(res[,'p_val_adj']),]
        gs = readRDS(paste0('/home/wucheng/imputation/deg/293T_Jurkat/bulk1/', sub('.rds','',f),'_diffgene.rds'))
        tmp <- mean(sapply(c(1:100)*10,function(i) {
          mean(rownames(res)[1:i] %in% gs)  ## discuss
        }))  
      } else {
        return(NA)
      }
  })
})
ove = t(ove)
colnames(ove) = sub('.rds','', colnames(ove))
saveRDS(ove, '/home/wucheng/imputation/deg/293T_Jurkat/wilcox_bulk_sc_overlaps1.rds')


allmtd = list.files('/home/wucheng/imputation/deg/293T_Jurkat/Method/')
mtd = allmtd[2]
allf = list.files(paste0('/home/wucheng/imputation/deg/293T_Jurkat/Method/',mtd,'/MAST/'))
ove <- sapply(allmtd, function(mtd){
  sapply(allf,function(f) {
      print(f)
      if (file.exists(paste0('/home/wucheng/imputation/deg/293T_Jurkat/Method/',mtd,'/MAST/',f))){
        res = readRDS(paste0('/home/wucheng/imputation/deg/293T_Jurkat/Method/',mtd,'/MAST/',f))
        res = res[order(res[,'p_val_adj']),]
        gs = readRDS(paste0('/home/wucheng/imputation/deg/293T_Jurkat/bulk1/', sub('.rds','',f),'_diffgene.rds'))
        tmp <- mean(sapply(c(1:100)*10,function(i) {
          mean(rownames(res)[1:i] %in% gs)  ## discuss
        }))  
      } else {
        return(NA)
      }
  })
})
ove = t(ove)
colnames(ove) = sub('.rds','', colnames(ove))
saveRDS(ove, '/home/wucheng/imputation/deg/293T_Jurkat/MAST_bulk_sc_overlaps1.rds')

###
library(reshape2)
library(ggplot2)
library(ggrepel)
ove = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/wilcox_bulk_sc_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o1 <- t(matrix(ove))
colnames(o1) <-c("DCA","deepImpute","DISC","MAGIC","Raw","scImpute","scScope","scVI","VIPER")
ove = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/MAST_bulk_sc_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o2 <- t(matrix(ove))
colnames(o2) <-c("DCA","deepImpute","DISC","MAGIC","Raw","scImpute","scScope","scVI","VIPER")
int <- intersect(colnames(o1),colnames(o2))
pd <- data.frame(MAST=o2[,int],Wilcox=o1[,int],mtd=int)
pdf('/home/wucheng/imputation/deg/293T_Jurkat/wilcox_mast_compare.pdf',width=4,height=4)
ggplot(pd,aes(x=MAST,y=Wilcox,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')
dev.off()
##
ove = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/wilcox_bulk_sc_overlaps1.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o1 <- t(matrix(ove))
colnames(o1) <-c("DCA","deepImpute","DISC","MAGIC","Raw","scImpute","scScope","scVI","VIPER")
ove = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/MAST_bulk_sc_overlaps1.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o2 <- t(matrix(ove))
colnames(o2) <-c("DCA","deepImpute","DISC","MAGIC","Raw","scImpute","scScope","scVI","VIPER")
int <- intersect(colnames(o1),colnames(o2))
pd <- data.frame(MAST=o2[,int],Wilcox=o1[,int],mtd=int)
pdf('/home/wucheng/imputation/deg/293T_Jurkat/wilcox_mast_compare1.pdf',width=4,height=4)
ggplot(pd,aes(x=MAST,y=Wilcox,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')
dev.off()

#########################################################################NULLDE ##293T
pbmc.data = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_mc_10_mce_1.loom")
ct = sub('_.*','',colnames(pbmc.data))
imp <-pbmc.data[,which(ct=="293T")]
df = expand.grid(c(10,50,100,500),c(10,50,100,500))
colnames(df) =  c('n1','n2')
df = df[df[,'n1']<=df[,'n2'],]
for (i in 1:nrow(df)){
  print(i)
  cn1 = df[i,'n1']
  cn2 = df[i,'n2']
  set.seed(12345)
  id = sample(1:ncol(imp), cn1+cn2)
  expr = imp[,id]
  pbmc <- CreateSeuratObject(counts = as.data.frame(expr))
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc@meta.data$V1 <-c(rep("cn1",cn1),rep("cn2",cn2))
  pbmc@active.ident <-as.data.frame(as.matrix(pbmc@meta.data))[,4]
  cluster.markers <- FindMarkers(pbmc, ident.1 = "cn1", ident.2 = "cn2", min.pct = 0.1,logfc.threshold=0,test.use = "wilcox")
  saveRDS(cluster.markers, paste0('/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/Method/Raw/wilcox/',cn1,'_',cn2,'.rds'))

  cluster1.markers <- FindMarkers(pbmc, ident.1 = "cn1", ident.2 = "cn2", min.pct = 0.1,logfc.threshold=0,test.use = "MAST")  
  saveRDS(cluster1.markers, paste0('/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/Method/Raw/MAST/',cn1,'_',cn2,'.rds')) 
}

#####
## sc_293T
source('/home/wucheng/imputation/DEG/function.R')
allmtd = list.files('/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/Method/')
df <- sapply(allmtd, function(mtd){
  rdir = paste0('/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/Method/', mtd,'/wilcox/')
  af = list.files(rdir)
  sapply(af, function(f){
    res <- readRDS(paste0(rdir,f))
	up <-intersect(which(res$p_val<=0.01),which(res$avg_logFC >=0.25))
	down <-intersect(which(res$p_val<=0.01),which(res$avg_logFC <(-0.25)))
    sum(length(up),length(down))  
  })
})

if (is.list(df)){
  df = df[sapply(df,length)>0]
  pd = as.matrix(do.call(cbind, df))  
} else {
  pd = df
}
if (grepl('_0_0.rds.rds',rownames(pd)[1])){
  rownames(pd) = sub('_0_0.rds.rds','',rownames(pd))  
}  else if (grepl('_0_0.rds',rownames(pd)[1])){
  rownames(pd) = sub('_0_0.rds','',rownames(pd))  
} else {
  rownames(pd) = sub('.rds','',rownames(pd))
}

library(reshape2)
pd = melt(pd)
colnames(pd) = c('data','method','Num')

mtdorder = names(sort(tapply(pd[,'Num'],list(pd[,'method']), mean), decreasing = T))
stat = tapply(pd[,'Num'],list(pd[,'method']), mean)
saveRDS(stat,'/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/293T_wilcox.rds')

#######
source('/home/wucheng/imputation/DEG/function.R')
library(RColorBrewer)
allmtd = list.files('/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/Method/')
df <- sapply(allmtd, function(mtd){
  rdir = paste0('/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/Method/', mtd,'/MAST/')
  af = list.files(rdir)
  sapply(af, function(f){
    res <- readRDS(paste0(rdir,f))
	up <-intersect(which(res$p_val<=0.01),which(res$avg_logFC >=0.25))
	down <-intersect(which(res$p_val<=0.01),which(res$avg_logFC <(-0.25)))
    sum(length(up),length(down))  
  })
})
if (is.list(df)){
  df = df[sapply(df,length)>0]
  pd = as.matrix(do.call(cbind, df))  
} else {
  pd = df
}
if (grepl('_0_0.rds.rds',rownames(pd)[1])){
  rownames(pd) = sub('_0_0.rds.rds','',rownames(pd))  
}  else if (grepl('_0_0.rds',rownames(pd)[1])){
  rownames(pd) = sub('_0_0.rds','',rownames(pd))  
} else {
  rownames(pd) = sub('.rds','',rownames(pd))
}
library(reshape2)
pd = melt(pd)
colnames(pd) = c('data','method','Num')
mtdorder = names(sort(tapply(pd[,'Num'],list(pd[,'method']), mean), decreasing = T))
stat = tapply(pd[,'Num'],list(pd[,'method']), mean)
saveRDS(stat,'/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/293T_MAST.rds')

############

##############
source('/home/wucheng/imputation/DEG/function.R')
library(reshape2)
library(ggplot2)
library(ggrepel)
o1 = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/293T_MAST.rds')
o2 = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/293T_wilcox.rds')
int <- intersect(names(o1),names(o2))
pd <- data.frame(MAST=o1,Wilcoxon=o2,mtd=int)
xmin <- (0)
ymin <- (0)
pdf('/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/wilcox_mast_compare.pdf',width=4,height=4)
ggplot(pd,aes(x=MAST,y=Wilcoxon,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
  xlim(c(xmin,max((pd$MAST)+2))) + ylim(c(ymin,max(pd$Wilcoxon)+2))
dev.off()

#####JURKAT
## sc_JURKAT
source('/home/wucheng/imputation/DEG/function.R')
library(RColorBrewer)
allmtd = list.files('/home/wucheng/imputation/deg/293T_Jurkat/Ju_NULLDE/Method/')
df <- sapply(allmtd, function(mtd){
  rdir = paste0('/home/wucheng/imputation/deg/293T_Jurkat/Ju_NULLDE/Method/', mtd,'/wilcox/')
  af = list.files(rdir)
  sapply(af, function(f){
    res <- readRDS(paste0(rdir,f))
	up <-intersect(which(res$p_val<=0.01),which(res$avg_logFC >=0.25))
	down <-intersect(which(res$p_val<=0.01),which(res$avg_logFC <(-0.25)))
    sum(length(up),length(down))  
  })
})
if (is.list(df)){
  df = df[sapply(df,length)>0]
  pd = as.matrix(do.call(cbind, df))  
} else {
  pd = df
}
if (grepl('_0_0.rds.rds',rownames(pd)[1])){
  rownames(pd) = sub('_0_0.rds.rds','',rownames(pd))  
}  else if (grepl('_0_0.rds',rownames(pd)[1])){
  rownames(pd) = sub('_0_0.rds','',rownames(pd))  
} else {
  rownames(pd) = sub('.rds','',rownames(pd))
}
library(reshape2)
pd = melt(pd)
colnames(pd) = c('data','method','Num')
mtdorder = names(sort(tapply(pd[,'Num'],list(pd[,'method']), mean), decreasing = T))
stat = tapply(pd[,'Num'],list(pd[,'method']), mean)
saveRDS(stat,'/home/wucheng/imputation/deg/293T_Jurkat/Ju_NULLDE/JURKAT_wilcox.rds')

#######
source('/home/wucheng/imputation/DEG/function.R')
library(RColorBrewer)
allmtd = list.files('/home/wucheng/imputation/deg/293T_Jurkat/Ju_NULLDE/Method')
df <- sapply(allmtd, function(mtd){
  rdir = paste0('/home/wucheng/imputation/deg/293T_Jurkat/Ju_NULLDE/Method/', mtd,'/MAST/')
  af = list.files(rdir)
  sapply(af, function(f){
    res <- readRDS(paste0(rdir,f))
	up <-intersect(which(res$p_val<=0.01),which(res$avg_logFC >=0.25))
	down <-intersect(which(res$p_val<=0.01),which(res$avg_logFC <(-0.25)))
    sum(length(up),length(down))  
  })
})

if (is.list(df)){
  df = df[sapply(df,length)>0]
  pd = as.matrix(do.call(cbind, df))  
} else {
  pd = df
}
if (grepl('_0_0.rds.rds',rownames(pd)[1])){
  rownames(pd) = sub('_0_0.rds.rds','',rownames(pd))  
}  else if (grepl('_0_0.rds',rownames(pd)[1])){
  rownames(pd) = sub('_0_0.rds','',rownames(pd))  
} else {
  rownames(pd) = sub('.rds','',rownames(pd))
}

library(reshape2)
pd = melt(pd)
colnames(pd) = c('data','method','Num')

mtdorder = names(sort(tapply(pd[,'Num'],list(pd[,'method']), mean), decreasing = T))
stat = tapply(pd[,'Num'],list(pd[,'method']), mean)
saveRDS(stat,'/home/wucheng/imputation/deg/293T_Jurkat/Ju_NULLDE/JURKAT_MAST.rds')

############

##############
source('/home/wucheng/imputation/DEG/function.R')
library(reshape2)
library(ggplot2)
library(ggrepel)
o1 = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/Ju_NULLDE/JURKAT_MAST.rds')
o2 = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/Ju_NULLDE/JURKAT_wilcox.rds')
int <- intersect(names(o1),names(o2))
pd <- data.frame(MAST=o1,Wilcoxon=o2,mtd=int)
xmin <- (0)
ymin <- (0)
pdf('/home/wucheng/imputation/deg/293T_Jurkat/Ju_NULLDE/wilcox_mast_compare1.pdf',width=4,height=4)
ggplot(pd,aes(x=MAST,y=Wilcoxon,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
  xlim(c(xmin,max((pd$MAST)+2))) + ylim(c(ymin,max(pd$Wilcoxon)+2))
dev.off()


##########
###############plot
library(ggplot2)
ove = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/MAST_bulk_sc_overlaps1.rds')
o1 <-as.matrix(ove)
o2 = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/293T_MAST.rds')
int <- names(o2)
pd <- data.frame(MAST=(1-o1[1,]),Wilcoxon=o2,mtd=int)
pdf('/home/wucheng/imputation/deg/293T_Jurkat/mast_compare_293T.pdf',width=4,height=4)
ggplot(pd,aes(x=Wilcoxon,y=MAST,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
 xlab("diff_number") + ylab("1-overlap")
 dev.off()

library(ggplot2)
ove = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/wilcox_bulk_sc_overlaps.rds')
o1 <-as.matrix(ove)
o2 = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/293T_wilcox.rds')
int <- names(o2)
pd <- data.frame(MAST=(1-o1[1,]),Wilcoxon=o2,mtd=int)
pdf('/home/wucheng/imputation/deg/293T_Jurkat/wilcox_compare_393T.pdf',width=4,height=4)
ggplot(pd,aes(x=Wilcoxon,y=MAST,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
 xlab("diff_number") + ylab("1-overlap")
 dev.off()
 
ove = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/MAST_bulk_sc_overlaps1.rds')
o1 <-as.matrix(ove)
o2 = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/293T_MAST.rds')
int <- names(o2)
pd <- data.frame(MAST=(1-o1[1,]),Wilcoxon=o2,mtd=int)
pdf('/home/wucheng/imputation/deg/293T_Jurkat/mast_compare_293T1.pdf',width=4,height=4)
ggplot(pd,aes(x=Wilcoxon,y=MAST,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
 xlab("diff_number") + ylab("1-overlap")
 dev.off()

library(ggplot2)
ove = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/wilcox_bulk_sc_overlaps1.rds')
o1 <-as.matrix(ove)
o2 = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/NULLDE/293T_wilcox.rds')
int <- names(o2)
pd <- data.frame(MAST=(1-o1[1,]),Wilcoxon=o2,mtd=int)
pdf('/home/wucheng/imputation/deg/293T_Jurkat/wilcox_compare_393T1.pdf',width=4,height=4)
ggplot(pd,aes(x=Wilcoxon,y=MAST,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
 xlab("diff_number") + ylab("1-overlap")
 dev.off()
 
####
library(ggplot2)
ove = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/MAST_bulk_sc_overlaps.rds')
o1 <-as.matrix(ove)
o2 = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/Ju_NULLDE/JURKAT_MAST.rds')
int <- names(o2)
pd <- data.frame(MAST=(1-o1[1,]),Wilcoxon=o2,mtd=int)
pdf('/home/wucheng/imputation/deg/293T_Jurkat/mast_compare_Jurkat.pdf',width=4,height=4)
ggplot(pd,aes(x=Wilcoxon,y=MAST,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
 xlab("diff_number") + ylab("1-overlap")
 dev.off()

library(ggplot2)
ove = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/wilcox_bulk_sc_overlaps.rds')
o1 <-as.matrix(ove)
o2 = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/Ju_NULLDE/JURKAT_wilcox.rds')
int <- names(o2)
pd <- data.frame(MAST=(1-o1[1,]),Wilcoxon=o2,mtd=int)
pdf('/home/wucheng/imputation/deg/293T_Jurkat/wilcox_compare_Jurkat.pdf',width=4,height=4)
ggplot(pd,aes(x=Wilcoxon,y=MAST,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
 xlab("diff_number") + ylab("1-overlap")
 dev.off()
 
ove = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/MAST_bulk_sc_overlaps1.rds')
o1 <-as.matrix(ove)
o2 = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/Ju_NULLDE/JURKAT_MAST.rds')
int <- names(o2)
pd <- data.frame(MAST=(1-o1[1,]),Wilcoxon=o2,mtd=int)
pdf('/home/wucheng/imputation/deg/293T_Jurkat/mast_compare_Jurkat1.pdf',width=4,height=4)
ggplot(pd,aes(x=Wilcoxon,y=MAST,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
 xlab("diff_number") + ylab("1-overlap")
 dev.off()

library(ggplot2)
ove = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/wilcox_bulk_sc_overlaps1.rds')
o1 <-as.matrix(ove)
o2 = readRDS('/home/wucheng/imputation/deg/293T_Jurkat/Ju_NULLDE/JURKAT_wilcox.rds')
int <- names(o2)
pd <- data.frame(MAST=(1-o1[1,]),Wilcoxon=o2,mtd=int)
pdf('/home/wucheng/imputation/deg/293T_Jurkat/wilcox_compare_Jurkat1.pdf',width=4,height=4)
ggplot(pd,aes(x=Wilcoxon,y=MAST,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
 xlab("diff_number") + ylab("1-overlap")
 dev.off()


################








