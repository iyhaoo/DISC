##########
data =read.table("C:/Users/ADMIN/Desktop/imputation/bulk/FPKM.txt",header=T,row.names=1,sep="\t")
data <-data[rowSums(data) > 0,]
gen <-matrix(unlist(strsplit(rownames(data),"\\.")),2)[1,]
infer <-as.matrix(read.table("C:/Users/ADMIN/Desktop/imputation/bulk/fitter.GRCh38.93.gtf"))
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
rownames(data) <-gene[kid]
data <- data[!grepl('^MT-',row.names(data)),]
saveRDS(data,"C:/Users/ADMIN/Desktop/imputation/bulk/FPKM_gene1.rds")

###########293T
bulk <-readRDS("/home/wucheng/imputation/bulk/FPKM_gene1.rds")
bulk <-log2(bulk+1)
pbmc.data = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/imputation/raw_mc_10_mce_1.loom")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
gene <-intersect(rownames(scData),rownames(bulk))
aa <-rowMeans(bulk[gene,][,1:2])
bb <-rowMeans(scData[gene,])
library(ggplot2)
df <- data.frame(x = bb, y = aa)
library(psych)
coorda <-corr.test(aa,bb,method="spearman")
Raw <-coorda$r

pbmc.data = readloom("/home/yuanhao/DISC_imputation_result/293T/result/imputation.loom")
raw = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/imputation/raw_mc_10_mce_1.loom")
pbmc.data <-pbmc.data[rownames(raw),]
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
DISC <-coorda$r

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/imputation/raw_DCA_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
DCA <-coorda$r

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/imputation/raw_scScope_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
scScope <-coorda$r

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/imputation/raw_VIPER_gene_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
VIPER <-coorda$r

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/imputation/raw_scVI_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
scVI <-coorda$r

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/imputation/raw_deepImpute_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
deepImpute <-coorda$r

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/imputation/raw_scImpute_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
scImpute <-coorda$r

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/imputation/raw_MAGIC_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
MAGIC <-coorda$r

T293 <-c(Raw,DISC,DCA,scScope,VIPER,scVI,deepImpute,scImpute,MAGIC)

###########JURKAT
bulk <-readRDS("/home/wucheng/imputation/bulk/FPKM_gene1.rds")
bulk <-log2(bulk+1)
pbmc.data = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/imputation/raw_mc_10_mce_1.loom")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
gene <-intersect(rownames(scData),rownames(bulk))
aa <-rowMeans(bulk[gene,][,3:4])
bb <-rowMeans(scData[gene,])
library(ggplot2)
df <- data.frame(x = bb, y = aa)
library(psych)
coorda <-corr.test(aa,bb,method="spearman")
Raw <-coorda$r
#pdf(file="/home/wucheng/imputation/bulk/pcc_293T_Raw.pdf",width=6,height=4)	
#p <-ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(size = 1)+ labs(x="pseudobulk", y = "Log2(FPKM+1)")
#p <-p +theme_classic()+theme(axis.text.x = element_text(size = 12,angle = 0, hjust = 1, vjust = 1),axis.text.y = element_text(size = 12,angle = 0, hjust = 1, vjust = 1))   
#p +annotate("text", label = paste0('R==',coorda$r), parse = TRUE, x = 4, y = 6,size =6 )
#dev.off()	

pbmc.data = readloom("/home/yuanhao/DISC_imputation_result/JURKAT/result/imputation.loom")
raw = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/imputation/raw_mc_10_mce_1.loom")
pbmc.data <-pbmc.data[rownames(raw),]
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
DISC <-coorda$r

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/imputation/raw_DCA_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
DCA <-coorda$r

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/imputation/raw_scScope_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
scScope <-coorda$r

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/imputation/raw_VIPER_gene_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
VIPER <-coorda$r

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/imputation/raw_scVI_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
scVI <-coorda$r

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/imputation/raw_deepImpute_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
deepImpute <-coorda$r

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/imputation/raw_scImpute_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
scImpute <-coorda$r

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT/imputation/raw_MAGIC_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
scData <-as.matrix(pbmc[["RNA"]]@data)
bb <-rowMeans(scData[gene,])
coorda <-corr.test(aa,bb,method="spearman")
MAGIC <-coorda$r

JURKAT <-c(Raw,DISC,DCA,scScope,VIPER,scVI,deepImpute,scImpute,MAGIC)
data <-rbind(T293,JURKAT)
rownames(data) <-c("293T","JURKAT")
colnames(data) <-c("Raw","DISC","DCA","scScope","VIPER","scVI","deepImpute","scImpute","MAGIC")
saveRDS(data,"/home/wucheng/imputation/bulk/sc_bulk_scc.rds")

######plot
library(ggplot2)
library(reshape2)
#data <-readRDS("C:/Users/ADMIN/Desktop/imputation/DEG/sc_bulk_scc.rds")
data <-readRDS("/home/wucheng/imputation/bulk/sc_bulk_scc.rds")
df <-melt(data)
levels <-c("Raw","DISC","scImpute","VIPER","MAGIC","DCA","deepImpute","scScope","scVI")
ggplot(df, aes(x=factor(Var2,levels=levels), y=value,fill = Var1))+
geom_bar(position = "dodge",stat = "identity",width = 0.5)+
ylim(0,1)+theme_classic()+theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+theme(axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"))

##############
bulk <-readRDS("/home/wucheng/imputation/bulk/FPKM_gene1.rds")
bulk <-log2(bulk+1)
bulk = cbind(rowMeans(bulk[,1:2]),rowMeans(bulk[,3:4]))
colnames(bulk) = c('293T','JURKAT')
pbmc.data = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_mc_10_mce_1.loom")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
sexpr <-as.matrix(pbmc[["RNA"]]@data)
cl = sub('_.*','',colnames(sexpr))
intgene = intersect(rownames(bulk),rownames(sexpr))
bulk = bulk[intgene,]
sexpr = sexpr[intgene,]
imp1 = sexpr[, which(cl=='293T')]
imp2 = sexpr[, which(cl=='JURKAT')]
aa <-rowMeans(imp1)-rowMeans(imp2)
bb <-bulk[,1]-bulk[,2]
library(ggplot2)
df <- data.frame(x = aa, y = bb)
library(psych)
coorda <-cor(aa,bb,method="spearman")
Raw <-coorda
#pdf(file="/home/wucheng/imputation/bulk/pcc_293T_JURKAT.pdf",width=6,height=4)	
#p <-ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(size = 1)+ labs( x="pseudobulk", y = "Log2(FPKM+1)")
#p <-p +theme_classic()+theme(axis.text.x = element_text(size = 12,angle = 0, hjust = 1, vjust = 1),axis.text.y = element_text(size = 12,angle = 0, hjust = 1, vjust = 1))   
#p +annotate("text", label =paste0('R==',coorda), parse = TRUE, x = 2, y = 0,size =6 )
#dev.off()	

pbmc.data = readloom("/home/yuanhao/DISC_imputation_result/JURKAT_293T/result/imputation.loom")
raw = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_mc_10_mce_1.loom")
pbmc.data <-pbmc.data[rownames(raw),]
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
sexpr <-as.matrix(pbmc[["RNA"]]@data)
sexpr = sexpr[intgene,]
imp1 = sexpr[, which(cl=='293T')]
imp2 = sexpr[, which(cl=='JURKAT')]
aa <-rowMeans(imp1)-rowMeans(imp2)
df <- data.frame(x = aa, y = bb)
coorda <-cor(aa,bb,method="spearman")
DISC <-coorda


pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_MAGIC_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
sexpr <-as.matrix(pbmc[["RNA"]]@data)
sexpr = sexpr[intgene,]
imp1 = sexpr[, which(cl=='293T')]
imp2 = sexpr[, which(cl=='JURKAT')]
aa <-rowMeans(imp1)-rowMeans(imp2)
df <- data.frame(x = aa, y = bb)
coorda <-cor(aa,bb,method="spearman")
MAGIC <-coorda

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_DCA_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
sexpr <-as.matrix(pbmc[["RNA"]]@data)
sexpr = sexpr[intgene,]
imp1 = sexpr[, which(cl=='293T')]
imp2 = sexpr[, which(cl=='JURKAT')]
aa <-rowMeans(imp1)-rowMeans(imp2)
df <- data.frame(x = aa, y = bb)
coorda <-cor(aa,bb,method="spearman")
DCA <-coorda

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_scVI_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
sexpr <-as.matrix(pbmc[["RNA"]]@data)
sexpr = sexpr[intgene,]
imp1 = sexpr[, which(cl=='293T')]
imp2 = sexpr[, which(cl=='JURKAT')]
aa <-rowMeans(imp1)-rowMeans(imp2)
df <- data.frame(x = aa, y = bb)
coorda <-cor(aa,bb,method="spearman")
scVI <-coorda

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_scScope_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
sexpr <-as.matrix(pbmc[["RNA"]]@data)
sexpr = sexpr[intgene,]
imp1 = sexpr[, which(cl=='293T')]
imp2 = sexpr[, which(cl=='JURKAT')]
aa <-rowMeans(imp1)-rowMeans(imp2)
df <- data.frame(x = aa, y = bb)
coorda <-cor(aa,bb,method="spearman")
scScope <-coorda

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_scImpute_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
sexpr <-as.matrix(pbmc[["RNA"]]@data)
sexpr = sexpr[intgene,]
imp1 = sexpr[, which(cl=='293T')]
imp2 = sexpr[, which(cl=='JURKAT')]
aa <-rowMeans(imp1)-rowMeans(imp2)
df <- data.frame(x = aa, y = bb)
coorda <-cor(aa,bb,method="spearman")
scImpute <-coorda

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_deepImpute_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
sexpr <-as.matrix(pbmc[["RNA"]]@data)
sexpr = sexpr[intgene,]
imp1 = sexpr[, which(cl=='293T')]
imp2 = sexpr[, which(cl=='JURKAT')]
aa <-rowMeans(imp1)-rowMeans(imp2)
df <- data.frame(x = aa, y = bb)
coorda <-cor(aa,bb,method="spearman")
deepImpute <-coorda

pbmc.data = readh5("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_VIPER_gene_mc_10_mce_1.hdf5")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
sexpr <-as.matrix(pbmc[["RNA"]]@data)
sexpr = sexpr[intgene,]
imp1 = sexpr[, which(cl=='293T')]
imp2 = sexpr[, which(cl=='JURKAT')]
aa <-rowMeans(imp1)-rowMeans(imp2)
df <- data.frame(x = aa, y = bb)
coorda <-cor(aa,bb,method="spearman")
VIPER <-coorda

dif_bulk <-c(Raw,DISC,DCA,scScope,VIPER,scVI,deepImpute,scImpute,MAGIC)
saveRDS(dif_bulk,"/home/wucheng/imputation/bulk/sc_dif_bulk_scc.rds")

#value <-readRDS("C:/Users/ADMIN/Desktop/imputation/DEG/sc_dif_bulk_scc.rds")
value <-readRDS("/home/wucheng/imputation/bulk/sc_dif_bulk_scc.rds")
method <-c("Raw","DISC","DCA","scScope","VIPER","scVI","deepImpute","scImpute","MAGIC")
df <-data.frame(method,value)
levels <-c("Raw","DISC","scImpute","VIPER","MAGIC","DCA","deepImpute","scScope","scVI")
ggplot(df, aes(x=factor(method,levels=levels), y=value))+
geom_bar(stat = "identity",width=0.4) +
ylim(0,0.65)+theme_classic()+theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+theme(axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"))


#####(cell-cell VS bulk)-dif
#bulk <-readRDS("/home/wucheng/imputation/bulk/FPKM_gene.rds")
#bulk <-log2(bulk+1)
#pbmc.data = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/293T/imputation/raw_mc_10_mce_1.loom")
#pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
#pbmc
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#scData <-as.matrix(pbmc[["RNA"]]@data)
#gene <-intersect(rownames(scData),rownames(bulk))
#aa <-rowMeans(bulk[gene,][,1:2])
#bb <-scData[gene,]
#scc <-NULL
#for(i in 1:ncol(bb)){
#coorda <-corr.test(aa,bb[,i],method="spearman")
#scc <-c(scc,coorda$r)
#}

#bulk <-readRDS("/home/wucheng/imputation/bulk/FPKM_gene.rds")
#bulk <-log2(bulk+1)
#bulk = cbind(rowMeans(bulk[,1:2]),rowMeans(bulk[,3:4]))
#colnames(bulk) = c('293T','JURKAT')
#pbmc.data = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/JURKAT_293T/imputation/raw_mc_10_mce_1.loom")
#pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
#pbmc
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#sexpr <-as.matrix(pbmc[["RNA"]]@data)
#sexpr <-cbind(sexpr[,1001:1100],sexpr[,4001:4100])
#cl = sub('_.*','',colnames(sexpr))
#intgene = intersect(rownames(bulk),rownames(sexpr))
#bulk = bulk[intgene,]
#sexpr = sexpr[intgene,]
#get_scCellType_bulkCellType_cor <- function(ct1, ct2, bulkDiff){
  imp1 = sexpr[, which(cl==ct1)]
  imp2 = sexpr[, which(cl==ct2)]
  corvec = NULL
  corvec <- sapply(1:ncol(imp1),function(i) {
    sapply(1:ncol(imp2), function(j) {
      cor((imp1[,i] - imp2[,j]), bulkDiff,method='spearman')
    })
  })
  as.vector(corvec)
#}
#cn = NULL
#v = sapply(1:(ncol(bulk)-1), function(i){
#  sapply((i+1):ncol(bulk), function(j){
#    cn = c(cn,print(paste0(colnames(bulk)[i],'_',colnames(bulk)[j])))
#    get_scCellType_bulkCellType_cor(ct1=colnames(bulk)[i], ct2=colnames(bulk)[j], bulkDiff = bulk[,colnames(bulk)[i]]-bulk[,colnames(bulk)[j]])
#  })
#})
#colnames(v) = cn
#saveRDS(v,'/home/wucheng/imputation/bulk/sc_sc/raw1.rds')

