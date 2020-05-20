###########10X_5cl
process10x_rmDupGenes <- function(genebycellmat){
      tb = genebycellmat
      tb <- tb[rowSums(tb) > 0,]
      gn = rownames(tb)  
      rs <- rowSums(tb)
      kid <- sapply(unique(gn),function(sid) {
            tmp <- which(gn==sid)
            if (length(tmp)==1) {
                  tmp
            } else {
                  tmp[which.max(rs[tmp])]
            }
      })
      tb <- tb[kid,]
      row.names(tb) <- gn[kid]
      tb <- tb[!grepl('^MT-',row.names(tb)),]
      tb = round(tb)
}

##########
bk = read.csv('/home/wucheng/imputation/DEG/GSE86337_reverse.stranded.unfiltered.count.matrix.txt', as.is=T, sep='\t')
tb = as.matrix(read.csv('/home/wucheng/imputation/DEG/GSE86337_anno.txt',sep='\t'))
suppressMessages(library(org.Hs.eg.db))
genename <- select(org.Hs.eg.db, key=as.character(bk$Entrez.Gene.IDs),columns=c("SYMBOL"),keytype="ENTREZID")$SYMBOL
bk = bk[,-1]
bk = as.matrix(bk)
rownames(bk) = genename
bk <- bk[!is.na(row.names(bk)),]
mat = process10x_rmDupGenes(bk)
colnames(mat) = tb[match(sapply(colnames(mat), function(i) paste0(strsplit(i,'_')[[1]][1:2],collapse='_')), tb[,'Sample.name']), 'characteristics.cell.line']
colnames(mat) = paste0(colnames(mat),'_', 1:ncol(mat))
saveRDS(mat,'/home/wucheng/imputation/DEG/GSE86337_processed_count.rds')

#####bulk deg
cnt <- readRDS("/home/wucheng/imputation/DEG/GSE86337_processed_count.rds")
ct <- c("HCC827","HCC827","H2228","H2228","H838","H838","A549","A549","H1975","H1975")
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
#	  ind <-intersect(c(which(res[,1]>=2),which(res[,1]<=(-2))),which(res[,'adj.P.Val']<=0.05))
#     res <- res[res[,'adj.P.Val']<0.05,]
	  ind <-intersect(c(which(res[,1]>=2),which(res[,1]<=(-2))),which(res[,'adj.P.Val']<=0.05))
      res <- res[ind,]
      gs <- rownames(res)
      saveRDS(gs,paste0('/home/wucheng/imputation/deg/10x_5cl/bulk2/',i,'_',j,'_diffgene.rds'))
  }
}

#######sc deg
pbmc.data = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/10X_5CL/imputation/raw_mc_10_mce_1.loom")
pbmc <- CreateSeuratObject(counts = as.data.frame(pbmc.data))
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

metadata <-as.data.frame(as.matrix(read.table("/home/wucheng/imputation/DEG/sc_10x_5cl.metadata.csv",header=T,row.names=1,sep=",")))
pbmc@active.ident <-metadata[,29]

ct <-c("HCC827","H2228","H838","A549","H1975")
for (n1 in 1:(length(ct)-1)){
  i = unique(ct)[n1]
  for (n2 in ((n1+1):length(ct))){
      j = unique(ct)[n2]
      print(paste0(i,'_',j))
      cluster.markers <- FindMarkers(pbmc, ident.1 = i, ident.2 = j, min.pct = 0,logfc.threshold=0,test.use = "wilcox")
      saveRDS(cluster.markers,paste0('/home/wucheng/imputation/deg/10x_5cl/Method/Raw/wilcox/',i,'_',j,'_diffgene.rds'))
      clu.markers <- FindMarkers(pbmc, ident.1 = i, ident.2 = j, min.pct = 0,logfc.threshold=0,test.use = "MAST")
	  saveRDS(clu.markers,paste0('/home/wucheng/imputation/deg/10x_5cl/Method/Raw/MAST/',i,'_',j,'_diffgene.rds'))
}}

#########overlap
##############bulk
allmtd = list.files('/home/wucheng/imputation/deg/10x_5cl/Method/')
mtd = allmtd[1]
allf = list.files(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/wilcox/'))
ove <- sapply(allmtd, function(mtd){
  sapply(allf,function(f) {
      print(f)
      if (file.exists(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/wilcox/',f))){
        res = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/wilcox/',f))
        res = res[order(res[,'p_val']),]
        gs = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/bulk/', sub('.rds','',f),'.rds'))
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
saveRDS(ove, '/home/wucheng/imputation/deg/10x_5cl/wilcox_bulk_sc_overlaps.rds')

allmtd = list.files('/home/wucheng/imputation/deg/10x_5cl/Method/')
mtd = allmtd[5]
allf = list.files(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/MAST/'))
ove <- sapply(allmtd, function(mtd){
  sapply(allf,function(f) {
      print(f)
      if (file.exists(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/MAST/',f))){
        res = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/MAST/',f))
        res = res[order(res[,'p_val']),]
        gs = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/bulk/', sub('.rds','',f),'.rds'))
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
saveRDS(ove, '/home/wucheng/imputation/deg/10x_5cl/MAST_bulk_sc_overlaps.rds')

######bulk1
allmtd = list.files('/home/wucheng/imputation/deg/10x_5cl/Method/')
mtd = allmtd[1]
allf = list.files(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/wilcox/'))
ove <- sapply(allmtd, function(mtd){
  sapply(allf,function(f) {
      print(f)
      if (file.exists(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/wilcox/',f))){
        res = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/wilcox/',f))
        res = res[order(res[,'p_val']),]
        gs = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/bulk2/', sub('.rds','',f),'.rds'))
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
saveRDS(ove, '/home/wucheng/imputation/deg/10x_5cl/wilcox_bulk_sc_overlaps2.rds')
##
allmtd = list.files('/home/wucheng/imputation/deg/10x_5cl/Method/')
mtd = allmtd[5]
allf = list.files(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/MAST/'))
ove <- sapply(allmtd, function(mtd){
  sapply(allf,function(f) {
      print(f)
      if (file.exists(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/MAST/',f))){
        res = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/MAST/',f))
        res = res[order(res[,'p_val']),]
        gs = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/bulk2/', sub('.rds','',f),'.rds'))
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
saveRDS(ove, '/home/wucheng/imputation/deg/10x_5cl/MAST_bulk_sc_overlaps2.rds')

#######plot
library(reshape2)
library(ggplot2)
library(ggrepel)
ove = readRDS('/home/wucheng/imputation/deg/10x_5cl/wilcox_bulk_sc_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o1 <- rowMeans(ove)
ove = readRDS('/home/wucheng/imputation/deg/10x_5cl/MAST_bulk_sc_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o2 <- rowMeans(ove)
int <- intersect(names(o1),names(o2))
pd <- data.frame(MAST=o2,Wilcox=o1,mtd=int)
pdf('/home/wucheng/imputation/deg/10x_5cl/wilcox_mast_compare.pdf',width=4,height=4)
ggplot(pd,aes(x=MAST,y=Wilcox,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')
dev.off()
##
ove = readRDS('/home/wucheng/imputation/deg/10x_5cl/wilcox_bulk_sc_overlaps2.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o1 <- rowMeans(ove)
ove = readRDS('/home/wucheng/imputation/deg/10x_5cl/MAST_bulk_sc_overlaps2.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o2 <- rowMeans(ove)
int <- intersect(names(o1),names(o2))
pd <- data.frame(MAST=o2,Wilcox=o1,mtd=int)
pdf('/home/wucheng/imputation/deg/10x_5cl/wilcox_mast_compare2.pdf',width=4,height=4)
ggplot(pd,aes(x=MAST,y=Wilcox,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')
dev.off()

#################################################
pbmc.data = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/10X_5CL/imputation/raw_mc_10_mce_1.loom")
metadata <-as.data.frame(as.matrix(read.table("/home/wucheng/imputation/DEG/sc_10x_5cl.metadata.csv",header=T,row.names=1,sep=",")))
imp <-pbmc.data[,which(metadata[29]=="A549")]
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
  cluster1.markers <- FindMarkers(pbmc, ident.1 = "cn1", ident.2 = "cn2", min.pct = 0.1,logfc.threshold=0,test.use = "MAST")  
  saveRDS(cluster.markers, paste0('/home/wucheng/imputation/deg/10x_5cl/NULLDE/Method/Raw/wilcox/',cn1,'_',cn2,'.rds'))
  saveRDS(cluster1.markers, paste0('/home/wucheng/imputation/deg/10x_5cl/NULLDE/Method/Raw/MAST/',cn1,'_',cn2,'.rds'))
}
######
source('/home/wucheng/imputation/DEG/function.R')
allmtd = list.files('/home/wucheng/imputation/deg/10x_5cl/NULLDE/Method/')
df <- sapply(allmtd, function(mtd){
  rdir = paste0('/home/wucheng/imputation/deg/10x_5cl/NULLDE/Method/', mtd,'/wilcox/')
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
saveRDS(stat,paste0('/home/wucheng/imputation/deg/10x_5cl/NULLDE/10x_5cl_wilcox.rds'))

#######
source('/home/wucheng/imputation/DEG/function.R')
library(RColorBrewer)
allmtd = list.files('/home/wucheng/imputation/deg/10x_5cl/NULLDE/Method/')
df <- sapply(allmtd, function(mtd){
  rdir = paste0('/home/wucheng/imputation/deg/10x_5cl/NULLDE/Method/', mtd,'/MAST/')
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
saveRDS(stat,paste0('/home/wucheng/imputation/deg/10x_5cl/NULLDE/10x_5cl_MAST.rds'))

##############
source('/home/wucheng/imputation/DEG/function.R')
library(reshape2)
library(ggplot2)
library(ggrepel)
o1 = readRDS('/home/wucheng/imputation/deg/10x_5cl/NULLDE/10x_5cl_MAST.rds')
o2 = readRDS('/home/wucheng/imputation/deg/10x_5cl/NULLDE/10x_5cl_wilcox.rds')
int <- intersect(names(o1),names(o2))
pd <- data.frame(MAST=o1,Wilcoxon=o2,mtd=int)
xmin <- (0)
ymin <- (0)
pdf('/home/wucheng/imputation/deg/10x_5cl/NULLDE/wilcox_mast_compare.pdf',width=4,height=4)
ggplot(pd,aes(x=MAST,y=Wilcoxon,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
 xlim(c(xmin,max((pd$MAST)+2))) + ylim(c(ymin,max(pd$Wilcoxon)+2))
 dev.off()

###############
library(ggplot2)
ove = readRDS('/home/wucheng/imputation/deg/10x_5cl/MAST_bulk_sc_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o1 <- rowMeans(ove)
o2 = readRDS('/home/wucheng/imputation/deg/10x_5cl/NULLDE/10x_5cl_MAST.rds')
int <- intersect(names(o1),names(o2))
pd <- data.frame(MAST=(1-o1),Wilcoxon=o2,mtd=int)
pdf('/home/wucheng/imputation/deg/10x_5cl/mast_compare.pdf',width=4,height=4)
ggplot(pd,aes(x=Wilcoxon,y=MAST,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
xlab("diff_number") + ylab("1-overlap")
 dev.off()

ove = readRDS('/home/wucheng/imputation/deg/10x_5cl/wilcox_bulk_sc_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o1 <- rowMeans(ove)
o2 = readRDS('/home/wucheng/imputation/deg/10x_5cl/NULLDE/10x_5cl_wilcox.rds')
int <- intersect(names(o1),names(o2))
pd <- data.frame(MAST=(1-o1),Wilcoxon=o2,mtd=int)

pdf('/home/wucheng/imputation/deg/10x_5cl/wilcox_compare.pdf',width=4,height=4)
ggplot(pd,aes(x=Wilcoxon,y=MAST,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
xlab("diff_number") + ylab("1-overlap")
 dev.off()

###
ove = readRDS('/home/wucheng/imputation/deg/10x_5cl/MAST_bulk_sc_overlaps2.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o1 <- rowMeans(ove)
o2 = readRDS('/home/wucheng/imputation/deg/10x_5cl/NULLDE/10x_5cl_MAST.rds')
int <- intersect(names(o1),names(o2))
pd <- data.frame(MAST=(1-o1),Wilcoxon=o2,mtd=int)
pdf('/home/wucheng/imputation/deg/10x_5cl/mast_compare2.pdf',width=4,height=4)
ggplot(pd,aes(x=Wilcoxon,y=MAST,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
xlab("diff_number") + ylab("1-overlap")
 dev.off()

ove = readRDS('/home/wucheng/imputation/deg/10x_5cl/wilcox_bulk_sc_overlaps2.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o1 <- rowMeans(ove)
o2 = readRDS('/home/wucheng/imputation/deg/10x_5cl/NULLDE/10x_5cl_wilcox.rds')
int <- intersect(names(o1),names(o2))
pd <- data.frame(MAST=(1-o1),Wilcoxon=o2,mtd=int)

pdf('/home/wucheng/imputation/deg/10x_5cl/wilcox_compare2.pdf',width=4,height=4)
ggplot(pd,aes(x=Wilcoxon,y=MAST,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
xlab("diff_number") + ylab("1-overlap")
 dev.off()

############
#############bulk 1.5
allmtd = list.files('/home/wucheng/imputation/deg/10x_5cl/Method/')
mtd = allmtd[5]
allf = list.files(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/wilcox/'))
#raw = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/10X_5CL/imputation/raw_mc_10_mce_1.loom")
#metadata <-as.matrix(read.table("/home/wucheng/imputation/DEG/sc_10x_5cl.metadata.csv",header=T,row.names=1,sep=","))
#ct <-matrix(metadata[,29])[,1]
f=allf[1]
ove <- sapply(allmtd, function(mtd){
  print(mtd)
  t(sapply(allf,function(f) {
    print(f)
    if (file.exists(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/wilcox/',f))){
      res = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/wilcox/',f))
      res = res[order(res[,'p_val_adj']),]
      ct1 <- sub('_.*','',f)
      ct2 <- sub('.*_','',sub('_diffgene.rds','',f))
	  raw_res = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',"Raw",'/wilcox/',f))
      fc <- abs(as.matrix(raw_res)[,2])      
      highid = names(which(fc>quantile(fc,0.9)))## express in more cells
      lowid = names(which(fc<quantile(fc,0.1)))## less
      #highid = names(which(rowMeans(tmpmat>0)>quantile(rowMeans(tmpmat>0),0.9)))## express in more cells
      #lowid = names(which(rowMeans(tmpmat>0)<quantile(rowMeans(tmpmat>0),0.1)))## less
      #bfex <- list.files('/home/wucheng/imputation/deg/10x_5cl/bulk/')
      #bf <- c(paste0(ct1,'_',ct2,'_diffgene.rds'),paste0(ct2,'_',ct1,'_diffgene.rds'))
      #bf <- intersect(bfex,bf)
      gs = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/bulk1/',f))
      tmpres <- res[rownames(res) %in% highid,]
      tmp1 <- mean(sapply(c(1:100)*10,function(i) {
        mean(rownames(tmpres[1:i,]) %in% gs)  ## discuss
      }))
      tmpres <- res[rownames(res) %in% lowid,]
      tmp2 <- mean(sapply(c(1:100)*10,function(i) {
        mean(rownames(tmpres[1:i,]) %in% gs)  ## discuss
      }),na.rm=T)
      c(tmp1,tmp2)
    } else {
      return(c(NA,NA))
    }
  }))
}) ## first 10 high, last 10 low
ove_high = t(ove[1:10,])
ove_low = t(ove[11:20,])
colnames(ove_high) <- colnames(ove_low) <- sub('.rds','', allf)
rdir = '/home/wucheng/imputation/deg/10x_5cl/high_low1.5/wilcox/'
dir.create(rdir,showWarnings = F, recursive = T)
saveRDS(ove_high, paste0(rdir,'bulk_sc_diffgene_overlaps_moreExprGene.rds'))
saveRDS(ove_low, paste0(rdir,'bulk_sc_diffgene_overlaps_lessExprGene.rds'))

###
allmtd = list.files('/home/wucheng/imputation/deg/10x_5cl/Method/')
mtd = allmtd[5]
allf = list.files(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/MAST/'))
#raw = readloom("/home/yuanhao/github_repositories/DISC/reproducibility/data/10X_5CL/imputation/raw_mc_10_mce_1.loom")
#metadata <-as.matrix(read.table("/home/wucheng/imputation/DEG/sc_10x_5cl.metadata.csv",header=T,row.names=1,sep=","))
#ct <-matrix(metadata[,29])[,1]
f=allf[1]
ove <- sapply(allmtd, function(mtd){
  print(mtd)
  t(sapply(allf,function(f) {
    print(f)
    if (file.exists(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/MAST/',f))){
      res = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',mtd,'/MAST/',f))
      res = res[order(res[,'p_val_adj']),]
      ct1 <- sub('_.*','',f)
      ct2 <- sub('.*_','',sub('_diffgene.rds','',f))
	  raw_res = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/Method/',"Raw",'/MAST/',f))
      fc <- abs(as.matrix(raw_res)[,2])      
      highid = names(which(fc>quantile(fc,0.9)))## express in more cells
      lowid = names(which(fc<quantile(fc,0.1)))## less
      #highid = names(which(rowMeans(tmpmat>0)>quantile(rowMeans(tmpmat>0),0.9)))## express in more cells
      #lowid = names(which(rowMeans(tmpmat>0)<quantile(rowMeans(tmpmat>0),0.1)))## less
      #bfex <- list.files('/home/wucheng/imputation/deg/10x_5cl/bulk/')
      #bf <- c(paste0(ct1,'_',ct2,'_diffgene.rds'),paste0(ct2,'_',ct1,'_diffgene.rds'))
      #bf <- intersect(bfex,bf)
      gs = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/bulk1/',f))
      tmpres <- res[rownames(res) %in% highid,]
      tmp1 <- mean(sapply(c(1:100)*10,function(i) {
        mean(rownames(tmpres[1:i,]) %in% gs)  ## discuss
      }))
      tmpres <- res[rownames(res) %in% lowid,]
      tmp2 <- mean(sapply(c(1:100)*10,function(i) {
        mean(rownames(tmpres[1:i,]) %in% gs)  ## discuss
      }),na.rm=T)
      c(tmp1,tmp2)
    } else {
      return(c(NA,NA))
    }
  }))
}) ## first 10 high, last 10 low
ove_high = t(ove[1:10,])
ove_low = t(ove[11:20,])
colnames(ove_high) <- colnames(ove_low) <- sub('.rds','', allf)
rdir = '/home/wucheng/imputation/deg/10x_5cl/high_low1.5/MAST/'
dir.create(rdir,showWarnings = F, recursive = T)
saveRDS(ove_high, paste0(rdir,'bulk_sc_diffgene_overlaps_moreExprGene.rds'))
saveRDS(ove_low, paste0(rdir,'bulk_sc_diffgene_overlaps_lessExprGene.rds'))
##########

###########boxplot
ove_high <-readRDS("/home/wucheng/imputation/deg/10x_5cl/high_low1.5/MAST/bulk_sc_diffgene_overlaps_moreExprGene.rds")
c(min(ove_high),max(ove_high))
df <-melt(ove_high)
levels <-c("Raw","DISC","scImpute","VIPER","MAGIC","DCA","deepImpute","scScope","scVI")
p  <- ggplot(df, aes(x=factor(Var1,levels=levels), y=value, fill=factor(Var1,levels=levels))) + geom_boxplot()
p <-p+ ylim(0.30,0.75)+theme_classic()+theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+theme(axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"))
pdf(file="/home/wucheng/imputation/deg/10x_5cl/high_low1.5/MAST/high_boxplot.pdf",width=6,height=4)
p+theme(legend.title=element_blank())
dev.off()

ove_low <-readRDS("/home/wucheng/imputation/deg/10x_5cl/high_low1.5/MAST/bulk_sc_diffgene_overlaps_lessExprGene.rds")
c(min(ove_low),max(ove_low))
df <-melt(ove_low)
levels <-c("Raw","DISC","scImpute","VIPER","MAGIC","DCA","deepImpute","scScope","scVI")
p  <- ggplot(df, aes(x=factor(Var1,levels=levels), y=value, fill=factor(Var1,levels=levels))) + geom_boxplot()
p <-p+ ylim(0,0.16)+theme_classic()+theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+theme(axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"))
pdf(file="/home/wucheng/imputation/deg/10x_5cl/high_low1.5/MAST/low_boxplot.pdf",width=6,height=4)
p+theme(legend.title=element_blank())
dev.off()

####
cnt <- readRDS("/home/wucheng/imputation/DEG/GSE86337_processed_count.rds")
ct <- c("HCC827","HCC827","H2228","H2228","H838","H838","A549","A549","H1975","H1975")
for (i in unique(ct)){colnames(cnt)[ct == i] <- paste0(i,'_',1:sum(ct==i))
}
cutoff <-c(0.0,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0)
for(cut in 1:length(cutoff)){
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
#	  ind <-intersect(c(which(res[,1]>=2),which(res[,1]<=(-2))),which(res[,'adj.P.Val']<=0.05))
#     res <- res[res[,'adj.P.Val']<0.05,]
	  ind <-intersect(c(which(res[,1]>=(cutoff[cut])),which(res[,1]<=(-cutoff[cut]))),which(res[,'adj.P.Val']<=0.05))
      res <- res[ind,]
	  rdir = paste0('/home/wucheng/imputation/deg/10x_5cl/bulk_cutoff/','bulk',cutoff[cut])
      dir.create(rdir,showWarnings = F, recursive = T)
      saveRDS(res,paste0('/home/wucheng/imputation/deg/10x_5cl/bulk_cutoff/','bulk',cutoff[cut],'/',i,'_',j,'_diffgene.rds'))
  }
}
}

allmtd = list.files('/home/wucheng/imputation/deg/10x_5cl/bulk1/')
ove <- sapply(allmtd, function(mtd){
res = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/bulk1/',mtd))
mtd <-length(res)
})


allmtd = list.files('/home/wucheng/imputation/deg/10x_5cl/bulk_cutoff/')
allf <-list.files('/home/wucheng/imputation/deg/10x_5cl/bulk_cutoff/bulk0')

ove <- sapply(allmtd, function(mtd){
      print(mtd)
  sapply(allf,function(f) {
      print(f)
      res = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/bulk_cutoff/',mtd,'/',f))
	 tmp <- nrow(res)
  })
})
saveRDS(ove, '/home/wucheng/imputation/deg/10x_5cl/bulk_cutoff/cutoff.rds')




###########

library(reshape2)
high1 <-melt(ove_high1)
high1$V1 <-"top1"
high2 <-melt(ove_high2)
high2$V1 <-"top2"
high3 <-melt(ove_high3)
high3$V1 <-"top3"
high4 <-melt(ove_high4)
high4$V1 <-"top4"
high5 <-melt(ove_high5)
high5$V1 <-"top5"
high6 <-melt(ove_high6)
high6$V1 <-"top6"
high7 <-melt(ove_high7)
high7$V1 <-"top7"
high8 <-melt(ove_high8)
high8$V1 <-"top8"
high9 <-melt(ove_high9)
high9$V1 <-"top9"
high10 <-melt(ove_high10)
high10$V1 <-"top10"
high <-rbind(high1,high2,high3,high4,high5,high6,high7,high8,high9,high10)
saveRDS(high, "/home/wucheng/imputation/deg/10x_5cl/high_low1.5/MAST/high.rds")

high <-readRDS("C:/Users/ADMIN/Desktop/imputation/DEG/high.rds")
index <-which(high[,1]=="Raw")
index1 <-which(high[,1]=="DISC")
high <-high[c(index,index1),]
levels <-c("top1","top2","top3","top4","top5","top6","top7","top8","top9","top10")
dodge <- position_dodge(width = 0.4)
ggplot(high, aes(x=factor(V1,levels=levels), y=value,fill = Var1))  + stat_boxplot(geom="errorbar",width=0.15,position = dodge)+geom_boxplot(width=0.4)+
ylim(0,1)+theme(legend.title=element_blank())+theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+   
labs(x="foldchange interval", y = "overlap")+theme_classic()+ scale_fill_manual(values = brewer.pal(3,"Set1")[c(1,2)])


mast <-readRDS("C:/Users/ADMIN/Desktop/imputation/DEG/MAST_bulk_sc_diffgene_overlaps_morelessExprGene.rds")
library(reshape2)
library(ggplot2)
df <-melt(mast)
p <-ggplot(df, aes(x=factor(Var1), y=value, colour=Var2,group=Var2)) + geom_line(size=1)+geom_point(size=2)
p<- p + labs( x="cutoff", y = "precent")+theme_classic()
p+ ylim(0,0.85)+theme(legend.title=element_blank())   
p

mast <-readRDS("C:/Users/ADMIN/Desktop/imputation/DEG/wilcox_bulk_sc_diffgene_overlaps_morelessExprGene.rds")
library(reshape2)
library(ggplot2)
df <-melt(mast)
p <-ggplot(df, aes(x=factor(Var1), y=value, colour=Var2,group=Var2)) + geom_line(size=1)+geom_point(size=2)
p<- p + labs( x="cutoff", y = "precent")+theme_classic()
p+ ylim(0,0.85)+theme(legend.title=element_blank())   
p
































