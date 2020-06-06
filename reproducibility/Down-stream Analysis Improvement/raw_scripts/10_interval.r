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
      fc <- sort(abs(as.matrix(raw_res)[,2]),decreasing=T)      
      highid1= names(which(fc>quantile(fc,0.9)))## express in more cells
	  highid2 = intersect(names(which(fc<quantile(fc,0.9))),names(which(fc>quantile(fc,0.8))))## express in more cells
      highid3 = intersect(names(which(fc>quantile(fc,0.8))),names(which(fc>quantile(fc,0.7))))## express in more cells
      highid4 = intersect(names(which(fc>quantile(fc,0.7))),names(which(fc>quantile(fc,0.6))))## express in more cells
      highid5 = intersect(names(which(fc>quantile(fc,0.6))),names(which(fc>quantile(fc,0.5))))## express in more cells
      highid6 = intersect(names(which(fc>quantile(fc,0.5))),names(which(fc>quantile(fc,0.4))))## express in more cells
      highid7 = intersect(names(which(fc>quantile(fc,0.4))),names(which(fc>quantile(fc,0.3))))## express in more cells
      highid8 = intersect(names(which(fc>quantile(fc,0.3))),names(which(fc>quantile(fc,0.2))))## express in more cells
      highid9 = intersect(names(which(fc>quantile(fc,0.2))),names(which(fc>quantile(fc,0.1))))## express in more cells
      lowid = names(which(fc<quantile(fc,0.1)))## less
      gs = readRDS(paste0('/home/wucheng/imputation/deg/10x_5cl/bulk1/',f))
      tmpres <- res[rownames(res) %in% highid1,]
      tmp1 <- mean(sapply(c(1:10)*10,function(i) {
        mean(rownames(tmpres[1:i,]) %in% gs)  ## discuss
      }))
	  tmpres <- res[rownames(res) %in% highid2,]
      tmp2 <- mean(sapply(c(1:10)*10,function(i) {
        mean(rownames(tmpres[1:i,]) %in% gs)  ## discuss
      }))
      tmpres <- res[rownames(res) %in% highid3,]
      tmp3 <- mean(sapply(c(1:10)*10,function(i) {
        mean(rownames(tmpres[1:i,]) %in% gs)  ## discuss
      }))
      tmpres <- res[rownames(res) %in% highid4,]
      tmp4 <- mean(sapply(c(1:10)*10,function(i) {
        mean(rownames(tmpres[1:i,]) %in% gs)  ## discuss
      }))
      tmpres <- res[rownames(res) %in% highid5,]
      tmp5 <- mean(sapply(c(1:10)*10,function(i) {
        mean(rownames(tmpres[1:i,]) %in% gs)  ## discuss
      }))
      tmpres <- res[rownames(res) %in% highid6,]
      tmp6 <- mean(sapply(c(1:10)*10,function(i) {
        mean(rownames(tmpres[1:i,]) %in% gs)  ## discuss
      }))
      tmpres <- res[rownames(res) %in% highid7,]
      tmp7 <- mean(sapply(c(1:10)*10,function(i) {
        mean(rownames(tmpres[1:i,]) %in% gs)  ## discuss
      }))
      tmpres <- res[rownames(res) %in% highid8,]
      tmp8 <- mean(sapply(c(1:10)*10,function(i) {
        mean(rownames(tmpres[1:i,]) %in% gs)  ## discuss
      }))
	  tmpres <- res[rownames(res) %in% highid9,]
      tmp9 <- mean(sapply(c(1:10)*10,function(i) {
        mean(rownames(tmpres[1:i,]) %in% gs)  ## discuss
      }))
      tmpres <- res[rownames(res) %in% lowid,]
      tmp10 <- mean(sapply(c(1:10)*10,function(i) {
        mean(rownames(tmpres[1:i,]) %in% gs)  ## discuss
      }),na.rm=T)
      c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10)
    } else {
      return(c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
    }
  }))
}) ## first 10 high, last 10 low
ove_high1 = t(ove[1:10,])
ove_high2 = t(ove[11:20,])
ove_high3 = t(ove[21:30,])
ove_high4 = t(ove[31:40,])
ove_high5 = t(ove[41:50,])
ove_high6 = t(ove[51:60,])
ove_high7 = t(ove[61:70,])
ove_high8 = t(ove[71:80,])
ove_high9 = t(ove[81:90,])
ove_high10 = t(ove[91:100,])
colnames(ove_high1) <- colnames(ove_high2)<-colnames(ove_high3) <- colnames(ove_high4) <-colnames(ove_high5) <- colnames(ove_high6)<-colnames(ove_high7) <- colnames(ove_high8) <-colnames(ove_high9) <- colnames(ove_high10) <- sub('.rds','', allf)
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
