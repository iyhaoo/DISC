library(pheatmap)
library(reshape2)
library(ggplot2)
library(plyr)
df <-as.matrix(read.table("C:/Users/ADMIN/Desktop/imputation/figure/figure4/retina-merge/ds0.3_Jaccard.txt",row.names=1,header=T,sep="\t"))
aa <-read.table("C:/Users/ADMIN/Desktop/imputation/figure/figure4/retina-merge1/deep_Jaccard3.txt",row.names=1,header=T,sep="\t")
cc <-read.table("C:/Users/ADMIN/Desktop/imputation/figure/figure4/retina-merge1/VIPER_Jaccard3.txt",row.names=1,header=T,sep="\t")
aa <-rbind(aa[6,],aa[10,],aa[1,],aa[4,],aa[3,],aa[8,],aa[2,],aa[5,],aa[12,],aa[9,],aa[7,],aa[8,])
cc <-rbind(cc[6,],cc[10,],cc[1,],cc[4,],cc[3,],cc[8,],cc[2,],cc[5,],cc[12,],cc[9,],cc[7,],cc[8,])
df1 <-rbind(df[1:7,],aa[,1],rep(0,12),cc[,1])
df1 <-df1[-3,]
rownames(df1) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df1) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","Fibroblasts","VE","Pericytes","Microglia","Rods")

df2 <-rbind(df[8:14,],aa[,2],rep(0,12),cc[,2])
df2 <-df2[-3,]
rownames(df2) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df2) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","Fibroblasts","VE","Pericytes","Microglia","Rods")

df3 <-rbind(df[15:21,],aa[,4],rep(0,12),cc[,4])
df3 <-df3[-3,]
rownames(df3) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df3) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","Fibroblasts","VE","Pericytes","Microglia","Rods")

df4 <-rbind(df[22:28,],aa[,3],rep(0,12),cc[,3])
df4 <-df4[-3,]
rownames(df4) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df4) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","Fibroblasts","VE","Pericytes","Microglia","Rods")

df5 <-rbind(df[29:35,],aa[,5],rep(0,12),cc[,5])
df5 <-df5[-3,]
rownames(df5) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df5) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","Fibroblasts","VE","Pericytes","Microglia","Rods")
dff <-(df1+df2+df3+df4+df5)/5
dff <-dff[c(1,2,8,9,4,3,7,5,6),]
pdf(file="C:/Users/ADMIN/Desktop/imputation/figure/figure4/Ja_Heat_RETINA0.3.pdf",width=13.5,height=9)
pheatmap(dff,scale = "none",color = colorRampPalette(colors = c("yellow","white","red"))(100),cluster_rows = FALSE,cluster_cols = FALSE,display_numbers = TRUE,labels_col=1)
dev.off()
pdf(file="C:/Users/ADMIN/Desktop/imputation/figure/figure4/Ja_Heat_RETINA0.31.pdf",width=13.5,height=9)
pheatmap(dff,scale = "none",color = colorRampPalette(colors = c("white","yellow","red"))(100),cluster_rows = FALSE,cluster_cols = FALSE,display_numbers = TRUE,labels_col=1)
dev.off()
##########################Jaccard0.3-PBMC
df <-as.matrix(read.table("C:/Users/ADMIN/Desktop/imputation/figure/figure4/pbmc0.3-merge/ds_Jaccard.txt",row.names=1,header=T,sep="\t"))
aa <-read.table("C:/Users/ADMIN/Desktop/imputation/figure/figure4/pbmc-merge1/deep_Jaccard3.txt",row.names=1,header=T,sep="\t")
bb <-read.table("C:/Users/ADMIN/Desktop/imputation/figure/figure4/pbmc-merge1/ScI_Jaccard3.txt",row.names=1,header=T,sep="\t")
cc <-read.table("C:/Users/ADMIN/Desktop/imputation/figure/figure4/pbmc-merge1/VIPER_Jaccard3.txt",row.names=1,header=T,sep="\t")
aa <-rbind(aa[3,],aa[2,],aa[1,],aa[4,],aa[6,],aa[7,],aa[5,],aa[8,])
bb <-rbind(bb[3,],bb[2,],bb[1,],bb[4,],bb[6,],bb[7,],bb[5,],bb[8,])
cc <-rbind(cc[3,],cc[2,],cc[1,],cc[4,],cc[6,],cc[7,],cc[5,],cc[8,])

df1 <-rbind(df[1:7,],aa[,1],bb[,1],cc[,1])
df1 <-df1[-3,]
rownames(df1) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df1) <-c("CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")

df2 <-rbind(df[8:14,],aa[,2],bb[,2],cc[,2])
df2 <-df2[-3,]
rownames(df2) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df2) <-c("CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")

df3 <-rbind(df[15:21,],aa[,3],bb[,3],cc[,3])
df3 <-df3[-3,]
rownames(df3) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df3) <-c("CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")

df4 <-rbind(df[22:28,],aa[,4],bb[,4],cc[,4])
df4 <-df4[-3,]
rownames(df4) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df4) <-c("CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")

df5 <-rbind(df[29:35,],aa[,5],bb[,5],cc[,5])
df5 <-df5[-3,]
rownames(df5) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df5) <-c("CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
dff <-(df1+df2+df3+df4+df5)/5
dff <-dff[c(1,2,8,9,4,3,7,5,6),]
pdf(file="C:/Users/ADMIN/Desktop/imputation/figure/figure4/Ja_Heat_PBMC0.3.pdf",width=9.5,height=9)
pheatmap(dff,scale = "none",color = colorRampPalette(colors = c("yellow","white","red"))(100),cluster_rows = FALSE,cluster_cols = FALSE,display_numbers = TRUE,labels_col=1)
dev.off()

pdf(file="C:/Users/ADMIN/Desktop/imputation/figure/figure4/Ja_Heat_PBMC0.31.pdf",width=9.5,height=9)
pheatmap(dff,scale = "none",color = colorRampPalette(colors = c("white","yellow","red"))(100),cluster_rows = FALSE,cluster_cols = FALSE,display_numbers = TRUE,labels_col=1)
dev.off()


###########
df <-as.matrix(read.table("C:/Users/ADMIN/Desktop/imputation/figure/figure4/split-merge1/top50_merge_neurous_ds0.3_Jaccard.txt",row.names=1,header=T,sep="\t"))
df1 <-rbind(df[1:8,],rep(0,11),rep(0,11))
df1 <-df1[-3,]
rownames(df1) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df1) <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus","Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","MI")

df2 <-rbind(df[9:16,],rep(0,11),rep(0,11))
df2 <-df2[-3,]
rownames(df2) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df2) <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus","Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","MI")

df3 <-rbind(df[17:24,],rep(0,11),rep(0,11))
df3 <-df3[-3,]
rownames(df3) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df3) <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus","Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","MI")

df4 <-rbind(df[25:32,],rep(0,11),rep(0,11))
df4 <-df4[-3,]
rownames(df4) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df4) <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus","Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","MI")

df5 <-rbind(df[33:40,],rep(0,11),rep(0,11))
df5 <-df5[-3,]
rownames(df5) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df5) <-c("Olfactory Bulb","Striatum","Cortex","Rostral Midbrain","Thalamus","Cerebellum","Medulla","Basal Ganglia","Hippocampus","Spinalcord","MI")
dff <-(df1+df2+df3+df4+df5)/5
dff <-dff[c(1,2,8,9,4,3,7,5,6),]
pdf(file="C:/Users/ADMIN/Desktop/imputation/figure/figure4/Ja_Heat_neur0.3.pdf",width=12.5,height=9)
pheatmap(dff,scale = "none",color = colorRampPalette(colors = c("yellow","white","red"))(100),cluster_rows = FALSE,cluster_cols = FALSE,display_numbers = TRUE,labels_col=1)
dev.off()
dff <-rbind(df1,df2,df3,df4,df5)
write.table(dff,"C:/Users/ADMIN/Desktop/imputation/figure/figure4/neur_Jaccard0.3.txt",row.names=TRUE, col.names=TRUE,sep="\t")

pdf(file="C:/Users/ADMIN/Desktop/imputation/figure/figure4/Ja_Heat_neur0.31.pdf",width=12.5,height=9)
pheatmap(dff,scale = "none",color = colorRampPalette(colors = c("white","yellow","red"))(100),cluster_rows = FALSE,cluster_cols = FALSE,display_numbers = TRUE,labels_col=1)
dev.off()


########
df <-as.matrix(read.table("C:/Users/ADMIN/Desktop/imputation/figure/figure4/split-merge1/top50_merge_non_neurous_ds0.3_Jaccard.txt",row.names=1,header=T,sep="\t"))
df1 <-rbind(df[1:8,],rep(0,8),rep(0,8))
df1 <-df1[-3,]
rownames(df1) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df1) <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")

df2 <-rbind(df[9:16,],rep(0,8),rep(0,8))
df2 <-df2[-3,]
rownames(df2) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df2) <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")

df3 <-rbind(df[17:24,],rep(0,8),rep(0,8))
df3 <-df3[-3,]
rownames(df3) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df3) <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")

df4 <-rbind(df[25:32,],rep(0,8),rep(0,8))
df4 <-df4[-3,]
rownames(df4) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df4) <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")

df5 <-rbind(df[33:40,],rep(0,8),rep(0,8))
df5 <-df5[-3,]
rownames(df5) <-c("Observed", "DISC","DCA","MAGIC","scScope","scVI","deepImpute","scImpute","VIPER")
colnames(df5) <-c("Astrocyte","Epend","Immune","OEC","Oligo","OPC","Vasc","VLMC")
dff <-(df1+df2+df3+df4+df5)/5
dff <-dff[c(1,2,8,9,4,3,7,5,6),]
pdf(file="C:/Users/ADMIN/Desktop/imputation/figure/figure4/Ja_Heat_non_neur0.3.pdf",width=9.5,height=9)
pheatmap(dff,scale = "none",color = colorRampPalette(colors = c("yellow","white","red"))(100),cluster_rows = FALSE,cluster_cols = FALSE,display_numbers = TRUE,labels_col=1)
dev.off()
dff <-rbind(df1,df2,df3,df4,df5)
write.table(dff,"C:/Users/ADMIN/Desktop/imputation/figure/figure4/non-neur_Jaccard0.3.txt",row.names=TRUE, col.names=TRUE,sep="\t")

pdf(file="C:/Users/ADMIN/Desktop/imputation/figure/figure4/Ja_Heat_nneur0.31.pdf",width=9.5,height=9)
pheatmap(dff,scale = "none",color = colorRampPalette(colors = c("white","yellow","red"))(100),cluster_rows = FALSE,cluster_cols = FALSE,display_numbers = TRUE,labels_col=1)
dev.off()
