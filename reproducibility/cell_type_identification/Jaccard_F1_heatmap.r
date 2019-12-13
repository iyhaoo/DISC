
df<-as.matrix(read.table("C:/Users/Xie-lab/Desktop/imputation/retina_analysis/ds1/merge/ds0.3_Jaccard.txt",row.names=1,header=T,sep="\t"))
df <-df[,-8]
df1 <-df[1:7,]
rownames(df1) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df1) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df1 <-melt(df1)
df2 <-df[8:14,]
rownames(df2) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI") 
colnames(df2) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df2 <-melt(df2)
df3 <-df[15:21,]
rownames(df3) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df3) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df3 <-melt(df3)
df4 <-df[22:28,]
rownames(df4) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df4) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df4 <-melt(df4)
df5 <-df[29:35,]
rownames(df5) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df5) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df5 <-melt(df5)

dff <-rbind(df1,df2,df3,df4,df5)
dff1 <- data_summary(dff, varname="value", groupnames=c("Var1", "Var2"))

dff1$Var1 <-factor(dff1$Var1,levels =c("scVI","scScope","DCA","MAGIC","SAVER","DeSCI","RAW"))
p <-ggplot(data = dff1, aes(x=Var2, y=Var1, fill=value))+ geom_tile()+ labs(x="Method", y = "Cell Type")+theme_classic()
p <-p +theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+theme(axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"))
p +scale_fill_gradient(low='white',high='red')+theme(legend.title=element_blank()) +geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4)
#############
df<-as.matrix(read.table("C:/Users/Xie-lab/Desktop/imputation/retina_analysis/ds1/merge/ds0.5_Jaccard.txt",row.names=1,header=T,sep="\t"))
df <-df[,-8]
df1 <-df[1:7,]
rownames(df1) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df1) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df1 <-melt(df1)
df2 <-df[8:14,]
rownames(df2) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI") 
colnames(df2) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df2 <-melt(df2)
df3 <-df[15:21,]
rownames(df3) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df3) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df3 <-melt(df3)
df4 <-df[22:28,]
rownames(df4) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df4) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df4 <-melt(df4)
df5 <-df[29:35,]
rownames(df5) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df5) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df5 <-melt(df5)

dff <-rbind(df1,df2,df3,df4,df5)
dff1 <- data_summary(dff, varname="value", groupnames=c("Var1", "Var2"))

dff1$Var1 <-factor(dff1$Var1,levels =c("scVI","scScope","DCA","MAGIC","SAVER","DeSCI","RAW"))
p <-ggplot(data = dff1, aes(x=Var2, y=Var1, fill=value))+ geom_tile()+ labs(x="Method", y = "Cell Type")+theme_classic()
p <-p +theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+theme(axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"))
p +scale_fill_gradient(low='white',high='red')+theme(legend.title=element_blank()) +geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4)



###########

df<-as.matrix(read.table("C:/Users/Xie-lab/Desktop/imputation/retina_analysis/ds1/merge/ds0.3_F1.txt",row.names=1,header=T,sep="\t"))
df <-df[,-8]
df1 <-df[1:7,]
rownames(df1) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df1) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df1 <-melt(df1)
df2 <-df[8:14,]
rownames(df2) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI") 
colnames(df2) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df2 <-melt(df2)
df3 <-df[15:21,]
rownames(df3) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df3) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df3 <-melt(df3)
df4 <-df[22:28,]
rownames(df4) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df4) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df4 <-melt(df4)
df5 <-df[29:35,]
rownames(df5) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df5) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df5 <-melt(df5)

dff <-rbind(df1,df2,df3,df4,df5)
dff1 <- data_summary(dff, varname="value", groupnames=c("Var1", "Var2"))

dff1$Var1 <-factor(dff1$Var1,levels =c("scVI","scScope","DCA","MAGIC","SAVER","DeSCI","RAW"))
p <-ggplot(data = dff1, aes(x=Var2, y=Var1, fill=value))+ geom_tile()+ labs(x="Method", y = "Cell Type")+theme_classic()
p <-p +theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+theme(axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"))
p +scale_fill_gradient(low='white',high='red')+theme(legend.title=element_blank()) +geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4)
#############
df<-as.matrix(read.table("C:/Users/Xie-lab/Desktop/imputation/retina_analysis/ds1/merge/ds0.5_F1.txt",row.names=1,header=T,sep="\t"))
df <-df[,-8]
df1 <-df[1:7,]
rownames(df1) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df1) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df1 <-melt(df1)
df2 <-df[8:14,]
rownames(df2) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI") 
colnames(df2) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df2 <-melt(df2)
df3 <-df[15:21,]
rownames(df3) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df3) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df3 <-melt(df3)
df4 <-df[22:28,]
rownames(df4) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df4) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df4 <-melt(df4)
df5 <-df[29:35,]
rownames(df5) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df5) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df5 <-melt(df5)

dff <-rbind(df1,df2,df3,df4,df5)
dff1 <- data_summary(dff, varname="value", groupnames=c("Var1", "Var2"))

dff1$Var1 <-factor(dff1$Var1,levels =c("scVI","scScope","DCA","MAGIC","SAVER","DeSCI","RAW"))
p <-ggplot(data = dff1, aes(x=Var2, y=Var1, fill=value))+ geom_tile()+ labs(x="Method", y = "Cell Type")+theme_classic()
p <-p +theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+theme(axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"))
p +scale_fill_gradient(low='white',high='red')+theme(legend.title=element_blank()) +geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4)


######
df<-as.matrix(read.table("C:/Users/Xie-lab/Desktop/imputation/retina_analysis/ds1/merge_ds/ds0.7_Jaccard.txt",row.names=1,header=T,sep="\t"))
df <-df[,-8]
df1 <-df[1:7,]
rownames(df1) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df1) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df1 <-melt(df1)
df2 <-df[8:14,]
rownames(df2) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI") 
colnames(df2) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df2 <-melt(df2)
df3 <-df[15:21,]
rownames(df3) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df3) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df3 <-melt(df3)
df4 <-df[22:28,]
rownames(df4) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df4) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df4 <-melt(df4)
df5 <-df[29:35,]
rownames(df5) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df5) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df5 <-melt(df5)

dff <-rbind(df1,df2,df3,df4,df5)
dff1 <- data_summary(dff, varname="value", groupnames=c("Var1", "Var2"))

dff1$Var1 <-factor(dff1$Var1,levels =c("scVI","scScope","DCA","MAGIC","SAVER","DeSCI","RAW"))
p <-ggplot(data = dff1, aes(x=Var2, y=Var1, fill=value))+ geom_tile()+ labs(x="Method", y = "Cell Type")+theme_classic()
p <-p +theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+theme(axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"))
p +scale_fill_gradient(low='white',high='red')+theme(legend.title=element_blank()) +geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4)


df<-as.matrix(read.table("C:/Users/Xie-lab/Desktop/imputation/retina_analysis/ds1/merge_ds/ds0.7_F1.txt",row.names=1,header=T,sep="\t"))
df <-df[,-8]
df1 <-df[1:7,]
rownames(df1) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df1) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df1 <-melt(df1)
df2 <-df[8:14,]
rownames(df2) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI") 
colnames(df2) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df2 <-melt(df2)
df3 <-df[15:21,]
rownames(df3) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df3) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df3 <-melt(df3)
df4 <-df[22:28,]
rownames(df4) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df4) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df4 <-melt(df4)
df5 <-df[29:35,]
rownames(df5) <-c("RAW","DeSCI","SAVER","DCA","MAGIC","scScope","scVI")
colnames(df5) <-c("HC","RGC","Amacrine cells","Cones","Bipolar","Muller glia","Astrocytes","VE","Pericytes","Microglia","Rods")
df5 <-melt(df5)

dff <-rbind(df1,df2,df3,df4,df5)
dff1 <- data_summary(dff, varname="value", groupnames=c("Var1", "Var2"))

dff1$Var1 <-factor(dff1$Var1,levels =c("scVI","scScope","DCA","MAGIC","SAVER","DeSCI","RAW"))
p <-ggplot(data = dff1, aes(x=Var2, y=Var1, fill=value))+ geom_tile()+ labs(x="Method", y = "Cell Type")+theme_classic()
p <-p +theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+theme(axis.text.y = element_text(size = 12, hjust = 1, vjust = 1, face = "bold"))
p +scale_fill_gradient(low='white',high='red')+theme(legend.title=element_blank()) +geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4)
