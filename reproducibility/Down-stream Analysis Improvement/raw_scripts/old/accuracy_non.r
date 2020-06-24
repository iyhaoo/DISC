data <-read.table("/home/wucheng/imputation/split-seq/cluster/DeSCI_2.6.0.0.15/26/non_neurous/markers.txt",header=T,row.names=1)

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
a
b
c
d
e
f
g
h

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
meta <-read.table("/home/wucheng/imputation/split-seq/cluster/DeSCI_2.6.0.0.15/26/non_neurous/meta.data.txt",header=T,row.names=1)
clu <-length(unique(meta[,5]))
A <-c(0:(clu-1))
Others <-setdiff(A,hebing)
Oth <-setdiff(A,unique(c(a,b,c,d,e,f,g,h)))
Oligo1
OPC1
Immune1
VLMC1
Vasc1
Astrocyte1
Epend1
OEC1
Others





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
data <-read.table("/home/wucheng/imputation/split-seq/cluster/DeSCI_2.6.0.0.15/26/non_neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/cluster/DeSCI_2.6.0.0.15/26/non_neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
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
write.table(cc,"/home/wucheng/imputation/split-seq/cluster/DeSCI_2.6.0.0.15/26/non_neurous/Accuracy.txt",col.names=F,sep="\t")
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
write.table(res,"/home/wucheng/imputation/split-seq/cluster/DeSCI_2.6.0.0.15/26/non_neurous/Recall.txt",row.names=T,col.names=T,sep="\t")