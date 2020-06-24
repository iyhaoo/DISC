data <-read.table("/home/wucheng/imputation/split-seq/cluster/DeSCI_2.6.0.0.15/26/neurous/markers.txt",header=T,row.names=1)
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

a
b
c
d
e
f
g
h
i
j
k
l
m
n

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

Olfa
Stri
Cort
Rost
Thal
Cere
Medu
Basal
Hipp
Spin
Mirg

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
dat <-read.table("/home/wucheng/imputation/split-seq/cluster/DeSCI_2.6.0.0.15/26/neurous/meta.data.txt",header=T,row.names=1)
dat <-dat[,-5]
clu <-length(unique(dat[,5]))
A<-c(0:(clu-1))

Others <-setdiff(A,hebing)

Olfa1
Stri1
Cort1
Rost1
Thal1
Cere1
Medu1
Basal1
Hipp1
Spin1
Mirg1
Others
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
data <-read.table("/home/wucheng/imputation/split-seq/cluster/DeSCI_2.6.0.0.15/26/neurous/meta.data.txt",header=T,row.names=1)
data <-data[,-5]
clu <-length(unique(data[,5]))
y <-c(0:(clu-1))
y
data$V7 <-plyr::mapvalues(x = data[,5], from = y, to = x)
write.table(data,"/home/wucheng/imputation/split-seq/cluster/DeSCI_2.6.0.0.15/26/neurous/meta.dataxin.txt",row.names=TRUE,col.names=TRUE)
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
write.table(cc,"/home/wucheng/imputation/split-seq/cluster/DeSCI_2.6.0.0.15/26/neurous/Accuracy.txt",col.names=F,sep="\t")

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
write.table(res,"/home/wucheng/imputation/split-seq/cluster/DeSCI_2.6.0.0.15/26/neurous/Recall.txt",row.names=T,col.names=T,sep="\t")