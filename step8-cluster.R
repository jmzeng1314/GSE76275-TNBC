## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2019-03-09 19:56:32
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2019-03-09  First version
###
### ---------------

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
if(F){
  load(file = 'step1-output.Rdata')
  # 每次都要检测数据
  dat[1:4,1:4]  
  library(hgu133plus2.db)
  ids=toTable(hgu133plus2SYMBOL) #toTable这个函数：通过看hgu133plus2.db这个包的说明书知道提取probe_id（探针名）和symbol（基因名）的对应关系的表达矩阵的函数为toTable
  head(ids) #head为查看前六行
  dat=dat[ids$probe_id,] #ids提取出probe_id这列，这列的每行都为一个探针，接着在dat这个矩阵中，按照刚刚取出的探针所在的行，再取出来组成一个新的矩阵dat，此操纵为取出与注视ids相对于的dat
  dat[1:4,1:4] 
  ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
  ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
  dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
  rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
  dat[1:4,1:4]  #保留每个基因ID第一次出现的信息
  
  dat=dat[,group_list=='TNBC']
  ## 首先聚类，发现需要 剔除 GSM1974605 这个样本
  dat=dat[,colnames(dat)!='GSM1974605']
  save(dat,file = 'GSE76275-TNBC-for-cluster.Rdata')
}

if(F){
  x=1:10
  y=2*x
  z=rnorm(10) #产生10个随机数
  tmp=data.frame(x,y,z)#转换为数据框
  dist(tmp) #dist()函数计算变量间距离
  head(tmp)
  dist(t( tmp ))
  dist(t(scale(tmp)))
  #scale函数进行数据的中心化和标准化
  cor(tmp) #cov,协方差系数,可用来计算相关性
  
}

load(file = 'GSE76275-TNBC-for-cluster.Rdata')
dat[1:4,1:4] 
dim(dat)
M=cor(dat)
pheatmap::pheatmap(M)

d=dist(t(dat))
pheatmap::pheatmap(d)
pheatmap::pheatmap(t(scale(d)),show_rownames = F)
hc=hclust(d) #层次聚类
?plot.hclust
plot(hc,labels =F)
# 对 dat 进行 scale 
hc=hclust(dist(t(dat))) 
plot(hc,labels =F)
clus = cutree(hc, 4)#分类
group_list=as.factor(clus)#转换为因子
table(group_list) #统计频率

# BiocManager::install("ConsensusClusterPlus", version = "3.8")

d=dat
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:5000],] #取top5000的基因
#rev函数可以实现向量或矩阵的翻转；当进行矩阵的翻转时,rev函数将矩阵当做一个向量处理，矩阵转换为向量的原则为列优先
d = sweep(d,1, apply(d,1,median,na.rm=T))
title="./" #所有的图片以及数据都会输出到这里的
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(as.matrix(d),maxK=6,reps=50,
                               pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",
                               seed=1262118388.71279,plot="png")
# 可以看到这4组和hclust层次聚类的结果差不多。
# 参考：https://uc-r.github.io/kmeans_clustering

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
suppressPackageStartupMessages(library("sigclust2"))
data_pc <- prcomp(t(dat))
df=data_pc$x ## 主成分对样本进行分类。
# 后面需要继续更新。
distance <- get_dist(t(df))
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
k2 <- kmeans(df, centers = 2, nstart = 25)
str(k2)
fviz_cluster(k2, data = df)

k3 <- kmeans(df, centers = 3, nstart = 25)
k4 <- kmeans(df, centers = 4, nstart = 25)
k5 <- kmeans(df, centers = 5, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = df) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = df) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = df) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = df) + ggtitle("k = 5")

library(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2)


# install.packages("remotes")
# BiocManager::install("pkimes/sigclust2")
suppressPackageStartupMessages(library("sigclust2"))
data_pc <- prcomp(t(dat))
par(mfrow=c(1, 2))
plot(data_pc$x[, 2], data_pc$x[, 1], xlab="PC2", ylab="PC1")
plot(data_pc$x[, 3], data_pc$x[, 1], xlab="PC3", ylab="PC1") 
