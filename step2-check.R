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
### Update Log: 2020-01-13  second version
###
### ---------------



rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'step1-output.Rdata')
# 每次都要检测数据
dat[1:4,1:4]
## 下面是画PCA的必须操作，需要看说明书。
dat=t(dat)#画PCA图时要求是行名是样本名，列名时探针名，因此此时需要转换
dat=as.data.frame(dat)#将matrix转换为data.frame
dat=cbind(dat,group_list) #cbind横向追加，即将分组信息追加到最后一列
dat[1:4,54672:54676]
library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") 
# The variable group_list (index = 54676) is removed
# before PCA analysis
# 这一步稍微有点耗时。
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
dat.pca
## 建议自己耗费2个小时，详细了解 dat.pca 这个对象里面的数据，帮助你理解什么是PCA分析。
head(dat.pca$ind$coord)
plot(dat.pca$ind$coord[,1:2])
# 这样绘图实在是太丑，所以需要其R包自带高级函数。

#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list, # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
# 其实这里完全没有必要把 group_list 变量添加到 dat 这个表达矩阵里面。
ggsave('all_samples_PCA.png')

### PCA 分析的结果，可以使用kmeans来进行聚类，比如下面挑选前两个主成分的。
df=dat.pca$ind$coord
k3 <- kmeans(df[,1:2], centers = 3, nstart = 25)
fviz_cluster(k3, geom = "point", data = df) + ggtitle("k = 3")
# 前两个主成分很明显，可以把样本分成3组，其中TNBC是两组，另外的nonTNBC独立成为一个组别。
# 也可以看看不同样本如果只使用前5个主成分，它们的距离如何。
distance <- get_dist(df)
fviz_dist(distance, 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
plot(hclust(dist(df)))
# 这些聚类的知识点，到第8步骤再谈。


rm(list = ls())  ## 魔幻操作，一键清空~
load(file = 'step1-output.Rdata')
#此步为一个小插曲，即计算一下从第一行开始计算每一行的sd值，直到最后一行所需要的时间
dat[1:4,1:4]
dat[1,]
system.time(apply(dat,1,sd))
#system.time返回CPU的使用时间
system.time(for(i in 1:nrow(dat)){     
  #i是一个变量，i从第一行开始取sd值，不知道多少行，直到最后一行（nrow），所以是i in 1:nrow(dat)，
  sd(dat[i,])
})

#取出方差最大的1000个基因
cg=names(tail(sort(apply(dat,1,sd)),1000))
#apply按行（'1'是按行取，'2'是按列取）取每一行的方差，从小到大排序，取最大的1000个
#sort 排序，默认升序排序。tail从后向前显示，默认后6行
library(pheatmap)
pheatmap(dat[cg,],show_colnames =F,show_rownames = F) 
#把那些提取出来的方差最大的1000个基因所在的每一行取出，组合起来为一个新的表达矩阵
n=t(scale(t(dat[cg,]))) 
#通过“scale”对log-ratio数值进行归一化，现在的dat是行名为探针，列名为样本名，
#由于scale这个函数应用在不同组数据间存在差异时，需要行名为样本，因此需要用t(dat[cg,])来转换，最后再转换回来
#t 对矩阵进行行列转换
n[n>2]=2 #限定上限，使表达量大于2的等于2
n[n< -2]= -2 #限定下限，使表达量小于-2的等于-2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(g=group_list) #把character向量转换为数据框
rownames(ac)=colnames(n) 
#把ac的行名给到n的列名，即对每一个探针标记上分组信息（是'noTNBC'还是'TNBC'）
## 可以看到TNBC具有一定的异质性，拿它来区分乳腺癌亚型指导临床治疗还是略显粗糙。
## 首先TNBC本身就可以继续分组，其次那些不是TNBC的病人跟T属于NBC病人的表达量是无法区分的。
pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac,filename = 'heatmap_top1000_sd_all.png')

## 单独看 noTNBC 的 top 1000 基因的热图
#取出方差最大的1000个基因
table(group_list)
n_dat=dat[,group_list=='noTNBC']
cg=names(tail(sort(apply(n_dat,1,sd)),1000)) 
library(pheatmap)
n=t(scale(t(n_dat[cg,])))
n[n>2]=2 #限定上限，使表达量大于2的等于2
n[n< -2]= -2 #限定下限，使表达量小于-2的等于-2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F,
          filename = 'heatmap_top1000_sd_noTNBC.png')
# 可以看到 虽然都属于 noTNBC ，但是其实里面也是有异质性的。 
dev.off() 

## 单独看 TNBC 的 top 1000 基因的热图
#取出方差最大的1000个基因
table(group_list)
tnbc_dat=dat[,group_list=='TNBC']
cg=names(tail(sort(apply(tnbc_dat,1,sd)),1000)) 
library(pheatmap)
n=t(scale(t(tnbc_dat[cg,])))
n[n>2]=2 #限定上限，使表达量大于2的等于2
n[n< -2]= -2 #限定下限，使表达量小于-2的等于-2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
pheatmap(n,show_colnames =F,show_rownames = F,
         filename = 'heatmap_top1000_sd_TNBC.png')
# 可以看到 虽然都属于 TNBC，但是其实里面也是有异质性的。 
dev.off()





