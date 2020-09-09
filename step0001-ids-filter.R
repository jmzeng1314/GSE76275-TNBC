rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F) #避免character类型自动转化为factor类型
load(file = 'step1-output.Rdata')
# 每次都要检测数据
dat[1:4,1:4]  
library(hgu133plus2.db)
ids=toTable(hgu133plus2SYMBOL) #toTable这个函数：通过看hgu133plus2.db这个包的说明书知道提取probe_id（探针名）和symbol（基因名）的对应关系的表达矩阵的函数为toTable
head(ids) #head为查看前六行
dat=dat[ids$probe_id,] #ids提取出probe_id这列，这列的每行都为一个探针，接着在dat这个矩阵中，按照刚刚取出的探针所在的行，再取出来组成一个新的矩阵dat，此操纵为取出与注视ids相对于的dat
#保证ids矩阵和dat矩阵长度相等
dat[1:4,1:4] 
ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$symbol),]#将symbol这一列去除重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
##确保两个矩阵长度一致
dat[1:4,1:4]  #保留每个基因ID第一次出现的信息

boxplot(apply(dat,1,mad))#mad：绝对中位差
pheatmap::pheatmap(dat[order(apply(dat,1,mad), decreasing = T)[1:50],])#mad：绝对中位差
pheatmap::pheatmap(dat[order(apply(dat,1,mad), decreasing = F)[1:50],])


library(corrplot) 
M <- cor( dat ) #cor计算相关系数
ac=data.frame(g=group_list)
rownames(ac)=colnames(M) 
pheatmap::pheatmap(M,show_colnames =F,show_rownames = F,
                   annotation_col=ac)
pheatmap::pheatmap(M,show_colnames =F,show_rownames = F,
         annotation_col=ac,filename = 'heatmap_cor.png')
dev.off()

