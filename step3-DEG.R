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
table(group_list) #table函数，查看group_list中的分组个数
#通过为每个数据集绘制箱形图，比较数据集中的数据分布
boxplot(dat[1,]~group_list) #按照group_list分组画箱线图
t.test(dat[1,]~group_list)
bp=function(g){         #定义一个函数bp，函数为{}里的内容
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
bp(dat[1,]) ## 调用上面定义好的函数，避免同样的绘图代码重复多次敲。
bp(dat[2,])
bp(dat[3,])
bp(dat[4,])


library(limma)
design=model.matrix(~factor( group_list ))#创建一个因子矩阵
#model.matrix创建一个矩阵模型（设计）
fit=lmFit(dat,design)
#lmFit 给定阵列数据为每个基因拟合线性模型
fit=eBayes(fit)
#eBayes 给出了一个微阵列线性模型拟合，通过经验贝叶斯调整标准误差到一个共同的值来计算修正后的t统计量、修正后的f统计量和微分表达式的对数概率。
## 上面是limma包用法的一种方式 
options(digits = 4) #设置全局的数字有效位数为4
#topTable(fit,coef=2,adjust='BH') 
topTable(fit,coef=2,adjust='BH') 
#coef可以是column number，也可以是column name，这样就可以指定你所感兴趣的两两比较的结果
#topTable从线性模型拟合中提取排名靠前的基因表。
bp(dat['241662_x_at',])
bp(dat['207639_at',])
deg=topTable(fit,coef=2,adjust='BH',number = Inf)

head(deg)
deg$probe_id=rownames(deg) #deg重新取一列命名为probe_id,将deg的行名给新取的这一列作为每一行的内容
library(hgu133plus2.db)
ids=toTable(hgu133plus2SYMBOL)
deg=merge(deg,ids,by='probe_id')#通过merge函数，由于deg和ids都有probe_id这一列，因此通过'probe_id'合并为新的deg

save(deg,file = 'deg.Rdata')

## for volcano 火山图
if(T){
  nrDEG=deg
  head(nrDEG)
  attach(nrDEG)
  plot(logFC,-log10(P.Value))
  library(ggpubr)
  df=nrDEG
  df$v= -log10(P.Value) #df新增加一列'v',值为-log10(P.Value)
  ggscatter(df, x = "logFC", y = "v",size=0.5)
  
  df$g=ifelse(df$P.Value>0.01,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
              ifelse( df$logFC >1.5,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
                      ifelse( df$logFC < -1.5,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
  )
  table(df$g)
  df$name=rownames(df)
  ggscatter(df, x = "logFC", y = "v",size=0.5,color = 'g')
  ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
            label = "symbol", repel = T,
            #label.select = rownames(df)[df$g != 'stable'] ,
            label.select = c('PROM1','AGR3','AGR2'), #挑选一些基因在图中显示出来
            palette = c("#00AFBB", "#E7B800", "#FC4E07") )
  
  ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
  df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                  ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
  table(df$p_c )
  ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
            palette = c("green", "red", "black") )
  
  
}


## for heatmap 热图
if(T){ 
  load(file = 'step1-output.Rdata')
  # 每次都要检测数据
  dat[1:4,1:4]
  table(group_list)
  x=deg$logFC #deg取logFC这列并将其重新赋值给x
  names(x)=deg$probe_id #deg取probe_id这列，并将其作为名字给x
  cg=c(names(head(sort(x),100)),#对x进行从小到大排列，取前100及后100，并取其对应的探针名，作为向量赋值给cg
       names(tail(sort(x),100)))
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F) #对dat按照cg取行，所得到的矩阵来画热图
  n=t(scale(t(dat[cg,])))#通过“scale”对log-ratio数值进行归一化，现在的dat是行名为探针，列名为样本名，由于scale这个函数应用在不同组数据间存在差异时，需要行名为样本，因此需要用t(dat[cg,])来转换，最后再转换回来
  
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n) #将ac的行名也就分组信息（是‘no TNBC’还是‘TNBC’）给到n的列名，即热图中位于上方的分组信息
  pheatmap(n,show_colnames =F,
           show_rownames = F,
          cluster_cols = F, 
           annotation_col=ac) #列名注释信息为ac即分组信息
  
  
}

# 从 top 的差异基因看
# 部分TNBC样本，实际上跟nonTNBC组内更相似。

# 一般来说，差异基因的mad或者sd值要大于其它基因。
# 同学们可以自行探索。




