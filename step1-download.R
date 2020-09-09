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
options(stringsAsFactors = F)#在调用as.data.frame的时，将stringsAsFactors设置为FALSE可以避免character类型自动转化为factor类型
# 注意查看下载文件的大小，检查数据 
f='GSE76275_eSet.Rdata'

# trying URL 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76275/matrix/GSE76275_series_matrix.txt.gz'
# Content type 'application/x-gzip' length 78267951 bytes (74.6 MB)
# 如果你在中国大陆，可以考虑使用 GEOChina 镜像下载

library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
if(!file.exists(f)){
  gset <- getGEO('GSE76275', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE76275_eSet.Rdata')  ## 载入数据
class(gset)  #查看数据类型
length(gset)  #
class(gset[[1]])
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]] #
dat=exprs(a) #a现在是一个对象，取a这个对象通过看说明书知道要用exprs这个函数
dim(dat)#看一下dat这个矩阵的维度
# [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
dat[1:4,1:4] #查看dat这个矩阵的1至4行和1至4列，逗号前为行，逗号后为列
pd=pData(a) #通过查看说明书知道取对象a里的临床信息用pData
## 挑选一些感兴趣的临床表型。
trait=pd[,51:53]
head(trait) #从前到后显示矩阵信息（默认前六行）
trait$T=substring(trait[,2],2,2) #提取第二列的第二个字符
#substring 提取字符串的一部分，字符串处理函数
#substring(x,first,last) x是字符向量输入，first是第一个字符要被提取的位置，last是最后一个字符要被提取的位置
trait$N=substring(trait[,2],4,4)
trait$M=substring(trait[,2],6,6)
colnames(trait)=c('age','tmn','bmi','T','M','N') #给列重新命名
head(trait)
# 把挑选到的 感兴趣的临床表型存储成为R文件
save(trait,file='trait.Rdata')

group_list = ifelse(pd$characteristics_ch1.1=='triple-negative status: not TN',
   'noTNBC','TNBC')
# $符号:是按列取，即对pd取characteristics_ch1.1这一列，判断是否为'triple-negative status:not TN'，如果是，就改为'noTNBC',否则改为'TNBC'
table(group_list) #统计频率
save(dat,group_list,file = 'step1-output.Rdata')


## 下面的代码是更为详细的注释，不需要运行的。
######### 导入数据后查看数据情况，行、列名，数据维度、目标数据的大致情况，比如我们想看不同分组的各个基因的表达量，那么你就要查看分组在哪一行/列，怎么分组的，基因是什么情况，只要在掌握已有数据架构的情况下，才能随心所欲的调整表达矩阵（以下代码可以自行增减，达到查看数据情况的目的即可）

# b = gset[[1]]  ## 降级提取b
# exprSet=exprs(b)  ## 获取表达矩阵（如下图2），发现已是log处理后的数据。通常表达矩阵的原始数字从0但好几百万都不等，需要进行归一化处理。14488875 elements
# pdata = pData(b)  ## 使用函数?pData获取样本临床信息（如性别、年龄、肿瘤分期等等） 265 obs. of 69 
# colnames(pdata) ## 查看列名，即查看所包含的临床信息类型，找到triple negative status在第67列
# length(colnames(pdata)) ## 查看pdata列数 69
# pdata[,67]   ## 查看第67列triple negative status的数据情况,发现按照"TN"，"not TN"排列
# group_list = as.character(pdata[, 67])  ## 将67列改成字符
# dim(exprSet) ## 查看矩阵的维度 54675行 265列
# table(group_list) ## 查看67列信息
# 







