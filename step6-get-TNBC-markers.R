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
load(file = 'step1-output.Rdata')
dat[1:4,1:4]
library(hgu133plus2.db)# BiocManager::install("hgu133plus2.db")
p2s=toTable(hgu133plus2SYMBOL)
# ESR1	Estrogen Receptor 1	
# ESR2	Estrogen Receptor 2
# PGR	Progesterone Receptor
# ERBB2	Erb-B2 Receptor Tyrosine Kinase 2
k=p2s$symbol %in% c('ERBB2','ESR1','ESR2','PGR')
# %in% 取交集
np=p2s[k,1]
ng=p2s[k,2]
x=dat[np,]
rownames(x)=paste(ng,np,sep = ':')
#paste 粘贴字符，seq指定分割符
library(pheatmap)
tmp=data.frame(group=group_list)
rownames(tmp)=colnames(x)
pheatmap(x,annotation_col = tmp)

# Histograms show the distribution and frequency of tumors using 
# relative ER, PR, and HER2 GE levels (log2) and bimodal fit to identify TN tumor samples. 

# a 2-component Gaussian mixture distribution model and 
# parameters were estimated by maximum likelihood optimization, using optim function 

hist(x["ESR1:205225_at",],breaks = 1:15)
#hist画直方图
wdata=data.frame(v=as.numeric(x["ESR1:205225_at",]))
library(ggpubr)
gghistogram(wdata, x = 'v',  y = "..density..",
            add_density = T,
          add = "mean", rug = TRUE)
          #color = "sex", palette = c("#00AFBB", "#E7B800")



