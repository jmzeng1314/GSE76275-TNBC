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
# 每次都要检测数据
dat[1:4,1:4]  
dim(dat)
class(dat)
str(dat)
library(hgu133plus2.db)
ids=toTable(hgu133plus2SYMBOL)
head(ids)
dat=dat[ids$probe_id,]
dat[1:4,1:4] 

ids$median=apply(dat,1,median)
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]
dat=dat[ids$probe_id,]
rownames(dat)=ids$symbol
dat[1:4,1:4]  
dim(dat)
dat['GAPDH',]
boxplot(dat[,1])

if(T){
  ddata=t(dat)
  ddata[1:4,1:4]
  s=colnames(ddata);head(s)
  library(org.Hs.eg.db)
  s2g=toTable(org.Hs.egSYMBOL)
  g=s2g[match(s,s2g$symbol),1];head(g)
  #  probe Gene.symbol Gene.ID
  dannot=data.frame(probe=s,
                    "Gene.Symbol" =s, 
                    "EntrezGene.ID"=g)
  ddata=ddata[,!is.na(dannot$EntrezGene.ID)]#去除列的NA值
  dannot=dannot[!is.na(dannot$EntrezGene.ID),] #去除行的NA值
  head(dannot)
  library(genefu)
  # c("scmgene", "scmod1", "scmod2","pam50", "ssp2006", "ssp2003", "intClust", "AIMS","claudinLow")
  
  s<-molecular.subtyping(sbt.model = "pam50",data=ddata,
                         annot=dannot,do.mapping=TRUE)
  table(s$subtype)
  tmp=as.data.frame(s$subtype)
  subtypes=as.character(s$subtype)
}

library(genefu)
pam50genes=pam50$centroids.map[c(1,3)]
pam50genes[pam50genes$probe=='CDCA1',1]='NUF2'
pam50genes[pam50genes$probe=='KNTC2',1]='NDC80'
pam50genes[pam50genes$probe=='ORC6L',1]='ORC6'
x=dat
x=x[pam50genes$probe[pam50genes$probe %in% rownames(x)] ,]
table(group_list)
tmp=data.frame(group=group_list,
                  subtypes=subtypes)
rownames(tmp)=colnames(x)

library(pheatmap)

pheatmap(x,show_rownames = T,show_colnames = F,
         annotation_col = tmp)
pheatmap(x,show_rownames = T,show_colnames = F,
         annotation_col = tmp,
         filename = 'ht_by_pam50_raw.png') 
dev.off()

x=t(scale(t(x)))
x[x>1.6]=1.6
x[x< -1.6]= -1.6

pheatmap(x,show_rownames = T,show_colnames = F,
         annotation_col = tmp)
pheatmap(x,show_rownames = T,show_colnames = F,
         annotation_col = tmp,
         filename = 'ht_by_pam50_scale.png') 
dev.off()






