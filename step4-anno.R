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

# 这里面的代码非常复杂，涉及到3个内容，而且几乎是独立的内容。

rm(list = ls())  ## 魔幻操作，一键清空~
load(file = 'deg.Rdata')
head(deg)
## 不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就很有可能不一致哦。
logFC_t=1 
deg$g=ifelse(deg$P.Value>0.01,'stable',
            ifelse( deg$logFC > logFC_t,'UP',
                    ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
##if 判断：如果这一基因的P.Value>0.01，则为stable基因,否则再进行if判断（对于那些P.Value<0.01的基因），
#如果这一基因logFC > logFC_t(logFC_t=1),则为up基因，否则再进行if判断，如果logFC < -logFC_t，则为DOWN基因，否则为stable基因
table(deg$g)

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
#提取包中的标准基因名（SYMBOL）与ENTREZID的对应关系
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG=deg
head(DEG)

DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol') 
#通过DEG矩阵的symbol列和df矩阵的SYMBOL列的对应来合并两个矩阵
#merge函数，通过相同的列的对应关系来合并两个矩阵
head(DEG)
save(DEG,file = 'anno_DEG.Rdata')


gene_up= DEG[DEG$g == 'UP','ENTREZID'] #取出UP基因及其对应的ENTREZID
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] #取出DOWN基因及其对应的ENTREZID
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] ) #转换为字符向量
data(geneList, package="DOSE") 
#DOSE做DO语义相似性的包
#这个包实现了Resnik、Schlicker、Jiang、Lin和Wang分别提出的测量DO术语和基因产物之间语义相似性的五种方法。
#通过超几何模型和基因集富集分析等富集分析方法，发现高通量生物数据的疾病相关性。
head(geneList)
boxplot(geneList)
boxplot(DEG$logFC)

geneList=DEG$logFC
names(geneList)=DEG$ENTREZID #把DEG$ENTREZID作为geneList的名字
geneList=sort(geneList,decreasing = T) #排序，从大到小（降序）

boxplot(geneList)

## KEGG pathway analysis
### 做KEGG数据集超几何分布检验分析，重点在结果的可视化及生物学意义的理解。
if(T){
  ###   over-representation test
  # 下面把3个基因集分开做超几何分布检验
  # 首先是上调基因集。
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  kk=kk.up
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,'kk.up.csv')
  
  # 首先是下调基因集。
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  kk=kk.down
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,'kk.down.csv')
  
  # 最后是上下调合并后的基因集。
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  kk=kk.diff
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,'kk.diff.csv')
  
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1
  source('functions.R') #画图设置, 这个图很丑，大家可以自行修改。
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  
  ggsave(g_kegg,filename = 'kegg_up_down.png')
  
  ###  GSEA 
  ## GSEA算法跟上面的使用差异基因集做超几何分布检验不一样。
  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 20,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)
  head(kk_gse)[,1:6]
  gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
  gseaplot(kk_gse, 'hsa04110',title = 'Cell cycle') 
  kk=DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  tmp=kk@result
  write.csv(kk@result,'kk.gsea.csv')
  
  
  # 这里找不到显著下调的通路，可以选择调整阈值，或者其它。
  # down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
  # up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
  # 
  # g_kegg=kegg_plot(up_kegg,down_kegg)
  # print(g_kegg)
  # ggsave(g_kegg,filename = 'kegg_up_down_gsea.png')
  # 
  
}

### GO database analysis 
### 做GO数据集超几何分布检验分析，重点在结果的可视化及生物学意义的理解。
{
  
  g_list=list(gene_up=gene_up,
              gene_down=gene_down,
              gene_diff=gene_diff)
  # 因为go数据库非常多基因集，所以运行速度会很慢。
  if(F){
    go_enrich_results <- lapply( g_list , function(gene) {
      lapply( c('BP','MF','CC') , function(ont) {
        cat(paste('Now process ',ont ))
        ego <- enrichGO(gene          = gene,
                        universe      = gene_all,
                        OrgDb         = org.Hs.eg.db,
                        ont           = ont ,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.99,
                        qvalueCutoff  = 0.99,
                        readable      = TRUE)
        
        print( head(ego) )
        return(ego)
      })
    })
    save(go_enrich_results,file = 'go_enrich_results.Rdata')
    
  }
  # 只有第一次需要运行，就保存结果哈，下次需要探索结果，就载入即可。
  load(file = 'go_enrich_results.Rdata')
  
  n1= c('gene_up','gene_down','gene_diff')
  n2= c('BP','MF','CC') 
  for (i in 1:3){
    for (j in 1:3){
      fn=paste0('dotplot_',n1[i],'_',n2[j],'.png')
      cat(paste0(fn,'\n'))
      png(fn,res=150,width = 1080)
      print( dotplot(go_enrich_results[[i]][[j]] ))
      dev.off()
    }
  }
  
  
}


### 对 MigDB中的全部基因集 做GSEA分析。
# http://www.bio-info-trainee.com/2105.html
# http://www.bio-info-trainee.com/2102.html 
{
  load(file = 'anno_DEG.Rdata')
  geneList=DEG$logFC
  names(geneList)=DEG$symbol
  geneList=sort(geneList,decreasing = T)
  #选择gmt文件（MigDB中的全部基因集）
  d='./MsigDB/symbols'
  gmts <- list.files(d,pattern = 'all')
  gmts
  #GSEA分析
  library(GSEABase) # BiocManager::install('GSEABase')
  ## 下面使用lapply循环读取每个gmt文件，并且进行GSEA分析
  gsea_results <- lapply(gmts, function(gmtfile){
    # gmtfile=gmts[2]
    # 如果不确定循环里面的代码，就尝试赋值，测试下面的代码。
    geneset <- read.gmt(file.path(d,gmtfile)) 
    print(paste0('Now process the ',gmtfile))
    egmt <- GSEA(geneList, TERM2GENE=geneset, verbose=FALSE)
    head(egmt)
    # gseaplot(egmt, geneSetID = rownames(egmt[1,]))
    
    return(egmt)
  })
  # 上面的代码耗时，所以保存结果到本地文件
  save(gsea_results,file = 'gsea_results.Rdata')
  #提取gsea结果，熟悉这个对象
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_list)
  gseaplot(gsea_results[[2]],'KEGG_CELL_CYCLE',title = 'KEGG_CELL_CYCLE') 
  gseaplot(gsea_results[[2]],'FARMER_BREAST_CANCER_CLUSTER_6') 
  gseaplot(gsea_results[[5]],'GO_CONDENSED_CHROMOSOME_OUTER_KINETOCHORE') 
  
}

### 对 MigDB中的全部基因集 做GSVA分析。
## 还有ssGSEA, PGSEA
{
  load(file = 'step1-output.Rdata')
  # 每次都要检测数据
  dat[1:4,1:4]  
  library(hgu133plus2.db)
  ids=toTable(hgu133plus2SYMBOL)
  #通过看hgu133plus2.db这个包的说明书知道提取probe_id（探针名）和symbol（基因名）的对应关系的表达矩阵的函数为toTable
  head(ids)
  dat=dat[ids$probe_id,]
  dat[1:4,1:4] 
  ids$median=apply(dat,1,median)
  #对dat这个矩阵按行操作，取每一行的中位数，将结果添加到ids矩阵median列
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]
  #对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
  ids=ids[!duplicated(ids$symbol),]
  #将symbol这一列去除重复项
  dat=dat[ids$probe_id,]
  rownames(dat)=ids$symbol
  dat[1:4,1:4]  
  
  X=dat
  table(group_list)
  ## Molecular Signatures Database (MSigDb) 
  d='./MSigDB/symbols/'
  gmts=list.files(d,pattern = 'all')
  gmts
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GSVA) # BiocManager::install('GSVA')
  
  if(T){
    es_max <- lapply(gmts, function(gmtfile){ 
      #gmtfile=gmts[8];gmtfile
      geneset <- getGmt(file.path(d,gmtfile))  
      es.max <- gsva(X, geneset, 
                     mx.diff=FALSE, verbose=FALSE, 
                     parallel.sz=1)
      return(es.max)
    })
    adjPvalueCutoff <- 0.001
    logFCcutoff <- log2(2)
    es_deg <- lapply(es_max, function(es.max){
      table(group_list)
      dim(es.max)
      design <- model.matrix(~0+factor(group_list))
      colnames(design)=levels(factor(group_list))
      rownames(design)=colnames(es.max)
      design
      library(limma)
      contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),
                                     levels = design)
      contrast.matrix<-makeContrasts("TNBC-noTNBC",
                                     levels = design)
      
      contrast.matrix ##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
      
      deg = function(es.max,design,contrast.matrix){
        ##step1
        fit <- lmFit(es.max,design)
        ##step2
        fit2 <- contrasts.fit(fit, contrast.matrix) 
        ##这一步很重要，大家可以自行看看效果
        
        fit2 <- eBayes(fit2)  ## default no trend !!!
        ##eBayes() with trend=TRUE
        ##step3
        res <- decideTests(fit2, p.value=adjPvalueCutoff)
        summary(res)
        tempOutput = topTable(fit2, coef=1, n=Inf)
        nrDEG = na.omit(tempOutput) 
        #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
        head(nrDEG)
        return(nrDEG)
      }
      
      re = deg(es.max,design,contrast.matrix)
      nrDEG=re
      head(nrDEG) 
      return(nrDEG)
    })
  } 
  
  gmts
  
  save(es_max,es_deg,file='gsva_msigdb.Rdata')
  
  load(file='gsva_msigdb.Rdata')
  
  library(pheatmap)
  lapply(1:length(es_deg), function(i){
    # i=2
    print(i)
    dat=es_max[[i]]
    df=es_deg[[i]]
    df=df[df$P.Value<0.01 & abs(df$logFC) > 0.3,]
    print(dim(df))
    if(nrow(df)>5){
      n=rownames(df)
      dat=dat[match(n,rownames(dat)),]
      ac=data.frame(g=group_list)
      rownames(ac)=colnames(dat)
      rownames(dat)=substring(rownames(dat),1,50)
      pheatmap::pheatmap(dat, 
                         fontsize_row = 8,height = 11,
                         annotation_col = ac,show_colnames = F,
                         filename = paste0('gsva_',strsplit(gmts[i],'[.]')[[1]][1],'.pdf'))
      
    }
})
  
  adjPvalueCutoff <- 0.001
  logFCcutoff <- log2(2)
  df=do.call(rbind ,es_deg)
  es_matrix=do.call(rbind ,es_max)
  df=df[df$P.Value<0.01 & abs(df$logFC) > 0.5,]
  write.csv(df,file = 'GSVA_DEG.csv')
}







