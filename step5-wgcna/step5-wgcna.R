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

# 这个代码运行结果有点诡异，暂时放在这里。

rm(list = ls())  ## 魔幻操作，一键清空~
# 首先构建 WGCNA 的input 信息。
library(WGCNA)
## step1:input
if(F){
  load(file = '../step1-output.Rdata')
  # 每次都要检测数据
  dat[1:4,1:4]  
  library(hgu133plus2.db)
  ids=toTable(hgu133plus2SYMBOL)#toTable提取包中的probe_id（探针名）和symbol（基因名）的对应关系的表达矩阵
  head(ids)
  dat=dat[ids$probe_id,] #过滤，筛出不存在包中的探针
  dat[1:4,1:4] 
  ids$median=apply(dat,1,median)
  #对dat这个矩阵按行操作，取每一行的中位数，将结果添加到ids矩阵median列
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]
  #对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
  ids=ids[!duplicated(ids$symbol),]
  #去除ids矩阵symbol重复项
  dat=dat[ids$probe_id,]
  rownames(dat)=ids$symbol
  #把ids矩阵的symbol名作为dat矩阵的行名
  dat[1:4,1:4]  
  dim(dat)
  load(file='../trait.Rdata') 
 # dat=dat[,group_list=='TNBC']
 # trait=trait[group_list=='TNBC',]
  head(trait)
  
  WGCNA_matrix = t(dat[order(apply(dat,1,mad), decreasing = T)[1:10000],])
  #mad绝对中位差
  #对dat矩阵每行求绝对中位差，order排序，降序排序，取前5000
  # 这个地方有问题。
  datExpr0 <- WGCNA_matrix  ## top 10000 mad genes
  datExpr <- datExpr0 
  
  datTraits=trait
  datTraits$gsm=rownames(datTraits)#添加gsm的列
  
  ## 下面主要是为了防止临床表型与样本名字对不上
  sampleNames = rownames(datExpr);
  traitRows = match(sampleNames, datTraits$gsm) 
  #match函数，返回sampleNames向量中每个元素在datTraits$gsm向量中的位置
  rownames(datTraits) = datTraits[traitRows, 'gsm'] 
  # pheatmap::pheatmap(datExpr)
  save(datExpr,datTraits,group_list,file = 'input_for_WGCNA.Rdata')
}

## 值得注意的是，这个时候，我们简化了运算量
# 仅仅是选取了top5000的MAD基因

## step 2 
rm(list = ls())
library(WGCNA)
load(file = 'input_for_WGCNA.Rdata')
## 查看是否有离群样品，这个时候可以选择性去除它们
sampleTree = hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
datExpr[1:4,1:4]
if(T){
  powers = c(c(1:10), seq(from = 12, to=50, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  sft
  sft$powerEstimate
  #设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
  png("step2-beta-value.png",width = 800,height = 600)
  # Plot the results:
  ##sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
}

sft$powerEstimate 
## step3 构建加权共表达网络（Weight co-expression network)
## 首先是一步法完成网络构建
cor <- WGCNA::cor
if(T){
  net = blockwiseModules(
    datExpr,
    power = sft$powerEstimate,
    maxBlockSize = 5000,
    TOMType = "unsigned", minModuleSize = 30,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = F, 
    verbose = 3
  )
  table(net$colors) 
}
cor<-stats::cor
## 然后是分布法完成网络构建，仅供有探索精神的同学挑战。

## 构建加权共表达网络分为两步：
## 1. 计算邻近值，也是就是两个基因在不样品中表达量的表达相关系数(pearson correlation rho)，
## 参考 2.b.2 in https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf
## 2. 计算topology overlap similarity (TOM)。 WGCNA认为，只通过计算两个基因的表达相关系数构建共表达网络是不足够的。
## 于是他们用TOM表示两个基因在网络结构上的相似性，即两个基因如果具有相似的邻近基因，这两个基因更倾向于有相互作用。
## 参考 2.b.3 in https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf

# 并不需要做
if(F){
  #(1)网络构建 Co-expression similarity and adjacency 
  adjacency = adjacency(datExpr, power = sft$powerEstimate) 
  #(2) 邻近矩阵到拓扑矩阵的转换，Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency);
  dissTOM = 1-TOM
  # (3) 聚类拓扑矩阵 Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average");
  # Plot the resulting clustering tree (dendrogram)
  sizeGrWindow(12,9)
  ## 这个时候的geneTree与一步法的 net$dendrograms[[1]] 性质类似，但是还需要进行进一步处理
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04);
  #(4) 聚类分支的修整 dynamicTreeCut 
  # We like large modules, so we set the minimum module size relatively high:
  minModuleSize = 30;
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);
  table(dynamicMods)
  #4. 绘画结果展示
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  #sizeGrWindow(8,6)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  #5. 聚类结果相似模块的融合，Merging of modules whose expression profiles are very similar
  #在聚类树中每一leaf是一个短线，代表一个基因，
  #不同分之间靠的越近表示有高的共表达基因，将共表达极其相似的modules进行融合
  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  # Plot the result
  #sizeGrWindow(7, 6)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  #选择有75%相关性的进行融合
  MEDissThres = 0.25
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs
  
}

## step 4 ： 模块可视化
if(T){
  
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  table(mergedColors)
  moduleColors=mergedColors
  # Plot the dendrogram and the module colors underneath
  png("step4-genes-modules.png",width = 800,height = 600)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  ## assign all of the gene to their corresponding module 
  ## hclust for the genes.
}

if(F){
  #明确样本数和基因
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  #首先针对样本做个系统聚类
  datExpr_tree<-hclust(dist(datExpr), method = "average")
  par(mar = c(0,5,2,0))
  plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
       cex.axis = 1, cex.main = 1,cex.lab=1)
  ## 如果这个时候样本是有性状，或者临床表型的，可以加进去看看是否聚类合理
  #针对前面构造的样品矩阵添加对应颜色
  head(datTraits)
  sample_colors1 <- numbers2colors(as.numeric(factor(datTraits$group_list)), 
                                   colors = rainbow(length(levels(factor(datTraits$group_list)))),signed = FALSE)
  
  ssss=as.matrix(  data.frame(group=sample_colors1 ))
  par(mar = c(1,4,3,1),cex=0.8)
  
  png("sample-subtype-cluster.png",width = 800,height = 600)
  plotDendroAndColors(datExpr_tree, ssss,
                      groupLabels = colnames(sample),
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
}

## step 5 (最重要的) 模块和性状的关系
## 这一步主要是针对于连续变量，如果是分类变量，需要转换成连续变量方可使用
table(group_list)

if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  design1=model.matrix(~0+ as.factor(group_list))
  design=design1
  colnames(design) 
  colnames(design)=c(levels(as.factor(group_list)) ) 
  moduleColors <- labels2colors(net$colors)
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩 (样本vs模块)
  moduleTraitCor = cor(MEs, design , use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  png("step5-Module-trait-relationships.png",width = 800,height = 1200,res = 120)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  table( labels2colors(net$colors))
  # 除了上面的热图展现形状与基因模块的相关性外
  # 还可以是条形图,但是只能是指定某个形状
  # 或者自己循环一下批量出图。
  colnames(design)
  adolescent = as.data.frame(design[,1]);
  names(adolescent) = "adolescent"
  y=adolescent
  GS1=as.numeric(cor(y,datExpr, use="p"))
  GeneSignificance=abs(GS1)
  # Next module significance is defined as average gene significance.
  ModuleSignificance=tapply(GeneSignificance,
                            moduleColors, mean, na.rm=T)
  sizeGrWindow(8,7)
  par(mfrow = c(1,1))
  # 如果模块太多，下面的展示就不友好
  # 不过，我们可以自定义出图。
  plotModuleSignificance(GeneSignificance,moduleColors,las=2)
  
}


## step 8 
# 主要是关心具体某个模块内部的基因
table(moduleColors)
group_g=data.frame(gene=colnames(datExpr),
                   group=moduleColors)
save(group_g,file='wgcna_group_g.Rdata')
save(net,file='wgcna_results_top10000.Rdata')

rm(list = ls())
load(file='wgcna_group_g.Rdata')
library(clusterProfiler)
# Convert gene ID into entrez genes
head(group_g)
tmp <- bitr(group_g$gene, fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")

de_gene_clusters=merge(tmp,group_g,by.x='SYMBOL',by.y='gene')
table(de_gene_clusters$group)
head(de_gene_clusters)

list_de_gene_clusters <- split(de_gene_clusters$ENTREZID, 
                               de_gene_clusters$group)

library(ggplot2)
gcSample= list_de_gene_clusters
source('function.R')
# 下一步非常耗时，保守估计半个小时
# 主要是对我们的模块进行功能注释，就是GO/KEGG数据库
com_kegg_go(gcSample,'top10000')




