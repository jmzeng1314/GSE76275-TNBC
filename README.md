# 三阴性乳腺癌表达矩阵探索(公共数据库挖掘实战)

交流群在：https://mp.weixin.qq.com/s/N9YFEkh0TjZ4BzZvP5OT7g 

生信技能树联盟创始人jimmy手把手带你完成一个GEO数据库(`GSE76275`)挖掘实例，从必备R包安装，表达矩阵下载，PCA/boxplot/heatmp的数据检查，探针ID到基因ID的转换，根据生物学分组进行差异分析并且绘制火山图、热图，还有简单的超几何分布检验的KEGG等数据库注释结果。

表达矩阵进行`GSEA/GSVA`分析。

根据`TNBC(三阴性乳腺癌)`的生物学特征提取指定基因的表达量，使用`genefu`这个R包进行`PAM50`分类。

下载`TCGA数据库`的BRCA的芯片表达矩阵，同样进行PAM50分类，结合临床信息，并且对应比较GSE76275数据集。

- 代码在：https://github.com/jmzeng1314/TCGA_BRCA

下载`TCGA数据库`的BRCA的RNA-seq表达矩阵，同样进行PAM50分类，结合临床信息，并且对应比较GSE76275数据集。

- 代码在：https://github.com/jmzeng1314/TCGA_BRCA

下载`GTEx数据库`的RNA-seq表达矩阵，并且提取其中**属于breast的样本**，同样进行PAM50分类，结合临床信息，并且对应比较GSE76275数据集。

- 代码在： https://github.com/jmzeng1314/gtex_BRCA

下载`METEBRIC数据库`的RNA-seq表达矩阵，同样进行PAM50分类，结合临床信息，并且对应比较GSE76275数据集。

- 代码在：https://github.com/jmzeng1314/METABRIC

### 首先需要安装必备R包

需要自行下载学习R语言及熟练使用Rstudio编辑器，根据课程配套代码安装R包。

### 第一步：下载表达矩阵

首先需要了解GEO数据库，建议通读： [解读GEO数据存放规律及下载，一文就够](http://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247486063&idx=1&sn=156bee5397e979722b36b78284188538&chksm=9b484ad4ac3fc3c2d025b9e4bb1c3c8392839c08d84697754d7d95d041b539479a45f19cf5d5&scene=21#wechat_redirect)

然后根据课件代码使用GEOquery包进行下载，详细教程见：[从GEO数据库下载得到表达矩阵 一文就够](http://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247486087&idx=1&sn=1e775a1c3e215384e381953a9fa74ec3&chksm=9b484a3cac3fc32ac6492cc3faf1e152277f6861c973510b384fe9bd352886006a87ba6f44f1&scene=21#wechat_redirect)

示例数据集：GEO数据库(`GSE76275`) 对应的3篇SCI文章见课程附件。

### 第三步：检查表达矩阵

首先根据全局基因表达绘制PCA图如下：

![](http://www.bio-info-trainee.com/wp-content/uploads/2018/12/all_samples_PCA.png)

代码见课程附件，具体参数视频里有介绍，结果很容易理解，可以看到那些不属于TNBC的样本跟TNBC组有着很清晰的界限，而且可以看到TNBC组本身也有着界限分明的两列。

还有top 1000的sd的基因提取出来绘制热图：

![](http://www.bio-info-trainee.com/wp-content/uploads/2018/12/heatmap_top1000_sd.png)

这个结果跟PCA分析结果相呼应！

### 第4步：差异分析

这里走标准的limma流程，详细推文介绍：[根据分组信息做差异分析- 这个一文不够的](http://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247486112&idx=1&sn=67a2104c62222bcb139623699f874a6c&scene=21#wechat_redirect)

### 第5步：GO/KEGG数据库注释

得到的差异分析结果，也可以走标准的火山图，热图，GO/KEGG数据库注释，见推文：[差异分析得到的结果注释一文就够](https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247486120&idx=1&sn=14d7892c1beec2fb9cdfc0ec0aba3e4e&scene=21#wechat_redirect)

### 第6步：GSEA/GSVA分析

代码在：https://github.com/jmzeng1314/GSE76275-TNBC 

### TCGA数据库挖掘

代码在：https://github.com/jmzeng1314/TCGA_BRCA

### GTEx数据库挖掘

代码在：https://github.com/jmzeng1314/gtex_BRCA

### METABRIC数据库挖掘

代码在：https://github.com/jmzeng1314/METABRIC

### 增加WGCNA流程代码

对表达矩阵挑选top5000的MAD基因，以及top10000后，分别独立走WGCNA流程看结果。

背景知识参考：https://github.com/jmzeng1314/my_WGCNA
