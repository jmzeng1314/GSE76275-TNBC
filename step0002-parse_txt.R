rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F) #避免character类型自动转化为factor类型

# 下面的代码是演示，GEOquery包的getGEO函数的细节
# 因为一个 readr 升级问题，曾经GEOquery包的getGEO函数瘫痪

fname='GSE76275_series_matrix.txt.gz' 
AnnotGPL=FALSE #不获得注释文件
destdir=tempdir() #destdir设置下载目录。tempdir临时文件存放路径
getGPL=TRUE #下载平台文件
parseCharacteristics=TRUE #解析GSE矩阵文件的特征信息

library(readr )                       
dat <- read_lines(fname)
## get the number of !Series and !Sample lines
series_header_row_count <- sum(grepl("^!Series_", dat))
#grepl返回所有的查询结果，并用逻辑向量表示有没有找到匹配
sample_header_start <- grep("^!Sample_", dat)[1]
#grep返回匹配项的下标
samples_header_row_count <- sum(grepl("^!Sample_", dat))
series_table_begin_line = grep("^!series_matrix_table_begin", dat) 
##  colClasses <- c('character',rep('numeric',nrow(sampledat)))
datamat <- read_tsv(fname,quote='"',
                    na=c('NA','null','NULL','Null'), skip = series_table_begin_line,
                    comment = '!series_matrix_table_end')
datamat[1:4,1:4]
