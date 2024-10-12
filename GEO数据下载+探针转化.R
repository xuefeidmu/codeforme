##手工下载soft.gz文件
##
rm(list = ls())
options(stringsAsFactors = F)

#加载R包
library(GEOquery)
#手工读取geo文件 filename这里下载的是soft文件
gse <- getGEO(filename = "GSE76250_family.soft.gz",destdir = ".")
str(gse)
## 查看用几个平台
length(gse)
## 提取探针 基因及其位于染色体上的位置等信息
## 这里是手动提取
id_probe <- gse@gpls$GPL17586@dataTable@table
dim(id_probe)
head(id_probe)
View(head(id_probe))## you need to check this , which column do you need
#提取probeid和geneassignment（包括基因symbol）
probe2gene <- id_probe[,c(2,8)]

# 提取第8列 gene_assignment中的基因名称,并添加到probe2gene
library(stringr) 
#把gene_assignment用//分开，提取第2部分
symbol <- str_split(probe2gene$gene_assignment,'//',simplify = T)[,2]
#去掉空格，加到probe2gene
probe2gene$symbol=trimws(symbol)
plot(table(table(probe2gene$symbol)),xlim=c(1,50))
head(probe2gene)

dim(probe2gene)
View(head(probe2gene))
ids2 <- probe2gene[,c(1,3)]
View(head(ids2))
ids2[1:20,1:2]#含有缺失值
table(table(unique(ids2$symbol)))#30907 ,30906个基因，一个空字符
save(ids2,probe2gene,file='gse-probe2gene.Rdata')

## 表达矩阵和基因id的合并
## 下载gse数据

## 下面二选一
library(GEOmirror)
gset <- geoChina("GSE76250")
## or
library(GEOquery)
gset <- getGEO( "GSE76250",destdir = ".")
gset
a=exprs(gset[[1]])
a[1:4,1:4]
gset[[1]]@annotation


#过滤表达矩阵
exprSet <- a
library(dplyr)
##把在ids2中有的探针保留
exprSet <- exprSet[rownames(exprSet) %in% ids2$probeset_id,]
dim(exprSet)
exprSet[1:5,1:5]

#ids过滤探针
#rownames(exprSet) 中每个元素在 ids2$probeset_id 中的位置
ids <- ids2[match(rownames(exprSet),ids2$probeset_id),]
dim(ids)
ids[1:5,1:2]
ids2[1:5,1:2]

#合并表达矩阵和ids
#定义合并用的函数idcombine
idcombine <- function(exprSet, ids){
  #以ids$symbol为标准，对exprSet分组，同一个symbol保留表达量最高的一行数据
  tmp <- by(exprSet,
            ids$symbol,
            function(x) rownames(x)[which.max(rowMeans(x))])
  #每个符号基因的最佳代表探针集ID
  probes <- as.character(tmp)
  print(dim(exprSet))
  #把矩阵按probes顺序排
  exprSet <- exprSet[rownames(exprSet) %in% probes,]
  print(dim(exprSet))
  #根据ids中的对应关系，把probe id变成symbol
  rownames(exprSet) <- ids[match(rownames(exprSet), ids$probeset_id),2]
  return(exprSet)#返回处理后的exprSet
}
new_exprSet <- idcombine(exprSet,ids)
new_exprSet[1:4,1:6]
dim(new_exprSet)

rownames(new_exprSet)
save(new_exprSet,file = "new_exprSet.Rdata")








