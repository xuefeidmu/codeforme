## 查看TCGA的33个project
rm(list=ls())
library(TCGAbiolinks)
library(stringr)
library(SummarizedExperiment)
library(tibble)


## 列出TCGA里面的33中癌症
projs<-getGDCprojects()$project_id%>%
 str_subset("TCGA")
projs

## 选择想要的
cancer_type="TCGA-BRCA"
expquery=GDCquery(project=cancer_type,
 data.category="Transcriptome Profiling",
 data.type="Gene Expression Quantification",
 workflow.type="STAR - Counts")

GDCdownload(expquery,directory = "TCGA-BRCA") 

expquery=GDCprepare(expquery,directory="TCGA-BRCA") 

## 保存数据
save(expquery,file="BRCA.TCGA.rda") 
