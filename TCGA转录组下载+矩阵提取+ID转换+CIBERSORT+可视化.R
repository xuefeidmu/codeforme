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


## 加载数据
load("BRCA.TCGA.rda")

## 提取表达矩阵count
exp=assay(expquery) %>% 
  as.data.frame()

## 提取表达矩阵tpm
tpm_exp=assay(expquery,4) %>% 
  as.data.frame()


## COUNT表达矩阵行名ID转换
library(tinyarray)
exp=trans_exp_new(exp) %>% 
  as.data.frame()
exp[1:4,1:4]


## TPM表达矩阵行名ID转换
library(tinyarray)
tpm_exp=trans_exp(tpm_exp,mrna_only=T)
tpm_exp[1:4,1:4]


## 行名从样本长名变成样本短名
colnames(tpm_exp) <- substr(colnames(tpm_exp), 1, 16)
colnames(exp) <- substr(colnames(exp), 1, 16)



## 读取临床信息
clinical <- colData(expquery) %>% 
 as.data.frame()


## 选择三阴乳腺癌"病人"
table(clinical$paper_BRCA_Subtype_PAM50)
tnbc_patient <- clinical$patient[which(clinical$paper_BRCA_Subtype_PAM50=="Basal")]


## 提取三阴病人的tpm表达谱
tpm_exp <- tpm_exp[, substr(colnames(tpm_exp), 1, 12) %in% tnbc_patient]




## 把TPM整理成CIBERSORT需要的格式
exp2=as.data.frame(tpm_exp)
exp2=rownames_to_column(exp2)
write.table(exp2,file="exp.txt",row.names=F,quote=F,sep="\t")

## 加载cibersort，分析
source("CIBERSORT.R")
TME.results=CIBERSORT("LM22.txt",
 "exp.txt",
 perm=1000,
 QN=T)
TME.results[1:4,1:4]
re<-TME.results[,-(23:25)] 

## 保存
write.csv(re,"TME.results.txt")

## 每列（免疫细胞种类）中一半数值>0的，选中
library(pheatmap)
k<-apply(re,
         2, ## 按列进行函数处理
         function(x){sum(x==0)<nrow(TME.results)/2}) ##0数值<1/2行的列
table(k)
re2<-as.data.frame(t(re[,k]))

## 热图

pheatmap(re2,scale="row",
 show_colnames=F,
 color=colorRampPalette(c("navy","white","firebrick3"))(50))


## 直方图展示

# 加载必要的库
library(RColorBrewer) # 用于色彩调色板
library(tidyr) # 用于数据操作（gather 函数）
library(ggplot2) # 用于创建可视化

# 使用 "Set1" 调色板的 8 种颜色定义自定义调色板
mypalette <- colorRampPalette(brewer.pal(8, "Set1"))

# 将 're' 转换为数据框并对其进行重塑，以便于 ggplot 使用
dat <- re %>%
 as.data.frame() %>% # 将 're' 转换为数据框
 rownames_to_column("Sample") %>% # 将行名创建为名为 'Sample' 的列
 gather(key = Cell_type, value = Proportion, -Sample) # 将数据重塑为长格式

# 使用 ggplot2 创建柱状图
ggplot(dat, aes(Sample, Proportion, fill = Cell_type)) + # 设置美学映射：x 轴 = Sample，y 轴 = Proportion，填充颜色 = Cell_type
 geom_bar(stat = "identity") + # 使用堆叠柱状图展示比例
 labs(fill = "Cell Type", x = "", y = "Estimated Proportion") + # 添加坐标轴和图例的标签
 theme_bw() + # 使用白色背景主题
 theme(
 axis.text.x = element_blank(), # 移除 x 轴的文本
 axis.ticks.x = element_blank(), # 移除 x 轴的刻度
 legend.position = "bottom" # 将图例放置在底部
 ) + 
 scale_y_continuous(expand = c(0.01, 0)) + # 调整 y 轴，使其扩展最小化
 scale_fill_manual(values = mypalette(22)) # 使用自定义调色板中的 22 种颜色


## 箱线图
ggplot(dat,aes(Cell_type,Proportion,fill=Cell_type))+
 geom_boxplot(outlier.shape=21,color="black")+
 theme_bw()+
 labs(x="CellType",y="EstimatedProportion")+
 theme(axis.text.x=element_blank(),
 axis.ticks.x=element_blank(),
 legend.position="bottom")+
 scale_fill_manual(values=mypalette(22))



## 肿瘤正常样本对比
dat$Group=ifelse(as.numeric(str_sub(dat$Sample,14,15))<10,"tumor","normal")
table(dat$Group)
library(ggpubr)
ggplot(dat,aes(Cell_type,Proportion,fill=Group))+
 geom_boxplot(outlier.shape=21,color="black")+
 theme_bw()+
 labs(x="CellType",y="EstimatedProportion")+
 theme(legend.position="top")+
 theme(axis.text.x=element_text(angle=80,vjust=0.5))+
 scale_fill_manual(values=mypalette(22)[c(6,1)])+stat_compare_means(aes(group=Group,label=..p.signif..),method="kruskal.test")









# 计算各细胞类型的中位数并按降序排列，以调整箱线图顺序
library(dplyr)

a = dat %>%
 group_by(Cell_type) %>% # 按 Cell_type 分组
 summarise(m = median(Proportion)) %>% # 计算每组的中位数
 arrange(desc(m)) %>% # 按中位数降序排列
 pull(Cell_type) # 提取排列后的 Cell_type 列表

# 调整 Cell_type 的因子顺序，使其在箱线图中按照中位数降序排列
dat$Cell_type = factor(dat$Cell_type, levels = a)

# 绘制调整顺序后的箱线图
ggplot(dat, aes(Cell_type, Proportion, fill = Cell_type)) +
 geom_boxplot(outlier.shape = 21, color = "black") + # 绘制箱线图，设置异常值形状和边框颜色
 theme_bw() + # 使用白色背景主题
 labs(x = "Cell Type", y = "Estimated Proportion") + # 添加坐标轴标签
 theme(axis.text.x = element_blank(), # 移除 x 轴文本
 axis.ticks.x = element_blank(), # 移除 x 轴刻度
 legend.position = "bottom") + # 将图例放在底部
 scale_fill_manual(values = mypalette(22)) # 使用自定义调色板中的 22 种颜色



## 首先在re中选中三阴肿瘤样本（排除正常样本）
tumor_re <- re[grepl("^.{13}01", rownames(re)), ]



##  根据免疫评分tumor_re选中M0高浸润和低浸润

# 加载 dplyr 包
library(dplyr)
tumor_re <- as.data.frame(tumor_re)
# 按照某列的中值将数据分为高组和低组
tumor_re <- tumor_re %>%
  mutate(Group = ifelse(tumor_re$`Macrophages M0` > median(tumor_re$`Macrophages M0`, na.rm = TRUE), "High", "Low"))

# 查看分组后的数据框
print(tumor_re)



## sample group
sample_m0_group <- as.data.frame(row.names = rownames(tumor_re),tumor_re$Group)


## 提取count表达谱
tnbc_tumor_count <- exp[,rownames(tumor_re)]




