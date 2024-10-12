## 加载数据
load("BRCA.TCGA.rda")


## 提取表达矩阵count
exp=assay(expquery) %>% 
  as.data.frame()

## COUNT表达矩阵行名ID转换
library(tinyarray)
exp=trans_exp_new(exp) %>% 
  as.data.frame()
exp[1:4,1:4]


## 读取肿瘤表型
clinical <- colData(expquery) %>% 
  as.data.frame()



tnbc_tumor <- clinical$barcode[which(clinical$paper_BRCA_Subtype_PAM50=="Basal")]

## 提取count exp中三阴乳腺癌肿瘤样本表达谱
tnbc.tumor.count <- exp[,tnbc_tumor]
range(tnbc.tumor.count)


## 加载estimate包
library(utils)
rforge<-"http://r-forge.r-project.org"
if(!require("estimate"))install.packages("estimate",repos=rforge,dependencies=TRUE)
library(estimate)
#help(package="estimate")


## 把自己的count存入exprSet
exprSet <- tnbc.tumor.count



dat=log2(edgeR::cpm(exprSet)+1)

library(estimate)
estimate<-function(dat,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file=input.f,sep='\t',quote=F)
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f,
                    id="GeneSymbol")##注意这个id
  estimateScore(input.ds=output.f,
                output.ds=output.ds,
                platform="illumina")##注意platform，转录组数据illumina，芯片affymetrix
  scores=read.table(output.ds,skip=2,header=T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro='BRCA'
scores=estimate(dat,pro) %>% as.data.frame()


## 查看结果
head(scores)
view(scores)



## 根据公式计算TumorPurity
TumorPurity=cos(0.6049872018+0.0001467884*scores[,3])
head(TumorPurity)
view(TumorPurity)



## 根据immunescore分组
tnbc.immune.score.group <- ifelse(scores$ImmuneScore > median(scores$ImmuneScore), 'high', 'low')
tnbc.immune.score.group<-factor(tnbc.immune.score.group,levels=c("high","low"))
table(tnbc.immune.score.group)




## 使用deseq2进行差异分析

## 准备数据
expr <- tnbc.tumor.count
group_list <- tnbc.immune.score.group



# Load the DESeq2 library for RNA-seq analysis
library(DESeq2)
# Create a data frame for metadata (colData), setting row names as column names of the expression matrix
# `group_list` is used to define experimental groups
colData <- data.frame(row.names = colnames(expr), condition = group_list)
# Create a DESeqDataSet object, providing the count data, colData, and design formula
dds <- DESeqDataSetFromMatrix(
  countData = expr,    # Expression data matrix (count data)
  colData = colData,   # Metadata about the samples
  design = ~ condition # Experimental design formula
)
# Set reference level for the factor 'condition' to the control group
dds$condition <- relevel(dds$condition, ref = "low") # "low" is the reference condition
# Run the DESeq function to perform differential expression analysis
dds <- DESeq(dds)
# Perform pairwise comparisons between conditions (the reverse order to get "treatment vs control")
res <- results(dds, contrast = c("condition", rev(levels(group_list))))
# Order the results by p-value to prioritize most significant changes
resOrdered <- res[order(res$pvalue),] 
# Convert results to a data frame for easy manipulation
DEG <- as.data.frame(resOrdered)
# Display the first few rows of the results
head(DEG)
# Remove rows with NA values (missing data)
DEG <- na.omit(DEG)
# Add a column named `change` to categorize genes as 'UP', 'DOWN', or 'NOT' based on log2FoldChange and p-value
# Set cutoff for log2 fold change
logFC_cutoff <- 2 # Custom cutoff for filtering differentially expressed genes
# Create a new column in DEG to mark genes as upregulated ('UP'), downregulated ('DOWN'), or not significant ('NOT')
DEG$change <- as.factor(
  ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
         ifelse(DEG$log2FoldChange > logFC_cutoff, 'UP', 'DOWN'), 'NOT')
)
# Display the first few rows with the added 'change' column
head(DEG)

# Note:
# The script sets a log2 fold change threshold of 1, and marks genes as 'UP' if log2FoldChange > cutoff, 'DOWN' if < -cutoff.
# p-value < 0.05 is used to define significance. Genes that do not meet these criteria are marked as 'NOT'.

table(DEG$change)

## 保存DEG为xls
write.csv(DEG,"tcga.brca.tnbc.immune.high.vs.low.deg.csv")



## 火山图可视化
# 加载 ggplot2 库
library(ggplot2)

# 设置 logFC 和 p 值的阈值
logFC_cutoff <- 2
pvalue_cutoff <- 0.05

# 计算上调、下调和不显著基因的数量
deg_counts <- table(DEG$change)
up_count <- deg_counts["UP"]
down_count <- deg_counts["DOWN"]

# 绘制火山图并添加参考线和上调、下调基因的数量
ggplot(DEG, aes(x = log2FoldChange, y = -log10(pvalue), color = change)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NOT" = "grey")) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), color = "black", linetype = "dashed") + # logFC 阈值参考线
  geom_hline(yintercept = -log10(pvalue_cutoff), color = "black", linetype = "dashed") + # p 值阈值参考线
  annotate("text", x = max(DEG$log2FoldChange) - 1, y = max(-log10(DEG$pvalue)) - 1, 
           label = paste("Upregulated:", up_count), color = "red", size = 3, hjust = 1) +
  annotate("text", x = min(DEG$log2FoldChange) + 1, y = max(-log10(DEG$pvalue)) - 1, 
           label = paste("Downregulated:", down_count), color = "blue", size = 3, hjust = 0) +
  theme_minimal() +
  labs(title = "TCGA.TNBC.Immunescore high vs. Immunescore low", 
       x = "Log2 Fold Change", 
       y = "-Log10 P-value") +
  theme(legend.position = "top")






