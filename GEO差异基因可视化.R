### 0. 项目组织
proj = "GSE210164"

### 1. 读取和整理数据
#### 1.1 表达矩阵
# 先读取一个文件，检查能否成功读取
a = read.delim("GSE210164_RAW/GSM6422899_Lewis-Lewis-1.txt.gz")

# 读取文件时保留列名、设置第一列作为行名，并只保留count列
a = read.delim("GSE210164_RAW/GSM6422899_Lewis-Lewis-1.txt.gz", check.names = F, row.names = 1)[1]

# 循环读取文件夹内所有的文件
files <- dir("GSE210164_RAW/")
dat = list()
setwd("GSE210164_RAW/")
for(i in 1:length(files)){
  dat[[i]] = read.delim(files[i], check.names = F, row.names = 1)[1]
}
dat = do.call(cbind, dat)  # 按列拼接数据，形成完整的表达矩阵
setwd("../")

# 转换成矩阵
exp = as.matrix(dat)

#### 1.2 临床信息
library(GEOquery)
eSet = getGEO("GSE210164", destdir = '.', getGPL = F)
eSet = eSet[[1]]
clinical <- pData(eSet)  # 提取临床信息

# 简化临床信息
library(tinyarray)
geo = geo_download("GSE210164", colon_remove = T)  # 去除多余列
clinical <- geo$pd

#### 1.3 表达矩阵行名ID转换
library(tinyarray)
exp = trans_exp_new(exp, species = "rat")  # 基因ID转换（大鼠）

#### 1.4 基因过滤
# 保留一半以上样本中有表达的基因
exp = exp[apply(exp, 1, function(x) sum(x > 0) > 0.5 * ncol(exp)), ]

#### 1.5 分组信息获取
library(stringr)
Group = clinical$title %>% str_split_i("-", 2)  # 按“-”分割，保留第二部分
Group = factor(Group, levels = c("Lewis", "BN"))

# 保存数据
save(exp, Group, proj, clinical, file = paste0(proj, ".Rdata"))

### 2. 差异分析
rm(list = ls())  # 清空环境变量
load("GSE210164.Rdata")

#### 2.1 DESeq2 差异分析
library(DESeq2)
colData <- data.frame(row.names = colnames(exp), condition = Group)
if (!file.exists(paste0(proj, "_dd.Rdata"))) {
  dds <- DESeqDataSetFromMatrix(countData = exp, colData = colData, design = ~ condition)
  dds <- DESeq(dds)
  save(dds, file = paste0(proj, "_dd.Rdata"))
}
load(file = paste0(proj, "_dd.Rdata"))
res <- results(dds, contrast = c("condition", rev(levels(Group))))
DEG1 <- as.data.frame(res)

#### 2.2 edgeR 差异分析
library(edgeR)
dge <- DGEList(counts = exp, group = Group)
dge <- calcNormFactors(dge)
design <- model.matrix(~Group)
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
fit <- glmLRT(fit)
DEG2 = as.data.frame(topTags(fit, n = Inf))

#### 2.3 limma 差异分析
library(limma)
dge <- edgeR::DGEList(counts = exp)
v <- voom(dge, design, normalize = "quantile")
fit <- lmFit(v, design)
fit= eBayes(fit)
DEG3 = topTable(fit, coef = 2, n = Inf)

#### 2.4 差异基因筛选
# 设置logFC和p值阈值
logFC_t = 2
pvalue_t = 0.05

# 标记差异基因上调、下调
DEG1$change = ifelse((DEG1$pvalue < pvalue_t) & (DEG1$log2FoldChange > logFC_t), "UP",
                     ifelse((DEG1$pvalue < pvalue_t) & (DEG1$log2FoldChange < -logFC_t), "DOWN", "NOT"))

### 3. 差异基因可视化
library(ggplot2)
library(tinyarray)

# PCA图
dat = log2(cpm(exp) + 1)
pca.plot = draw_pca(dat, Group)

# 热图
h1 = draw_heatmap(dat[rownames(DEG1)[DEG1$change != "NOT"], ], Group)

### 4. 富集分析
library(clusterProfiler)
deg = DEG1
s2e = bitr(deg$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)
deg = inner_join(deg, s2e, by = c("symbol" = "SYMBOL"))

ekk <- enrichKEGG(gene = deg$ENTREZID, organism = 'rat')
ego <- enrichGO(gene = deg$ENTREZID, OrgDb = org.Rn.eg.db, ont = "ALL", readable = TRUE)

# 可视化
plot_GO <- barplot(ego, split = "ONTOLOGY")
plot_KEGG <- barplot(ekk)

# 保存可视化结果
ggsave("GSE210164_GO.png", width = 15, height = 10)
ggsave("GSE210164_KEGG.png", width = 15, height = 10)
