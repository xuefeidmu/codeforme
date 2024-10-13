## 多个GEO数据合并


## 手动设置工作文件夹


## 1. 加载所需R包
# Load necessary libraries (avoiding redundancy)
library(readxl)
library(tidyverse)
library(GEOquery)
library(limma)
library(affy)
library(stringr)


## 2. 获取GEO数据
# Download GEO data或者手动下载放在文件夹里
# 运行下面代码可以自动识别工作文件夹的相应文件
gset <- getGEO('GSE205185', destdir=".", AnnotGPL = TRUE, getGPL = TRUE)
gset2 <- getGEO('GSE29431', destdir=".", AnnotGPL = TRUE, getGPL = TRUE)
gset3 <- getGEO('GSE20711', destdir=".", AnnotGPL = TRUE, getGPL = TRUE)

# Check the class and structure
class(gset)
gset[[1]]
class(gset2)
gset2[[1]]
class(gset3)
gset3[[1]]


## 3. 读取平台文件
# Get platform annotations
plf1 <- gset[[1]]@annotation
plf2 <- gset2[[1]]@annotation
plf3 <- gset3[[1]]@annotation

# Read platform files (example filenames; adjust as necessary)
GPL_data_11 <- Table(getGEO(filename = "GPL21185.soft.gz", AnnotGPL = TRUE))
GPL_data_22 <- Table(getGEO(filename = "GPL570.annot.gz", AnnotGPL = TRUE))
GPL_data_33 <- Table(getGEO(filename = "GPL570.annot.gz", AnnotGPL = TRUE))






## 4. 提取基因表达矩阵
# Extract expression data for each dataset
exp <- exprs(gset[[1]])
exp2 <- exprs(gset2[[1]])
exp3 <- exprs(gset3[[1]])

# Store probe names
probe_name <- rownames(exp)
probe_name2 <- rownames(exp2)
probe_name3 <- rownames(exp3)




## 5. ID转换


## 常用的转换ID的方法
find_anno("GPL96")
ids <- AnnoProbe::idmap('GPL96')
rownames(exp3)=ids$symbol[match(rownames(exp3),ids$probe_id)]



# For gset1
loc1 <- match(GPL_data_11[,1], probe_name)
probe_exp1 <- exp[loc1,]
raw_geneid1 <- as.matrix(GPL_data_11[,"GENE_SYMBOL"])
index1 <- which(!is.na(raw_geneid1))
geneid1 <- raw_geneid1[index1]
exp_matrix1 <- probe_exp1[index1,]
geneidfactor1 <- factor(geneid1)
gene_exp_matrix1 <- apply(exp_matrix1, 2, function(x) tapply(x, geneidfactor1, mean))
rownames(gene_exp_matrix1) <- levels(geneidfactor1)

# For gset2
loc2 <- match(GPL_data_22[,1], probe_name2)
probe_exp2 <- exp2[loc2,]
raw_geneid2 <- as.matrix(GPL_data_22[,"Gene symbol"])
index2 <- which(!is.na(raw_geneid2))
geneid2 <- raw_geneid2[index2]
exp_matrix2 <- probe_exp2[index2,]
geneidfactor2 <- factor(geneid2)
gene_exp_matrix2 <- apply(exp_matrix2, 2, function(x) tapply(x, geneidfactor2, mean))
rownames(gene_exp_matrix2) <- levels(geneidfactor2)

# For gset3
loc3 <- match(GPL_data_33[,1], probe_name3)
probe_exp3 <- exp3[loc3,]
raw_geneid3 <- as.matrix(GPL_data_33[,"Gene symbol"])
index3 <- which(!is.na(raw_geneid3))
geneid3 <- raw_geneid3[index3]
exp_matrix3 <- probe_exp3[index3,]
geneidfactor3 <- factor(geneid3)
gene_exp_matrix3 <- apply(exp_matrix3, 2, function(x) tapply(x, geneidfactor3, mean))
rownames(gene_exp_matrix3) <- levels(geneidfactor3)




## 6. 多矩阵合并
# Convert to data frames and merge datasets
geo_exp_1 <- as.data.frame(gene_exp_matrix)
geo_exp_2 <- as.data.frame(gene_exp_matrix2)
geo_exp_3 <- as.data.frame(gene_exp_matrix3)

# Find common samples across datasets
sameSample <- intersect(rownames(geo_exp_1), rownames(geo_exp_2))
sameSample <- intersect(rownames(geo_exp_3), sameSample)

# Filter datasets to common samples
gene_exp1 <- geo_exp_1[sameSample, , drop = FALSE]
gene_exp2 <- geo_exp_2[sameSample, , drop = FALSE]
gene_exp3 <- geo_exp_3[sameSample, , drop = FALSE]

# Combine the data
bindgeo <- cbind(gene_exp1, gene_exp2, gene_exp3)





## 7. 读取分组信息 这里需要自行修改
# Read grouping information (modify based on your actual pdata structure)
pdata <- pData(gset[[1]])
group1 <- as.matrix(pdata[,"source_name_ch1"])

                          
pdata2 <- pData(gset2[[1]])
group2 <- as.matrix(pdata2[,"group"]) ##修改"group"

                        
pdata3 <- pData(gset3[[1]])
group3 <- as.matrix(pdata3[,"group"]) ##修改"group"

# Combine group data
talgroup <- as.data.frame(rbind(group1, group2, group3))
talgroup_list <- factor(talgroup$group, levels = c("N", "T"))

# Save group information
write.csv(talgroup, file = "group.csv")


## 8. 数据矫正
# Boxplot before normalization
boxplot(bindgeo, outline = TRUE, notch = TRUE, col = talgroup_list, las = 2)
dev.off()

# Normalize between arrays
bindgeo_normal <- normalizeBetweenArrays(bindgeo)
bindgeo_normal <- log2(bindgeo_normal + 1)

# Save normalized data
bindgeo_normal <- as.data.frame(bindgeo_normal)
bindgeo_normal <- na.omit(bindgeo_normal)
write.csv(bindgeo_normal, file = "bindgeo_exp.csv")









## 9. 差异分析
# Differential expression analysis
design <- model.matrix(~talgroup_list)
fit <- lmFit(bindgeo_normal, design)
fit <- eBayes(fit)
deg <- topTable(fit, coef = 2, number = Inf)

# Save DEGs
write.table(deg, file = "deg_all.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# Identify up and down regulated genes
logFC <- 1
adj.P.Val <- 0.05
k1 <- (deg$adj.P.Val < adj.P.Val) & (deg$logFC < -logFC)
k2 <- (deg$adj.P.Val < adj.P.Val) & (deg$logFC > logFC)
deg$change <- ifelse(k1, "down", ifelse(k2, "up", "stable"))

# Save up and down regulated genes
write.csv(deg, file = "upanddown.csv")










