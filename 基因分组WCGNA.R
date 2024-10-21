# 加载必要的库并导入数据
library(WGCNA)
options(stringsAsFactors = FALSE)

# 读取你的表达矩阵数据
data <- read.csv("你的表达矩阵文件.csv", row.names = 1)

# 根据某个基因按中位值分组
# 假设你要用基因geneA来将样本分成高表达组和低表达组
gene_name <- "geneA"
median_value <- median(data[gene_name, ])
group <- ifelse(data[gene_name, ] > median_value, "High", "Low")

# 将数据按分组进行分开
data_high <- data[, group == "High"]
data_low <- data[, group == "Low"]

# 转置表达矩阵，因为WGCNA要求行是样本，列是基因
datExpr <- t(data)

# 检查数据是否适合进行WGCNA分析
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 样本聚类以检测异常
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")

# 选择软阈值
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 可视化软阈值选择
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")

# 选择合适的软阈值
softPower <- 6  # 假设选择了6作为软阈值
adjacency <- adjacency(datExpr, power = softPower)

# 将邻接矩阵转换为拓扑重叠矩阵 (TOM)，并计算相似度
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# 基于TOM对基因进行聚类
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "Gene clustering on TOM-based dissimilarity",
     xlab = "", sub = "")

# 模块识别：使用动态树切割方法识别模块
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

table(dynamicMods)

# 将模块用颜色表示
dynamicColors <- labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# 计算模块的特征基因 (module eigengenes)
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

# 计算每个模块与高低表达分组的相关性
group_numeric <- as.numeric(factor(group, levels = c("Low", "High")))
moduleTraitCor <- cor(MEs, group_numeric, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = ncol(datExpr))

# 可视化模块与性状的相关性
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "Group",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))
