##根据一个基因高低表达分组后，对两组进行WCGNA分析
# 加载必要的库并导入数据
library(WGCNA)
options(stringsAsFactors = FALSE)

# 读取你的表达矩阵数据（注意：表达矩阵数据要求为TPM值，并在WGCNA分析之前进行log2(TPM+1)处理）
data <- read.csv("你的表达矩阵文件.csv", row.names = 1)

# 对TPM表达矩阵进行log2(TPM+1)处理
data <- log2(data + 1)

# 根据某个基因按中位值分组
# 假设你要用基因geneA来将样本分成高表达组和低表达组
gene_name <- "geneA"
median_value <- median(as.numeric(data[gene_name, ]))
group <- ifelse(data[gene_name, ] > median_value, "High", "Low")

# 将数据按分组进行分开
data_high <- data[, group == "High"]
data_low <- data[, group == "Low"]

# 转置表达矩阵，因为WGCNA要求行是样本，列是基因
datExpr <- t(data)

# 检查数据是否适合进行WGCNA分析
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  # 排除有问题的样本和基因
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 样本聚类以检测异常
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", labels = FALSE)

# 排除异常样本（通过树状图手动选择要排除的样本）
# 假设我们手动选择了一些异常样本并将它们排除
clust <- cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)
keepSamples <- (clust == 1)
datExpr <- datExpr[keepSamples, ]

# 选择软阈值
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 自动选择合适的软阈值
softPower <- sft$powerEstimate

# 检查softPower是否为NA并重新计算
if (is.na(softPower)) {
  softPower <- 6  # 如果未找到合适的软阈值，使用默认值6
}

# 可视化软阈值选择
if (!is.na(softPower)) {
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")

# 在最佳软阈值处画一条横线
  abline(v = softPower, col = "blue", lty = 2)
}

# 选择合适的软阈值
softPower <- sft$powerEstimate  # 自动选择的软阈值
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
