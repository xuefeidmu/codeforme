# 1. 加载必要的包
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  install.packages("WGCNA")
}
library(WGCNA)

# 设置参数，防止分析结果受数据集内的行名影响
options(stringsAsFactors = FALSE)

# 2. 输入文件路径
input_file <- "path/to/your/input_file.csv"  # 替换为实际文件路径

# 3. 读取数据并进行详细解析
#    WGCNA分析需要每一行为一个样本，每一列为一个基因的表达量。表达量是TPM。
#    第一列通常是样本ID。
data <- read.csv(input_file, header = TRUE, row.names = 1)

# 检查输入数据的基本信息
cat("数据集的维度: ", dim(data), "\n")
cat("前几行数据: \n")
print(head(data))

# 4. 数据检查和清洗
#    检查是否有缺失值并进行处理
if (any(is.na(data))) {
  cat("数据包含缺失值。开始处理缺失值...\n")
  # 使用替代方法处理缺失值，例如使用列的平均值替代
  data[is.na(data)] <- apply(data, 2, function(x) mean(x, na.rm = TRUE))
} else {
  cat("数据没有缺失值。\n")
}

# 确保数据足够的样本和基因数目
if (ncol(data) < 20) {
  stop("基因数目少于20个，可能无法进行有效的网络分析。请检查输入数据。")
}
if (nrow(data) < 15) {
  stop("样本数少于15个，可能无法进行有效的网络分析。请检查输入数据。")
}

# 5. 数据标准化
#    WGCNA要求基因表达数据通常进行标准化（如log2转换），TPM数据也应log2处理
cat("对数据进行log2标准化...\n")
data <- log2(data + 1)

# 6. 样本聚类，检查是否有样本的质量不佳
sampleTree <- hclust(dist(data), method = "average")
plot(sampleTree, main = "样本聚类", sub = "", xlab = "", cex.main = 1.2, cex.lab = 1.2)

# 设定距离阈值来切割树状图，自动划分样本组
cutHeight <- 100  # 根据你的数据调整该值

# 使用cutree函数根据距离阈值将样本划分为不同的组
clust <- cutree(sampleTree, h = cutHeight)

# 找到离群样本（属于单独组的样本）
table(clust)  # 查看每个组的样本数量
outliers <- names(clust[clust %in% which(table(clust) == 1)])  # 找出异常样本的名字

# 批量去除这些异常样本
data_clean <- data[!rownames(data) %in% outliers, ]

# 查看处理后的数据
print(paste("移除了", length(outliers), "个异常样本：", paste(outliers, collapse = ", ")))

# 7. 软阈值功率选择
#    找到合适的软阈值功率来确保网络的无标度属性
powers <- c(1:20)
sft <- pickSoftThreshold(data_clean, powerVector = powers, verbose = 5)

# 绘制软阈值选择图
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     type = "n", xlab = "软阈值功率", ylab = "无标度拓扑模型适合度 (R^2)",
     main = "软阈值功率选择")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.9, col = "red")  # R^2 = 0.9 为参考线

# 8. 网络构建和模块识别
#    根据软阈值功率选择适合的值
softPower <- 6  # 你可以根据上面选出的软阈值调整这个值
adjacency <- adjacency(data_clean, power = softPower)

# 将邻接矩阵转换为拓扑重叠矩阵（TOM），并计算模块
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# 层次聚类和模块切割
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "基因聚类树", xlab = "", sub = "")

# 动态树切割，检测基因模块
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
moduleColors <- labels2colors(dynamicMods)

# 9. 模块可视化
plotDendroAndColors(geneTree, moduleColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# 10. 输出模块和基因的对应关系
module_gene_relationship <- data.frame(Gene = colnames(data_clean), Module = moduleColors)
write.csv(module_gene_relationship, "module_gene_relationship.csv", row.names = FALSE)

cat("WGCNA分析完成。模块和基因的对应关系已保存为module_gene_relationship.csv\n")
