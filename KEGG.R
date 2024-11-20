# 安装必要的R包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "pathview"))

# 加载必要的R包
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(ggplot2)

# 准备基因列表（Gene Symbol 格式）
gene_symbols <- c("TP53", "EGFR", "AKT1", "BRCA1", "MTOR")

# 将Gene Symbol转换为Entrez ID
gene_list <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# 移除NA值
gene_list <- na.omit(gene_list)

# 执行KEGG富集分析
kegg_result <- enrichKEGG(gene = gene_list,
                          organism = 'hsa', # 'hsa'表示人类
                          pvalueCutoff = 0.05)

# 查看分析结果
print(head(kegg_result))

# 点图可视化
dotplot(kegg_result, showCategory = 10)

# 条形图可视化
barplot(kegg_result, showCategory = 10, title = "KEGG Pathway Enrichment")

# 保存KEGG结果为CSV文件
write.csv(as.data.frame(kegg_result), file = "kegg_results.csv", row.names = FALSE)

# 可选：生成KEGG通路图
# 选择一个感兴趣的通路
pathway_id <- "hsa05200" # 替换为感兴趣的KEGG通路ID
pathview(gene.data = gene_list,
         pathway.id = pathway_id,
         species = "hsa")
