# 加载所需的R包
library(readxl) # 用于读取Excel文件
library(tidyverse) # 包含ggplot2等数据处理和可视化包
library(GEOquery) # 用于从GEO数据库下载数据
library(limma) # 用于差异表达分析
library(affy) # 用于Affymetrix芯片数据处理
library(stringr) # 用于字符串处理

# 从GEO数据库下载GSE14520数据集的soft文件
# 并存储在当前目录下
gse <- getGEO(filename = "GSE14520_family.soft.gz", destdir = ".")

# 查看gse对象的结构
str(gse)

# 查看gse对象中包含的平台数量
length(gse)

# 提取第一个样本的表达数据作为示例
# (GSM362948只是一个例子，你可以根据需要选择其他样本)
sample_data <- gse@gsms[["GSM362948"]]@dataTable@table

# 查看示例样本数据的前几行
head(sample_data)

# 获取所有样本的名称
sample_names <- names(gse@gsms)

# 创建一个数据框，用于存储样本和平台的对应关系
sample_platform_df <- data.frame(Sample = character(), Platform = character())

# 遍历每个样本，提取其所属平台信息
for (sample in sample_names) {
  platform <- gse@gsms[[sample]]@header$platform_id
  sample_platform_df <- rbind(sample_platform_df, data.frame(Sample = sample, Platform = platform))
}

# 打印样本与平台的对应关系
print(sample_platform_df)

# 根据平台信息，将样本分为两组
samples_GPL571 <- sample_platform_df$Sample[sample_platform_df$Platform == "GPL571"]
samples_GPL3921 <- sample_platform_df$Sample[sample_platform_df$Platform == "GPL3921"]

# 创建两个列表，用于存储不同平台的表达数据
expression_list_GPL571 <- list()
expression_list_GPL3921 <- list()

# 提取属于GPL571平台的样本的表达数据
for (sample in samples_GPL571) {
  sample_data <- gse@gsms[[sample]]@dataTable@table
  expression_list_GPL571[[sample]] <- sample_data
}

# 提取属于GPL3921平台的样本的表达数据
for (sample in samples_GPL3921) {
  sample_data <- gse@gsms[[sample]]@dataTable@table
  expression_list_GPL3921[[sample]] <- sample_data
}

# 使用Reduce函数将同一平台的样本表达数据合并成一个矩阵
# 不同的样本数据根据"ID_REF"列进行合并
expression_matrix_GPL571 <- Reduce(function(x, y) merge(x, y, by = "ID_REF", all = TRUE), expression_list_GPL571)
expression_matrix_GPL3921 <- Reduce(function(x, y) merge(x, y, by = "ID_REF", all = TRUE), expression_list_GPL3921)

# 查看合并后的表达矩阵的前几行
head(expression_matrix_GPL571)
head(expression_matrix_GPL3921)

# 从GEO平台注释信息中提取探针ID和基因Symbol
id_probe_GPL571 <- gse@gpls[["GPL571"]]@dataTable@table
id_probe_GPL3921 <- gse@gpls[["GPL3921"]]@dataTable@table

# 只保留"ID"和"Gene Symbol"两列
probe2gene_GPL571 <- id_probe_GPL571[, c("ID", "Gene Symbol")]
probe2gene_GPL3921 <- id_probe_GPL3921[, c("ID", "Gene Symbol")]

# 定义一个函数，用于从Gene Symbol中提取标准的基因Symbol
# (去除多余的///等符号)
extract_gene_symbol <- function(gene_assignment) {
  symbols <- str_extract(gene_assignment, "^[^///]+")
  return(symbols)
}

# 应用函数提取基因Symbol
probe2gene_GPL571$Symbol <- extract_gene_symbol(probe2gene_GPL571$`Gene Symbol`)
probe2gene_GPL3921$Symbol <- extract_gene_symbol(probe2gene_GPL3921$`Gene Symbol`)

# 查看转换后的探针信息
head(probe2gene_GPL571)

# 将探针ID转换为基因Symbol
# 并将表达矩阵与探针信息合并
expression_matrix_GPL571 <- merge(probe2gene_GPL571[, c("ID", "Symbol")], expression_matrix_GPL571, by.x = "ID", by.y = "ID_REF", all.y = TRUE)
expression_matrix_GPL3921 <- merge(probe2gene_GPL3921[, c("ID", "Symbol")], expression_matrix_GPL3921, by.x = "ID", by.y = "ID_REF", all.y = TRUE)

# 定义一个函数，用于去除重复的基因Symbol
# (如果一个基因对应多个探针，则取表达量的平均值)
aggregate_expression <- function(expr_matrix) {
  expr_matrix <- expr_matrix %>%
    group_by(Symbol) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) # 对所有列计算平均值
  return(expr_matrix)
}

# 应用函数去除重复的基因Symbol
expression_matrix_GPL571 <- aggregate_expression(expression_matrix_GPL571)
expression_matrix_GPL3921 <- aggregate_expression(expression_matrix_GPL3921)

# 删除"ID"列(探针ID已经不再需要)
expression_matrix_GPL571 <- expression_matrix_GPL571[,-1]
expression_matrix_GPL3921 <- expression_matrix_GPL3921[,-1]

# 定义一个函数来处理单个表达矩阵
process_expression_matrix <- function(expression_matrix, matrix_name) {
  # 1. 找出重复的基因符号
  duplicated_symbols <- expression_matrix$Symbol[duplicated(expression_matrix$Symbol)]

  # 打印重复的基因符号，方便查看
  print(paste("矩阵", matrix_name, "中重复的基因符号："))
  print(duplicated_symbols)

  # 2. 对重复的基因符号取平均值
  averaged_expression_matrix <- aggregate(. ~ Symbol, data = expression_matrix, FUN = mean)

  # 3. 将基因符号设置为行名
  rownames(averaged_expression_matrix) <- averaged_expression_matrix$Symbol
  averaged_expression_matrix$Symbol <- NULL

  # 4. 打印处理后的表达矩阵信息
  print(paste("处理后的矩阵", matrix_name, "信息："))
  print(paste("行数 (基因数量):", nrow(averaged_expression_matrix)))
  print(paste("列数 (样本数量):", ncol(averaged_expression_matrix)))

  return(averaged_expression_matrix)
}

# 处理 expression_matrix_GPL571
expression_matrix_GPL571_processed <- process_expression_matrix(expression_matrix_GPL571, "expression_matrix_GPL571")

# 处理 expression_matrix_GPL3921
expression_matrix_GPL3921_processed <- process_expression_matrix(expression_matrix_GPL3921, "expression_matrix_GPL3921")

# 合并两个表达矩阵
merged_expression_matrix <- merge(expression_matrix_GPL571_processed,
                                 expression_matrix_GPL3921_processed,
                                 by = "row.names",
                                 all = FALSE)

# 将行名设置回数据框的行名
rownames(merged_expression_matrix) <- merged_expression_matrix$Row.names
merged_expression_matrix$Row.names <- NULL

# 打印合并后的表达矩阵信息
print("合并后的表达矩阵信息：")
print(paste("行数 (基因数量):", nrow(merged_expression_matrix)))
print(paste("列数 (样本数量):", ncol(merged_expression_matrix)))

# 提取GPD2基因的表达数据
gpd2_expression <- merged_expression_matrix["GPD2", ]

# 转换为适合绘图的数据格式
