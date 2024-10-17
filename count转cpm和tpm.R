
## 安装IOBR package
if (!requireNamespace("IOBR", quietly = TRUE))
  devtools::install_github("IOBR/IOBR")

## 加载IOBR
library(IOBR) 

# 更新 IOBR
devtools::install_github("IOBR/IOBR", ref = "master")

# 加载 IOBR 包
library("IOBR")

# 查看帮助文档
help("count2tpm")

# 使用 count2tpm 函数
count2tpm(
  countMat,               # 表达矩阵，行为基因
  idType = "Ensembl",     # "Ensembl", "ENTREZ", "SYMBOL"
  org = "hsa",            # Organism, hsa 或者 mmu, 只对通过 biomaRt 注释的流程有效
  source = "web",         # 默认是通过 biomaRt 来获取长度, 如果设为“default”且 idType = "Ensembl",
                          # 提供了 length_ensembl 对象，通过 gencode v22 计算获得长度
  effLength = NULL,       # 用户可自己提供基因长度的数据，如果提供了 effLength 需要明确
  id = "id",              # 对应表达矩阵行名
  gene_symbol = "symbol", # 哪一列是基因 symbol
  length = "eff_length",  # 哪一列是基因长度
  remove_redundancy = "mean" # 使用什么方法来去掉重复的基因：默认 mean, 还可以选择 sd 或者 median
)



##或者
# 使用TCGA的count数据作为示例
data(eset_stad, package = "IOBR")  # 从 IOBR 包中加载名为 "eset_stad" 的数据，作为表达矩阵示例

# 进行转换，并伴随基因注释
eset <- count2tpm(                 # 使用 count2tpm 函数，将 count 数据转换为 TPM（Transcripts Per Million）
  countMat = eset_stad,            # 输入的表达矩阵是 "eset_stad" 
  source = "local",                # 设置 source 为 "local"，表示基因长度数据是本地已有的，不通过biomaRt获取
  idType = "ensembl"               # 基因的ID类型为 Ensembl ID
)

# 显示转换后的结果的前几行
head(eset)                         # 使用 head 函数查看转换后的表达矩阵的前几行






