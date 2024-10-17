# --------------------------
# 将 count 矩阵转换为 CPM 和 TPM 的完整流程
# --------------------------

# Step 1: 转换 count 表达矩阵为 CPM
# 使用 edgeR 包中的 cpm() 函数，将 count 数据转换为 CPM 表达量，并通过 log2 变换
exprSet <- log2(edgeR::cpm(exprSet) + 1)

# --------------------------
# Step 2: 转换 count 矩阵为 TPM
# 需要处理基因组注释文件，并计算每个基因的外显子长度

# 安装并加载 GenomicFeatures 包，如果没有安装则进行安装
if (!require(GenomicFeatures)) BiocManager::install("GenomicFeatures")
library(GenomicFeatures)

# 读取 gencode 注释文件，构建转录本数据库
txdb <- makeTxDbFromGFF("gencode.v36.annotation.gtf.gz", format = "gtf")

# 获取每个基因的外显子数据
exons.list.per.gene <- exonsBy(txdb, by = "gene")

# 对每个基因，将所有外显子减少为一组非重叠的外显子区域，计算这些区域的长度（宽度）并求和
exonic.gene.sizes <- sum(width(GenomicRanges::reduce(exons.list.per.gene)))

# 将基因ID和对应的外显子总长度整理为一个数据框
gfe <- data.frame(gene_id = names(exonic.gene.sizes), length = exonic.gene.sizes)
head(gfe)[1:5, 1:2]

# 输出前5行数据：
# gene_id                length
# ENSG00000000003.15     4536
# ENSG00000000005.6      1476
# ENSG00000000419.13     1207
# ENSG00000000457.14     6883
# ENSG00000000460.17     5970

# 保存计算的基因长度数据，以便后续使用
save(gfe, file = "gfe.Rdata")

# --------------------------
# Step 3: 使用基因长度将 count 数据转换为 TPM
# 加载基因长度数据
load("gfe.Rdata")

# 提取基因长度列
effLen <- gfe$length

# 定义将 count 数据转换为 TPM 的函数
Counts2TPM <- function(counts, effLen) {
  rate <- log(counts) - log(effLen)     # 计算表达速率
  denom <- log(sum(exp(rate)))           # 计算归一化因子
  exp(rate - denom + log(1e6))           # 计算 TPM
}

# 使用 Counts2TPM 函数，将 count 矩阵 hnsc 转换为 TPM
hnsc_tpm_raw <- apply(hnsc, 2, Counts2TPM, effLen = effLen)

# 输出转换后的 TPM 数据的前3行和前3列
head(hnsc_tpm_raw)[1:3, 1:3]

# 示例输出：
#                      TCGA-UF-A7JF-01A-11R-A34R-07  TCGA-CN-4725-01A-01R-1436-07  TCGA-D6-6827-01A-11R-1915-07
# ENSG00000000003.15              31.13713                        15.74248                      72.37614
# ENSG00000000005.6                0.00000                         0.00000                       0.00000
# ENSG00000000419.13             117.46366                       136.35307                     112.14231

# --------------------------
