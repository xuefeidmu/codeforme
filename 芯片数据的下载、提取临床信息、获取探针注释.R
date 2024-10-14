# 生信技能树-day8 GEO数据挖掘-芯片数据的下载、提取临床信息、获取探针注释

# 第一步：安装并加载R包
# 安装所需的R包
install.packages(c("GEOquery", "stringr", "limma", "tinyarray", "hgu133plus2.db"))

# 加载R包
library(GEOquery)
library(stringr)
library(limma)
library(tinyarray)
library(hgu133plus2.db)

# 第二步：下载数据，提取表达矩阵和临床信息
rm(list = ls())
options(timeout = 100000)  # 设置下载时间限制
options(scipen = 20)  # 不使用科学计数法显示
eSet = getGEO("GSE7305", destdir = '.', getGPL = F)  # 下载数据

# 检查下载结果
class(eSet)
eSet = eSet[[1]]  # 提取表达矩阵
exp <- exprs(eSet)  # 提取表达矩阵
dim(exp)  # 查看表达矩阵的维度
range(exp)  # 查看数据范围，判断是否需要log转换
exp = log2(exp + 1)  # log转换数据
boxplot(exp, las = 2)  # 绘制箱线图，检查异常样本

# 标准化数据
exp = normalizeBetweenArrays(exp)

# 第三步：提取临床信息
pd <- pData(eSet)
# 确保表达矩阵和临床数据列名一致
p = identical(rownames(pd), colnames(exp))
if(!p) {
  s = intersect(rownames(pd), colnames(exp))
  exp = exp[, s]
  pd = pd[s, ]
}

# 保存结果
save(pd, exp, file = "step1output.Rdata")

# 第四步：设置实验分组
k = str_detect(pd$title, "Normal")  # 根据关键词对分组进行处理
Group = ifelse(k, "Normal", "Disease")  # 分配分组信息
Group = factor(Group, levels = c("Normal", "Disease"))  # 设置因子并设置对照组为参考水平

# 第五步：获取探针注释
gpl_number <- eSet@annotation  # 获取GPL编号
find_anno(gpl_number)  # 使用tinyarray获取注释
ids <- toTable(hgu133plus2SYMBOL)  # 获取探针注释

# 保存结果
save(exp, Group, ids, file = "step2output.Rdata")

# 代码流程结束
