
library(TCGAbiolinks)
library(stringr)
library(SummarizedExperiment)
library(tibble)
##加载数据
load("BRCA.TCGA.rda")
## 提取表达矩阵count
exp=assay(expquery) %>% 
  as.data.frame()
## 读取肿瘤表型
clinical <- colData(expquery) %>% 
  as.data.frame()
library(IOBR)
# Using the TCGA count data as an example
data(eset_stad, package = "IOBR")
# Transformation is accompanied by gene annotation
exp_tpm <- count2tpm(countMat = exp, source = "local", idType = "Ensembl")
head(exp_tpm)
exp_tpm[1:2,1:2]

## 提取三阴乳腺癌样本名
sample <- rownames(clinical)[which(clinical$paper_BRCA_Subtype_PAM50=="Basal")]


## 提取生存信息,去除生存时长≤0的
surv <- clinical[sample,c("vital_status","paper_days_to_last_followup")];surv$status <- ifelse(surv$vital_status=="Alive",0,1);surv$time <- as.numeric(surv$paper_days_to_last_followup)
surv$sample <- sample
surv <- na.omit(surv);surv <- surv[surv$time > 0, ]


## 文件夹放入jVenn.csv
jVenn <- read.csv("jVenn.csv")
jVenn <- as.data.frame(jVenn)
gene <- jVenn$List.1.List.2
gene <- gene[gene != ""]



# 加载必要的包
library(glmnet)
library(survival)

# 准备生存对象
surv_object <- Surv(time = surv$time, event = surv$status)

# 将表达矩阵转换为矩阵格式##行名是样本，列名是基因symbol
x <- t(as.matrix(exp_tpm[,rownames(surv)]))
match <- match(gene,colnames(x))
# 去除NA
match <- match[!is.na(match)]

# 输出去除 NA 后的向量
print(match)

x <- x[,match]

# 使用 LASSO（alpha = 0.5）进行 Cox 回归，使用交叉验证选择最佳参数
fit <- cv.glmnet(x, surv_object, family = "cox", alpha = 1)

# 参数详解：
# - x: 解释变量矩阵（基因表达矩阵），行表示样本，列表示特征（基因）
# - surv_object: 生存对象，包含生存时间和生存状态，通过 Surv() 函数创建
# - family = "cox": 指定模型类型为 Cox 回归（适合生存分析）
# - alpha = 0.5: 正则化混合参数，0.5 表示结合 L1 和 L2 惩罚（弹性网络）

# 获取具有最优 λ 的模型系数，并转换为矩阵格式
coef_matrix <- as.matrix(coef(fit, s = "lambda.min"))

# 参数详解：
# - coef(fit, s = 0.01): 从交叉验证得到的模型中提取系数，s = 0.01 表示指定的 λ 值
# - as.matrix(): 将提取的系数对象（S4 类）转换为标准矩阵格式，方便后续处理

# 提取具有非零系数的基因名称
selected_genes <- rownames(coef_matrix)[coef_matrix != 0]

# 输出具有预测意义的基因
print(selected_genes)


# 可视化部分

# 1. 绘制交叉验证曲线
plot(fit)
title("Cross-Validation Curve for LASSO Cox Regression", line = 2.5)

# 参数详解：
# - plot(fit): 绘制交叉验证曲线，横轴为 Log(λ)，纵轴为偏似然偏差（Partial Likelihood Deviance）
# - title(): 为图添加标题，描述图的意义

# 2. 绘制 LASSO 路径图（系数路径图）
fit_glmnet <- glmnet(x, surv_object, family = "cox", alpha = 1)
plot(fit_glmnet, xvar = "lambda", label = TRUE)
title("Coefficient Path for LASSO Cox Regression", line = 2.5)

# 参数详解：
# - glmnet(x, surv_object, family = "cox", alpha = 1): 直接拟合 LASSO Cox 模型，alpha = 1 表示完全使用 LASSO 正则化
# - plot(fit_glmnet, xvar = "lambda", label = TRUE): 绘制系数路径图，横轴为 Log(λ)，展示不同基因系数随着 λ 的变化
# - label = TRUE: 为路径图中的系数添加标签，标识出不同的基因名称
# - title(): 为图添加标题，描述图的意义

# 通过以上步骤，LASSO 选择了一些具有预测意义的基因，并通过可视化帮助理解正则化过程中的特征选择与 λ 值的影响。

