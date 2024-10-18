
library(TCGAbiolinks)
library(stringr)
library(SummarizedExperiment)
library(tibble)

load("BRCA.TCGA.rda")## 这个数据需要放在工作文件夹里

## 提取表达矩阵count
exp=assay(expquery) %>% 
  as.data.frame()

## COUNT表达矩阵行名ID转换(可以暂时不改，用count2tpm函数直接变)
library(tinyarray)
exp=trans_exp_new(exp) %>% 
  as.data.frame()
exp[1:4,1:4]


## 读取肿瘤表型
clinical <- colData(expquery) %>% 
  as.data.frame()


library(IOBR)## 这个包用来count转tpm，同时ensembl转gene symbol
# Transformation is accompanied by gene annotation
exp_tpm <- count2tpm(countMat = exp, source = "local", idType = "Ensembl")
head(exp_tpm)
exp_tpm[1:2,1:2]


## 提取basal的样本名
sample <- rownames(clinical)[which(clinical$paper_BRCA_Subtype_PAM50=="Basal")]


##设置要研究的基因
gene <- "COL1A1"


## 分组
exp1 <- as.data.frame(t(exp_tpm[gene,sample]))
group <- ifelse(exp1[,1]>median(exp1[,1]),"high","low")
  

## 提取生存信息
surv <- clinical[sample,c("vital_status","paper_days_to_last_followup")]
surv$status <- ifelse(surv$vital_status=="Alive",0,1)
surv$group <- group
surv$time <- as.numeric(surv$paper_days_to_last_followup)
my_data <- surv


# 加载生存分析和可视化包
library(survival)
library(survminer)
# 构造生存对象
surv_object <- Surv(time = my_data$time, event = my_data$status)

# 根据分组变量 group 进行生存分析
fit <- survfit(surv_object ~ group, data = my_data)


# 打印生存分析结果
print(fit)

# 可视化生存曲线，区分不同的分组
ggsurvplot(
  fit,                      # Kaplan-Meier拟合对象
  data = my_data,            # 数据集
  pval = TRUE,               # 显示p值，判断组间生存差异是否显著
  conf.int = TRUE,           # 显示置信区间
  risk.table = TRUE,         # 显示风险表
  legend.title = "Group",    # 图例标题
  legend.labs = c("Group 1", "Group 2"), # 根据您的分组修改图例标签
  xlab = "Time (day)",    # X轴标签
  ylab = "Survival Probability", # Y轴标签
  palette = c("#E7B800", "#2E9FDF")  # 设置不同组的曲线颜色
)






