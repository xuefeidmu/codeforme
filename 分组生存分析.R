
# 加载生存分析和可视化包
library(survival)
library(survminer)

# 假设您的数据集命名为 my_data
# my_data 数据集包含以下列：
# time: 生存时间
# status: 生存状态 (0 = 存活, 1 = 死亡)
# group: 分组变量 (如 group1, group2)

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
  xlab = "Time (months)",    # X轴标签
  ylab = "Survival Probability", # Y轴标签
  palette = c("#E7B800", "#2E9FDF")  # 设置不同组的曲线颜色
)
