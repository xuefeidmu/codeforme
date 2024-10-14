# 清空环境
rm(list = ls())  
load(file = "step2output.Rdata")
# 输入数据：exp 和 Group

# PCA分析
# 转置表达矩阵，并转化为数据框
dat=as.data.frame(t(exp))

# 进行主成分分析（PCA）
library(FactoMineR)
library(factoextra) 
dat.pca <- PCA(dat, graph = FALSE)

# 可视化PCA结果
fviz_pca_ind(dat.pca,
             geom.ind = "point", # 显示点
             col.ind = Group, # 按照Group着色
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # 添加浓度椭圆
             legend.title = "Groups"
)

# 绘制top 1000 sd基因的热图
g = names(tail(sort(apply(exp,1,sd)),1000)) 
# 获取标准差最大的1000个基因
n = exp[g,]

# 加载绘图包并设置注释
library(pheatmap)
annotation_col = data.frame(row.names = colnames(n),
                            Group = Group)

# 绘制热图
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         scale = "row", # 按行标准化
         breaks = seq(-3,3,length.out = 100) # 设置色带范围
)

# 差异分析
rm(list = ls()) 
load(file = "step2output.Rdata")

# 加载差异分析所需的包
library(limma)
design = model.matrix(~Group)
fit = lmFit(exp,design)
fit = eBayes(fit)

# 获取差异分析结果
deg = topTable(fit,coef = 2,number = Inf)

# 给deg添加probe_id列
library(dplyr)
deg = mutate(deg,probe_id = rownames(deg))

# 去重并连接探针注释
ids = distinct(ids,symbol,.keep_all = T)
deg = inner_join(deg,ids,by="probe_id")

# 设置logFC和pValue的阈值，标记上下调基因
logFC_t = 1
p_t = 0.05
k1 = (deg$P.Value < p_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < p_t)&(deg$logFC > logFC_t)
deg = mutate(deg,change = ifelse(k1,"down",ifelse(k2,"up","stable")))

# 打印上下调基因的数量
table(deg$change)

# 火山图
library(ggplot2)
ggplot(data = deg, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, aes(color=change)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",linewidth=0.8) +
  geom_hline(yintercept = -log10(p_t),lty=4,col="black",linewidth=0.8) +
  theme_bw()

# 差异基因热图
exp = exp[deg$probe_id,]
rownames(exp) = deg$symbol
diff_gene = deg$symbol[deg$change !="stable"]
n = exp[diff_gene,]

annotation_col = data.frame(group = Group)
rownames(annotation_col) = colnames(n) 

pheatmap(n,show_colnames =F,
         show_rownames = F,
         scale = "row",
         annotation_col=annotation_col,
         breaks = seq(-3,3,length.out = 100)
)

# 添加ENTREZID列，用于后续富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
s2e = bitr(deg$symbol, 
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = org.Hs.eg.db)

deg = inner_join(deg,s2e,by=c("symbol"="SYMBOL"))

# 保存结果
save(exp,Group,deg,logFC_t,p_t,file = "step4output.Rdata")
