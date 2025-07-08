# 加载必要的库
library(DESeq2)

# 读取数据
exprSet <- read.csv("counts.csv", header = TRUE)
exprSet <- data.frame(exprSet)
rownames(exprSet) <- exprSet$X
exprSet <- exprSet[which(!(rownames(exprSet) %in% "-")), ]
exprSet <- exprSet[, -1]

# 创建样本信息数据框
data <- data.frame(row.names = colnames(exprSet),
                   group = c("SU_WT", "SU_WT", "SU_WT", "SU_KO", "SU_KO", "SU_KO"),
                   sample = colnames(exprSet))

# 过滤数据
exprSet <- exprSet[rowSums(exprSet) > 10 * ncol(exprSet), ]

# 创建分组向量
group <- factor(c(rep("SU_WT", 3), rep("SU_KO", 3)))
Data <- data.frame(row.names = colnames(exprSet), group = group)

# 创建 DESeq2 数据集
dds <- DESeqDataSetFromMatrix(countData = exprSet, colData = Data, design = ~ group)
dds2 <- DESeq(dds)

# 进行差异表达分析
tmp <- results(dds2, contrast = c("group", "SU_KO", "SU_WT"))
DEG_DESeq2 <- as.data.frame(tmp)
DEG_DESeq2 <- na.omit(DEG_DESeq2)
FC <- 1  # 设置 fold change 阈值
padj <- 0.05  # 设置 p-value 阈值
DEG_DESeq2$Significant <- "no"
up <- intersect(which(DEG_DESeq2$log2FoldChange > FC), which(DEG_DESeq2$padj < padj))
down <- intersect(which(DEG_DESeq2$log2FoldChange < -FC), which(DEG_DESeq2$padj < padj))
DEG_DESeq2$Significant[up] <- "up"
DEG_DESeq2$Significant[down] <- "down"
table(DEG_DESeq2$Significant)

# 保存结果
write.csv(DEG_DESeq2, file = "SUWT_vs_SUKO_DEGs.csv")








# 读取数据
exprSet <- read.csv("counts.csv", header = TRUE)
exprSet <- data.frame(exprSet)
rownames(exprSet) <- exprSet$X
exprSet <- exprSet[which(!(rownames(exprSet) %in% "-")), ]
exprSet <- exprSet[, -1]

# 将非整数值转换为整数
exprSet <- round(exprSet)

# 检查转换后的数据
head(exprSet)


# 将负数和NA替换为0
exprSet[is.na(exprSet) | exprSet < 0] <- 0

# 再次检查数据
head(exprSet)

# 创建样本信息数据框
Data <- data.frame(row.names = colnames(exprSet),
                   group = c("SU_WT", "SU_WT", "SU_WT", "SU_KO", "SU_KO", "SU_KO"))

# 查看样本信息数据框
head(Data)

# 创建 DESeq2 数据集
dds <- DESeqDataSetFromMatrix(countData = exprSet, colData = Data, design = ~ group)

# 查看创建好的对象
dds


# 进行差异表达分析
dds2 <- DESeq(dds)
tmp <- results(dds2, contrast = c("group", "SU_KO", "SU_WT"))
DEG_DESeq2 <- as.data.frame(tmp)
DEG_DESeq2 <- na.omit(DEG_DESeq2)

# 设置 fold change 和 p-value 阈值
FC <- 1  # 设置 fold change 阈值
padj <- 0.01  # 设置 p-value 阈值

# 标记显著差异基因
DEG_DESeq2$Significant <- "no"
up <- intersect(which(DEG_DESeq2$log2FoldChange > FC), which(DEG_DESeq2$padj < padj))
down <- intersect(which(DEG_DESeq2$log2FoldChange < -FC), which(DEG_DESeq2$padj < padj))
DEG_DESeq2$Significant[up] <- "up"
DEG_DESeq2$Significant[down] <- "down"

# 查看显著差异基因的分布
table(DEG_DESeq2$Significant)

# 保存结果
write.csv(DEG_DESeq2, file = "SUWT_vs_SUKO_DEGs.csv")




library(ggplot2)

# 设置 fold change 和 p-value 阈值
FC <- 1  # 设置 fold change 阈值
padj <- 0.01  # 设置 p-value 阈值

# 计算 -log10(padj)
DEG_DESeq2$neg_log10_padj <- -log10(DEG_DESeq2$padj)


# 筛选显著上调和下调的基因
up_genes <- DEG_DESeq2[DEG_DESeq2$Significant == "up" & DEG_DESeq2$log2FoldChange > FC & DEG_DESeq2$padj < padj, ]
down_genes <- DEG_DESeq2[DEG_DESeq2$Significant == "down" & DEG_DESeq2$log2FoldChange < -FC & DEG_DESeq2$padj < padj, ]

# 获取前十个显著上调和下调的基因
top_up_genes <- head(up_genes, 10)
top_down_genes <- head(down_genes, 10)



# 绘制火山图
p <- ggplot(DEG_DESeq2, aes(x = log2FoldChange, y = neg_log10_padj, color = Significant)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = c(-FC, FC), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(padj), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot of DEGs (SU_WT vs SU_KO)",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted P-value)") +
  scale_color_manual(values = c("no" = "gray", "up" = "red", "down" = "blue")) +
  theme_minimal() +
  theme(legend.position = "right")

# 添加显著上调基因的标注
p <- p + geom_text(data = top_up_genes, aes(label = rownames(top_up_genes)), vjust = -0.5, hjust = 0.5, color = "red", size = 3)

# 添加显著下调基因的标注
p <- p + geom_text(data = top_down_genes, aes(label = rownames(top_down_genes)), vjust = 0.5, hjust = 0.5, color = "blue", size = 3)

# 显示图表
print(p)


# 定义要标注的基因
genes_to_label <- c("Cyp1a1", "Cyp1b1")

# 筛选出要标注的基因
genes_to_label_data <- DEG_DESeq2[rownames(DEG_DESeq2) %in% genes_to_label, ]

# 添加特定基因的标注
p <- p + geom_text(data = genes_to_label_data, aes(label = rownames(genes_to_label_data)), vjust = -0.5, hjust = 0.5, 

# 筛选出显著下调的基因
down_genes <- DEG_DESeq2[DEG_DESeq2$Significant == "down" & DEG_DESeq2$log2FoldChange < -FC & DEG_DESeq2$padj < padj, ]

# 添加显著下调基因的标注
p <- p + geom_text(data = down_genes, aes(label = rownames(down_genes)), vjust = -0.5, hjust = 0.5, color = "blue", size = 2.5, check_overlap = TRUE)
                   
p                  
                   
# 保存火山图
ggsave("volcano_plot_with_labels.png", width = 10, height = 8, dpi = 300)














# 筛选出显著下调的基因
down_genes <- DEG_DESeq2[DEG_DESeq2$Significant == "down" & DEG_DESeq2$log2FoldChange < -FC & DEG_DESeq2$padj < padj, ]

# 找到左上角最显著的点
top_down_gene <- down_genes[which.max(down_genes$neg_log10_padj), ]
gene_name <- rownames(top_down_gene)


# 绘制火山图
p <- ggplot(DEG_DESeq2, aes(x = log2FoldChange, y = neg_log10_padj, color = Significant)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = c(-FC, FC), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(padj), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot of DEGs (SU_WT vs SU_KO)",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted P-value)") +
  scale_color_manual(values = c("no" = "gray", "up" = "red", "down" = "blue")) +
  theme_minimal() +
  theme(legend.position = "right")

# 添加显著下调基因的标注
p <- p + geom_text(data = down_genes, aes(label = rownames(down_genes)), vjust = -0.5, hjust = 0.5, color = "blue", size = 2.5, check_overlap = TRUE)

# 添加左上角最显著基因的标注
p <- p + geom_text(data = top_down_gene, aes(x = log2FoldChange, y = neg_log10_padj, label = gene_name), vjust = -1, hjust = -0.5, color = "blue", size = 3, fontface = "bold")

# 显示图表
print(p)

# 保存火山图
ggsave("volcano_plot_with_top_down_label.png", width = 10, height = 8, dpi = 300)
