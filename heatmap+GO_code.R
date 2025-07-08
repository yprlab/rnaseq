# 载入必要的库
library(dplyr)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tibble)

# 设置工作目录
setwd("/home/dongxiaotong/R/lianxi/zl")

# 载入DEG数据，设置行名为第一列
DEG_DESeq2 <- read.csv("SUWT_vs_SUKO_DEGs.csv", row.names = 1)

# 检查数据结构
print(head(DEG_DESeq2))
print(names(DEG_DESeq2))  # 确保列名正确

# 将行名转换为一个列
DEG_DESeq2 <- DEG_DESeq2 %>%
  rownames_to_column(var = "gene")

# 提取显著差异表达基因
significant_genes <- DEG_DESeq2 %>%
  filter(Significant %in% c("up", "down")) %>%
  dplyr::select(gene, log2FoldChange, pvalue, padj, Significant)

# 检查 significant_genes 是否成功创建
if (nrow(significant_genes) == 0) {
  stop("No significant genes found. Check the filtering criteria.")
}

# 进行GO富集分析
gene_list <- significant_genes$gene

# 确保 clusterProfiler 包已加载
library(clusterProfiler)

# 使用 bitr 函数进行基因标识符转换
eg <- bitr(gene_list, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)

# 获取每个基因对应的GO term
go <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

# 获取每个基因对应的GO term
go_terms <- as.data.frame(go@result)
print(head(go_terms))  # 确保列名正确

# 重命名列
go_terms <- go_terms %>%
  dplyr::select(ID, Description) %>%
  dplyr::rename(GeneID = ID, GO_Description = Description)

# 将基因名转换为ENTREZID
gene_to_entrez <- eg %>%
  dplyr::select(SYMBOL, ENTREZID) %>%
  dplyr::rename(Gene = SYMBOL, GeneID = ENTREZID)

# 合并基因名和GO term
gene_go <- left_join(gene_to_entrez, go_terms, by = "GeneID")

# 确保每个基因只保留一个GO term
gene_go <- gene_go %>%
  group_by(Gene) %>%
  dplyr::slice(1) %>%
  ungroup()

# 将GO term添加到significant_genes中
significant_genes <- left_join(significant_genes, gene_go, by = c("gene" = "Gene"))

# 提取这些基因的表达数据
# 假设 count_data 已经加载
count_data <- read.csv("count1.csv", row.names = 1)
top_gene_expression <- count_data[result, ]

# 确保annotation_row是data.frame类型
annotation_row <- result %>%
  dplyr::select(gene, Significant, GO_Description) %>%
  dplyr::rename(significance = Significant)

# 确保annotation_row的行名与top_gene_expression的行名匹配
rownames(annotation_row) <- result$gene


pdf("pdf19.pdf", width = 6, height = 9)

# 绘制热图，并对每一行进行 z-score 标准化
pheatmap(top_gene_expression, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
       
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = F,
         #show_colnames = TRUE,
         color = colorRampPalette(c("#6238ff", "#fffffF", "#ff220e"))(50),
         main = "Top Differential Genes Heatmap",
         scale = "row")  # 添加这一行进行 z-score 标准化
dev.off()

library(pheatmap)



# 生成热图并保存结果
heatmap_result <-pheatmap(top_gene_expression, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("#6238ff", "#fffffF", "#ff220e"))(50),
         main = "Top Differential Genes Heatmap",
         scale = "row")  # 添加这一行进行 z-score 标准化


# 聚类后的数据矩阵
heatmap_result$tree_row$order  # 行的顺序
heatmap_result$tree_col$order  # 列的顺序

# 重新排列原始数据矩阵以匹配聚类顺序
clustered_data <- top_gene_expression[heatmap_result$tree_row$order, heatmap_result$tree_col$order]
print(clustered_data)
clustered_data$gene <- rownames(clustered_data)




# 提取前50行和后50行的完整数据
top_50_rows <- head(clustered_data, 50)  # 前50行
bottom_50_rows <- tail(clustered_data, 50)  # 后50行

# 如果只需要 gene 列
top_50_gene <- top_50_rows$gene
bottom_50_gene <- bottom_50_rows$gene

# 创建结果数据框
result <- c(
  top_50_gene,
  bottom_50_gene
)

# 查看结果
print(result)

result <- c(result , "Cyp1a1")












# 示例数据框（假设 significant_genes 已经定义）
# significant_genes <- data.frame(log2FoldChange = rnorm(200, mean = 0, sd = 2))

# 提取 log2FoldChange 列
log2FoldChange <- significant_genes$log2FoldChange

# 分别提取正值和负值
positive_values <- log2FoldChange[log2FoldChange > 0]
negative_values <- log2FoldChange[log2FoldChange < 0]

# 提取正值最大的前50个值
top_50_positive_indices <- order(positive_values, decreasing = TRUE)[1:min(50, length(positive_values))]
top_50_positive_rows <- which(log2FoldChange %in% positive_values[top_50_positive_indices])

# 提取负值最小的前50个值（绝对值最大的负值）
top_50_negative_indices <- order(negative_values, decreasing = FALSE)[1:min(50, length(negative_values))]
top_50_negative_rows <- which(log2FoldChange %in% negative_values[top_50_negative_indices])

# 合并结果
top_50_rows <- unique(c(top_50_positive_rows, top_50_negative_rows))

# 提取完整的行数据
top_50_significant_genes <- significant_genes[top_50_rows, ]

# 查看结果
print(top_50_significant_genes)


top_50_significant_genes$gene


write.csv(top_50_significant_genes, file = "top_50_significant_genes.csv", row.names = TRUE)



# 保存 top_gene_expression 数据
write.csv(top_gene_expression, file = "top_gene_expression.csv", row.names = TRUE)







# 筛选上升和下降前十的基因
top_up_genes <- significant_genes %>%
  filter(Significant == "up") %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

top_down_genes <- significant_genes %>%
  filter(Significant == "down") %>%
  arrange(log2FoldChange) %>%
  head(10)

# 合并上升和下降的基因
top_genes <- bind_rows(top_up_genes, top_down_genes)

# 提取这些基因的表达数据
count_data <- read.csv("counts.csv", row.names = 1)
top_gene_expression <- count_data[top_genes$gene, ]

# 检查 top_gene_expression 是否为空
if (nrow(top_gene_expression) == 0) {
  stop("No expression data found for significant genes.")
}

# 保存 top_gene_expression 数据
write.csv(top_gene_expression, file = "top_gene_expression.csv", row.names = TRUE)

# 确保annotation_row是data.frame类型
annotation_row <- top_genes %>%
  dplyr::select(gene, Significant) %>%
  dplyr::rename(significance = Significant)

# 确保annotation_row的行名与top_gene_expression的行名匹配
rownames(annotation_row) <- top_genes$gene
