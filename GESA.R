library(clusterProfiler)  # GSEA 和富集分析主力包
library(org.Hs.eg.db)     # 人类注释数据库（ENTREZID 与 SYMBOL 等 ID 转换）
library(enrichplot)       # 用于 GSEA 排名图等可视化
library(ggplot2)          # 用于气泡图
library(dplyr)            # 数据处理辅助

# Step0 工作路径设置 #
setwd("C:/Users/Lamarck/Desktop")

# Step1 读取数据 (必须含ENTREZID和log2FoldChange两列)
deg <- read.csv("genes_ENSEMBL_ENTREZID.csv")
head(deg)

# Step2 构建 geneList
geneList <- deg$log2FoldChange
names(geneList) <- as.character(deg$ENTREZID)
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)

# Step3 运行 GSEA 分析

## KEGG 富集分析（gseKEGG）
kegg_gsea <- gseKEGG(
  geneList      = geneList,
  organism      = "hsa",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)
kegg_gsea <- setReadable(kegg_gsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

## GO 富集分析（gseGO）
go_gsea <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)
go_gsea <- setReadable(go_gsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Step4 保存富集分析结果
kegg_result <- kegg_gsea@result
go_result   <- go_gsea@result

write.csv(kegg_result, file = "GSEA_result_KEGG_human.csv", row.names = FALSE)
write.csv(go_result, file = "GSEA_result_GO_human.csv",   row.names = FALSE)

# Step5 可视化输出

# 画KEGG分析第一个通路的GSEA曲线图并保存为PDF
pdf(file = "GSEA_KEGG_Enrichment_Score_Curve.pdf", width = 8, height = 6)
gseaplot2(
  kegg_gsea, 
  geneSetID    = 5,
  pvalue_table = TRUE,
  title        = kegg_result$Description[1]
)
dev.off()

# 画GO分析第一个通路的GSEA曲线图并保存为PDF
pdf(file = "GSEA_GO_Enrichment_Score_Curve.pdf", width = 8, height = 6)
gseaplot2(
  go_gsea, 
  geneSetID    = 20,
  pvalue_table = TRUE,
  title        = go_result$Description[1]
)
dev.off()
