library(clusterProfiler)  # GSEA 和富集分析主力包
library(org.Hs.eg.db)     # 人类注释数据库（ENTREZID 与 SYMBOL 等 ID 转换）
library(enrichplot)       # 用于 GSEA 排名图等可视化
library(ggplot2)          # 用于气泡图
library(dplyr)            # 数据处理辅助

setwd('C:/Users/Lamarck/Desktop')

# Step1 读取数据

# CSV文件包含两列ENTREZID 和 log2FoldChange
deg <- read.csv("genes_ENSEMBL_ENTREZID.csv")
head(deg)

# Step2 构建geneList

# geneList是GSEA分析的核心输入格式是 names(geneList) = ENTREZID
geneList <- deg$log2FoldChange
names(geneList) <- as.character(deg$ENTREZID)

# 按照log2FoldChange从大到小排序（GSEA 要求降序）
geneList <- sort(geneList, decreasing = TRUE)

# 查看排序后的 geneList
head(geneList)

# Step3 运行GSEA分析

# KEGG 富集分析（gseKEGG）
kegg_gsea <- gseKEGG(
  geneList = geneList,
  organism = "hsa",               # hsa 表示人类（human）
  minGSSize = 10,                 # 最小基因集大小
  maxGSSize = 500,                # 最大基因集大小
  pvalueCutoff = 0.05,            # p 值阈值
  pAdjustMethod = "BH",           # 多重检验方法（Benjamini-Hochberg）
  verbose = FALSE
)

# 将 KEGG 分析结果转为可读（基因名）形式
kegg_gsea <- setReadable(kegg_gsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# GO 富集分析（gseGO）
go_gsea <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",                    # 可选："BP"、"MF"、"CC"、"ALL"
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# 转换 GO 结果为可读基因名
go_gsea <- setReadable(go_gsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Step4 保存富集分析结果

# 提取结果数据框
kegg_result <- kegg_gsea@result
go_result <- go_gsea@result

# 写入 CSV 文件保存
write.csv(kegg_result, file = "GSEA_result_KEGG_human.csv", row.names = FALSE)
write.csv(go_result, file = "GSEA_result_GO_human.csv", row.names = FALSE)

# Step5 可视化

# ---  GSEA 排名图（enrichment score curve） ---
# 画出 KEGG 分析中第一个通路的 GSEA 曲线图
gseaplot2(kegg_gsea, geneSetID = 1, 
          pvalue_table = TRUE,
          title = kegg_result$Description[1])

# 画出 GO 分析中第一个通路的 GSEA 曲线图
gseaplot2(go_gsea, geneSetID = 1, 
          pvalue_table = TRUE,
          title = go_result$Description[1])
