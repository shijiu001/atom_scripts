#### RNA-seq样本GO富集分析 ####

# 获取差异基因列表
genename <- row.names(diffgene_Tet3KO)


# genename 从symbol转换到ID
library(org.Mm.eg.db)
geneID <-  mapIds(org.Mm.eg.db, genename, 'ENTREZID', 'SYMBOL')


# 富集分析
library(clusterProfiler)

BP.params <- enrichGO(   gene   = geneID,    
                         OrgDb  = org.Mm.eg.db,    
                         ont   = "BP"  ,    
                         pAdjustMethod = "BH",    
                         pvalueCutoff  = 0.05,    
                         qvalueCutoff  = 0.2)    #BP, Biological Process

BP.list <- setReadable(BP.params, org.Mm.eg.db, keyType = "ENTREZID")

# dotplot
dotplot(BP.list, showCategory=30)
barplot(BP.list,showCategory=20,drop=T)
