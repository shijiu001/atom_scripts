# 读入数据
setwd("~/project/20211109totalRNA-seqdata/morula/TEcount")

file_names <- list.files()
TEcounts <- read.table(file = file_names[1], col.names = file_names[1])
for(i in 2:length(file_names)){
  newdata <- read.table(file = file_names[i], col.names = file_names[i])
  if (rownames(TEcounts) == rownames(newdata)){TEcounts <- cbind(TEcounts,newdata)}
}

remove(file_names,i,newdata)

colnames(TEcounts) <- c("IVF-3","IVF-4","IVF-5","NT-1","NT-2","NT-3","Farrerol-1","Farrerol-2","Farrerol-3","SCR7-1","SCR7-2","SCR7-3","IVF-1","IVF-2")

TEcounts <- TEcounts[,c(4:9,13:14)] #提取其中IVF, NT, Farrerol的数据

# 构建dds对象
library('DESeq2')
pheno <- read.table('../proper/pheno-morula.txt',header = TRUE) # ‘phenodata.txt’是一个对实验样本的描述性文件
pheno <- pheno[c(4:9,13:14),]
dds <- DESeqDataSetFromMatrix(TEcounts, colData=pheno, design= ~ treatment)

#执行DESeq
dds <- DESeq(dds)

# 提取标准化结果
# vsd <- vst(dds)
# ex_vsd <- assay(vsd)
# write.csv(ex_vsd, '../result/VST.csv')
#
# rld <- rlog(dds)
# ex_rld <- assay(dds)
#
# rld <- rlogTransformation(dds)
# exprSet_new=assay(rld)

# res <- results(dds)
# resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=TRUE)
# resultsNames(dds)

# PCA
plotPCA(vsd, intgroup = 'treatment')
plotPCA(vsd, intgroup = 'id')
plotPCA(vsd, intgroup = 'treatment', returnData = TRUE)

# pca.info <- prcomp(ex_vsd)
# head(pca.info$rotation)
# pca.data <- data.frame(sample = rownames(pca.info$rotation),Type = c('IVF','IVF','Farrerol(50nM)','Farrerol(50nM)','Farrerol(50nM)','SCR7(10μM)','SCR7(10μM)','SCR7(10μM)','NT','NT','NT'),pca.info$rotation)
# ggscatter(pca.data,x = "PC1",y = "PC2",color = "Type") + theme_base()

# 计算差异基因
resFarr_VS_NT <- results(dds, contrast=c("treatment", "Farrerol(50nM)", "NT"))  # 把实验组分别和control比较
# resSCR_VS_NT <- results(dds, contrast=c("treatment", "SCR7(10μM)", "NT"))
resFarr_VS_IVF <- results(dds, contrast=c("treatment", "Farrerol(50nM)", "IVF"))  # 把实验组分别和control比较
# resSCR_VS_IVF <- results(dds, contrast=c("treatment", "SCR7(10μM)", "IVF"))
resIVF_VS_NT <- results(dds, contrast=c("treatment", "IVF", "NT"))  # 把实验组分别和control比较


resFarr_VS_NT <- na.omit(resFarr_VS_NT)
# resSCR_VS_NT <- na.omit(resSCR_VS_NT)
resFarr_VS_IVF <- na.omit(resFarr_VS_IVF)
# resSCR_VS_IVF <- na.omit(resSCR_VS_IVF)
resIVF_VS_NT <- na.omit(resIVF_VS_NT)


Farr_VS_NT <- as.data.frame(resFarr_VS_NT)
# SCR_VS_NT <- as.data.frame(resSCR_VS_NT)
Farr_VS_IVF <- as.data.frame(resFarr_VS_IVF)
# SCR_VS_IVF <- as.data.frame(resSCR_VS_IVF)
IVF_VS_NT <- as.data.frame(resIVF_VS_NT)

# 筛选差异基因
diffTE_Farr_VS_NT <- Farr_VS_NT[which(abs(Farr_VS_NT$log2FoldChange) >1.3 & Farr_VS_NT$padj < 0.001),]
# diffTE_SCR_VS_NT <- SCR_VS_NT[which(abs(SCR_VS_NT$log2FoldChange) >1.5 & SCR_VS_NT$padj < 0.01),]
diffTE_Farr_VS_IVF <- Farr_VS_IVF[which(abs(Farr_VS_IVF$log2FoldChange) >2 & Farr_VS_IVF$padj < 0.001),]
# diffTE_SCR_VS_IVF <- SCR_VS_IVF[which(abs(SCR_VS_IVF$log2FoldChange) >1.5 & SCR_VS_IVF$padj < 0.01),]
diffTE_IVF_VS_NT <-IVF_VS_NT[which(abs(IVF_VS_NT$log2FoldChange) >2 & IVF_VS_NT$padj < 0.001),]


# write.csv(diffTE_Farr_VS_IVF, '../result/diff_Farr_VS_IVF.csv')
# write.csv(diffTE_Farr_VS_NT, '../result/diff_Farr_VS_NT.csv')
# write.csv(diffTE_SCR_VS_IVF, '../result/diff_SCR_VS_IVF.csv')
# write.csv(diffTE_SCR_VS_NT, '../result/diff_SCR_VS_NT.csv')
# write.csv(diffTE_IVF_VS_NT, '../result/diff_IVF_VS_NT.csv')

# 提取IVF和Farr与NT差异基因的并集
diff_list <- union(rownames(diffTE_Farr_VS_NT),rownames(diffTE_IVF_VS_NT))

# 将结果中的差异基因找到并标记‘up’ or 'down'
Farr_VS_NT$threshold = factor(ifelse(Farr_VS_NT$padj < 0.001 & abs(Farr_VS_NT$log2FoldChange) >= 1.3, ifelse(Farr_VS_NT$log2FoldChange>= 1.3 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
# SCR_VS_NT$threshold = factor(ifelse(SCR_VS_NT$padj < 0.01 & abs(SCR_VS_NT$log2FoldChange) >= 1, ifelse(SCR_VS_NT$log2FoldChange>= 1 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
Farr_VS_IVF$threshold = factor(ifelse(Farr_VS_IVF$padj < 0.001 & abs(Farr_VS_IVF$log2FoldChange) >= 2, ifelse(Farr_VS_IVF$log2FoldChange>= 2 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
# SCR_VS_IVF$threshold = factor(ifelse(SCR_VS_IVF$padj < 0.01 & abs(SCR_VS_IVF$log2FoldChange) >= 2, ifelse(SCR_VS_IVF$log2FoldChange>= 2 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
IVF_VS_NT$threshold = factor(ifelse(IVF_VS_NT$padj < 0.001 & abs(IVF_VS_NT$log2FoldChange) >= 2, ifelse(IVF_VS_NT$log2FoldChange>= 2 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))


#绘制火山图
library(ggplot2)
library(ggrepel)

Farr_VS_NT$gene <- rownames(Farr_VS_NT) # 把行名中的基因名专门赋值给一列，用于之后的调用
# SCR_VS_NT$gene <- rownames(SCR_VS_NT)
Farr_VS_IVF$gene <- rownames(Farr_VS_IVF)
# SCR_VS_IVF$gene <- rownames(SCR_VS_IVF)
IVF_VS_NT$gene <- rownames(IVF_VS_NT)

ggplot(Farr_VS_NT,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
  ggtitle('Farr_VS_NT log2FC=1.3')+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  geom_text_repel(
    data = Farr_VS_NT[Farr_VS_NT$padj<0.001&abs(Farr_VS_NT$log2FoldChange)>1.3,],
    aes(label = gene),
    size = 3,
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_vline(xintercept=c(-1.3,1.3),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.001),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05


# ggplot(SCR_VS_NT,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
#   geom_point()+
#   scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
#   geom_text_repel(
#     data = SCR_VS_NT[SCR_VS_NT$padj<0.01&abs(SCR_VS_NT$log2FoldChange)>1,],
#     aes(label = gene),
#     size = 3,
#     segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
#   theme_bw()+#修改图片背景
#   theme(
#     legend.title = element_blank()#不显示图例标题
#   )+
#   ylab('-log10 (p-adj)')+#修改y轴名称
#   xlab('log2 (FoldChange)')+#修改x轴名称
#   geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
#   geom_hline(yintercept = -log10(0.01),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
#



ggplot(Farr_VS_IVF,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
  ggtitle('Farr_VS_IVF log2FC=2')+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  geom_text_repel(
    data = Farr_VS_IVF[Farr_VS_IVF$padj<0.001&abs(Farr_VS_IVF$log2FoldChange)>2,],
    aes(label = gene),
    size = 3,
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_vline(xintercept=c(-2,2),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.001),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05


ggplot(IVF_VS_NT,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
  ggtitle('IVF_VS_NT log2FC=2')+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  geom_text_repel(
    data = IVF_VS_NT[IVF_VS_NT$padj<0.001&abs(IVF_VS_NT$log2FoldChange)>2,],
    aes(label = gene),
    size = 3,
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_vline(xintercept=c(-2,2),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.001),lty=3,col="black",lwd=0.5)


# ggplot(SCR_VS_IVF,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
#   geom_point()+
#   scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
#   geom_text_repel(
#     data = SCR_VS_IVF[SCR_VS_IVF$padj<0.01&abs(SCR_VS_IVF$log2FoldChange)>2,],
#     aes(label = gene),
#     size = 3,
#     segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
#   theme_bw()+#修改图片背景
#   theme(
#     legend.title = element_blank()#不显示图例标题
#   )+
#   ylab('-log10 (p-adj)')+#修改y轴名称
#   xlab('log2 (FoldChange)')+#修改x轴名称
#   geom_vline(xintercept=c(-2,2),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
#   geom_hline(yintercept = -log10(0.01),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
