# RNA-seq worrkflow in R

# 读入数据
setwd("~/project/Tet_qzb/featurecounts")
data <- read.table("AllSampleCounts.txt", header=TRUE, quote="\t", skip=1)
names(data)[7:16] <- c("r1","Tet1-2","Tet1-4","Tet2-1","Tet2-8","Tet2-9","Tet3-13","Tet3-2","Tet3-6","Tet3-7")
countData <- as.matrix(data[7:16])
rownames(countData) <- data$Geneid

# 构建dds对象
library('DESeq2')
phenodata <- read.table('phenodata.txt') # ‘phenodata.txt’是一个对实验样本的描述性文件
dds <- DESeqDataSetFromMatrix(countData, colData=phenodata, design= ~ treatment)

#执行DESeq
dds <- DESeq(dds)
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=TRUE)
resultsNames(dds)
resTet3KO<-results(dds, contrast=c("treatment", "Tet3KO", "control"))  # 把实验组分别和control比较
resTet2KO<-results(dds, contrast=c("treatment", "Tet2KO", "control"))
resTet1KO<-results(dds, contrast=c("treatment", "Tet1KO", "control"))

resdata <- na.omit(resdata)
resTet3KO <- na.omit(resTet3KO)

Tet3KO <- as.data.frame(resTet3KO)

# 筛选差异基因
diff_gene_deseq2 <- resdata[which(abs(resdata$log2FoldChange) >2 & resdata$padj < 0.05),]
diffgene_Tet3KO <- Tet3KO[which(abs(Tet3KO$log2FoldChange) >2 & Tet3KO$padj < 0.05),]

# 将结果中的差异基因找到并标记‘up’ or 'down'
resdata$threshold = factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) >= 1, ifelse(resdata$log2FoldChange>= 1 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
Tet3KO$threshold = factor(ifelse(Tet3KO$padj < 0.05 & abs(Tet3KO$log2FoldChange) >= 1, ifelse(Tet3KO$log2FoldChange>= 1 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
Tet3KO$gene<-rownames(Tet3KO)

#绘制火山图
ggplot(Tet3KO,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  geom_text_repel(
    data = Tet3KO[Tet3KO$padj<0.05&abs(Tet3KO$log2FoldChange)>1,],
    aes(label = gene),
    size = 3,
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05