#### RNA-seq样本PCA分析####
#加载的R包
install.packages(c("ggpubr","ggthemes","gmodels"))
library(ggpubr)
#加载差异基因表达矩阵
library(gmodels)
library(ggpubr)
library(ggplot2)
library(ggthemes)
data<-read.csv("C:/Users/Administrator/Desktop/I.csv",header = T,row.names = 1)
head(data)#每一列为一个样本，每一行为一个基因
#计算PCA
pca.info <- fast.prcomp(Tet3KO_PCAdata)
??fast.prcomp()
#显示PCA计算结果
head(pca.info$rotation)
#计算Y1与Y28之间的差异
pca.data <- data.frame(sample = rownames(pca.info$rotation),Type = c(rep("r1",1),rep("Tet3",4)),pca.info$rotation)
#绘图
ggscatter(pca.data,x = "PC1",y = "PC2",color = "Type") + theme_base()
###其它图形修饰参数自己摸索吧，哈哈~