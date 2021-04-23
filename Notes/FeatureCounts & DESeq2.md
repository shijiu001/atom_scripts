# FeatureCounts & DESeq2

DEseq2的输入中需要一个计数矩阵，我本以为之前FeatureCounts的结果可以直接进行输入，后来发现之前的结果是每个样本输入了一个文件，需要先将文件进行合并，虽然这个操作并不困难，但我还是用FeatureCounts对全部样本重新进行了计数。

`featureCounts -T 8 -p -t exon -g gene_id -a ../chrX_data/genes/chrX.gtf -o AllSampleCounts.txt *.bam`

这样得到的结果中，第七列及之后会分别记录各个**样本的名称**及计数结果。

将结果读入R studio，忽略第一行

`data <- read.table("all_feature.txt", header=TRUE, quote="\t", skip=1)`

将pheno_data的第一列，即样本名称写入一个新的列表

`SampleNames <- pheno_data[,1]`

>这里不能用pheno_data[1]直接赋值，否则SampleNames仍是个表格，或许表格也可以进行下去（修改data中的列名），但我不知道如何实现。

得到的SmapleNames的数据为factor类，将其转换为char类，并将data中7-18列的列名改为样本标签

`SampleNames <- as.character(SampleNames)`

`names(data)[7:18]<- SampleNames`

>这里如果不进行数据类型转换，改变列名之后就是“1”,“2”, ... ,"12"，我也不清楚为什么

从data中提取出计数结果countData

`countData <- as.matrix(data[7:18])`
`rownames(countData) <- data$Geneid`

>写到这里我突然会了如何从列表中直接提取样本标签用来修改列名

构建dds矩阵，这里我用sex作为因变量

`dds <- DESeqDataSetFromMatrix(countData, colData=pheno_data, design= ~ sex)`

DEseq

```
dds <- DESeq(dds)
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
```

>DEseq的过程中，会进行一次标准化，我觉得这个应该就是算一下FPKM，但是我比较了DESeq的结果和之前的FPKM，发现两者的区别还是挺大的。

另，我试着画了个火山图来看看差异，结果只有一个基因显示出了差异/xjj