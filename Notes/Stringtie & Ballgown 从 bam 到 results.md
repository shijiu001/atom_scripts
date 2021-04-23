# Stringtie 从 bam 到 gtf 
# Ballgown 从 gtf 到 results

### 组装成transcripts，并产生Ballgown的输入

获得 bam 文件后，利用 stringtie 将其中的测序片段组装成transcripts：
>这里需要提供一个参考gtf文件

`for i in *.bam; do i=${i%.bam*}; stringtie -p 8 -G ../chrX_data/genes/chrX.gtf -o ${i}.gtf -l ${i} ${i}.bam; done`

利用 stringtie --merge 将共12个样本中的transcripts合并为一个跨全部样本的转录本参考：
>这里同样需要提供一个参考gtf文件（但我觉得这是不必要的，我认为此处简单合并就可以了，并不明白这个参考gtf文件有什么用）

`stringtie --merge -p 8 -G ../chrX_data/genes/chrX.gtf -o stringtie_merged.gtf ../chrX_data/mergelist.txt`

以这个合并的 gtf 文件为参考，重新组装transcripts，并作为 Ballgown 的输入：

`for i in *.bam;do i=${i%_chrX.bam*}; stringtie -e -B -p 8 -G ../mid_gtf/stringtie_merged.gtf -o ../ballgown/${i}/${i}_chrX.gtf ${i}_chrX.bam; done`

### 读入数据，计算差异，绘制图片

接下来的操作在R studio中进行

配置环境：

```
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
```

读入样本的基本信息：

`pheno_data = read.csv("geuvadis_phenodata.csv")`

读入前文中产生的转录本数据：

`bg_chrX = ballgown(dataDir = "ballgown", samplePattern = "ERR", pData=pheno_data)`
>此处我并不理解 ballgown 读入了怎样的数据，我以为它最重要的输入就是 gtf 文件，而实际上它的输入是一个类似dataframe的东西，而且我还并不很清楚它的结构。

过滤低丰度的基因：

`bg_chrX_filt = subset(bg_chrX,"rowVars(texpr(bg_chrX)) >1",genomesubset=TRUE)`

>为什么过滤掉样本间差异少于1的转录本可以被认为是在过滤低丰度的基因？还是说这只是为了方便寻找差异较大的基因？
>**texpr**这个函数（？）是什么意思？

确认组间有差异的transcripts，此处比较male和female之间的基因差异：

`results_transcripts = stattest(bg_chrX_filt, feature="transcript",covariate="sex",adjustvars = c("population"), getFC=TRUE, meas="FPKM")`

>**getFC**允许结果中显示flodchange

同上一步，不过这里针对gene进行分析：

`results_genes = stattest(bg_chrX_filt, feature="gene", covariate="sex", adjustvars = c("population"), getFC=TRUE, meas="FPKM")`

向针对transcripts分析的结果中加入gene name和gene id：

`results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_chrX_filt), geneIDs=ballgown::geneIDs(bg_chrX_filt), results_transcripts)`

根据p值对结果排序：
```
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)
```
>这里我不太理解q value是什么。

## 关于 FeatureCounts
我利用 FeatureCounts 计算了 bam 文件中的 transcripts，以exon为一个feature，一个gene为一个meta-transcript（或者是一个transcript为一个meta-feature，我更倾向于认为一个transcript是一个meta-feature）

` for i in *.bam; do i=${i%.bam*}; featureCounts -T 8 -p -t exon -g gene_id -a ../chrX_data/genes/chrX.gtf -o ${i}.txt ${i}.bam;done`

我下载到PC，并用excel打开了这个结果文件（.txt文件）我觉得我大致可以理解。

但这个过程中，大概只有50%的测序数据被成功匹配到了参考gtf文件，这个比例是正常的吗？为什么会产生如此多不能成功匹配的数据呢？

## 其他
我试图在markdown中引用服务器中保存的图片，我以为这应该是可行的，但我并没有成功，我想问下这是否可行，如果可行应该如何实现？