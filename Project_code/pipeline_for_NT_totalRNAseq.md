# code used in NT project

## total RNA-seq 20211109totalRNA-seqdata

### morula part

#### Procession in shell

```sh
# 1.Fastqc
wd=~/project/20211109totalRNA-seqdata/morula/rawdata
cd $wd
mkdir fastqc
for i in *fastq.gz ;do
nohup fastqc $i -o ./fastqc &
done

multiqc -n "morula_rna_raw" ./ &


# 2.Cutadapt
wd=~/project/20211109totalRNA-seqdata/morula/rawdata
nwd=~/project/20211109totalRNA-seqdata/morula/cutdata
mkdir $nwd
cd $wd
for i in *R1.fastq.gz
do
t=${i/R1/R2}
echo $i
echo $t
nohup cutadapt -j 3 -a AGATCGGAAGAGC  -A AGATCGGAAGAGC --trim-n -m 80 -q 20,20 -o $nwd/${i%.fastq*}.cut.fq -p $nwd/${t%.fastq*}.cut.fq $i $t > $nwd/${i%%.*}.cut.log &
done

wd=~/project/20211109totalRNA-seqdata/morula/cutdata
cd $wd
mkdir fastqc
for i in *fq ;do nohup fastqc $i -o ./fastqc & done


# 3.1.align
wd=~/project/20211109totalRNA-seqdata/morula/cutdata
nwd=~/project/20211109totalRNA-seqdata/morula/align
mkdir $nwd
cd $wd
for i in *R1.cut.fq
do
echo $i
t=${i/R1.cut.fq/R2.cut.fq}
echo $t
nohup hisat2 -p 8 --dta-cufflinks --no-discordant  -t -x /home1/share/hisat2_index/mm10/genome -1 $i -2 $t -S $nwd/${i%%.*}.sam  > $nwd/${i%%.*}.log &
echo ${i%%.*}.sam
done


# 3.2.SNPsplit
# 3.2.1 align
wd=~/project/20211109totalRNA-seqdata/morula/cutdata
nwd=~/project/20211109totalRNA-seqdata/morula/SNPsplit/aligned
cd $wd
for i in ZK*R1_cut.fq
do t=${i/R1_cut.fq/R2_cut.fq}
y=${i/R1.cut.fq/pwk-Nmasked.sam}
echo $i,$t $y
nohup bowtie2 -p 8 -x ~/ref/snp_mouse/PWK/PWK_PhJ_N-masked/bowtie2_index/GRCm38_PWK_PhJ_N-masked -1 $i -2 $t -S $y &
done
mkdir $nwd
mv *sam $nwd

# 3.2.2.sam2bam (sort by name)
wd=~/project/20211109totalRNA-seqdata/morula/SNPsplit/aligned
cd $wd
for i in *sam; do t=${i/sam/sortbyname.bam}; nohup samtools sort -n -@ 8 $i -o $t & done

# 3.2.3.snpsplit
wd=~/project/20211109totalRNA-seqdata/morula/SNPsplit/aligned
nwd=~/project/20211109totalRNA-seqdata/morula/SNPsplit/snpsplit
cd $wd
mkdir $nwd
for i in *bam; do
  nohup SNPsplit --snp_file ~/ref/snp_mouse/PWK/all_SNPs_PWK_PhJ_GRCm38.txt.gz --paired -o ../snpsplit/ --no_sort $i &
done

# 3.2.4.featureCounts
wd=~/project/20211109totalRNA-seqdata/morula/SNPsplit/snpsplit
nwd=~/project/20211109totalRNA-seqdata/morula/SNPsplit/snpsplit/featurecounts
cd $wd
featureCounts -T 8 -p -t exon -g gene_id -a /home1/share/gtf/mm10.gtf -o featurecounts/SNP_morula_featureCounts.txt *.genome[12].bam


# 4.flagstat
wd=~/project/20211109totalRNA-seqdata/morula/align
cd $wd
mkdir flagstat
for i in *.sam
do
echo $i
nohup samtools flagstat -@ 1 $i > ./flagstat/${i/.sam/.flagstat} &
done


# 5.sam2bam
wd=~/project/20211109totalRNA-seqdata/morula/align
nwd=~/project/20211109totalRNA-seqdata/morula/proper
mkdir $nwd
cd $wd
for i in *sam
do
nohup samtools sort -@ 2 -o ${i%.*}.sorted.bam $i &
done

for i in *sorted.bam
do
nohup samtools view -bf 0x2 -q 20 -@ 4 $i -o $nwd/${i%sort*}proper.bam  &
done


# 6.1.stringtie
wd=~/project/20211109totalRNA-seqdata/morula/proper
nwd=~/project/20211109totalRNA-seqdata/morula/stringtie
cd $wd
for i in *proper.bam; do
   outfolder=${i%.proper*}
   mkdir -p $nwd/$outfolder
   nohup stringtie -p 8 $i -b $nwd/$outfolder -e -G /home1/share/gtf/mm10.gtf -o $nwd/$outfolder/$outfolder_stringtie.gtf > $nwd/$outfolder.log &
done


# 6.2.featurecounts
wd=~/project/20211109totalRNA-seqdata/morula/proper
nwd=~~/project/20211109totalRNA-seqdata/morula/featureCounts
mkdir $nwd
cd $wd
nohup featureCounts -T 8 -p -t exon -g gene_id -a /home1/share/gtf/mm10.gtf -o $nwd/morulaCounts.txt *.bam &


#6.3.TEcount
wd=~/project/20211109totalRNA-seqdata/morula/align
nwd=~/project/20211109totalRNA-seqdata/morula/TEcount
cd $wd
for i in *sorted.bam
do echo $i
t=${i/sorted.bam/TE}
echo $t
TEcount --sortByPos --format BAM --mode multi -b $i --GTF /home1/share/gtf/mm10.gtf --TE ~/ref/gtf/mm10_rmsk_TE.gtf --project $t &
done
mkdir $nwd
mv *.cntTable $nwd

#6.4.SQuIRE
# fill the part after learn and try

#7.RRR region
  #7.0.RRR region gtf file prepare
  #bedToGenePred whole_rrr_class_mm10_lftover_tran.bed /dev/stdout | genePredToGtf file /dev/stdin whole_rrr_class_mm10_lftover_tran.gtf
  # awk -F "\t" '{OFS="\t"} $7="?"' whole_rrr_class_mm10_lftover_tran.gtf > whole_rrr_class_mm10_lftover_tran.gtf

wd=~/project/20211109totalRNA-seqdata/2C/proper
nohup featureCounts -T 8 -F GTF -p -t CDS -g gene_id -a ~/project/20211109totalRNA-seqdata/2C/ref_publised_data/RRR_region/whole_rrr_class_mm10_lftover_tran.gtf -o ~/project/20211109totalRNA-seqdata/2C/RRR/2c_RRRregion_counts.txt *bam &
```

#### Procession in R

```R
# Ballgown for FPKM and normalized FPKM
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

setwd("~/project/20211109totalRNA-seqdata/morula/stringtie/")
bgmorula=ballgown(dataDir='IFN',samplePattern='wmz',meas='FPKM')
fpkm_morula=gexpr(bgmorula)

tpm_morula <- t(t(fpkm_morula)/colSums(fpkm_morula))*10^6
log2tpm_morula <- log2(tpm_morula+1)


# PCA
library(ggpubr)
library(gmodels)
library(ggpubr)
library(ggplot2)
library(ggthemes)

tlog2tpm_morula <- t(log2tpm_morula)
tlog2tpm_morula <- tlog2tpm_morula[,which(colSums(tlog2tpm_morula)>0)]
pca.morula <- prcomp(tlog2tpm_morula, scale. = TRUE)
summary(pca.morula)
pca.morula.data <- data.frame(sample = rownames(pca.morula$x),Type = c(rep('IVF',2),rep('Farrerol(50nM)',3),rep('NT',3)),pca.morula$x)
ggscatter(pca.morula.data,x = "PC1",y = "PC2",color = "Type", xlab = "PC1 (32.15%)", ylab = "PC2 (19.87%)") + theme_base()


# Ballgown for different gene expression
# which not used in this project
# this part will be commented
  # pData(bg) = data.frame(id = sampleNames(bg), treatment = c(rep('IVF',2),rep('Farrerol(50nM)',3),rep('NT',3)))
  # bg_IVFvsNT = subset(bg, "treatment != 'Farrerol(50nM)' ", genomesubset = FALSE)
  # ff = gexpr(bg_IVFvsNT)
  # stat_results = stattest(bg_IVFvsNT, feature='transcript', meas='FPKM', covariate='treatment', getFC = TRUE)
  # stat_results$log2fc = log2(stat_results$fc)
  # diff_gene = stat_results[which(stat_results$qval<0.05 & abs(stat_results$log2fc)>1.5),]

  # code for plot from wangmt
  # ggscatter(scatterdata,x="log_ctrl",y="log_ko",color="change",palette=c("blue","grey","red"),
  #         xlab="log2(FPKM+1) in ctrl AF",
  #         ylab="log2(FPKM+1) in ko AF",
  #         ggtheme = theme_bw())+
  # geom_abline(linetype="dashed",intercept = 1,alpha=0.8,size=0.8,color = "black")+
  # geom_abline(linetype="dashed",intercept = -1,alpha=0.8,size=0.8,color = "black")+
  # xlim(-1,12) + ylim(-1,12)+coord_cartesian(expand = F)+theme(aspect.ratio=1)+
  # ggtitle(this_tile)+
  # theme(panel.grid = element_blank(),
  #       panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
  #       axis.text=element_text(size=12,face="bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))


# corelation of samples based on FPKM
library(pheatmap)

corelation_morula = cor(log2(fpkm_morula+1))
colnames(corelation_morula)=c(rep('IVF',2),rep('Farrerol(50nM)',3),rep('NT',3))
rownames(corelation_morula)=c(rep('IVF',2),rep('Farrerol(50nM)',3),rep('NT',3))
pheatmap(corelation_morula, cellwidth = 37, cellheight = 35, color = colorRampPalette(colors = c("blue","white","red"))(100))


# heatmap of union diffgene (IVFvsNT & FarrvsNT)
diff_list = intersect(rownames(log2tpm_morula),diff_list)
M = t(scale(t(log2tpm_morula[diff_list,]+1)))
M = na.omit(M)
colnames(M)=c(rep('IVF',2),rep('Farrerol(50nM)',3),rep('NT',3))

MNT = apply(M[,6:8], 1, function(x) max(x)-min(x))
MIVF = apply(M[,1:2], 1, function(x) max(x)-min(x))
MFarr = apply(M[,3:5], 1, function(x) max(x)-min(x))
M = cbind(M,MIVF,MFarr,MNT)
M = M[which(M[,9]<1.5),]
M = M[which(M[,10]<1.5),]
M = M[which(M[,11]<1.5),]
  # the part above removes gene which differ too much within the group


# pheatmap(M, method = 'average', show_rownames = FALSE, cellwidth = 50, cellheight = 0.44, color = colorRampPalette(colors = c("blue","white","red"))(100))
dd = pheatmap(M[,1:8], method = 'average', show_rownames = FALSE, cutree_rows = 8, color = colorRampPalette(colors = c("blue","white","red"))(100), cluster_cols = FALSE)
row_cluster = cutree(dd$tree_row, k=8)
newOrder = M[dd$tree_row$order,1:8]
newOrder = cbind(newOrder, row_cluster[match(rownames(newOrder),names(row_cluster))])
colnames(newOrder)[9]='cluster'
c1 = newOrder[newOrder[,"cluster"]==1,]
c2 = newOrder[newOrder[,"cluster"]==2,]
c3 = newOrder[newOrder[,"cluster"]==3,]
c4 = newOrder[newOrder[,"cluster"]==4,]
c5 = newOrder[newOrder[,"cluster"]==5,]
c6 = newOrder[newOrder[,"cluster"]==6,]
c7 = newOrder[newOrder[,"cluster"]==7,]
c8 = newOrder[newOrder[,"cluster"]==8,]
heat = rbind(c7,c8,c1,c6,c2,c4,c5,c3)
pheatmap(heat[,1:8], cluster_row = FALSE, cluster_col = FALSE, show_rownames = FALSE, cellwidth = 50, cellheight = 0.8, gaps_row = c(34,508,597,634,665), color = colorRampPalette(colors = c("blue","white","red"))(100))
  # the par above reorder the clusters and draw the heatmap in the proper order

write.csv(heat[1:34,1:8], file = 'morula_heatmap_part1.csv')
write.csv(heat[35:508,1:8], file = 'morula_heatmap_part2.csv')
write.csv(heat[509:597,1:8], file = 'morula_heatmap_part3.csv')
write.csv(heat[598:634,1:8], file = 'morula_heatmap_part4.csv')
write.csv(heat[635:665,1:8], file = 'morula_heatmap_part5.csv')
write.csv(heat[666:744,1:8], file = 'morula_heatmap_part6.csv')

###############################下次从这开始

ICMTEheat = t(scale(t(log2tpm_morula[ICMTE$ICMTEgene,])))
ICMTEheat = na.omit(ICMTEheat)
colnames(ICMTEheat)=c(rep('IVF',2),rep('Farrerol(50nM)',3),rep('NT',3))
pheatmap(ICMTEheat, method = 'average', show_rownames = FALSE)

DNAdamage <- read.table('DNAdamage.txt',header = TRUE)
DNAdamageheat = t(scale(t(log2tpm_morula[DNAdamage$DNAdamage,])))
DNAdamageheat = na.omit(DNAdamageheat)
DNAdamageheat = unique(DNAdamageheat)
colnames(DNAdamageheat)=c(rep('IVF',2),rep('Farrerol(50nM)',3),rep('NT',3))
pheatmap(DNAdamageheat, clustering_method = 'ward.D2', show_rownames = FALSE)

NHEJ <- read.table('~/project/20211109totalRNA-seqdata/morula/ref_publised_data/NHEJ.txt', header = TRUE)
NHEJ <- unique(NHEJ)
NHEJheat = t(scale(t(log2tpm_morula[NHEJ$Symbol,])))
NHEJheat = na.omit(NHEJheat)
colnames(NHEJheat)=c(rep('IVF',2),rep('Farrerol(50nM)',3),rep('NT',3))
n = pheatmap(NHEJheat, clustering_method = 'ward', cluster_cols = FALSE, show_rownames = FALSE)


HR <- read.table('HR.txt', header = TRUE)
HR <- unique(HR)
HRheat = t(scale(t(log2tpm_morula[HR$Symbol,])))
HRheat = na.omit(HRheat)
colnames(HRheat)=c(rep('IVF',2),rep('Farrerol(50nM)',3),rep('NT',3))
n = pheatmap(HRheat, clustering_method = 'ward', cluster_cols = FALSE, show_rownames = FALSE)


# 绘制散点图

fpkm_mean = data.frame(IVF = apply(fpkm[,1:2], 1, mean), Farr = apply(fpkm[,3:5], 1, mean), NT = apply(fpkm[,6:8], 1, mean))
fpkm_mean = log2(fpkm_mean+1)


fpkm_mean$IVFvsNT = diffTE_IVF_VS_NT$threshold[match(rownames(fpkm_mean),diffTE_IVF_VS_NT$gene)]
fpkm_mean$IVFvsNT[is.na(fpkm_mean$IVFvsNT)] = 'NoSignifi'

ggplot(fpkm_mean,aes(x=NT,y=IVF,color = IVFvsNT))+
  scale_color_manual(values=c(c("red","blue","grey"))) +
  geom_point() +
  theme_bw()


fpkm_mean$FarrvsNT = diffTE_Farr_VS_NT$threshold[match(rownames(fpkm_mean),diffTE_Farr_VS_NT$gene)]
fpkm_mean$FarrvsNT[is.na(fpkm_mean$FarrvsNT)] = 'NoSignifi'

ggplot(fpkm_mean,aes(x=NT,y=Farr,color = FarrvsNT))+
  scale_color_manual(values=c(c("red","blue","grey"))) +
  geom_point() +
  theme_bw()
```
