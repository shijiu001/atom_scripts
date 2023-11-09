# CUT&TAG pipeline(coded during the first analysis of CUT&TAG data, 2023.5.17-2023.5.28)

## WT and Dux-KO embryos at late 2-cell stage, targeting H3K9me3

### 0. enviroment

```sh
conda activate chip
```

### 1. Fastqc

```sh
wd=~/project/DUX/batch1_CUTandTAG/1.rawdata
cd $wd
mkdir fastqc
for i in *fastq.gz; do
nohup fastqc $i -o ./fastqc &
done

multiqc -n 'batch1_H3K9me3_CUTandTAG_raw' ./fastqc &
```

### 2. cutadapt

```sh
wd=~/project/DUX/batch1_CUTandTAG/1.rawdata
nwd=~/project/DUX/batch1_CUTandTAG/2.cutdata
mkdir $nwd
cd $wd
for i in *R1_*.fastq.gz; do 
t=${i/_R1/_R2}
echo $i
echo $t
nohup cutadapt -j 3 -a AGATCGGAAGAGC  -A AGATCGGAAGAGC --trim-n -m 25 -q 20,20 -o $nwd/${i%_S*}_R1.cut.fq -p $nwd/${t%_S*}_R2.cut.fq $i $t > $nwd/${i%%_*}.cut.log & 
done

cd $nwd
mkdir log
mv *log log/
```

#### 2ex1. rename files

#### 2ex2. fastqc after cutadapt

```sh
wd=~/project/DUX/batch1_CUTandTAG/2.cutdata
cd $wd
mkdir fastqc
for i in *fq ;do nohup fastqc $i -o ./fastqc & done

multiqc -n 'batch1_H3K9me3_CUTandTAG_cutadapt' ./fastqc &
```

### 3. align

#### 3.1 align to mm10

```sh
wd=~/project/DUX/batch1_CUTandTAG/2.cutdata
nwd=~/project/DUX/batch1_CUTandTAG/3.align
mkdir $nwd
cd $wd
for i in *.R1.cut.fq
do
echo $i
t=${i/R1.cut.fq/R2.cut.fq}
echo $t
nohup bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 8 -x ~/ref/bowtie_index/mm10 -1 $i -2 $t -S $nwd/${i%%.R1*}.bowtie2.sam > $nwd/${i%%.R1*}.bowtie2.log &
done
```

#### 3.2 align to spike-in genome for spike-in calibration

```sh
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 8 -x ${spikeInRef} -1 ${projPath}/fastq/${histName}_R1.fastq.gz -2 ${projPath}/fastq/${histName}_R2.fastq.gz -S $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam & > $projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.log

seqDepthDouble=`samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam | wc -l`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth >$projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDepth
```

#### 3.3 alignment QC

```R
## Collect the alignment results from the bowtie2 alignment summary files
setwd("~/project/CUTandTAG_20230424/3.align")

library(dplyr)
library(magrittr)
file_names <- list.files(pattern = 'log')
alignResult = c()
for(file in file_names){
  alignRes = read.table(file, header = FALSE, fill = TRUE)
  alignRate = substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  Info = strsplit(file, c("_","."))[[1]]
  alignResult = data.frame(Genotype = Info[1], Replicate = substr(Info[4],1,4), 
                           SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric, 
                           MappedFragNum_mm10 = alignRes$V1[4] %>% as.character %>% as.numeric + alignRes$V1[5] %>% as.character %>% as.numeric, 
                           AlignmentRate_mm10 = alignRate %>% as.numeric)  %>% rbind(alignResult, .)
}
alignResult$Genotype = as.factor(alignResult$Genotype)
alignResult %>% mutate(AlignmentRate_mm10 = paste0(AlignmentRate_mm10, "%"))

# ##Spike-in
# spikeAlign = c()
# for(hist in sampleList){
#   spikeRes = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2_spikeIn.txt"), header = FALSE, fill = TRUE)
#   alignRate = substr(spikeRes$V1[6], 1, nchar(as.character(spikeRes$V1[6]))-1)
#   histInfo = strsplit(hist, "_")[[1]]
#   spikeAlign = data.frame(Histone = histInfo[1], Replicate = histInfo[2], 
#                           SequencingDepth = spikeRes$V1[1] %>% as.character %>% as.numeric, 
#                           MappedFragNum_spikeIn = spikeRes$V1[4] %>% as.character %>% as.numeric + spikeRes$V1[5] %>% as.character %>% as.numeric, 
#                           AlignmentRate_spikeIn = alignRate %>% as.numeric)  %>% rbind(spikeAlign, .)
# }
# spikeAlign$Histone = factor(spikeAlign$Histone, levels = histList)
# spikeAlign %>% mutate(AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))

# alignSummary = left_join(alignResult, spikeAlign, by = c("Histone", "Replicate", "SequencingDepth")) %>%
#   mutate(AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"), 
#          AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))
# alignSummary

## Generate sequencing depth boxplot
library(ggplot2)
library(viridis)
# library(ggpubr) <which can only be quited for the failure of installing ggpubr>
fig3A = alignResult %>% ggplot(aes(x = Genotype, y = SequencingDepth/1000000, fill = Genotype)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Sequencing Depth per Million") +
  xlab("") + 
  ggtitle("A. Sequencing Depth")

fig3B = alignResult %>% ggplot(aes(x = Genotype, y = MappedFragNum_mm10/1000000, fill = Genotype)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Mapped Fragments per Million") +
  xlab("") +
  ggtitle("B. Alignable Fragment (mm10)")

fig3C = alignResult %>% ggplot(aes(x = Genotype, y = AlignmentRate_mm10, fill = Genotype)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Mapped Fragments") +
  xlab("") +
  ggtitle("C. Alignment Rate (mm10)")

# fig3D = spikeAlign %>% ggplot(aes(x = Histone, y = AlignmentRate_spikeIn, fill = Histone)) +
#     geom_boxplot() +
#     geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
#     scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
#     scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
#     theme_bw(base_size = 18) +
#     ylab("Spike-in Alignment Rate") +
#     xlab("") +
#     ggtitle("D. Alignment Rate (E.coli)")

# ggarrange(fig3A, fig3B, fig3C, fig3D, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom") <which can only be quited for the failure of installing ggpubr>


```

### 4. remove dup

```sh
## Sort by coordinate
for i in *bowtie2.sam; do
t=${i%%..bowtie2.sam}_bowtie2.sorted.sam
echo $i $t
nohup picard SortSam I=$i O=$t SORT_ORDER=coordinate &
done

## mark duplicates
picard MarkDuplicates -I WT_l2C_H3K9me3_rep1_bowtie2.sorted.sam -O WT_l2C_H3K9me3_rep1_bowtie2.DupMarked.sam -METRICS_FILE WT_l2C_H3K9me3_rep1_picard.DupMark.txt

## remove duplicates
picardCMD MarkDuplicates I=$projPath/alignment/sam/${histName}_bowtie2.sorted.sam O=$projPath/alignment/removeDuplicate/${histName}_bowtie2.sorted.rmDup.sam REMOVE_DUPLICATES=true METRICS_FILE=$projPath/alignment/removeDuplicate/picard_summary/${histName}_picard.rmDup.txt

for i in *sorted.sam; do
nohup picard MarkDuplicates -I $i -O ${i/sorted/rmDup} --REMOVE_DUPLICATES true -METRICS_FILE ${i/bowtie2.sorted.sam/picard.rmDup.txt} &
done
```

### 5. mapped fragment size distribution

```sh
wd=~/project/2C_related/Dux/CUTandTAG_20230424/3.align
mkdir fragmentLen

for i in *sam; do
samtools view -F 0x04 $i | awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > ${i%.bowtie2.sam}.fragmentLen.txt
done
```

```R
setwd("~/project/2C_related/Dux/CUTandTAG_20230424/3.align")
file_names <- list.files(pattern = 'fragmentLen.txt')
fragLen = c()
for(file in file_names){
  Info = strsplit(file, c("_","."))[[1]]
  fragLen = read.table(file, header = FALSE) %>% mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)), Genotype = Info[1], Histone = Info[3], Sample = paste(Info[1],gsub('.fragmentLen.txt','',Info[4],),sep = '-')) %>% rbind(fragLen, .)
}
fragLen$Genotype = factor(fragLen$Genotype)
fragLen$Histone = factor(fragLen$Histone)
fragLen$Sample = factor(fragLen$Sample)

fig5A = fragLen %>% ggplot(aes(x = Sample, y = fragLen, weight = Weight, fill = Genotype)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  # ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("")

fig5B = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = Genotype, group = Sample)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

# ggarrange(fig5A, fig5B, ncol = 2)
```

### 6. Filtering mapped reads by the mapping quality filtering

```sh
wd=~/project/2C_related/Dux/CUTandTAG_20230424/4.markdup/rmDup
cd $wd

for i in *sam; do
nohup samtools view -h -q 2 $i > ${i/sam/qualityScore2.sam} &
done
```

### 7. File format conversion

```sh
wd=~/project/2C_related/Dux/CUTandTAG_20230424/4.markdup/rmDup
cd $wd

# Filter and keep the mapped read pairs
for i in *qualityScore2.sam; do 
nohup samtools view -bS -F 0x04 $i > ${i/sam/bam} &
done

# Convert into bed file format
for i in *qualityScore2.bam; do
nohup bedtools bamtobed -i $i -bedpe > ${i%.*}.bed &
done 

## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
for i in *bed; do
awk '$1==$4 && $6-$2 < 5000 {print $0}' $i > ${i/bed/clean.bed};
done

## Only extract the fragment related columns
for i in *clean.bed; do
cut -f 1,2,6 $i | sort -k1,1 -k2,2n -k3,3n  > ${i%.*}.fragments.bed;
done

for i in *fragments.bed; do
nohup bedtools genomecov -bg -i $i -g ~/ref/genome_fasta/mm10/mm10.chrom.size > ${i/bed/bedgraph} &
done
```

### 8. call peak
```sh
for i in *bedgraph; do
bash ~/miniconda3/envs/chip/bin/SEACR_1.3.sh $i 0.05 norm relaxed ${i%%.*}.SEACR.peaks;
done
```