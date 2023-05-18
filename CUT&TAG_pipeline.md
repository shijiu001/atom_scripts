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
fig3A = alignResult %>% ggplot(aes(x = Histone, y = SequencingDepth/1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Sequencing Depth per Million") +
    xlab("") + 
    ggtitle("A. Sequencing Depth")

fig3B = alignResult %>% ggplot(aes(x = Histone, y = MappedFragNum_hg38/1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Mapped Fragments per Million") +
    xlab("") +
    ggtitle("B. Alignable Fragment (hg38)")

fig3C = alignResult %>% ggplot(aes(x = Histone, y = AlignmentRate_hg38, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("% of Mapped Fragments") +
    xlab("") +
    ggtitle("C. Alignment Rate (hg38)")

fig3D = spikeAlign %>% ggplot(aes(x = Histone, y = AlignmentRate_spikeIn, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Spike-in Alignment Rate") +
    xlab("") +
    ggtitle("D. Alignment Rate (E.coli)")

ggarrange(fig3A, fig3B, fig3C, fig3D, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")


```


### 4. remove dup
```sh
## depending on how you load picard and your server environment, the picardCMD can be different. Adjust accordingly.
picardCMD="java -jar picard.jar"
mkdir -p $projPath/alignment/removeDuplicate/picard_summary

## Sort by coordinate
$picardCMD SortSam I=$projPath/alignment/sam/${histName}_bowtie2.sam O=$projPath/alignment/sam/${histName}_bowtie2.sorted.sam SORT_ORDER=coordinate

## mark duplicates
$picardCMD MarkDuplicates I=$projPath/alignment/sam/${histName}_bowtie2.sorted.sam O=$projPath/alignment/removeDuplicate/${histName}_bowtie2.sorted.dupMarked.sam METRICS_FILE=$projPath/alignment/removeDuplicate/picard_summary/${histName}_picard.dupMark.txt

## remove duplicates
picardCMD MarkDuplicates I=$projPath/alignment/sam/${histName}_bowtie2.sorted.sam O=$projPath/alignment/removeDuplicate/${histName}_bowtie2.sorted.rmDup.sam REMOVE_DUPLICATES=true METRICS_FILE=$projPath/alignment/removeDuplicate/picard_summary/${histName}_picard.rmDup.txt
```

### 5. mapped fragment size distribution