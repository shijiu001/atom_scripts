# Pipeline for pol2 stacc-seq

## data from XieWei

### 0.data download

```sh
# sra-tools一定要最新版的，不然报错报一堆
nohup prefetch --option-file SRR_Acc_List.txt &
for i in SRR9*; do nohup fasterq-dump --split-3 $i -O fastq & done
```

### ex0.1 data combine

```sh

cat SRR9910789_1.fastq SRR9910790_1.fastq SRR9910791_1.fastq SRR9910792_1.fastq > PN5_Pol2_rep1_R1.fastq
cat SRR9910789_2.fastq SRR9910790_2.fastq SRR9910791_2.fastq SRR9910792_2.fastq > PN5_Pol2_rep1_R2.fastq

cat SRR9910793_1.fastq SRR9910794_1.fastq > PN5_Pol2_rep2_R1.fastq
cat SRR9910793_2.fastq SRR9910794_2.fastq > PN5_Pol2_rep2_R2.fastq

cat SRR9910795_1.fastq SRR9910796_1.fastq SRR9910797_1.fastq > e2C_Pol2_rep1_R1.fastq
cat SRR9910795_2.fastq SRR9910796_2.fastq SRR9910797_2.fastq > e2C_Pol2_rep1_R2.fastq

cat SRR9910798_1.fastq SRR9910799_1.fastq SRR9910800_1.fastq > e2C_Pol2_rep2_R1.fastq
cat SRR9910798_2.fastq SRR9910799_2.fastq SRR9910800_2.fastq > e2C_Pol2_rep2_R2.fastq

cat SRR9910801_1.fastq SRR9910802_1.fastq SRR9910803_1.fastq > l2C_Pol2_rep1_R1.fastq
cat SRR9910801_2.fastq SRR9910802_2.fastq SRR9910803_2.fastq > l2C_Pol2_rep1_R2.fastq

cat SRR9910804_1.fastq SRR9910805_1.fastq SRR9910806_1.fastq > l2C_Pol2_rep2_R1.fastq
cat SRR9910804_2.fastq SRR9910805_2.fastq SRR9910806_2.fastq > l2C_Pol2_rep2_R2.fastq

cat SRR9910807_1.fastq SRR9910808_1.fastq SRR9910809_1.fastq > l2C_Pol2_rep3_R1.fastq
cat SRR9910807_2.fastq SRR9910808_2.fastq SRR9910809_2.fastq > l2C_Pol2_rep3_R2.fastq

```

### 1.fastqc

```sh
wd=~/project/DNAdamage/Mirin/IVF/Pol2_XieWei/1.rawdata
cd $wd
mkdir fastqc
for i in *fastq.gz; do
nohup fastqc $i -o ./fastqc &
done

cd ./fastqc/
multiqc -n 'IVF_mirin_ATAC_raw' ./ &

```

### 2.cutdata

```sh

wd=~/project/DNAdamage/Mirin/IVF/Pol2_XieWei/1.rawdata
cd $wd
nwd=~/project/DNAdamage/Mirin/IVF/Pol2_XieWei/2.cutdata
mkdir $nwd



```

### 2ex1.Fastqc after cutadapt

```sh
wd=~/project/DNAdamage/Mirin/IVF/Pol2_XieWei/2.cutdata
cd $wd
mkdir fastqc
for i in *cut.fq; do
nohup fastqc $i -o ./fastqc &
done

cd ./fastqc/
multiqc -n 'IVF_mirin_ATAC_cut' ./ &
```

### 3.align

```sh
wd=~/project/DNAdamage/Mirin/IVF/Pol2_XieWei/2.cutdata
nwd=~/project/DNAdamage/Mirin/IVF/Pol2_XieWei/3.align
mkdir $nwd
cd $wd

for i in *R1.cut.fq
do
echo $i
t=${i/R1.cut.fq/R2.cut.fq}
echo $t
nohup bowtie2 --local --very-sensitive-local --no-mixed --no-discordant --no-unal -p 8 -X 1000 -x ~/ref/bowtie_index/mm10 -1 $i -2 $t -S $nwd/${i%%_*}.sam > $nwd/${i%%_*}.log &
done
```

for i in *bam; do nohup bamCoverage -p 3 -b $i -o ${i/bam/bw} --binSize 25 --normalizeUsing RPKM --ignoreForNormalization chrM --centerReads & done