# Pipeline for mirin treated embryos

## take IVF embryos for instant

### 1.Fastqc

```sh
wd=~/project/DNAdamage/Mirin/IVF/1.rawdata
cd $wd
mkdir fastqc
for i in *fastq.gz; do
nohup fastqc $i -o ./fastqc &
done

cd ./fastqc/
multiqc -n 'IVF_mirin_raw' ./ &
```

### 2.Cutadapt

```sh
wd=~/project/DNAdamage/Mirin/IVF/1.rawdata
nwd=~/project/DNAdamage/Mirin/IVF/2.cutdata
mkdir $nwd
for i in *R1*.fastq.gz; do
t=${i/R1/R2}
echo $i
echo $t
nohup cutadapt -j 3 -a AGATCGGAAGAGC  -A AGATCGGAAGAGC --trim-n -m 80 -q 20,20 -o $nwd/${i%_S*}_R1.cut.fq -p $nwd/${t%_S*}_R2.cut.fq $i $t > $nwd/${i%%_*}.cut.log &
done
```

### 2ex1. Rename files

> sample_info.txt  
> sample_name     sample_type     experiment      treatment       time    stage repeat
> WMZ-0404-1      Total-RNAseq    IVF     ctrl    5hpf    PN3   rep1
> WMZ-0404-2      Total-RNAseq    IVF     ctrl    5hpf    PN3   rep2
> WMZ-0404-3      Total-RNAseq    IVF     ctrl    10hpf   PN5   rep1
> WMZ-0404-4      Total-RNAseq    IVF     ctrl    10hpf   PN5   rep2
> WMZ-0404-5      Total-RNAseq    IVF     mirin   10hpf   PN5   rep1
> ......

```sh
awk 'NR>1 {k=$1; i=$3"-"$4"-"$5"-"$6"-"$7; a=i"_R1.cut.fq"; b=i"_R2.cut.fq"; c=k"_R1.cut.fq"; d=k"_R2.cut.fq"; print c,a; print d,b }' ../1.rawdata/sample_info.txt | xargs -n2 mv
```

### 2ex2. Fastqc after cutadapt

```sh
wd=~/project/DNAdamage/Mirin/IVF/2.cutdata
cd $wd
mkdir fastqc
for i in *cut.fq; do
nohup fastqc $i -o ./fastqc &
done

cd ./fastqc/
multiqc -n 'cutdata' ./ &
```

### 3.align

```sh
wd=~/project/DNAdamage/Mirin/IVF/2.cutdata
nwd=~/project/DNAdamage/Mirin/IVF/3.align
mkdir $nwd
cd $wd
for i in *R1.cut.fq
do
echo $i
t=${i/R1.cut.fq/R2.cut.fq}
echo $t
nohup hisat2 -p 5 --dta-cufflinks --no-discordant  -t -x ~/ref/hisat2_index/mm10/genome -1 $i -2 $t -S $nwd/${i%%.*}.sam  > $nwd/${i%%.*}.log &
echo ${i%%.*}.sam
done
```

### 3ex1.flagstat

```sh
wd=~/project/DNAdamage/Mirin/IVF/3.align
cd $wd
mkdir flagstat
for i in *.sam; do
echo $i
nohup samtools flagstat -@ 4 $i > ./flagstat/${i/.sam/.flagstat} &
done
```

### 4.sam2bam

```sh
wd=~/project/DNAdamage/Mirin/IVF/3.align
cd $wd
for i in *sam
do
nohup samtools sort -@ 4 -o ${i%.*}.sorted.bam $i &
done

nwd=~/project/DNAdamage/Mirin/IVF/4.proper
mkdir $nwd
for i in *sorted.bam
do
nohup samtools view -bf 0x2 -q 20 -@ 4 $i -o $nwd/${i%sort*}proper.bam  &
done
```

### 5.stringtie

```sh
wd=~/project/DNAdamage/Mirin/IVF/4.proper
nwd=~/project/DNAdamage/Mirin/IVF/5.stringtie
cd $wd
for i in *proper.bam; do
outfolder=${i%.proper*}
mkdir -p $nwd/$outfolder
nohup stringtie -p 8 $i -b $nwd/$outfolder -e -G ~/ref/gtf/mm10.gtf -o $nwd/$outfolder/$outfolder_stringtie.gtf > $nwd/$outfolder.log &
done
```

### 6.TEtranscripts

```sh
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ~/ref/STAR_index --genomeFastaFiles ~/ref/genome_fasta/mm10/mm10.fa --sjdbGTFfile ~/ref/gtf/mm10.gtf --sjdbOverhang 149
```

```sh
nwd=~/project/DNAdamage/Mirin/IVF/6.TEtranscripts/6.3.align_STAR
mkdir -p $nwd
wd=~/project/DNAdamage/Mirin/IVF/2.cutdata
cd $wd

for i in *10hpf*R1.cut.fq; do
t=${i/R1/R2};
pf=${i%%_R*}.;
nohup STAR --runThreadN 6 --genomeDir ~/ref/STAR_index --readFilesIn ${i} ${t} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $nwd/${pf} --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 &
done
```

```sh
nohup TEtranscripts -t siSycp3-l2C-rep1.Aligned.sortedByCoord.out.bam  siSycp3-l2C-rep2.Aligned.sortedByCoord.out.bam -c ctrl-l2C-rep1.Aligned.sortedByCoord.out.bam  ctrl-l2C-rep2.Aligned.sortedByCoord.out.bam --GTF ~/ref/gtf/mm10.gtf --TE ~/ref/gtf/mm10_rmsk_TE.gtf --sortByPos &
```


### 7.RRR
```sh
wd=~/project/DNAdamage/Mirin/partheno/4.proper
cd $wd
nwd=~/project/DNAdamage/Mirin/partheno/7.RRR/7.5.featurecount_RRR
mkdir $nwd

nohup featureCounts -T 8 -F GTF -p -t CDS -g gene_id -a ~/project/NT/2C/ref_publised_data/RRR_region/whole_rrr_class_mm10_lftover_tran.gtf -o $nwd/partheno_RRRregion_counts.txt *bam &

nohup featureCounts -T 8 -F GTF -p -t CDS -g gene_id -a ~/ref/gtf/mm10.gtf -o $nwd/partheno_allgene_counts.txt *bam &