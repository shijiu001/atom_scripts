# Pipeline for Mirin treated embryos ATAC

## aaa

### 1.Fastqc

```sh
wd=~/project/2C_related/POL2_inhibit/H3K4me3_CUTandTAG_20231006/1.rawdata/seq_combine
cd $wd
mkdir fastqc
for i in *fastq.gz; do
nohup fastqc $i -o ./fastqc &
done

cd ./fastqc/
multiqc -n 'rawdata_fastqc' ./ &
```

### 2.cutdata

```sh
wd=~/project/2C_related/POL2_inhibit/H3K4me3_CUTandTAG_20231006/1.rawdata/seq_combine
cd $wd
nwd=~/project/2C_related/POL2_inhibit/H3K4me3_CUTandTAG_20231006/2.cutdata
mkdir $nwd

for i in *R1*.fastq.gz; do
t=${i/R1/R2}
echo $i
echo $t
nohup cutadapt -j 3 -a CTGTCTCTTATA  -A CTGTCTCTTATA --trim-n -m 80 -q 20,20 -o $nwd/${i%_*}_R1.cut.fq -p $nwd/${t%_*}_R2.cut.fq $i $t > $nwd/${i%%_*}.cut.log &
done
```

### 2ex1.Rename files

```sh
awk 'NR>1 {k=$1; i=$5"-"$7"-"$3"-"$8; a=i"_R1.cut.fq"; b=i"_R2.cut.fq"; c=k"_R1.cut.fq"; d=k"_R2.cut.fq"; print c,a; print d,b }' ../1.rawdata/sampleinfo.txt | xargs -n2 mv
```

### 2ex2.Fastqc after cutadapt

```sh
wd=~/project/2C_related/POL2_inhibit/H3K4me3_CUTandTAG_20231006/2.cutdata
cd $wd
mkdir fastqc
for i in *cut.fq; do
nohup fastqc $i -o ./fastqc &
done

cd ./fastqc/
multiqc -n 'cutdata_fastqc' ./ &
```

### 3.align

```sh
wd=~/project/2C_related/POL2_inhibit/H3K4me3_CUTandTAG_20231006/2.cutdata
nwd=~/project/2C_related/POL2_inhibit/H3K4me3_CUTandTAG_20231006/3.align
mkdir $nwd
cd $wd

for i in *R1.cut.fq
do
echo $i
t=${i/R1.cut.fq/R2.cut.fq}
echo $t
nohup bowtie2 --no-mixed --no-discordant --no-unal -p 8 -X 1000 -x ~/ref/bowtie_index/mm10 -1 $i -2 $t -S $nwd/${i%%_*}.sam > $nwd/${i%%_*}.log &
done
```

### 3ex1.flagstat

```sh
wd=~/project/2C_related/POL2_inhibit/H3K4me3_CUTandTAG_20231006/3.align
cd $wd
mkdir flagstat
for i in *.sam; do
echo $i
nohup samtools flagstat -@ 4 $i > ./flagstat/${i/.sam/.flagstat} &
done
```

### 4.sam2bam

```sh
wd=~/project/2C_related/POL2_inhibit/H3K4me3_CUTandTAG_20231006/3.align
cd $wd
for i in *sam
do
nohup samtools sort -@ 4 -o ${i%.*}.sorted.bam $i &
done
```

### 5.rmdup

```sh
wd=~/project/2C_related/POL2_inhibit/H3K4me3_CUTandTAG_20231006/3.align
cd $wd
for i in *bam
do 
samtools index $i
done

for i in *bam
do 
nohup sambamba markdup -t 4 -r $i ${i/bam/rmdup.bam} &
done

nwd=~/project/2C_related/POL2_inhibit/H3K4me3_CUTandTAG_20231006/4.rmdup
mkdir $nwd
mv *rmdup* $nwd
```

### 6.bamCoverage

```sh
for i in *sorted.bam; do bamCoverage -p 8 --normalizeUsing RPKM -b $i -o ${i%%.*}.bw; done
```

### 7.call_peaks

```sh

```

### 

computeMatrix reference-point --referencePoint center -b 4000 -a 4000 -S ctrl-PN3-H3K4me3-rep1.bw ctrl-e2C-H3K4me3-rep1.bw DRB-e2C-H3K4me3-rep1.bw -R ../6.peaks/all_summit.bed --skipZeros --sortRe
gions keep -o matrix_all_peak_summit.gz

plotHeatmap -m matrix_all_peak_summit.gz -o all_peaks_summit.png --zMax 15 --whatToShow 'plot, heatmap and colorbar'
