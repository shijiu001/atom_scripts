# Pipeline for Mirin treated embryos ATAC

## take IVF embryos for instant

### 0.md5check and rawdata extract

```sh
conda activate chip

wd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720
cd $wd
grep 'FQZ-0712' MD5.txt | md5sum -c > md5check.txt

nwd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720/1.rawdata
mkdir $nwd
find $wd -name '*gz' | xargs -i mv {} $nwd
```

### 0ex1.seq combine

```sh

wd1=~/project/2C_related/Otx1/20240103_ATAC/1.rawdata/seqbatch1/1.rawdata
wd2=~/project/2C_related/Otx1/20240103_ATAC/1.rawdata/seqbatch2/1.rawdata
nwd=~/project/2C_related/Otx1/20240103_ATAC/1.rawdata/seq_combined
cd $wd1

for i in *R1*fastq.gz; do
t=${i%%_*}; echo $t; 
cat ${wd1}/${t}*R1*fastq.gz ${wd2}/${t}*R1*fastq.gz > ${nwd}/${t}_R1.fastq.gz; 
cat ${wd1}/${t}*R2*fastq.gz ${wd2}/${t}*R2*fastq.gz > ${nwd}/${t}_R2.fastq.gz; 
done

```

### 1.Fastqc

```sh
wd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720/1.rawdata
cd $wd
mkdir fastqc
for i in *fastq.gz; do
nohup fastqc $i -o ./fastqc &
done

cd ./fastqc/
multiqc -n 'rawdata' ./ &
```

### 2.cutadapt

```sh
wd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720/1.rawdata
nwd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720/2.cutdata
mkdir $nwd
for i in *R1*.fastq.gz; do
t=${i/R1/R2}
echo $i
echo $t
nohup cutadapt -j 3 -a CTGTCTCTTATA  -A CTGTCTCTTATA --trim-n -m 80 -q 20,20 -o $nwd/${i%_S*}_R1.cut.fq -p $nwd/${t%_S*}_R2.cut.fq $i $t > $nwd/${i%%_*}.cut.log &
done
```

### 2ex1.Rename files

> sample_name   sample_type expriment   treatment   time    stage   repeat  batch   strain  sample_input_size   index_i5    index_i7    amplify    concentration(12μL)
> FQZ-0712-1    ATAC-seq    IVF    ctrl    22hpf    e2C    rep1    20230712    BDF1    ~150cells    N101    N931    19cycles    176
> FQZ-0712-2    ATAC-seq    IVF    ctrl    22hpf    e2C    rep2    20230712    BDF1    ~150cells    N101    N932    19cycles    182
> FQZ-0712-3    ATAC-seq    IVF    mirin    22hpf    e2C    rep1    20230712    BDF1    ~45cells    N101    N933    21cycles    78
> FQZ-0712-4    ATAC-seq    IVF    mirin    22hpf    e2C    rep2    20230712    BDF1    ~45cells    N101    N934    21cycles    102

```sh
awk 'NR>1 {k=$1; i=$3"-"$4"-"$5"-"$6"-"$7; a=i"_R1.cut.fq"; b=i"_R2.cut.fq"; c=k"_R1.cut.fq"; d=k"_R2.cut.fq"; print c,a; print d,b }' ../1.rawdata/sampleinfo.txt | xargs -n2 mv
```

### 2ex2.Fastqc after cutadapt

```sh
wd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720/2.cutdata
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
wd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720/2.cutdata
nwd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720/3.align
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

### 3ex1.flagstat

```sh
wd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720/3.align
cd $wd
mkdir flagstat
for i in *.sam; do
echo $i
nohup samtools flagstat -@ 4 $i > ./flagstat/${i/.sam/.flagstat} &
done
```

### 4.sam2bam

```sh
wd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720/3.align
cd $wd
for i in *sam
do
nohup samtools sort -@ 4 -o ${i%.*}.sorted.bam $i &
done
```

### 3ex1.ChrM rate & ChrM remove

```sh
wd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720/3.align
cd $wd

for i in *bam
do 
samtools index $i
done

for i in *bam
do echo $i
samtools view -c $i
samtools view -c $i chrM
done

nwd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720/4.rmdup
mkdir $nwd
for i in *bam
do
samtools view -h $i | grep -v 'chrM' | samtools view -b -o $nwd/${i/sorted.bam/chrM-rm.bam}
done
```

### 5.rmdup

```sh
for i in *bam
do 
nohup sambamba markdup -t 4 -r $i ${i/bam/rmdup.bam} &
done
```

### 6.ATAC_shift

```sh
nwd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720/5.ATAC_shift
mkdir $nwd

for i in *rmdup.bam; do 
alignmentSieve -b $i -o ${nwd}/${i/.bam/.shift.bam} -p 8 --ATACshift
done

for i in *bam; do 
nohup samtools sort -@ 8 -o ${i%.*}.sorted.bam $i &
done

for i in *shift.sorted.bam; do samtools index $i; done
```

### 7.bigwig

```sh
nwd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720/6.bigwig
mkidr $nwd

for i in *sorted.bam; do 
bamCoverage -p 8 --normalizeUsing RPKM -b $i -o ${nwd}/${i%%.*}.bw
done
```

### 8.call peaks

```sh
nwd=~/project/DNAdamage/Mirin/IVF/ATAC_20230720/7.peaks
for i in *.bam; do 
nohup macs3 callpeak -t $i -f BAMPE -g mm -n ${i%%.*} -B -q 0.001 --keep-dup 1 --outdir $nwd & 
done
```

### 9.DiffBind

```sh
macs3 bdgdiff --t1 ctrl-PN3-H3K4me3-rep1_treat_pileup.bdg --c1 ctrl-PN3-H3K4me3-rep1_control_lambda.bdg --t2 ctrl-e2C-H3K4me3-rep1_treat_pileup.bdg --c2 ctrl-e2C-H3K4me3-rep1_control_lambda.bdg -C 2 -g 200 --d1 17 --d2 22 --outdir ./macs3_bdgdiff/ --o-prefix ctrl_PN3-vs-e2C
```

computeMatrix scale-regions -b 3000 -a 3000 -S IVF-ctrl-22hpf-e2C-rep1.bw IVF-ctrl-22hpf-e2C-rep2.bw IVF-mirin-22hpf-e2C-rep1.bw IVF-mirin-22hpf-e2C-rep2.bw -R ~/project/DNAdamage/Mirin/resource_set/zga_genelist_position.bed --regionBodyLength 6000 --skipZeros -o matrix_zgagene_ATAC.gz

computeMatrix scale-regions -S IVF-ctrl-22hpf-e2C-rep1.bw IVF-ctrl-22hpf-e2C-rep2.bw IVF-mirin-24hpf-e2C-rep1.bw IVF-mirin-24hpf-e2C-rep2.bw -R ~/project/NT/2C/ref_publised_data/RRR_region/whole_rrr_class_mm10_lftover_tran.bed --regionBodyLength 10000 --sortRegions keep --skipZeros -o new_matrix_e2C_RRR_ATAC.gz
 1798  plotHeatmap -m new_matrix_e2C_RRR_ATAC.gz -o new_e2C_RRR_ATAC.png --zMax 5 --whatToShow 'plot, heatmap and colorbar'
 1799  plotHeatmap -m mwe_matrix_zgagene_ATAC.gz -o new_e2C_ZGAgene_ATAC.png --zMax 8 --whatToShow 'plot, heatmap and colorbar'
 1800  computeMatrix scale-regions -b 3000 -a 3000 -S IVF-ctrl-22hpf-e2C-rep1.bw IVF-ctrl-22hpf-e2C-rep2.bw IVF-mirin-24hpf-e2C-rep1.bw IVF-mirin-24hpf-e2C-rep2.bw IVF-mirin-22hpf-e2C-rep1.bw IVF-mirin-22hpf-e2C-rep2.bw -R ~/project/DNAdamage/Mirin/resource_set/zga_genelist_position.bed --regionBodyLength 6000 --skipZeros -o all_matrix_zgagene_ATAC.gz
 1801  plotHeatmap -m all_matrix_zgagene_ATAC.gz -o all_e2C_ZGAgene_ATAC.png --zMax 8 --whatToShow 'plot, heatmap and colorbar'
 1802  plotHeatmap -m all_matrix_zgagene_ATAC.gz -o all_e2C_ZGAgene_ATAC.png --zMax 12 --whatToShow 'plot, heatmap and colorbar'
 1803  plotHeatmap -m all_matrix_zgagene_ATAC.gz -o all_e2C_ZGAgene_ATAC.png --zMax 15 --whatToShow 'plot, heatmap and colorbar'

 computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -S IVF-ctrl-22hpf-e2C-rep1.bw IVF-ctrl-22hpf-e2C-rep2.bw IVF-mirin-22hpf-e2C-rep1.bw IVF-mirin-22hpf-e2C-rep2.bw -R ~/project/DNAdamage/Mirin/IVF/Pol2_XieWei/resorce_data/Pol2_bindsite_1Cspecific.bed ~/project/DNAdamage/Mirin/IVF/Pol2_XieWei/resorce_data/Pol2_bindsite_l2Cspecific.bed ~/project/DNAdamage/Mirin/IVF/Pol2_XieWei/resorce_data/Pol2_bindsite_shared.bed --regionBodyLength 5000 --skipZeros --sortRegions keep -o matrix_Pol2_bindsite.gz
 1023  plotHeatmap -m matrix_Pol2_bindsite.gz -o Pol2_bindsite.png --zMax 5 --whatToShow 'plot, heatmap and colorbar'

 findMotifsGenome.pl WT_enrich.HomerInput.txt mm10 ~/project/2C_related/Dux/ATAC_20230720/8.motif_homer/ -size 300 -p 10
