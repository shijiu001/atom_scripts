for i in *R1_001.fastq.gz; do t=${i/R1/R2}; echo $i; echo $t; nohup trim_galore --phred33 --paired --fastqc --clip_R1 5 --three_prime_clip_R1 5 --clip_R2 5 --three_prime_clip_R2 5 --length 80 -q 25 --trim-n -o ../cutdata $i $t & done

for i in *R1*; do t=${i/R1_001_val_1/R2_001_val_2}; echo $i $t; nohup bowtie2 -p 3 -x ~/ref/bowtie_index/mm10 -1 $i -2 $t | samtools sort -O bam -@ 3 -o - > ${i%%R1*}.bam & done

for i in WMZ-1016*R1*; do t=${i/R1_001_val_1/R2_001_val_2}; echo $i $t; nohup bowtie2 -p 10 -x ~/ref/bowtie_index/mm10 -1 $i -2 $t -S ${i%%_R1*}.sam & done


for i in *sam; do samtools sort -O bam -@ 10 -o ${i%%.*}.sorted.bam $i & done

for i in *.bam; do echo $i; nohup samtools flagstat $i > ./flagstat/${i/.bam/.flagstat} & done


for i in *sorted.bam
do
nohup samtools view -bh -f 0x2 -q 30 -@ 4 $i -o $nwd/${i%sorted*}proper.bam  &
done

for i in *;
do nohup picard MarkDuplicates -I $i -O ${i/proper/rmdup} -M ${i/proper.bam/dup_metrics.txt} --REMOVE_DUPLICATES &
done
# picard MarkDuplicates -I WMZ-1016-1_S2_L003
# samtools view -c -F 0x400 WMZ-1016-1_S2_L003.markdup.bam

ls *.bam | xargs -i samtools  index {}
for i in *bam; do nohup bamCoverage -p 3 -b $i -o ${i/bam/bw} --binSize 25 --normalizeUsing RPKM --ignoreForNormalization chrM --centerReads & done

nohup macs2 callpeak -t WMZ-1016-1_S2_L003.rmdup.bam -c WMZ-1016-3_S13_L004.rmdup.bam  -f BAMPE -g mm -n WMZ-1016-1 -B -q 0.05 --keep-dup 1 --outdir ~/project/NT/E7.5/ChIP/ExE/regular_peaks/ > ~/project/NT/E7.5/ChIP/ExE/regular_peaks/WMZ-1016-1.regular_peaks.log &

nohup macs2 callpeak -t WMZ-1016-1_S2_L003.rmdup.bam -c WMZ-1016-3_S13_L004.rmdup.bam -f BAMPE -g mm -n WMZ-1016-1 -B --broad --broad-cutoff 0.05 --keep-dup 1 --outdir ~/project/NT/E7.5/ChIP/ExE/broad_peaks > ~/project/NT/E7.5/ChIP/ExE/broad_peaks/WMZ-1016-1.broad_peaks.log &

for i in WMZ-1030-[56]_*bam; do nohup macs2 callpeak -t $i -c WMZ-1030-7_S4_L003.rmdup.bam -f BAMPE -g mm -n ${i%%_*} -B --broad --broad-cutoff 0.05 --keep-dup 1 --outdir ~/project/NT/E7.5/ChIP/ExE/broad_peaks > ~/project/NT/E7.5/ChIP/ExE/broad_peaks/${i%%_*}.broad_peaks.log & done

for i in FQZ-1029-[12]_*.bam; do nohup macs2 callpeak -t $i -c FQZ-1029-3_S50_L003.rmdup.bam -f BAMPE -g mm -n ${i%%_*} -B -q 0.05 --keep-dup 1 --outdir ~/project/NT/E7.5/ChIP_H3K27me3/VE/regular_peaks/ > ~/project/NT/E7.5/ChIP_H3K27me3/VE/regular_peaks/${i%%_*}.regular_peaks.log & done




for i in WMZ-1022-1[056]_*R1.cut.fq; do t=${i/R1.cut.fq/R2.cut.fq}; y=${i/R1.cut.fq/pwk-Nmasked.sam}; echo $i,$t $y; nohup bowtie2 -p 8 -x ~/ref/snp_mouse/PWK/PWK_PhJ_N-masked/bowtie2_index/GRCm38_PWK_PhJ_N-masked -1 $i -2 $t -S $y & done


awk -F "\t" 'BEGIN {OFS="\t"}  {print "POL2_ctrl_e2C_B1target_slop200_"NR,$1,$2,$3,"."}' POL2_ctrl_e2C_B1target_slop200.bed > POL2_ctrl_e2C_B1target_slop200.HomerInput.txt

nohup findMotifsGenome.pl POL2_ctrl_e2C_B1target_slop200.HomerInput.txt mm10 ./B1_target_motif/ -size 300 -p 10 &