0.check md5
nohup find ~/project/DDX5/nchip_H3K4me3/data -name '*.fastq.gz' -type f -print0 | xargs -0 md5sum > ./nchip——inputbu.md5 &


mv WMT-221018-1_S5_L003_R1_001.fastq.gz input_3w_ko_rep1_bu.R1.fastq.gz
mv WMT-221018-1_S5_L003_R2_001.fastq.gz input_3w_ko_rep1_bu.R2.fastq.gz
mv WMT-221018-2_S53_L004_R1_001.fastq.gz input_3w_ko_rep2_bu.R1.fastq.gz
mv WMT-221018-2_S53_L004_R2_001.fastq.gz input_3w_ko_rep2_bu.R2.fastq.gz
mv WMT-221018-3_S54_L004_R1_001.fastq.gz input_3w_ctrl_rep1_bu.R1.fastq.gz
mv WMT-221018-3_S54_L004_R2_001.fastq.gz input_3w_ctrl_rep1_bu.R2.fastq.gz
mv WMT-221018-4_S55_L004_R1_001.fastq.gz input_3w_ctrl_rep2_bu.R1.fastq.gz
mv WMT-221018-4_S55_L004_R2_001.fastq.gz input_3w_ctrl_rep2_bu.R2.fastq.gz

nohup cat WMT-0122-4_L1_704D01.R1.fastq.gz WMT-0122-4.R1.fastq.gz > 3w_gv_ko_rep1.R1.fastq.gz &
nohup cat WMT-0122-4_L1_704D01.R2.fastq.gz WMT-0122-4.R2.fastq.gz > 3w_gv_ko_rep1.R2.fastq.gz &
nohup cat WMT-0122-6_L1_706D01.R1.fastq.gz WMT-0122-6.R1.fastq.gz > 3w_gv_ko_rep3.R1.fastq.gz &
nohup cat WMT-0122-6_L1_706D01.R2.fastq.gz  WMT-0122-6.R2.fastq.gz > 3w_gv_ko_rep3.R2.fastq.gz &
1.rawdata
```
wd=~/project/DDX5/nchip_H3K4me3/data
cd $wd
mkdir fastqc
for i in *fastq.gz
do
nohup fastqc $i -o ./fastqc &
done

multiqc -n "raw_wgbs_fastqc.html" ./

```
2.clean data （clean完会直接对cleandata进行fastq)

wd=~/project/DDX5/nchip_H3K4me3/data
nwd=~/project/DDX5/nchip_H3K4me3/data/2.cutdata
mkdir $nwd
cd $wd
for i in *.R1.fastq.gz
do
t=${i/R1/R2}
echo $i
echo $t
nohup trim_galore -q 25  --phred33 -fastqc --stringency 3 --length 50 -j 4 --path_to_cutadapt /home1/wangmt/miniconda3/envs/python3.7/bin/cutadapt --gzip --paired $wd/$i $wd/$t -o $nwd > $nwd/${i%%.*}_cut.log &
done

multiqc -n "clean_fastqc.html" ./
```

3. align to mm10

```
conda activate python3.9
wd=~/project/DDX5/nchip_H3K4me3/data/2.cutdata
nwd=~/project/DDX5/nchip_H3K4me3/data/3.align
mkdir $nwd
cd $wd
for i in *.R1_val_1.fq.gz
do
t=${i/R1_val_1/R2_val_2}
echo $i
echo $t
nohup bowtie2 -p 10 -x /home1/share/bowtie2_index/mm10  --local --very-sensitive-local -t --no-mixed --no-discordant --no-unal -1 $i -2 $t -S $nwd/${i%%.*}.sam > $nwd/${i%%.*}.log &
done

echo -e 'sample\ttotal_paired(million)\tproperly_paired(million)\tratio\tmore_1_time(million)\tratio\talign_ratio_all' >  bowtie2_align_summary.tab
for i in *.log; do
echo $i
key=${i%%.*}
awk -v k=$key 'BEGIN{FS="[ \t();%]+";OFS="\t"} NR==5 {t=$1} NR==8{a=$2;b=$3} NR==9{c=$2;d=$3} NR==10{f=$1} END{print k,t/1e6,a/1e6,b"%",c/1e6,d"%",f"%"}' $i >> bowtie2_align_summary.tab
done
```



4. sam2bam

```
conda deactivate
wd=~/project/DDX5/nchip_H3K4me3/data/3.align
cd $wd
for i in *sam; do
nohup samtools sort -@ 5 -o ${i%.*}.sorted.bam $i &
done

5.filter proper reads
wd=~/project/DDX5/nchip_H3K4me3/data/3.align
nwd=~/project/DDX5/nchip_H3K4me3/data/4.proper
mkdir $nwd
cd $wd
for i in *sorted.bam; do
nohup samtools view -@ 5 -q 20 $i -o $nwd/${i%sort*}proper.bam &
done

```

6. (optional) Remove dup

```
wd=~/project/DDX5/nchip_H3K4me3/data/4.proper
nwd=~/project/DDX5/nchip_H3K4me3/data/5.unique
cd $wd
mkdir $nwd
for i in *.proper.bam;do
nohup sambamba markdup -r -t 4 -p $i $nwd/${i%proper*}duprm.bam > $nwd/${i%proper*}duprm.log &
done
```
#nohup samtools flagstat -@ 5 8Cell_Klf5_rep1.sam > 8Cell_Klf5_rep1.sam.flagstat &
nohup samtools flagstat -@ 5 3w_gv_H3K4me3_ctrl_rep1.duprm.bam > 3w_gv_H3K4me3_ctrl_rep1.duprm.flagstat &
nohup samtools flagstat -@ 5 3w_gv_H3K4me3_ko_rep1.duprm.bam > 3w_gv_H3K4me3_ko_rep1.duprm.flagstat &
nohup samtools flagstat -@ 5 3w_gv_H3K4me3_ctrl_input_rep1.duprm.bam > 3w_gv_H3K4me3_ctrl_input_rep1.duprm.flagstat &
nohup samtools flagstat -@ 5 3w_gv_H3K4me3_ctrl_input_rep2.duprm.bam > 3w_gv_H3K4me3_ctrl_input_rep2.duprm.flagstat &
nohup samtools flagstat -@ 5 3w_gv_H3K4me3_ko_input_rep1.duprm.bam > 3w_gv_H3K4me3_ko_input_rep1.duprm.flagstat &
nohup samtools flagstat -@ 5 3w_gv_H3K4me3_ko_input_rep2.duprm.bam > 3w_gv_H3K4me3_ko_input_rep2.duprm.flagstat &


6.1 bw
```
wd=~/project/DDX5/nchip_H3K4me3/data/5.unique
nwd=~/project/DDX5/nchip_H3K4me3/data/6.bw
mkdir $nwd
cd $wd
for i in *.bam
do
name=${i%%.*}
echo $name
nohup bamCoverage  -p 5 -b $i -o $nwd/${name}.bw \
--normalizeUsing RPKM  --centerReads --binSize 25 -ignore chrM > $nwd/${name}.bw.log &
done
```
nohup bamCoverage  -p 5 -b 3w_gv_H3K4me3_ctrl_rep2.duprm.bam -o ~/project/DDX5/3w_8w_H3K4me3_cutrun/6.bw/3w_gv_H3K4me3_ctrl_rep2.bw \
--normalizeUsing RPKM  --centerReads --binSize 25 -ignore chrM > ~/project/DDX5/3w_8w_H3K4me3_cutrun/6.bw/3w_gv_H3K4me3_ctrl_rep2.log &
nohup bamCoverage  -p 5 -b 3w_gv_H3K4me3_ko_rep1.duprm.bam -o ~/project/DDX5/3w_8w_H3K4me3_cutrun/6.bw/3w_gv_H3K4me3_ko_rep1.bw \
--normalizeUsing RPKM  --centerReads --binSize 25 -ignore chrM > ~/project/DDX5/3w_8w_H3K4me3_cutrun/6.bw/3w_gv_H3K4me3_ko_rep1.log &
```
7.callpeak (macs)
conda activate python3.9
wd=~/project/DDX5/H3K4me3_3w_GV_CQQ/5.unique
nwd=~/project/DDX5/H3K4me3_3w_GV_CQQ/8.narrowpeak
mkdir $nwd
cd $wd
for i in *bam
do
#没有input时
#xie wei call broadpeak时：–broad –nomodel –nolambda
##narrowpeak
nohup macs2 callpeak -t $i  -f BAMPE -g mm -n ${i%%.*} -B -q 0.05 --keep-dup 1 --outdir $nwd > $nwd/${i%%.*}.log &
##broadpeak
nohup macs2 callpeak -t $i  -f BAMPE -g mm -n ${i%%.*} -B --nomodel --nolambda --broad --broad-cutoff 0.05 --keep-dup all --outdir $nwd > $nwd/${i%%.*}.log &
done



#有input时,narrowPeak
nohup macs2 callpeak -t test_gv.duprm.bam -c test_input.duprm.bam  -f BAMPE -g mm -n test_gv -B -q 0.05 --keep-dup 1 --outdir ~/project/DDX5/new_gv_h3k4me3_nchip/7.narrowpeak > ~/project/DDX5/new_gv_h3k4me3_nchip/7.narrowpeak/test_gv_narrow.log &
nohup macs2 callpeak -t ctrl_gv.duprm.bam -c ctrl_input.duprm.bam  -f BAMPE -g mm -n ctrl_gv -B -q 0.05 --keep-dup 1 --outdir ~/project/DDX5/new_gv_h3k4me3_nchip/7.narrowpeak > ~/project/DDX5/new_gv_h3k4me3_nchip/7.narrowpeak/ctrl_gv_narrow.log &
nohup macs2 callpeak -t ko_gv.duprm.bam -c ko_input.duprm.bam  -f BAMPE -g mm -n ko_gv -B -q 0.05 --keep-dup 1 --outdir ~/project/DDX5/new_gv_h3k4me3_nchip/7.narrowpeak > ~/project/DDX5/new_gv_h3k4me3_nchip/7.narrowpeak/ko_gv_narrow.log &

nohup macs2 callpeak -t GV-H3K4me3-Cre1.duprm.bam -c ~/project/DDX5/new_gv_h3k4me3_nchip/5.unique/ko_input.duprm.bam  -f BAMPE -g mm -n GV-H3K4me3-Cre1 -B -q 0.05 --keep-dup 1 --outdir ~/project/DDX5/H3K4me3_3w_GV_CQQ/7.narrowpeak > ~/project/DDX5/H3K4me3_3w_GV_CQQ/7.narrowpeak/ko_gv_rep1_narrow.log &
nohup macs2 callpeak -t GV-H3K4me3-Cre2.duprm.bam -c ~/project/DDX5/new_gv_h3k4me3_nchip/5.unique/ko_input.duprm.bam  -f BAMPE -g mm -n GV-H3K4me3-Cre2 -B -q 0.05 --keep-dup 1 --outdir ~/project/DDX5/H3K4me3_3w_GV_CQQ/7.narrowpeak > ~/project/DDX5/H3K4me3_3w_GV_CQQ/7.narrowpeak/ko_gv_rep2_narrow.log &
nohup macs2 callpeak -t GV-H3K4me3-Cre-1.duprm.bam -c ~/project/DDX5/new_gv_h3k4me3_nchip/5.unique/ctrl_input.duprm.bam  -f BAMPE -g mm -n GV-H3K4me3-Cre-1 -B -q 0.05 --keep-dup 1 --outdir ~/project/DDX5/H3K4me3_3w_GV_CQQ/7.narrowpeak > ~/project/DDX5/H3K4me3_3w_GV_CQQ/7.narrowpeak/ctrl_gv_rep1_narrow.log &
nohup macs2 callpeak -t GV-H3K4me3-Cre-2.duprm.bam -c ~/project/DDX5/new_gv_h3k4me3_nchip/5.unique/ctrl_input.duprm.bam  -f BAMPE -g mm -n GV-H3K4me3-Cre-2 -B -q 0.05 --keep-dup 1 --outdir ~/project/DDX5/H3K4me3_3w_GV_CQQ/7.narrowpeak > ~/project/DDX5/H3K4me3_3w_GV_CQQ/7.narrowpeak/ctrl_gv_rep2_narrow.log &

#有input时,broadPeak
nohup macs2 callpeak -t test_gv.duprm.bam -c test_input.duprm.bam  -f BAMPE -g mm -n test_gv -B --nomodel --nolambda --broad --broad-cutoff 0.05 --keep-dup all --outdir ~/project/DDX5/new_gv_h3k4me3_nchip/8.broadpeak > ~/project/DDX5/new_gv_h3k4me3_nchip/8.broadpeak/test_gv_broad.log &
nohup macs2 callpeak -t ctrl_gv.duprm.bam -c ctrl_input.duprm.bam  -f BAMPE -g mm -n ctrl_gv -B --nomodel --nolambda --broad --broad-cutoff 0.05 --keep-dup all --outdir ~/project/DDX5/new_gv_h3k4me3_nchip/8.broadpeak > ~/project/DDX5/new_gv_h3k4me3_nchip/8.broadpeak/ctrl_gv_broad.log &
nohup macs2 callpeak -t ko_gv.duprm.bam -c ko_input.duprm.bam  -f BAMPE -g mm -n ko_gv -B --nomodel --nolambda --broad --broad-cutoff 0.05 --keep-dup all --outdir ~/project/DDX5/new_gv_h3k4me3_nchip/8.broadpeak > ~/project/DDX5/new_gv_h3k4me3_nchip/8.broadpeak/ko_gv_broad.log &

nohup macs2 callpeak -t GV-H3K4me3-Cre1.duprm.bam -c ~/project/DDX5/new_gv_h3k4me3_nchip/5.unique/ko_input.duprm.bam  -f BAMPE -g mm -n GV-H3K4me3-b-Cre1 -B --nomodel --nolambda --broad --broad-cutoff 0.05 --keep-dup 1 --outdir ~/project/DDX5/H3K4me3_3w_GV_CQQ/8.broadpeak > ~/project/DDX5/H3K4me3_3w_GV_CQQ/8.broadpeak/ko_gv_rep1_broad.log &
nohup macs2 callpeak -t GV-H3K4me3-Cre2.duprm.bam -c ~/project/DDX5/new_gv_h3k4me3_nchip/5.unique/ko_input.duprm.bam  -f BAMPE -g mm -n GV-H3K4me3-b-Cre2 -B --nomodel --nolambda --broad --broad-cutoff 0.05 --keep-dup 1 --outdir ~/project/DDX5/H3K4me3_3w_GV_CQQ/8.broadpeak > ~/project/DDX5/H3K4me3_3w_GV_CQQ/8.broadpeak/ko_gv_rep2_broad.log &
nohup macs2 callpeak -t GV-H3K4me3-Cre-1.duprm.bam -c ~/project/DDX5/new_gv_h3k4me3_nchip/5.unique/ctrl_input.duprm.bam  -f BAMPE -g mm -n GV-H3K4me3-b-Cre-1 -B --nomodel --nolambda --broad --broad-cutoff 0.05 --keep-dup 1 --outdir ~/project/DDX5/H3K4me3_3w_GV_CQQ/8.broadpeak > ~/project/DDX5/H3K4me3_3w_GV_CQQ/8.broadpeak/ctrl_gv_rep1_broad.log &
nohup macs2 callpeak -t GV-H3K4me3-Cre-2.duprm.bam -c ~/project/DDX5/new_gv_h3k4me3_nchip/5.unique/ctrl_input.duprm.bam  -f BAMPE -g mm -n GV-H3K4me3-b-Cre-2 -B --nomodel --nolambda --broad --broad-cutoff 0.05 --keep-dup 1 --outdir ~/project/DDX5/H3K4me3_3w_GV_CQQ/8.broadpeak > ~/project/DDX5/H3K4me3_3w_GV_CQQ/8.broadpeak/ctrl_gv_rep2_broad.log &

8.样本重复性检验
#利用各自callpeak的bed坐标
cat *.narrowPeak | cut -f1-3 | sort -k1,1 -k2n,2 | bedtools merge > merge.bed
nohup multiBigwigSummary BED-file --BED ~/project/DDX5/cutruntest/6.peak/merge.bed -b *.bw -out 4sample.npz &
plotCorrelation -in 4sample.npz --corMethod pearson --whatToPlot heatmap --skipZeros --plotNumbers --removeOutliers -o cor_heatmap.png
plotPCA -in 4sample.npz -o pca.png

nohup multiBamSummary bins -b *.bam -o input_compare.npz &
plotCorrelation -in input_compare.npz --corMethod pearson --whatToPlot heatmap --skipZeros --plotNumbers --removeOutliers -o cor_input_heatmap.png
plotPCA -in input_compare.npz -o pca.png

nohup plotFingerprint -b ctrl_gv.duprm.bam ctrl_input.duprm.bam ko_gv.duprm.bam ko_input.duprm.bam test_gv.duprm.bam test_input.duprm.bam ~/project/DDX5/nchip_H3K4me3/4.H3K4me3_3w_GV_CQQ/5.unique/GV-H3K4me3-Cre1.duprm.bam ~/project/DDX5/nchip_H3K4me3/4.H3K4me3_3w_GV_CQQ/5.unique/GV-H3K4me3-Cre2.duprm.bam ~/project/DDX5/nchip_H3K4me3/4.H3K4me3_3w_GV_CQQ/5.unique/GV-H3K4me3-Cre-1.duprm.bam ~/project/DDX5/nchip_H3K4me3/4.H3K4me3_3w_GV_CQQ/5.unique/GV-H3K4me3-Cre-2.duprm.bam  -l ctrl_gv ctrl_input ko_gv ko_input test_gv test_input ko_rep1 ko_rep2 ctrl_rep1 ctrl_rep2 -plot test_with_cqq.png &

#利用promoter区的坐标
重复性好的话合并bam文件callpeak
cat mm10.refFlat |cut -f1-5 |awk '{if($4=="+") {print $3"\t"$5-2000"\t"$5+500"\t"$1"\t"$2"\t"$4} else {print $3"\t"$6-500"\t"$6+2000"\t"$1"\t"$2"\t"$4}}' > mm10_promoter.bed
nohup multiBigwigSummary BED-file --BED ~/project/anotation/mm10_promoter.bed -b *.bw --outRawCounts allsample.txt -p 10 -out allsample.npz &
plotCorrelation -in allsample.npz --corMethod pearson --whatToPlot heatmap --skipZeros --plotNumbers --removeOutliers -o cor_heatmap_bw.png
plotPCA -in allsample.npz -o pca.png
#soothscatter

8.1合并bam文件
nohup samtools merge 3w_gv_H3K4me3_ctrl.bam 3w_gv_H3K4me3_ctrl_rep1.duprm.bam 3w_gv_H3K4me3_ctrl_rep2.duprm.bam &
nohup samtools merge 3w_gv_H3K4me3_ko.bam 3w_gv_H3K4me3_ko_rep1.duprm.bam 3w_gv_H3K4me3_ko_rep2.duprm.bam &

9.找差异peak（合并后的peak文件）
nohup bedtools intersect -a 8Cell_Klf5_rep1_peaks.narrowPeak -b 8Cell_siTfap2c_Klf5_rep1_peaks.narrowPeak -v > deg.peak &
plotCorrelation -in 4sample_bw.npz --corMethod pearson --whatToPlot heatmap --skipZeros --plotNumbers --removeOutliers -o cor_heatmap.png

10.Peak anno(注释差异peak或各自peak)
wd=~/project/DDX5/cutruntest/6.peak
nwd=~/project/DDX5/cutruntest/7.Peak_anno
mkdir $nwd
cd $wd
for i in *narrowPeak
do
nohup annotatePeaks.pl $i mm10 1> $nwd/${i%_*}.peakAnno.xls 2>$nwd/${i%_*}_annLog.txt &
done


11.motif(可针对差异peak或各自peak)
wd=~/project/DDX5/cutruntest/6.peak
nwd=~/project/DDX5/cutruntest/8.motif
mkdir $nwd
cd $wd
for i in *narrowPeak
do
nohup findMotifsGenome.pl $i mm10 $nwd/${i%_*}_peakAnalysis -size given -p 7 -len 8 -preparsedDir $nwd > $nwd/${i%_*}_peakAnalysis.log &
done

#准备bed文件，提取tss,tts位置信息
awk '{print $3"\t"$5"\t"$6"\t"$1"\t"$2"\t"$4}' mm10.refFlat >./mm10.TSS_TTS.bed
nohup computeMatrix scale-regions -p 3 -S ~/project/DDX5/cutruntest/bw/8Cell_Klf5_rep1.bw -R ~/project/anotation/mm10.TSS_TTS.bed \
-o 8Cell_Klf5_rep1_region.gz -m 5000 -b 2000 -a 2000 --skipZeros --binSize 25 &

nohup computeMatrix reference-point -p 8 -S ~/project/DDX5/cutruntest/bw/8Cell*.bw -R ~/project/anotation/mm10.TSS_TTS.bed \
--referencePoint TSS -a 2000 -b 500 --skipZeros --binSize 25 -out 8Cell_Klf5_tss_all.gz &

plotHeatmap -m 8Cell_Klf5_tss_all.gz --perGroup -out 8Cell_Klf5_tss_all_2.png
#--perGroup 多个样本画在一起

plotProfile -m 8Cell_Klf5_tss_all.gz --perGroup --plotType heatmap -out heatmap.png




#比较peak数目，统一深度
cd ~/work_space/4.ProjectET/analysis/fig1/build/1.Tfap2c_bam
mkdir flagstat
for i in *.bam
do
echo $i
nohup samtools flagstat -@ 4 $i > ./flagstat/${i/.bam/.flagstat} &
done

rm Samtools_Flagstat_summary.tab

echo -e 'sample\tmapped_pairs' >  Samtools_Flagstat_summary.tab
for i in *flagstat; do
echo $i
key=${i/_stat.log/}
cat $i | awk -v k=$key 'BEGIN{FS="[ \t();%]+";OFS="\t"} NR==5{a=$1} END{print k,a/2}' >>  Samtools_Flagstat_summary.tab
done

####1. downsample
# sample mapped_pairs min ds
# 8Cell_Klf5 32843017 32843017 1
# 8Cell_siTfap2c_Klf5 46166617 32843017 0.711401856
# Morula_Klf5 45216350 32843017 0.72635268
wd=~/work_space/4.ProjectET/analysis/fig1/build/1.Klf5_bam
cd $wd
nohup samtools view  -@ 10 -bs 1 8Cell_Klf5.duprm.bam > $wd/8Cell_Klf5_ds.bam  &
nohup samtools view  -@ 10 -bs 0.71 8Cell_siTfap2c_Klf5.duprm.bam > $wd/8Cell_siTfap2c_Klf5_ds.bam &
nohup samtools view  -@ 10 -bs 0.73 Morula_Klf5.duprm.bam > $wd/Morula_Klf5_ds.bam &
