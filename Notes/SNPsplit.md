# SNPsplit

### 下载vcf文件

`wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz`
`md5sum mgp.v5.merged.snps_all.dbSNP142.vcf.gz      # md5检查`

> 这个可能是因为网络问题，我下载了几次都是失败的，md5值都不匹配，后来我就直接下载只有CAST_EiJ的那个文件了

### 建立一个N-masked的参考基因组

```
SNPsplit_genome_preparation --vcf_file CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz --reference_genome ../fasta/ --strain CAST_EiJ --nmasking
cat *.fa > mm10_N-masked.fa
```

### N-masked的bowtie2 index

> 我看到说明文档里提到需要端到端的对比，提到了bowtie2，我就用了这个软件

`bowtie2-build mm10_N-masked.fa mm10_CASTEiJ_N-masked`

### 利用N-masked的index做mapping

> 我只用文章中的一个zygote数据来试

`bowtie2 -p 8 -x ~/ref/snp_mouse/CAST_EiJ/CAST_EiJ_N-masked/bowtie2_index/mm10_CASTEiJ_N-masked -U ./sra/SRR805449.fastq.gz -S zy1.sam`

### 利用snp_file区分品系

`SNPsplit --snp_file ~/ref/snp_mouse/CAST_EiJ/SNPs_CAST_EiJ/CASTEiJ_SNP.txt zy1.sam`


但是我最终得到的结果是这样的，显然我是我预期的结果，我又看了一遍，并没能发现问题出在何处

```
Allele-tagging report
=====================
Processed 50675822 read alignments in total
Reads were unaligned and hence skipped: 12229868 (24.13%)
38445954 reads were unassignable (75.87%)
0 reads were specific for genome 1 (0.00%)
0 reads were specific for genome 2 (0.00%)
39 reads did not contain one of the expected bases at known SNP positions (0.00%)
0 contained conflicting allele-specific SNPs (0.00%)


SNP coverage report
===================
SNP annotation file:	/home1/chenjiayu/ref/snp_mouse/CAST_EiJ/SNPs_CAST_EiJ/CASTEiJ_SNP.txt
SNPs stored in total:	20668296
N-containing reads:	39
non-N:			38445915
total:			50675822
Reads had a deletion of the N-masked position (and were thus dropped):	0 (0.00%)
Of which had multiple deletions of N-masked positions within the same read:	0

Of valid N containing reads,
N was present in the list of known SNPs:	0 (0.00%)
N was not present in the list of SNPs:		165 (100.00%)

Input file:						'zy1.allele_flagged.bam'
Writing unassigned reads to:				'zy1.unassigned.bam'
Writing genome 1-specific reads to:			'zy1.genome1.bam'
Writing genome 2-specific reads to:			'zy1.genome2.bam'


Allele-specific single-end sorting report
=========================================
Read alignments processed in total:		38445954
Reads were unassignable:			38445954 (100.00%)
Reads were specific for genome 1:		0 (0.00%)
Reads were specific for genome 2:		0 (0.00%)
Reads contained conflicting SNP information:	0 (0.00)
```

> 我本以为这是因为合子里大多是母源的RNA所以没有CAST/EiJ的SNP，但我又去看了文章，发现CAST/EiJ就是母本