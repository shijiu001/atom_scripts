下载工具：aspera

ENA数据库的数据存储格式：

单端测序：`/vol1/fastq/SRRxxx/00x/SRRxxxxxxx/SRRxxxxxxx.fastq.gz`

双端测序：`/vol1/fastq/SRRxxx/00x/SRRxxxxxxx/SRRxxxxxxx_1.fastq.gz`
`/vol1/fastq/SRRxxx/00x/SRRxxxxxxx/SRRxxxxxxx_2.fastq.gz`

>第三级目录`SRRxxx/`中‘xxx’是前三位数字
>第四级目录`00x`不一定存在
>单端测序结果仅一个文件，双端结果有两个

```
ftp='era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq'

cat SRR_Acc_List.txt |  while read i
 do
       SRR=$(echo ${i:0:6})
       a=$(echo ${i: -1})  #这里:(冒号)后面有个空格，不能漏掉
       ascp -QT -l 300m -P 33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh ${ftp}/${SRR}/${i}/${i}.fastq.gz ~/testANDtry/sra/ #单端测序的话只有
这一个文件
       # ascp -QT -l 300m -P 33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh ${ftp}/${SRR}/${i}/${i}_1.fastq.gz ./
       # ascp -QT -l 300m -P 33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh ${ftp}/${SRR}/${i}/${i}_2.fastq.gz ./
 done

>SRR_Acc_List.txt中写明了需要下载的SRR标号，每行写一个。