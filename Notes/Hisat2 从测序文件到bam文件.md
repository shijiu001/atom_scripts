# Hisat2 从测序文件到bam文件

> 应该说，并不是从一个真正的测序文件开始进行质量控制的

>文章：Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown

>中文教程：http://www.360doc.com/content/19/1224/14/68068867_881789468.shtml
http://www.360doc.com/content/18/0715/20/19913717_770622175.shtml

### 下载测序结果及软件安装

`wget ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol/chrX_data.tar.gz` 下载数据

软件的下载有一部分用PC下载之后上传到服务器

`scp` 用scp上传到服务器

`wget` 或者直接下载到服务器

对于已经编译好的软件直接解压即可使用，下载源码的软件需经过编译

`tar -x`用于解压，`-z`用于解压'gzip'文件，`-j`用于解压'bzip2'文件，`-v`显示操作过程，`-f`指定解压文件（虽然我不知道不指定解压文件要怎么实现解压）

`make`/`make install`用于编译

配置环境变量

将所需的软件用`mv`命令移动到`$/home/bin`目录下，并将此目录配置到环境变量

`export PATH=$HOME/bin:$PATH`

---

### 回贴

首先，关于回贴算法，自我按照说明书做到'bam'文件之后我就想搞懂回贴是怎样实现的，但是直到目前，我也只明白了回贴和比对之间的区别，以及建立index的原因

按照protocol，用hisat2进行回贴

`hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran\`

` -1chrX_data/samples/ERR188044_chrX_1.fastq.gz\` 

` -2chrX_data/samples/ERR188044_chrX_2.fastq.gz -S ERR188044_chrX.sam`

写成循环便是

`for i in {188044,188104,188234,188245,188257,188273,188337,188383,188401,188428,188454,204916};\`

`do hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR${i}_chrX_1.fastq.gz -2 chrX_data/samples/ERR${i}_chrX_2.fastq.gz -S ERR${i}_chrX.sam; done`

如此就得到了'sam'文件

利用samtools对得到的sam文件进行排序并将其转化为bam文件

`for i in {188044,188104,188234,188245,188257,188273,188337,188383,188401,188428,188454,204916};do samtools sort -@ 8 -o ERR${i}_chrX.bam ERR${i}_chrX.sam; done`

> sam文件开头的若干行以@开始，会对这个sam文件做一些说明
>
> sam文件的介绍https://en.wikipedia.org/wiki/SAM_(file_format)

---

### 一些我想要请教的问题

下载测序文件和索引文件时，下载速度极慢，有什么方法可以提高下载速度？

关于如何将一个任务挂到后台，在之后如何将他停止

关于索引文件的结构

推荐的目录结构

---

#### 我做的一些其他事情

更改配色方案

安装了conda简化之后的软件安装