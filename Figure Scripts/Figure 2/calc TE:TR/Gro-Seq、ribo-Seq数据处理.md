<!--下面是全局的风格-->
<style>
	h1 { color:purple }
	h2 { background:yellow; color:green }
	a { color:#333 } 
	a:hover{ color:#F00; }
	h4 {color: darkblue;}
	a:link{
		text-decoration:none;
	}
</style>

<link rel="stylesheet" href="http://olddriver.website/bootstrap/css/button.css">

<center>
# Gro-Seq、ribo-Seq数据处理
</center>


## 1. Gro-Seq数据处理

<a href="http://science.sciencemag.org/content/sci/suppl/2008/12/04/1162228.DC1/Core.SOM.pdf"><font size=5><b>Science 2008 Supplementary</b></font></a>

<center>
<img src="https://ws1.sinaimg.cn/large/006tNbRwly1fgm8v9htnsj30ys0p4q6o.jpg">
</center>


<a href="https://elifesciences.org/lookup/doi/10.7554/eLife.02407"><font size=5 color='green'><b>Elife 2014</b></font></a>
<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48895"><font size=5 color='blue'><b> GSE48895 </b></font></a>
<center>

<table>
<tr>
	<td><img src="https://ws2.sinaimg.cn/large/006tNc79ly1fgoe4zaozcj307907o3yu.jpg"></td>
	<td><img src="https://ws3.sinaimg.cn/large/006tNc79ly1fgodykydzsj30jj09ngn4.jpg"></td>
</tr>
</table>
</center>


Gro-Seq的数据是用抗体把带有Br-U的RNA抓出测序，从而达到测量转录速率的目的。

对测序数据的处理分为下面的三步即可，可以适当的切除5'端的序列，序列可以map回基因组上，也可以map回转录组上，具体看情况而定

```bash
function readCollapse
{
    ~/lipan/icSHAPE/icSHAPE/scripts/readCollapse.pl \
        -U $1 \
        -o $2
}

function trimmomatic
{
    java -jar ~/lipan/icSHAPE/icSHAPE/bin/trimmomatic-0.30.jar SE \
    -threads 8 \
    -phred33 \
    $1 \
    $2 \
    ILLUMINACLIP:/Share/home/zhangqf8/lipan/icSHAPE/icSHAPE/data/adapter/TruSeq2-PE.fa:2:30:7 \
    HEADCROP:5 \
    SLIDINGWINDOW:4:15 \
    MINLEN:30
}

function MAP
{
    STAR --runMode alignReads \
        --genomeDir ~/lipan/INDEX/smart_transcriptome/STAR/mm10 \
        --readFilesIn $1 \
        --outFileNamePrefix $2 \
        --outReadsUnmapped Fastx \
        --genomeLoad LoadAndKeep \
        --outSAMattributes All \
        --runThreadN 8 \
        --outSAMtype SAM \
        --quantMode TranscriptomeSAM GeneCounts
}
```

* <b>Map回基因组</b>：如果测到了大量的intron，希望看polII的速度；
* <b>Map回转录组</b>：看单位时间mature RNA的产量。

> 这里的转录组最好是每个基因只包含单个最长的转录本：这里推荐构建一个这样的转录组数据库，对于每一个mRNA，选CDS最长的转录本，如果CDS一样长，选UTR最长的转录本。这个数据库已经构建在目录`~/lipan/INDEX/smart_transcriptome/`下，其中还加上了rRNA序列。

<b>计算FPKM</b>：map结束以后，如果想看mature RNA的产量，下一步就是计算每一个基因的表达量：

```bash
~/lipan/python_utils/calcFPKM.py \
    -i SRR315594.Aligned.out.sam \
    -o SRR315594.fpkm \
    -m 1/0 \ #考虑multi-map ??
    -g 1/0 \ #考虑gapped-map reads ??
    -s 1/0 \ #考虑reversed complementary mapped reads ??
    -r 20 #只输出map上的reads数大于20的基因的信息
```

## 2. ribo-Seq数据处理

<center>
<img src="https://ws4.sinaimg.cn/large/006tNc79ly1fgof234b1rj30gq0c0dh2.jpg">
</center>

<table>
<tr>
	<th><a href="http://www.sciencedirect.com/science/article/pii/S0092867417303586"><font size=5 color='green'><b>Cell2017</b></font></a><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96643"><font size=5 color='blue'><b> GSE96643 </b></font></a></th>
	<th><a href="http://www.sciencedirect.com/science/article/pii/S0092867411011925"><font size=5 color='green'><b>Cell2011</b></font></a><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30839"><font size=5 color='blue'><b> GSE30839 </b></font></a>
</th>
</tr>
<tr>
	<td><img src="https://ws1.sinaimg.cn/large/006tNc79ly1fgof66q4ahj30an0aqq3v.jpg"></td>
	<td><img src="https://ws2.sinaimg.cn/large/006tNc79ly1fgof6j5yroj30af0afgm9.jpg"></td>
</tr>
</table>
</center>

Ribo-Seq可以看Ribosome翻译的速度，看Ribosome聚集的位置，Ribo-Seq为了避免RNA表达量的影响，需要做一个平行的RNA-Seq来做对照。Cell2017建议如下处理得到Translation Efficiency:

<b>RNA-Seq and Ribo-Seq Analysis</b>

Sequenced reads were aligned to a reference set of human curated protein-coding transcripts (plus the five human rRNA transcripts) using Bowtie (Langmead et al., 2009). This reference set of transcripts was based on Ensembl’s gene annotations (release 69). ==For genes with multiple isoforms, the one with longest coding DNA sequence (CDS) region and, in case not unique, the one with longest UTRs among the ones with the longest CDS, was selected to represent the gene==. Only uniquely mapped reads were used in subsequent analyses. RNA expression levels and ribosome occupancy were estimated by calculating reads per kilobase of mRNA per million reads (RPKM) per transcript, taking into account either all reads that map to the transcript (for estimation of RNA levels using RNA-seq data) or only those mapping to its CDS (for estimation of ribosome occupancy). ==In estimation of ribosome occupancy in CDS, 5′ ends of reads were offset 12 nucleotides to the 3′ direction to match the P-site location of ribosome (Ingolia et al., 2009).== Translation efficiency (TE) was estimated by the (log2) ratio between ribosome occupancy and mRNA level. Only genes with at least 20 reads in RNA-seq and Ribo-seq samples were included in analyses.


```bash
function readCollapse
{
    ~/lipan/icSHAPE/icSHAPE/scripts/readCollapse.pl \
        -U $1 \
        -o $2
}

function MAP
{
    STAR \
        --runMode alignReads \
        --genomeDir ~/lipan/INDEX/smart_transcriptome/STAR/mm10 \
        --readFilesIn $1 \
        --outFileNamePrefix $2 \
        --outReadsUnmapped Fastx \
        --genomeLoad NoSharedMemory \
        --outSAMattributes All \
        --runThreadN $3 \
        --outSAMtype SAM \
        --outFilterMismatchNmax 2 \
        --clip3pNbases $4
}
```

<b>Map</b>: Ribo-Seq的数据应该map到转录组上，5'端除非有barcode，一般不切除，注意ribo-Seq reads的3'端会有长短不一的adapter，在无法知道adapter的情况下，应该直接切除一段碱基，或者用--clip3pNbases参数挂起一段（~20），这样可以提高Map率；而对应的RNA-Seq样本一般不需要或者小(~5)。

<b>FPKM与Occupancy</b>：

```bash
#计算FPKM
~/lipan/python_utils/calcFPKM.py \
    -i SRR315594.Aligned.out.sam \
    -o SRR315594.fpkm \
    -m 1/0 \ #考虑multi-map ??
    -g 1/0 \ #考虑gapped-map reads ??
    -s 1/0 \ #考虑reversed complementary mapped reads ??
    -r 20 #只输出map上的reads数大于20的基因的信息

# 计算Occupancy
~/lipan/python_utils/calcFPKM.py \
    -i SRR315594.Aligned.out.sam \
    -o SRR315594.occu.fpkm \
    --ribo_occup scripts/mouse_rna_cds_info.txt \
    --offset 12 \ # reads向3'端偏移多少碱基
    -m 1/0 \ #考虑multi-map ??
    -g 1/0 \ #考虑gapped-map reads ??
    -s 1/0 \ #考虑reversed complementary mapped reads ??
    -r 20 #只输出map上的reads数大于20的基因的信息
```

Occupancy只统计CDS区域的reads，所以需要一个文件说明每一个基因的CDS区间(mouse\_rna\_cds\_info.txt)，这个文件可以用下面的脚本制作：

```python
from ParseTrans import *
import os

genomeCoorBedFileName = os.environ.get('HOME') + '/lipan/DYNAMIC/GTF/mm10.genomeCoor.bed'
parseTrans = ParseTransClass(genomeCoorBedFile = genomeCoorBedFileName)

for transID in parseTrans.TransIParser:
    trans = parseTrans.getTransFeature(transID)
    geneID = trans['gene_id']
    length = trans['trans_len']
    cds_start = trans['cds_start']
    cds_end = trans['cds_end']
    print "%s\t%d\t%d\t%d" % (transID+'_'+geneID, length, cds_start, cds_end)
```




