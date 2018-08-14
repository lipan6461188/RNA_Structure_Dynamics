
#########
##### BUild a index contains all human mRNA, Each gene has a isoform with longest CDS
#########


input_fasta=smart_hg38_transcriptome.fa
nuc_num=$(grep -v "^>" $input_fasta | wc -m)
rna_num=$(grep "^>" $input_fasta | wc -l)
genomeSAindexNbases=$(echo $nuc_num | awk '{cur_v=log($1)/log(2)/2-1; if(cur_v>14) print 14; else print cur_v}')
genomeChrBinNbits=$(echo -e $nuc_num"\t"$rna_num | awk '{cur_v=log($1/$2)/log(2); if(cur_v>18) print 18; else print cur_v}')
CMD="STAR \
    --runThreadN 20 \
    --runMode genomeGenerate \
    --genomeDir ./ \
    --genomeFastaFiles $input_fasta  \
    --genomeSAindexNbases $genomeSAindexNbases \
    --genomeChrBinNbits $genomeChrBinNbits"
echo $CMD
$CMD



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
    MINLEN:15
}

function MAP
{
    input_fasta=$1
    output_dir=$2

    bsub -q Z-ZQF -n 20 -e $output_dir/error -o $output_dir/log "STAR \
        --runMode alignReads \
        --genomeDir ~/lipan/INDEX/smart_transcriptome/STAR/hg38 \
        --readFilesIn $input_fasta \
        --outFileNamePrefix $output_dir \
        --outReadsUnmapped Fastx \
        --genomeLoad LoadAndKeep \
        --outSAMattributes All \
        --runThreadN 20 \
        --outSAMtype SAM \
        --alignEndsType Local"
}

# readCollapse

in_dir=/Share/home/zhangqf8/sunlei/data/otherdata/transcript_elongation/cell_2014
out_dir=/Share/home/zhangqf8/lipan/DYNAMIC/transcription/cell_2014

fq1=$in_dir/293_40min_18d_rep1_SRR942449.fastq
fq2=$in_dir/293_40min_18d_rep2_SRR942450.fastq
fq3=$in_dir/293_40min_18d_rep3_SRR942451.fastq

readCollapse $fq1 $out_dir/18d_rep1.fastq
readCollapse $fq2 $out_dir/18d_rep2.fastq
readCollapse $fq3 $out_dir/18d_rep3.fastq

# trimmomatic

cd /Share/home/zhangqf8/lipan/DYNAMIC/transcription/cell_2014

trimmomatic 18d_rep1.fastq 18d_rep1.trim.fastq
trimmomatic 18d_rep2.fastq 18d_rep2.trim.fastq
trimmomatic 18d_rep3.fastq 18d_rep3.trim.fastq

# mapping

MAP 18d_rep1.fastq 18S_rep1/
MAP 18d_rep2.fastq 18S_rep2/
MAP 18d_rep3.fastq 18S_rep3/

# CalcFPKM

calcFPKM.py \
    -i Aligned.out.sam \
    -o 18S_rep1.fpkm \
    -m 0 \
    -g 0 \
    -s 0 \
    -r 20

