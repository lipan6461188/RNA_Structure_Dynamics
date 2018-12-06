
##### 
##### buils a index with GAS5 and map all reads to it
#####

# load sequence
# import tools
# hg38_seq = tools.readSeq("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/transcriptome.fa")
# tid = 'ENST00000430245'
# print >>open("/Share/home/zhangqf8/lipan/DYNAMIC/human_GAS5/index/Gas5.fa", 'w'), ">%s\n%s" % (tid, hg38_seq[tid])


##### 1. Mapping

function icSHAPE_map()
{
    output_sam=$1
    input_fastq=$2
    genome_index=$3

    THREADS=2

    nohup $(bowtie2 -U $input_fastq -x $genome_index --non-deterministic --time --norc -p $THREADS | awk 'substr($0,1,1)=="@"||$2!=4{print $0}' | samtools view -bh | samtools sort | samtools view -h > $output_sam) &
}

ROOT=/150T/zhangqf/GenomeAnnotation/icSHAPE/hek293
INDEX=/Share/home/zhangqf8/lipan/DYNAMIC/human_GAS5/index/Gas5

cd ~/lipan/DYNAMIC/human_GAS5/mapping/hek_ch

icSHAPE_map hek_ch_D1.sam $ROOT/ch/293chD1.trim.fastq $INDEX
icSHAPE_map hek_ch_D2.sam $ROOT/ch/293chD2.trim.fastq $INDEX
icSHAPE_map hek_ch_N1.sam $ROOT/ch/293chN1.trim.fastq $INDEX
icSHAPE_map hek_ch_N2.sam $ROOT/ch/293chN2.trim.fastq $INDEX

cd ~/lipan/DYNAMIC/human_GAS5/mapping/hek_np

icSHAPE_map hek_np_D1.sam $ROOT/np/293npD1.trim.fastq $INDEX
icSHAPE_map hek_np_D2.sam $ROOT/np/293npD2.trim.fastq $INDEX
icSHAPE_map hek_np_N1.sam $ROOT/np/293npN1.trim.fastq $INDEX
icSHAPE_map hek_np_N2.sam $ROOT/np/293npN2.trim.fastq $INDEX

cd ~/lipan/DYNAMIC/human_GAS5/mapping/hek_cy

icSHAPE_map hek_cy_D1.sam $ROOT/cy/293cyD1.trim.fastq $INDEX
icSHAPE_map hek_cy_D2.sam $ROOT/cy/293cyD2.trim.fastq $INDEX
icSHAPE_map hek_cy_N1.sam $ROOT/cy/293cyN1.trim.fastq $INDEX
icSHAPE_map hek_cy_N2.sam $ROOT/cy/293cyN2.trim.fastq $INDEX


##### 2. Calculate SHAPE score

function calc_icSHAPE()
{
    # root-file-name
    OUT_DIR=$1

    # sam-file-name
    SAM_D1=$2
    SAM_D2=$3
    SAM_N1=$4
    SAM_N2=$5

    # Check true parameters
    if [ $# < 5 ]; then
        echo "Usage: calc_icSHAPE out_dir DMSO_1 DMSO_2 NAI_1 NAI_2";
        return;
    fi
    
    # remove all files in $OUT_DIR
    mkdir -p $OUT_DIR
    rm -rfi $OUT_DIR/*

    #  ===========>  Reset this path  <==========
    # icSHAPE scripts directory
    BIN=~/lipan/icSHAPE/icSHAPE/scripts

    #  ===========>  icSHAPE parameters  <==========
    # ( all parameters are from icSHAPE pipeline default parameters )
    MIN_RPKM=5                      # filter all transcript with RPKM lower than 5
    NORM_METHOD=mean:vigintile2
    HEAD_SKIP=32
    TAIL_SKIP=32
    SCALING_FORM=100
    DIV_FACTOR=10
    SUB_FACTOR=0.25
    WINSOR_METHOD=factor5:scaling1
    INPUT_COVERAGE=100              # a lower cutoff => 100
    TARGET_HIT=0.1                    # a lower cutoff => 1
    HEAD_NULL=5
    TAIL_NULL=30

    # Cluster parameters
    
    # input files
    RT_D1=$OUT_DIR/rt_D1.txt
    RT_D2=$OUT_DIR/rt_D2.txt
    RT_N1=$OUT_DIR/rt_N1.txt
    RT_N2=$OUT_DIR/rt_N2.txt

    # output files
    RPKM_D1=$OUT_DIR/rpkm_D1.txt
    RPKM_D2=$OUT_DIR/rpkm_D2.txt
    RPKM_N1=$OUT_DIR/rpkm_N1.txt
    RPKM_N2=$OUT_DIR/rpkm_N2.txt

    RT_D=$OUT_DIR/rt_D.txt
    RT_N=$OUT_DIR/rt_N.txt

    RT_N_D=$OUT_DIR/rt_norm_D.txt
    RT_N_N=$OUT_DIR/rt_norm_N.txt

    ENRICH=$OUT_DIR/shape.enrich
    SHAPE=$OUT_DIR/shape.out

    ## 1. estimate RPKM

    perl $BIN/estimateRPKM.pl -i $SAM_D1 -o $RPKM_D1

    perl $BIN/estimateRPKM.pl -i $SAM_D2 -o $RPKM_D2

    perl $BIN/estimateRPKM.pl -i $SAM_N1 -o $RPKM_N1

    perl ~/lipan/icSHAPE/icSHAPE/scripts/estimateRPKM.pl -i $SAM_N2 -o $RPKM_N2

    ## 2. calcRT

    perl $BIN/calcRT.pl -i $SAM_D1 -o $RT_D1 -r $RPKM_D1 -c $MIN_RPKM

    perl $BIN/calcRT.pl -i $SAM_D2 -o $RT_D2 -r $RPKM_D2 -c $MIN_RPKM

    perl $BIN/calcRT.pl -i $SAM_N1 -o $RT_N1 -r $RPKM_N1 -c $MIN_RPKM

    perl $BIN/calcRT.pl -i $SAM_N2 -o $RT_N2 -r $RPKM_N2 -c $MIN_RPKM

    ## 3. combineRTreplicates

    perl $BIN/combineRTreplicates.pl -i $RT_D1:$RT_D2 -o $RT_D

    perl $BIN/combineRTreplicates.pl -i $RT_N1:$RT_N2 -o $RT_N

    ## 4. normalizeRTfile

    perl $BIN/normalizeRTfile.pl -i $RT_D -o $RT_N_D -m $NORM_METHOD -d $HEAD_SKIP -l $TAIL_SKIP -f $SCALING_FORM

    perl $BIN/normalizeRTfile.pl -i $RT_N -o $RT_N_N -m $NORM_METHOD -d $HEAD_SKIP -l $TAIL_SKIP -f $SCALING_FORM

    ## 5. calcEnrich

    perl $BIN/calcEnrich.pl -f $RT_N_N -b $RT_N_D -o $ENRICH -w $WINSOR_METHOD -y $DIV_FACTOR -x $SUB_FACTOR -e complex

    ## 6. filterEnrich

    perl $BIN/filterEnrich.pl -i $ENRICH -o $SHAPE -t $INPUT_COVERAGE -T $TARGET_HIT -s $HEAD_NULL -e $TAIL_NULL
}

ROOT=/Share/home/zhangqf8/lipan/DYNAMIC/human_GAS5/mapping

calc_icSHAPE hek_ch $ROOT/hek_ch/hek_ch_D1.sam $ROOT/hek_ch/hek_ch_D2.sam $ROOT/hek_ch/hek_ch_N1.sam $ROOT/hek_ch/hek_ch_N2.sam

calc_icSHAPE hek_np $ROOT/hek_np/hek_np_D1.sam $ROOT/hek_np/hek_np_D2.sam $ROOT/hek_np/hek_np_N1.sam $ROOT/hek_np/hek_np_N2.sam

calc_icSHAPE hek_cy $ROOT/hek_cy/hek_cy_D1.sam $ROOT/hek_cy/hek_cy_D2.sam $ROOT/hek_cy/hek_cy_N1.sam $ROOT/hek_cy/hek_cy_N2.sam


