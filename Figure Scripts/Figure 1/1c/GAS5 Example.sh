

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



##### 3. Visualization

## show SHAPE map with VARNA

import visual, structure, tools

seq_gas5 = tools.readSeq("/Share/home/zhangqf8/lipan/DYNAMIC/human_GAS5/index/Gas5.fa")['ENST00000430245']
shape_gas5_ch = tools.loadicSHAPE("/Share/home/zhangqf8/lipan/DYNAMIC/human_GAS5/shape_score/hek_ch/shape.out")['ENST00000430245']
shape_gas5_np = tools.loadicSHAPE("/Share/home/zhangqf8/lipan/DYNAMIC/human_GAS5/shape_score/hek_np/shape.out")['ENST00000430245']
shape_gas5_cy = tools.loadicSHAPE("/Share/home/zhangqf8/lipan/DYNAMIC/human_GAS5/shape_score/hek_cy/shape.out")['ENST00000430245']

s = seq_gas5.find("CCAGTGGTC")
e = s + 25

shape_ch = shape_gas5_ch[s:e]
shape_np = shape_gas5_np[s:e]
shape_cy = shape_gas5_cy[s:e]

ss = structure.predictStructure(seq_gas5[s:e])

visual.Plot_RNAStructure_Shape(seq_gas5[s:e], ss, shape_ch, mode='heatmap')
visual.Plot_RNAStructure_Shape(seq_gas5[s:e], ss, shape_np, mode='heatmap')
visual.Plot_RNAStructure_Shape(seq_gas5[s:e], ss, shape_cy, mode='heatmap')

ROC_ch = tools.calc_shape_ROC(ss, shape_ch, step=0.05)
ROC_np = tools.calc_shape_ROC(ss, shape_np, step=0.05)
ROC_cy = tools.calc_shape_ROC(ss, shape_cy, step=0.05)

tools.plot_ROC(ROC_ch); tools.plt.show()
tools.plot_ROC(ROC_np)
tools.plot_ROC(ROC_cy); tools.plt.show()

tools.calc_AUC(ROC_ch)
tools.calc_AUC(ROC_np)
tools.calc_AUC(ROC_cy)


##### 3. Visualization with NAI RT

## show SHAPE map with VARNA

def read_normalized_RT(normRTFile):
    normRT = {}
    for line in open(normRTFile):
        if line[0] == '#':
            continue
        if 'RTstop' in line:
            data = line.strip().split()
            assert len(data[5:]) == int(data[1])
            normRT[ data[0] ] = [ float(it) for it in data[5:] ]
    return normRT

rt_ch = read_normalized_RT("/Share/home/zhangqf8/lipan/DYNAMIC/human_GAS5/shape_score/hek_ch/rt_norm_N.txt")['ENST00000430245']
rt_np = read_normalized_RT("/Share/home/zhangqf8/lipan/DYNAMIC/human_GAS5/shape_score/hek_np/rt_norm_N.txt")['ENST00000430245']
rt_cy = read_normalized_RT("/Share/home/zhangqf8/lipan/DYNAMIC/human_GAS5/shape_score/hek_cy/rt_norm_N.txt")['ENST00000430245']

shape_ch = rt_ch[s:e]
shape_np = rt_np[s:e]
shape_cy = rt_cy[s:e]

visual.Plot_RNAStructure_Shape(seq_gas5[s:e], ss, shape_ch, mode='heatmap')
visual.Plot_RNAStructure_Shape(seq_gas5[s:e], ss, shape_np, mode='heatmap')
visual.Plot_RNAStructure_Shape(seq_gas5[s:e], ss, shape_cy, mode='heatmap')




