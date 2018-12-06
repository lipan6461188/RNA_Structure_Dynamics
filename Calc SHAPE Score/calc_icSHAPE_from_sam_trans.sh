
function job_id()
{
    echo $2 | gawk 'match($0, /<([0-9]+)>/, id){print id[1]}'
}

function calc_icSHAPE()
{
    # root-file-name
    OUT_DIR=$1

    # sam-file-name
    SAM_D1=$2
    SAM_D2=$3
    SAM_N1=$4
    SAM_N2=$5

    # remove all files in $OUT_DIR
    mkdir -p $OUT_DIR
    rm -rf $OUT_DIR/*

    # icSHAPE scripts directory
    BIN=~/lipan/icSHAPE/icSHAPE/scripts

    # icSHAPE parameters ( all parameters are from icSHAPE pipeline default parameters )
    MIN_RPKM=0.01              # filter all transcript with RPKM lower than 5
    NORM_METHOD=mean:vigintile2
    HEAD_SKIP=32
    TAIL_SKIP=32
    SCALING_FORM=100
    DIV_FACTOR=10
    SUB_FACTOR=0.25
    WINSOR_METHOD=factor5:scaling1
    INPUT_COVERAGE=100      # a lower cutoff => 100
    TARGET_HIT=1            # a lower cutoff => 1
    HEAD_NULL=5
    TAIL_NULL=30



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

    BSUB_1=$(bsub -n 5 -J Pan_estimateRPKM_D1 -q Z-ZQF -e $OUT_DIR/error -o $OUT_DIR/log \
        "perl $BIN/estimateRPKM.pl -i $SAM_D1 -o $RPKM_D1")

    BSUB_2=$(bsub -n 5 -J Pan_estimateRPKM_D2 -q Z-ZQF -e $OUT_DIR/error -o $OUT_DIR/log \
        "perl $BIN/estimateRPKM.pl -i $SAM_D2 -o $RPKM_D2")

    BSUB_3=$(bsub -n 5 -J Pan_estimateRPKM_N1 -q Z-ZQF -e $OUT_DIR/error -o $OUT_DIR/log \
        "perl $BIN/estimateRPKM.pl -i $SAM_N1 -o $RPKM_N1")

    BSUB_4=$(bsub -n 5 -J Pan_estimateRPKM_N2 -q Z-ZQF -e $OUT_DIR/error -o $OUT_DIR/log \
        "perl ~/lipan/icSHAPE/icSHAPE/scripts/estimateRPKM.pl -i $SAM_N2 -o $RPKM_N2")

    ## 2. calcRT

    BSUB_5=$(bsub -n 10 -J Pan_calcRT_D1 -w "done($(job_id $BSUB_1))" -q Z-ZQF -e $OUT_DIR/error -o $OUT_DIR/log \
        "perl $BIN/calcRT.pl -i $SAM_D1 -o $RT_D1 -r $RPKM_D1 -c $MIN_RPKM")

    BSUB_6=$(bsub -n 10 -J Pan_calcRT_D2 -w "done($(job_id $BSUB_2))" -q Z-ZQF -e $OUT_DIR/error -o $OUT_DIR/log \
        "perl $BIN/calcRT.pl -i $SAM_D2 -o $RT_D2 -r $RPKM_D2 -c $MIN_RPKM")

    BSUB_7=$(bsub -n 10 -J Pan_calcRT_N1 -w "done($(job_id $BSUB_3))" -q Z-ZQF -e $OUT_DIR/error -o $OUT_DIR/log \
        "perl $BIN/calcRT.pl -i $SAM_N1 -o $RT_N1 -r $RPKM_N1 -c $MIN_RPKM")

    BSUB_8=$(bsub -n 10 -J Pan_calcRT_N2 -w "done($(job_id $BSUB_4))" -q Z-ZQF -e $OUT_DIR/error -o $OUT_DIR/log \
        "perl $BIN/calcRT.pl -i $SAM_N2 -o $RT_N2 -r $RPKM_N2 -c $MIN_RPKM")

    ## 3. combineRTreplicates

    BSUB_9=$(bsub -n 10 -J Pan_combine_D -w "done($(job_id $BSUB_5))&&done($(job_id $BSUB_6))" -q Z-ZQF -e $OUT_DIR/error -o $OUT_DIR/log \
        "perl $BIN/combineRTreplicates.pl -i $RT_D1:$RT_D2 -o $RT_D")

    BSUB_10=$(bsub -n 10 -J Pan_combine_N -w "done($(job_id $BSUB_7))&&done($(job_id $BSUB_8))" -q Z-ZQF -e $OUT_DIR/error -o $OUT_DIR/log \
        "perl $BIN/combineRTreplicates.pl -i $RT_N1:$RT_N2 -o $RT_N")

    ## 4. normalizeRTfile

    BSUB_11=$(bsub -n 10 -J Pan_norm_D -w "done($(job_id $BSUB_9))" -q Z-ZQF -e $OUT_DIR/error -o $OUT_DIR/log \
        "perl $BIN/normalizeRTfile.pl -i $RT_D -o $RT_N_D -m $NORM_METHOD -d $HEAD_SKIP -l $TAIL_SKIP -f $SCALING_FORM")

    BSUB_12=$(bsub -n 10 -J Pan_norm_N -w "done($(job_id $BSUB_10))" -q Z-ZQF -e $OUT_DIR/error -o $OUT_DIR/log \
        "perl $BIN/normalizeRTfile.pl -i $RT_N -o $RT_N_N -m $NORM_METHOD -d $HEAD_SKIP -l $TAIL_SKIP -f $SCALING_FORM")

    ## 5. calcEnrich

    BSUB_13=$(bsub -n 10 -J Pan_enrich_N -w "done($(job_id $BSUB_11))&&done($(job_id $BSUB_12))" -q Z-ZQF -e $OUT_DIR/error -o $OUT_DIR/log \
        "perl $BIN/calcEnrich.pl -f $RT_N_N -b $RT_N_D -o $ENRICH -w $WINSOR_METHOD -y $DIV_FACTOR -x $SUB_FACTOR -e complex")

    ## 6. filterEnrich

    BSUB_14=$(bsub -n 10 -J Pan_filter_N -w "done($(job_id $BSUB_13))" -q Z-ZQF -e $OUT_DIR/error -o $OUT_DIR/log \
        "perl $BIN/filterEnrich.pl -i $ENRICH -o $SHAPE -t $INPUT_COVERAGE -T $TARGET_HIT -s $HEAD_NULL -e $TAIL_NULL")
}

OUT_ROOT_Dir=/Share/home/zhangqf8/lipan/DYNAMIC/shape_score


# ************************************************
# ************************************************
#                       mES
# ************************************************
# ************************************************


# ================================================
#   mes ch vitro
# ================================================

OUT_DIR=$OUT_ROOT_Dir/mes_ch_vitro

SAM_D1=/Share/home/zhangqf8/sunlei/data/icshape/mes/ch_vitro_transcript_low/Chr_D1.sam
SAM_D2=/Share/home/zhangqf8/sunlei/data/icshape/mes/ch_vitro_transcript_low/Chr_D2.sam
SAM_N1=/Share/home/zhangqf8/sunlei/data/icshape/mes/ch_vitro_transcript_low/Chr_vitro_N1.sam
SAM_N2=/Share/home/zhangqf8/sunlei/data/icshape/mes/ch_vitro_transcript_low/Chr_vitro_N2.sam

calc_icSHAPE $OUT_DIR $SAM_D1 $SAM_D2 $SAM_N1 $SAM_N2


# ================================================
#   mes np vitro
# ================================================

OUT_DIR=$OUT_ROOT_Dir/mes_np_vitro

SAM_D1=/Share/home/zhangqf8/sunlei/data/icshape/mes/np_vitro_transcript_low/Nuc_D1.sam 
SAM_D2=/Share/home/zhangqf8/sunlei/data/icshape/mes/np_vitro_transcript_low/Nuc_D2.sam 
SAM_N1=/Share/home/zhangqf8/sunlei/data/icshape/mes/np_vitro_transcript_low/Nuc_vitro_N1.sam
SAM_N2=/Share/home/zhangqf8/sunlei/data/icshape/mes/np_vitro_transcript_low/Nuc_vitro_N2.sam

calc_icSHAPE $OUT_DIR $SAM_D1 $SAM_D2 $SAM_N1 $SAM_N2


# ================================================
#   mes cy vitro
# ================================================

OUT_DIR=$OUT_ROOT_Dir/mes_cy_vitro

SAM_D1=/Share/home/zhangqf8/sunlei/data/icshape/mes/cy_vitro_transcript_low/Cyt_D1.sam
SAM_D2=/Share/home/zhangqf8/sunlei/data/icshape/mes/cy_vitro_transcript_low/Cyt_D2.sam
SAM_N1=/Share/home/zhangqf8/sunlei/data/icshape/mes/cy_vitro_transcript_low/Cyt_vitro_N1.sam
SAM_N2=/Share/home/zhangqf8/sunlei/data/icshape/mes/cy_vitro_transcript_low/Cyt_vitro_N2.sam

calc_icSHAPE $OUT_DIR $SAM_D1 $SAM_D2 $SAM_N1 $SAM_N2



# ================================================
#   mes ch vivo
# ================================================

OUT_DIR=$OUT_ROOT_Dir/mes_ch_vivo

SAM_D1=/Share/home/zhangqf8/sunlei/data/icshape/mes/ch_vivo_transcript_low/Chr_D1.sam
SAM_D2=/Share/home/zhangqf8/sunlei/data/icshape/mes/ch_vivo_transcript_low/Chr_D2.sam
SAM_N1=/Share/home/zhangqf8/sunlei/data/icshape/mes/ch_vivo_transcript_low/Chr_vivo_N1.sam
SAM_N2=/Share/home/zhangqf8/sunlei/data/icshape/mes/ch_vivo_transcript_low/Chr_vivo_N2.sam

calc_icSHAPE $OUT_DIR $SAM_D1 $SAM_D2 $SAM_N1 $SAM_N2

bsub -q Z-ZQF perl ~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl -i $OUT_DIR/shape.enrich -o $OUT_DIR/low_shape.out -t 40 -T 0.1 -s 5 -e 30

# ================================================
#   mes np vivo
# ================================================

OUT_DIR=$OUT_ROOT_Dir/mes_np_vivo

SAM_D1=/Share/home/zhangqf8/sunlei/data/icshape/mes/np_vivo_transcript_low/Nuc_D1.sam
SAM_D2=/Share/home/zhangqf8/sunlei/data/icshape/mes/np_vivo_transcript_low/Nuc_D2.sam
SAM_N1=/Share/home/zhangqf8/sunlei/data/icshape/mes/np_vivo_transcript_low/Nuc_vivo_N1.sam
SAM_N2=/Share/home/zhangqf8/sunlei/data/icshape/mes/np_vivo_transcript_low/Nuc_vivo_N2.sam

calc_icSHAPE $OUT_DIR $SAM_D1 $SAM_D2 $SAM_N1 $SAM_N2

bsub -q Z-ZQF perl ~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl -i $OUT_DIR/shape.enrich -o $OUT_DIR/low_shape.out -t 40 -T 0.1 -s 5 -e 30

# ================================================
#   mes cy vivo
# ================================================

OUT_DIR=$OUT_ROOT_Dir/mes_cy_vivo

SAM_D1=/Share/home/zhangqf8/sunlei/data/icshape/mes/cy_vivo_transcript_low/Cyt_D1.sam
SAM_D2=/Share/home/zhangqf8/sunlei/data/icshape/mes/cy_vivo_transcript_low/Cyt_D2.sam
SAM_N1=/Share/home/zhangqf8/sunlei/data/icshape/mes/cy_vivo_transcript_low/Cyt_vivo_N1.sam
SAM_N2=/Share/home/zhangqf8/sunlei/data/icshape/mes/cy_vivo_transcript_low/Cyt_vivo_N2.sam

calc_icSHAPE $OUT_DIR $SAM_D1 $SAM_D2 $SAM_N1 $SAM_N2

bsub -q Z-ZQF perl ~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl -i $OUT_DIR/shape.enrich -o $OUT_DIR/low_shape.out -t 40 -T 0.1 -s 5 -e 30



# ************************************************
# ************************************************
#                  HEK293
# ************************************************
# ************************************************


# ================================================
#   hek ch vitro
# ================================================

OUT_DIR=$OUT_ROOT_Dir/hek_ch_vitro

SAM_D1=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/ch/vitro_transcript/293chD1.sam
SAM_D2=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/ch/vitro_transcript/293chD2.sam
SAM_N1=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/ch/vitro_transcript/293chT1.sam
SAM_N2=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/ch/vitro_transcript/293chT2.sam

calc_icSHAPE $OUT_DIR $SAM_D1 $SAM_D2 $SAM_N1 $SAM_N2


# ================================================
#   hek np vitro
# ================================================

OUT_DIR=$OUT_ROOT_Dir/hek_np_vitro

SAM_D1=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/np/vitro_transcript/293NPD1.sam
SAM_D2=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/np/vitro_transcript/293NPD2.sam
SAM_N1=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/np/vitro_transcript/293NPT.sam
SAM_N2=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/np/vitro_transcript/293NPT3.sam

calc_icSHAPE $OUT_DIR $SAM_D1 $SAM_D2 $SAM_N1 $SAM_N2


# ================================================
#   hek cy vitro
# ================================================

OUT_DIR=$OUT_ROOT_Dir/hek_cy_vitro

SAM_D1=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/cy/vitro_transcript_low/293cyD.sam
SAM_D2=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/cy/vitro_transcript_low/293cyD3.sam
SAM_N1=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/cy/vitro_transcript_low/293cyT1.sam
SAM_N2=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/cy/vitro_transcript_low/293cyT2.sam

calc_icSHAPE $OUT_DIR $SAM_D1 $SAM_D2 $SAM_N1 $SAM_N2


# ================================================
#   hek ch vivo
# ================================================

OUT_DIR=$OUT_ROOT_Dir/hek_ch_vivo

SAM_D1=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/ch/vivo_transcript/293chD1.sam
SAM_D2=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/ch/vivo_transcript/293chD2.sam
SAM_N1=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/ch/vivo_transcript/293chN.sam
SAM_N2=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/ch/vivo_transcript/293chN3.sam

calc_icSHAPE $OUT_DIR $SAM_D1 $SAM_D2 $SAM_N1 $SAM_N2

bsub -q Z-ZQF perl ~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl -i $OUT_DIR/shape.enrich -o $OUT_DIR/low_shape.out -t 40 -T 0.1 -s 5 -e 30


# ================================================
#   hek np vivo
# ================================================

OUT_DIR=$OUT_ROOT_Dir/hek_np_vivo

SAM_D1=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/np/vivo_transcript/293NPD1.sam
SAM_D2=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/np/vivo_transcript/293NPD2.sam
SAM_N1=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/np/vivo_transcript/293NPN1.sam
SAM_N2=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/np/vivo_transcript/293NPN2.sam

calc_icSHAPE $OUT_DIR $SAM_D1 $SAM_D2 $SAM_N1 $SAM_N2

bsub -q Z-ZQF perl ~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl -i $OUT_DIR/shape.enrich -o $OUT_DIR/low_shape.out -t 40 -T 0.1 -s 5 -e 30


# ================================================
#   hek cy vivo
# ================================================

OUT_DIR=$OUT_ROOT_Dir/hek_cy_vivo

SAM_D1=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/cy/vivo_transcript/293cyD.sam
SAM_D2=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/cy/vivo_transcript/293cyD3.sam
SAM_N1=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/cy/vivo_transcript/293cyN.sam
SAM_N2=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/cy/vivo_transcript/293cyN3.sam

calc_icSHAPE $OUT_DIR $SAM_D1 $SAM_D2 $SAM_N1 $SAM_N2

bsub -q Z-ZQF perl ~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl -i $OUT_DIR/shape.enrich -o $OUT_DIR/low_shape.out -t 40 -T 0.1 -s 5 -e 30


# ================================================
#   High cutoff
# ================================================


perl ~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl -i mes_ch_vivo/shape.enrich -o /dev/stdout -t 200 -T 2 -s 5 -e 30 | awk '$3>5{print $0}' > mes_ch_vivo/high_shape.out
perl ~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl -i mes_np_vivo/shape.enrich -o /dev/stdout -t 200 -T 2 -s 5 -e 30 | awk '$3>5{print $0}' > mes_np_vivo/high_shape.out
perl ~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl -i mes_cy_vivo/shape.enrich -o /dev/stdout -t 200 -T 2 -s 5 -e 30 | awk '$3>5{print $0}' > mes_cy_vivo/high_shape.out
perl ~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl -i mes_wc_vivo/shape.enrich -o /dev/stdout -t 200 -T 2 -s 5 -e 30 | awk '$3>5{print $0}' > mes_wc_vivo/high_shape.out

perl ~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl -i hek_ch_vivo/shape.enrich -o /dev/stdout -t 200 -T 2 -s 5 -e 30 | awk '$3>5{print $0}' > hek_ch_vivo/high_shape.out
perl ~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl -i hek_np_vivo/shape.enrich -o /dev/stdout -t 200 -T 2 -s 5 -e 30 | awk '$3>5{print $0}' > hek_np_vivo/high_shape.out
perl ~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl -i hek_cy_vivo/shape.enrich -o /dev/stdout -t 200 -T 2 -s 5 -e 30 | awk '$3>5{print $0}' > hek_cy_vivo/high_shape.out
perl ~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl -i hek_wc_vivo/shape.enrich -o /dev/stdout -t 200 -T 2 -s 5 -e 30 | awk '$3>5{print $0}' > hek_wc_vivo/high_shape.out









