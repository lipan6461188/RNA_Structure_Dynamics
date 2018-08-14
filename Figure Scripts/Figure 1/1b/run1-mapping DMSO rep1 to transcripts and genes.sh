


hek_trans_db=/Share/home/zhangqf/sunlei/data/human/index/transcriptomeindex/transcriptome
mes_trans_db=/Share/home/zhangqf/sunlei/data/mouse/index/transcriptomeindex/transcriptome

hek_gene_db=/Share/home/zhangqf/sunlei/data/human/index/onlygene_ENindex/genes
mes_gene_db=/Share/home/zhangqf/sunlei/data/mouse/index/onlygene_ENindex/onlygene


hek_ch_dir=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/ch/vitro_transcript_low/293chD1.trimmed.fastq
hek_np_dir=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/np/vivo_gene/293NPD2.trimmed.fastq
hek_cy_dir=/Share/home/zhangqf8/sunlei/data/icshape/final_icSHAPE/293/cy/vivo_gene/293cyD.trimmed.fastq

mes_ch_dir=/Share/home/zhangqf8/sunlei/data/icshape/mes/ch_vivo_gene_low/Chr_D1.trimmed.fastq
mes_np_dir=/Share/home/zhangqf8/sunlei/data/icshape/mes/np_vivo_gene_low/Nuc_D1.trimmed.fastq
mes_cy_dir=/Share/home/zhangqf8/sunlei/data/icshape/mes/cy_vivo_gene_low/Cyt_D1.trimmed.fastq

out_dir=/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution


########  HEK293  ########

#### Ch

bsub -q Z-ZQF -e $out_dir/error -o $out_dir/log -n 20 -J hek_ch_trans \
    bowtie2 \
        -U $hek_ch_dir \
        -S $out_dir/hek_ch_trans.sam \
        -x $hek_trans_db \
        --non-deterministic \
        --time \
        --norc \
        -p 20

bsub -q Z-ZQF -e $out_dir/error -o $out_dir/log -n 20 -J hek_ch_gene \
    bowtie2 \
        -U $hek_ch_dir \
        -S $out_dir/hek_ch_gene.sam \
        -x $hek_gene_db \
        --non-deterministic \
        --time \
        --norc \
        -p 20

#### Np

bsub -q Z-ZQF -e $out_dir/error -o $out_dir/log -n 20 -J hek_np_trans \
    bowtie2 \
        -U $hek_np_dir \
        -S $out_dir/hek_np_trans.sam \
        -x $hek_trans_db \
        --non-deterministic \
        --time \
        --norc \
        -p 20

bsub -q Z-ZQF -e $out_dir/error -o $out_dir/log -n 20 -J hek_np_gene \
    bowtie2 \
        -U $hek_np_dir \
        -S $out_dir/hek_np_gene.sam \
        -x $hek_gene_db \
        --non-deterministic \
        --time \
        --norc \
        -p 20

#### Cy

bsub -q Z-ZQF -e $out_dir/error -o $out_dir/log -n 20 -J hek_cy_trans \
    bowtie2 \
        -U $hek_cy_dir \
        -S $out_dir/hek_cy_trans.sam \
        -x $hek_trans_db \
        --non-deterministic \
        --time \
        --norc \
        -p 20

bsub -q Z-ZQF -e $out_dir/error -o $out_dir/log -n 20 -J hek_cy_gene \
    bowtie2 \
        -U $hek_cy_dir \
        -S $out_dir/hek_cy_gene.sam \
        -x $hek_gene_db \
        --non-deterministic \
        --time \
        --norc \
        -p 20



########  mES  ########


#### Ch

bsub -q Z-ZQF -e $out_dir/error -o $out_dir/log -n 20 -J mes_ch_trans \
    bowtie2 \
        -U $mes_ch_dir \
        -S $out_dir/mes_ch_trans.sam \
        -x $mes_trans_db \
        --non-deterministic \
        --time \
        --norc \
        -p 20

bsub -q Z-ZQF -e $out_dir/error -o $out_dir/log -n 20 -J mes_ch_gene \
    bowtie2 \
        -U $mes_ch_dir \
        -S $out_dir/mes_ch_gene.sam \
        -x $mes_gene_db \
        --non-deterministic \
        --time \
        --norc \
        -p 20

#### Np

bsub -q Z-ZQF -e $out_dir/error -o $out_dir/log -n 20 -J mes_np_trans \
    bowtie2 \
        -U $mes_np_dir \
        -S $out_dir/mes_np_trans.sam \
        -x $mes_trans_db \
        --non-deterministic \
        --time \
        --norc \
        -p 20

bsub -q Z-ZQF -e $out_dir/error -o $out_dir/log -n 20 -J mes_np_gene \
    bowtie2 \
        -U $mes_np_dir \
        -S $out_dir/mes_np_gene.sam \
        -x $mes_gene_db \
        --non-deterministic \
        --time \
        --norc \
        -p 20

#### Cy

bsub -q Z-ZQF -e $out_dir/error -o $out_dir/log -n 20 -J mes_cy_trans \
    bowtie2 \
        -U $mes_cy_dir \
        -S $out_dir/mes_cy_trans.sam \
        -x $mes_trans_db \
        --non-deterministic \
        --time \
        --norc \
        -p 20

bsub -q Z-ZQF -e $out_dir/error -o $out_dir/log -n 20 -J mes_cy_gene \
    bowtie2 \
        -U $mes_cy_dir \
        -S $out_dir/mes_cy_gene.sam \
        -x $mes_gene_db \
        --non-deterministic \
        --time \
        --norc \
        -p 20











