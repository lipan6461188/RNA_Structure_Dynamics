
### KO

samtools view -h KO_1.bam > KO_1.sam
samtools view -h KO_2.bam > KO_2.sam
samtools view -h KO_3.bam > KO_3.sam

sam2dg -in KO_1.sam -out /dev/stdout -mode single | awk '{print $2"\t"$6"\t"$7"\t"$1"\t"1"\t"$3}' > KO_1_raw.bed
sam2dg -in KO_2.sam -out /dev/stdout -mode single | awk '{print $2"\t"$6"\t"$7"\t"$1"\t"1"\t"$3}' > KO_2_raw.bed
sam2dg -in KO_3.sam -out /dev/stdout -mode single | awk '{print $2"\t"$6"\t"$7"\t"$1"\t"1"\t"$3}' > KO_3_raw.bed

cat KO_1_raw.bed KO_2_raw.bed KO_3_raw.bed > KO_all_raw.bed

### WT

samtools view -h WT_1.bam > WT_1.sam
samtools view -h WT_2.bam > WT_2.sam
samtools view -h WT_3.bam > WT_3.sam

sam2dg -in WT_1.sam -out /dev/stdout -mode single | awk '{print $2"\t"$6"\t"$7"\t"$1"\t"1"\t"$3}' > WT_1_raw.bed
sam2dg -in WT_2.sam -out /dev/stdout -mode single | awk '{print $2"\t"$6"\t"$7"\t"$1"\t"1"\t"$3}' > WT_2_raw.bed
sam2dg -in WT_3.sam -out /dev/stdout -mode single | awk '{print $2"\t"$6"\t"$7"\t"$1"\t"1"\t"$3}' > WT_3_raw.bed

#sam2dg -in WT_all.sam -out /dev/stdout -mode single | awk '{print $2"\t"$6"\t"$7"\t"$1"\t"1"\t"$3}' > WT_all.bed

cat WT_1_raw.bed WT_2_raw.bed WT_3_raw.bed > WT_all_raw.bed


### Combine and sample

ROOT_KO = "/Share/home/zhangqf8/lipan/DYNAMIC/LIN28A/KO/"
ROOT_WT = "/Share/home/zhangqf8/lipan/DYNAMIC/LIN28A/WT/"
tools.sample_file(ROOT_KO+"KO_all_raw.bed", ROOT_KO+"KO_all_sample.bed", ratio=0.7572236, skip='')

Piranha -b 50 -s KO_all_sample.bed > KO_all_pir.bed
Piranha -b 50 -s WT_all_raw.bed > WT_all_pir.bed


### Visualization

cp WT_1.sam WT_all.sam
samtools view WT_2.sam >> WT_all.sam
samtools view WT_3.sam >> WT_all.sam
samtools view WT_all.sam | wc -l
samtools view -bh WT_all.sam | samtools sort > WT_all.bam
samtools index WT_all.bam
bam2wig.py -i WT_all.bam -o WT_all -s /150T/zhangqf/GenomeAnnotation/size/mm10.genome.size

cp KO_1.sam KO_all.sam
samtools view KO_2.sam >> KO_all.sam
samtools view KO_3.sam >> KO_all.sam
samtools view KO_all.sam | wc -l
samtools view -bh -s 0.757 KO_all.sam | samtools sort > KO_all.bam
samtools view KO_all.bam | wc -l
samtools index KO_all.bam
bam2wig.py -i KO_all.bam -o KO_all -s /150T/zhangqf/GenomeAnnotation/size/mm10.genome.size


#########################
### RNA-Seq
#########################

INDEX=/150T/zhangqf/GenomeAnnotation/INDEX/Bowtie2/mm10/genome/mm10

bsub -q Z-BNODE -n 16 -e error "bowtie2 -p 16 -x $INDEX -U raw_data/SRR1207291.fastq \
    | awk 'substr(\$0,1,1)==\"@\"||\$2!=4{print \$0}' | samtools view -bh | samtools sort | samtools view -h > WT_1.sam"
bsub -q Z-BNODE -n 16 -e error "bowtie2 -p 16 -x $INDEX -U raw_data/SRR1207293.fastq \
    | awk 'substr(\$0,1,1)==\"@\"||\$2!=4{print \$0}' | samtools view -bh | samtools sort | samtools view -h > WT_2.sam"
bsub -q Z-BNODE -n 16 -e error "bowtie2 -p 16 -x $INDEX -U raw_data/SRR1207295.fastq \
    | awk 'substr(\$0,1,1)==\"@\"||\$2!=4{print \$0}' | samtools view -bh | samtools sort | samtools view -h > KO_1.sam"
bsub -q Z-BNODE -n 16 -e error "bowtie2 -p 16 -x $INDEX -U raw_data/SRR1207297.fastq \
    | awk 'substr(\$0,1,1)==\"@\"||\$2!=4{print \$0}' | samtools view -bh | samtools sort | samtools view -h > KO_2.sam"

samtools view -bh -s 0.139512 KO_1.sam > KO_sample.bam; samtools index KO_sample.bam
samtools view -bh -s 0.158688 WT_1.sam > WT_sample.bam; samtools index WT_sample.bam

bam2wig.py -i KO_sample.bam -o KO_sample -s /150T/zhangqf/GenomeAnnotation/size/mm10.genome.size
bam2wig.py -i WT_sample.bam -o WT_sample -s /150T/zhangqf/GenomeAnnotation/size/mm10.genome.size

samtools view -h WT_sample.bam > WT_sample.sam
samtools view -h KO_sample.bam > KO_sample.sam

sam2dg -in WT_sample.sam -out /dev/stdout -mode single | awk '{print $2"\t"$6"\t"$7"\t"$1"\t"1"\t"$3}' > WT_sample.bed
sam2dg -in KO_sample.sam -out /dev/stdout -mode single | awk '{print $2"\t"$6"\t"$7"\t"$1"\t"1"\t"$3}' > KO_sample.bed

~/usr/piranha-1.2.1/bin/Piranha -b 50 -i 50 -s KO_all_sample.bed KO_RNASeq.bed > KO_peaks.bed
~/usr/piranha-1.2.1/bin/Piranha -b 50 -i 50 -s WT_all_raw.bed WT_RNASeq_1.bed WT_RNASeq_2.bed > WT_peaks.bed

bsub -q Z-ZQF -n 20 -e KO_1_RPKM/error "cufflinks -o KO_1_RPKM/ -p 20 --GTF /150T/zhangqf/GenomeAnnotation/Gencode/GRCm38.gtf KO_1.sam"
bsub -q Z-ZQF -n 20 -e WT_1_RPKM/error "cufflinks -o WT_1_RPKM/ -p 20 --GTF /150T/zhangqf/GenomeAnnotation/Gencode/GRCm38.gtf WT_1.sam"

bsub -q Z-ZQF -n 20 -e KO_2_RPKM/error "cufflinks -o KO_2_RPKM/ -p 20 --GTF /150T/zhangqf/GenomeAnnotation/Gencode/GRCm38.gtf KO_2.sam"
bsub -q Z-ZQF -n 20 -e WT_2_RPKM/error "cufflinks -o WT_2_RPKM/ -p 20 --GTF /150T/zhangqf/GenomeAnnotation/Gencode/GRCm38.gtf WT_2.sam"
