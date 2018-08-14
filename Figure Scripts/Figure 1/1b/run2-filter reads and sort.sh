
############# mES CH/NP/CY

bsub -q Z-ZQF -e error "samtools view mes_ch_gene.sam | awk '\$2==0{print \$1\"\t\"\$3\"\t\"\$4}\' | sort -k 1 --parallel=10 > mes_ch_gene.tab"
bsub -q Z-ZQF -e error "samtools view mes_ch_trans.sam | awk '\$2==0{print \$1\"\t\"\$3\"\t\"\$4}' | sort -k 1 --parallel=10 > mes_ch_trans.tab"
bsub -q Z-ZQF -e error "samtools view mes_np_gene.sam | awk '\$2==0{print \$1\"\t\"\$3\"\t\"\$4}' | sort -k 1 --parallel=10 > mes_np_gene.tab"
bsub -q Z-ZQF -e error "samtools view mes_np_trans.sam | awk '\$2==0{print \$1\"\t\"\$3\"\t\"\$4}' | sort -k 1 --parallel=10 > mes_np_trans.tab"
bsub -q Z-ZQF -e error "samtools view mes_cy_gene.sam | awk '\$2==0{print \$1\"\t\"\$3\"\t\"\$4}' | sort -k 1 --parallel=10 > mes_cy_gene.tab"
bsub -q Z-ZQF -e error "samtools view mes_cy_trans.sam | awk '\$2==0{print \$1\"\t\"\$3\"\t\"\$4}' | sort -k 1 --parallel=10 > mes_cy_trans.tab"


############# hek CH/NP/CY

bsub -q Z-ZQF -e error "samtools view hek_ch_gene.sam | awk '\$2==0{print \$1\"\t\"\$3\"\t\"\$4}\' | sort -k 1 --parallel=10 > hek_ch_gene.tab"
bsub -q Z-ZQF -e error "samtools view hek_ch_trans.sam | awk '\$2==0{print \$1\"\t\"\$3\"\t\"\$4}' | sort -k 1 --parallel=10 > hek_ch_trans.tab"
bsub -q Z-ZQF -e error "samtools view hek_np_gene.sam | awk '\$2==0{print \$1\"\t\"\$3\"\t\"\$4}' | sort -k 1 --parallel=10 > hek_np_gene.tab"
bsub -q Z-ZQF -e error "samtools view hek_np_trans.sam | awk '\$2==0{print \$1\"\t\"\$3\"\t\"\$4}' | sort -k 1 --parallel=10 > hek_np_trans.tab"
bsub -q Z-ZQF -e error "samtools view hek_cy_gene.sam | awk '\$2==0{print \$1\"\t\"\$3\"\t\"\$4}' | sort -k 1 --parallel=10 > hek_cy_gene.tab"
bsub -q Z-ZQF -e error "samtools view hek_cy_trans.sam | awk '\$2==0{print \$1\"\t\"\$3\"\t\"\$4}' | sort -k 1 --parallel=10 > hek_cy_trans.tab"




