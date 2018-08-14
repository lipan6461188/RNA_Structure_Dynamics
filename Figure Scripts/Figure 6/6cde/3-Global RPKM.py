
from icSHAPE import *
import tools, ParseTrans


def read_FPKM(inFile, min_fpkm):
    FPKM = {}
    for line in open(inFile):
        if not line.startswith("ENS"):
            continue
        data = line.strip().split()
        if data[12] == 'OK':
            fpkm = float(data[9])
            if fpkm > min_fpkm:
                FPKM[data[0]] = fpkm
    return FPKM

ROOT = "/Share/home/zhangqf8/lipan/DYNAMIC/LIN28A/RNA-Seq/"
RPKM_WT_1 = read_FPKM(ROOT+"WT_1_RPKM/genes.fpkm_tracking", 5); print len(WT_1)
RPKM_WT_2 = read_FPKM(ROOT+"WT_2_RPKM/genes.fpkm_tracking", 5); print len(WT_2)
RPKM_KO_1 = read_FPKM(ROOT+"KO_1_RPKM/genes.fpkm_tracking", 5); print len(KO_1)
RPKM_KO_2 = read_FPKM(ROOT+"KO_2_RPKM/genes.fpkm_tracking", 5); print len(KO_2)

def generate_RPKM(WT_1, WT_2):
    WT = []
    for gid in set(WT_1) & set(WT_2):
        WT.append((np.log2(WT_1[gid]), np.log2(WT_2[gid])))
    KO_df = pd.DataFrame(WT, columns=['log2_rep1', 'log2_rep2'])
    sns.jointplot(data=KO_df, x='log2_rep1', y='log2_rep2')
    plt.tight_layout()
    plt.show()


generate_RPKM(RPKM_WT_1, RPKM_WT_2)
generate_RPKM(RPKM_KO_1, RPKM_KO_2)

def WTKO_RPKM(WT, KO):
    COM = []
    for gid in set(WT) & set(KO):
        COM.append((np.log2(WT[gid]), np.log2(KO[gid])))
    COM_df = pd.DataFrame(COM, columns=['log2_WT', 'log2_KO'])
    sns.jointplot(data=COM_df, x='log2_WT', y='log2_KO')
    plt.tight_layout()
    plt.show()

WTKO_RPKM(RPKM_WT_1, RPKM_KO_1)


Gene_RPKM = []
for locus in m6A_KO_peaks:
    chr_id, strand = locus
    for region in m6A_KO_peaks[locus]:
        gene_list = mm10_parser.genomeCoor2geneCoor(chr_id, region[0], region[1], strand)
        for gene_it in gene_list:
            gene = gene_it[3]
            try:
                WT = np.log2(np.mean( (RPKM_WT_1[gene], RPKM_WT_2[gene]) ))
                KO = np.log2(np.mean( (RPKM_KO_1[gene], RPKM_KO_2[gene]) ))
                Gene_RPKM.append( (WT, KO) )
            except KeyError, TypeError:
                continue

Gene_RPKM_df = pd.DataFrame(data=Gene_RPKM, columns=['log2_WT', 'log2_KO'])
sns.jointplot(data=Gene_RPKM_df, x='log2_WT', y='log2_KO')
plt.show()

