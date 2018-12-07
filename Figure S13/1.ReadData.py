
"""

#### Prepare .genomeCoor files

def rem_genomeCoor_version(raw_genomeCoor, new_genomeCoor):
    OUT = open(new_genomeCoor, 'w')
    for line in open(raw_genomeCoor):
        data = line.strip().split()
        data[4] = ".".join(data[4].split('.')[:-1])
        data[5] = ".".join(data[5].split('.')[:-1])
        print >>OUT, "\t".join(data)
    OUT.close()

raw_genomeCoor = "/150T/zhangqf/GenomeAnnotation/Gencode/hg38.genomeCoor.bed"
new_genomeCoor = "/tmp/hg38.genomeCoor.bed"
rem_genomeCoor_version(raw_genomeCoor, new_genomeCoor)

raw_genomeCoor = "/150T/zhangqf/GenomeAnnotation/Gencode/mm10.genomeCoor.bed"
new_genomeCoor = "/tmp/mm10.genomeCoor.bed"
rem_genomeCoor_version(raw_genomeCoor, new_genomeCoor)

"""

import tools
from icSHAPE import *
from ParseTrans import *

####################
### Font
####################

font = {'family': 'normal', 'weight': 'bold','size': 10}
matplotlib.rc('font', **font)

def read_half_life_to_Trans(inFile, mouse_geneInfo):
    half_life = pd.read_csv(inFile, sep="\t", header=None)
    half_life_dict = dict( zip(half_life.iloc[:,0], half_life.iloc[:,1]) )
    new_half_life_dict = dict()
    for gene_id in half_life_dict:
        try:
            trans_id_list = mouse_geneInfo[gene_id]['transcript']
        except:
            continue
        for trans_id in trans_id_list:
            new_half_life_dict[trans_id] = half_life_dict[gene_id]
    #return new_half_life_dict
    return pd.DataFrame(new_half_life_dict.items(), columns=['trans_id', 'HL'])

def read_half_life(inFile):
    half_life = pd.read_csv(inFile, sep="\t", header=None)
    return dict( zip(half_life.iloc[:,0], half_life.iloc[:,1]) )

def filter_common_dataSet(set_1, set_2):
    common_set = []
    for trans_id in set(set_1.keys()) & set(set_2.keys()):
        common_set.append( (trans_id, set_1[trans_id], set_2[trans_id]) )
    return common_set

def divide(list_1, list_2):
    list_1 = list(list_1)
    list_2 = list(list_2)
    assert len(list_1) == len(list_2)
    
    new_list = []
    for it1, it2 in zip(list_1, list_2):
        new_list.append( 1.0*it1/it2 )
    
    return new_list

def mean_list(list_1, list_2, list_3):
    list_1,list_2,list_3 = list(list_1),list(list_2),list(list_3)
    
    assert len(list_1) == len(list_2) == len(list_3)
    mean_list = []
    
    for v1,v2,v3 in zip(list_1, list_2, list_3):
        mean_list.append( 1.0*(v1+v2+v3)/3 )
    
    return mean_list

def calc_GINI(icSHAPE):
    import numpy as np
    
    gini = {'ch':{}, 'np':{}, 'cy':{}, 'wc':{}}
    
    for comp in ('ch', 'np', 'cy', 'wc'):
        for tid in icSHAPE[comp]:
            shape = icSHAPE[comp][tid]
            #gini_list = computeGINI_splice(shape, valid_cutoff=15, splice=50)
            gini_index = computeGINI(shape, valid_cutoff=15)
            #if len(gini_list) >= 1:
            if gini_index != -1:
                gini[comp][tid] = gini_index
    
    return gini

def calc_start_codon_gini(icSHAPE, Parser):
    import numpy as np
    
    gini = {'ch':{}, 'np':{}, 'cy':{}, 'wc':{}}
    
    for comp in ('ch', 'np', 'cy', 'wc'):
        for tid in icSHAPE[comp]:
            try:
                feature = Parser.getTransFeature(tid, verbose=False)
            except KeyError:
                continue
            gene_type = feature['gene_type']
            cds_start = feature['cds_start']
            if gene_type in('mRNA', 'protein_coding'):
                shape = icSHAPE[comp][tid][cds_start-10:cds_start+10]
                gini_index = computeGINI(shape, valid_cutoff=15)
                if gini_index != -1:
                    gini[comp][tid] = gini_index
    
    return gini

def calc_5UTR_gini(icSHAPE, Parser):
    import numpy as np
    
    gini = {'ch':{}, 'np':{}, 'cy':{}, 'wc':{}}
    
    for comp in ('ch', 'np', 'cy', 'wc'):
        for tid in icSHAPE[comp]:
            try:
                feature = Parser.getTransFeature(tid, verbose=False)
            except KeyError:
                continue
            gene_type = feature['gene_type']
            cds_start = feature['cds_start']
            if gene_type in('mRNA', 'protein_coding'):
                shape = icSHAPE[comp][tid][:cds_start]
                gini_index = computeGINI(shape, valid_cutoff=15)
                if gini_index != -1:
                    gini[comp][tid] = gini_index
    
    return gini

def calc_3UTR_gini(icSHAPE, Parser):
    import numpy as np
    
    gini = {'ch':{}, 'np':{}, 'cy':{}, 'wc':{}}
    
    for comp in ('ch', 'np', 'cy', 'wc'):
        for tid in icSHAPE[comp]:
            try:
                feature = Parser.getTransFeature(tid, verbose=False)
            except KeyError:
                continue
            gene_type = feature['gene_type']
            cds_end = feature['cds_end']
            if gene_type in('mRNA', 'protein_coding'):
                shape = icSHAPE[comp][tid][cds_end:]
                gini_index = computeGINI(shape, valid_cutoff=15)
                if gini_index != -1:
                    gini[comp][tid] = gini_index
    
    return gini

def tid_2_geneNames(tid_list, Parser):
    geneNames = []
    for tid in tid_list:
        geneName = Parser.getTransFeature(tid)['gene_name']
        geneNames.append(geneName)
    return geneNames

def remove_hist_genes(tid_dict, Parser):
    tid_list = tid_dict.keys()
    raw_len = len(tid_list)
    
    ## Count
    key_not_found = 0
    hist_gene = 0
    
    for tid in tid_list:
        try:
            gene_name = Parser.getTransFeature(tid, verbose=False)['gene_name']
        except KeyError:
            del tid_dict[tid]
            key_not_found += 1
            continue
        if gene_name[:4].upper() == 'HIST':
            del tid_dict[tid]
            hist_gene += 1
            continue
    
    print "total: %s key_not_found: %s, hist_gene: %s" % (raw_len, key_not_found, hist_gene)

def scatter_plot(dict_x, dict_y, color="", xlabel="", ylabel="", title="", filt_list=None, xlim=None, ylim=None, Parser=None, xlog=False, ylog=False, kind='reg'):
    common_set = set(dict_x) & set(dict_y)
    if Parser: common_set = common_set & set(Parser.TransIParser)
    if filt_list: common_set = common_set & set(filt_list)
    if not color: color = "black"
    
    df = []
    for tid in common_set:
        if Parser:
            gene_type = Parser.getTransFeature(tid)['gene_type']
            if gene_type not in ('mRNA', 'protein_coding'):
                continue
        x = np.log2(dict_x[tid]) if xlog else dict_x[tid]
        y = np.log2(dict_y[tid]) if ylog else dict_y[tid]
        df.append( (tid, x, y) )
    
    df = pd.DataFrame(df, columns=['tid', 'x', 'y'])
    print "data size: ", len(df)
    ax = sns.jointplot(data=df, x='x', y='y', kind=kind, color=color, xlim=xlim, ylim=ylim)
    
    ax.set_axis_labels(xlabel, ylabel)
    if title:
        plt.title(title)
    
    return ax



####################
###### Mouse Load Functions
####################


def load_mouse_parser():
    mm10_parseTrans = ParseTransClass(genomeCoorBedFile = "/Share/home/zhangqf8/lipan/DYNAMIC/GTF/mm10.genomeCoor.bed")
    mouse_geneInfo = mm10_parseTrans.getGeneInfo()
    return mm10_parseTrans, mouse_geneInfo

def load_mouse_shape(mm10_parseTrans, rev_hist=False):
    icSHAPE_ROOT_DIR = '/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/%s/shape.out'
    ch_icshape_file = icSHAPE_ROOT_DIR % ("mes_ch_vivo", )
    np_icshape_file = icSHAPE_ROOT_DIR % ("mes_np_vivo", )
    cy_icshape_file = icSHAPE_ROOT_DIR % ("mes_cy_vivo", )
    wc_icshape_file = icSHAPE_ROOT_DIR % ("mes_wc_vivo", )
    
    ch_icshape = tools.loadicSHAPE(ch_icshape_file)
    np_icshape = tools.loadicSHAPE(np_icshape_file)
    cy_icshape = tools.loadicSHAPE(cy_icshape_file)
    wc_icshape = tools.loadicSHAPE(wc_icshape_file)
    
    if rev_hist:
        remove_hist_genes(ch_icshape, mm10_parseTrans)
        remove_hist_genes(np_icshape, mm10_parseTrans)
        remove_hist_genes(cy_icshape, mm10_parseTrans)
        remove_hist_genes(wc_icshape, mm10_parseTrans)
    
    mouse_shape = {'ch':ch_icshape, 'np':np_icshape, 'cy':cy_icshape, 'wc':wc_icshape}
    
    return mouse_shape


def load_mouse_TS(mm10_parseTrans, rev_hist=False):
    rna_seq_1 = pd.read_csv("/Share/home/zhangqf8/lipan/DYNAMIC/transcription/elifeMouse/fastq/50minDMSO_1/SRR935117.FPKM", sep="\t", skiprows=1)
    rna_seq_2 = pd.read_csv("/Share/home/zhangqf8/lipan/DYNAMIC/transcription/elifeMouse/fastq/50minDMSO_2/SRR935118.FPKM", sep="\t", skiprows=1)
    
    rna_seq_1 = dict(rna_seq_1.loc[:, ['#id_name', 'fpkm']].values)
    rna_seq_2 = dict(rna_seq_2.loc[:, ['#id_name', 'fpkm']].values)
    
    TS = []
    for geneID in rna_seq_1:
        if geneID in rna_seq_2:
            TS.append([geneID, rna_seq_1[geneID], rna_seq_2[geneID]])
    
    TS.sort(key=lambda x: x[0])
    TS = pd.DataFrame(TS, columns=['geneID', 'gro_seq_rpkm_1', 'gro_seq_rpkm_2'])
    TS = TS.iloc[4:,]
    
    TS['TS'] = numpy.log2((TS.gro_seq_rpkm_1+TS.gro_seq_rpkm_2)/2)
    TS['transID'] = TS['geneID'].apply( func=lambda x: x.split('_')[0] )
    TS = dict( zip(list(TS['transID']), list(TS['TS'])) )
    
    if rev_hist:
        remove_hist_genes(TS, mm10_parseTrans)
    
    return TS

def load_mouse_TE(mm10_parseTrans, rev_hist=False):
    rna_seq = pd.read_csv("/Share/home/zhangqf8/lipan/DYNAMIC/translation/Cell2011/fastq/SRR315594.occu.fpkm", sep="\t", skiprows=1)
    ribo_seq = pd.read_csv("/Share/home/zhangqf8/lipan/DYNAMIC/translation/Cell2011/fastq/SRR315623.occu.fpkm", sep="\t", skiprows=1)
    
    rna_seq = dict(rna_seq.loc[:, ['#id_name', 'fpkm']].values)
    ribo_seq = dict(ribo_seq.loc[:, ['#id_name', 'fpkm']].values)
    
    TE = []
    for geneID in rna_seq:
        if geneID in ribo_seq:
            TE.append([geneID, rna_seq[geneID], ribo_seq[geneID]])
    
    TE.sort(key=lambda x: x[0])
    TE = pd.DataFrame(TE, columns=['geneID', 'rna_seq_rpkm', 'ribo_seq_occ'])
    
    TE['TE'] = numpy.log2( divide(TE.ribo_seq_occ, TE.rna_seq_rpkm) )
    TE['transID'] = TE['geneID'].apply( func=lambda x:x.split('_')[0] )
    TE = dict( zip(list(TE['transID']), list(TE['TE'])) )
    
    if rev_hist:
        remove_hist_genes(TE, mm10_parseTrans)
    
    return TE


####################
### Load Mouse GINI
####################

def load_mouse_gini(mouse_shape, mm10_parseTrans, rev_hist=False):
    mouse_gini = calc_GINI(mouse_shape)
    mouse_start_coding_gini = calc_start_codon_gini(mouse_shape, mm10_parseTrans)
    mouse_5UTR_gini = calc_5UTR_gini(mouse_shape, mm10_parseTrans)
    mouse_3UTR_gini = calc_3UTR_gini(mouse_shape, mm10_parseTrans)
    
    if rev_hist:
        for comp in ('ch', 'np', 'cy', 'wc'):
            remove_hist_genes(mouse_gini[comp], mm10_parseTrans)
            remove_hist_genes(mouse_start_coding_gini[comp], mm10_parseTrans)
            remove_hist_genes(mouse_5UTR_gini[comp], mm10_parseTrans)
            remove_hist_genes(mouse_3UTR_gini[comp], mm10_parseTrans)
    
    return mouse_gini, mouse_start_coding_gini, mouse_5UTR_gini, mouse_3UTR_gini



####################
###  Human Load Functions
####################

def load_human_parser():
    hg38_parseTrans = ParseTransClass(genomeCoorBedFile = "/Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed")
    human_geneInfo = hg38_parseTrans.getGeneInfo()
    return hg38_parseTrans, human_geneInfo


def load_human_degaradation(hg38_parseTrans, rev_hist=False):
    human_halflife_file_1 = "/Share/home/zhangqf8/lipan/DYNAMIC/half-life/human/human_ID_halfLife.txt"
    
    hl_t_1 = read_half_life_to_Trans(human_halflife_file_1, human_geneInfo); print len(hl_t_1)
    hl_t_1 = dict( zip(hl_t_1.trans_id, hl_t_1.HL) )
    
    if rev_hist:
        remove_hist_genes(hl_t_1, hg38_parseTrans)
    
    return hl_t_1

def load_human_shape(hg38_parseTrans, shape_name='low_shape.out', rev_hist=False):
    icSHAPE_ROOT_DIR = '/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/%s/'+shape_name
    ch_icshape_file = icSHAPE_ROOT_DIR % ("hek_ch_vivo", )
    np_icshape_file = icSHAPE_ROOT_DIR % ("hek_np_vivo", )
    cy_icshape_file = icSHAPE_ROOT_DIR % ("hek_cy_vivo", )
    wc_icshape_file = icSHAPE_ROOT_DIR % ("hek_wc_vivo", )
    
    ch_icshape = tools.loadicSHAPE(ch_icshape_file)
    np_icshape = tools.loadicSHAPE(np_icshape_file)
    cy_icshape = tools.loadicSHAPE(cy_icshape_file)
    wc_icshape = tools.loadicSHAPE(wc_icshape_file)
    
    if rev_hist:
        remove_hist_genes(ch_icshape, hg38_parseTrans)
        remove_hist_genes(np_icshape, hg38_parseTrans)
        remove_hist_genes(cy_icshape, hg38_parseTrans)
        remove_hist_genes(wc_icshape, hg38_parseTrans)
    
    human_shape = {'ch':ch_icshape, 'np':np_icshape, 'cy':cy_icshape, 'wc':wc_icshape}
    
    return human_shape

def load_human_gini(human_shape, hg38_parseTrans, rev_hist=False):
    human_gini = calc_GINI(human_shape)
    human_start_coding_gini = calc_start_codon_gini(human_shape, hg38_parseTrans)
    human_5UTR_gini = calc_5UTR_gini(human_shape, hg38_parseTrans)
    human_3UTR_gini = calc_3UTR_gini(human_shape, hg38_parseTrans)
    
    if rev_hist:
        for comp in ('ch', 'np', 'cy', 'wc'):
            remove_hist_genes(human_gini[comp], hg38_parseTrans)
            remove_hist_genes(human_start_coding_gini[comp], hg38_parseTrans)
            remove_hist_genes(human_5UTR_gini[comp], hg38_parseTrans)
            remove_hist_genes(human_3UTR_gini[comp], hg38_parseTrans)
    
    return human_gini, human_start_coding_gini, human_5UTR_gini, human_3UTR_gini

def partial_cor(XY, XZ, YZ):
    """ 
       X <--> Z <--> Y
    """
    return (XY-XZ*YZ)/(np.sqrt(1-XZ**2) * np.sqrt(1-YZ**2))



####################
### Load Data
####################

mm10_parseTrans, mouse_geneInfo = load_mouse_parser()
mouse_shape = load_mouse_shape(mm10_parseTrans, rev_hist=False)
mouse_TS = load_mouse_TS(mm10_parseTrans, rev_hist=False)
mouse_TE = load_mouse_TE(mm10_parseTrans, rev_hist=False)
mouse_gini, mouse_start_coding_gini, mouse_5UTR_gini, mouse_3UTR_gini = load_mouse_gini(mouse_shape, mm10_parseTrans, rev_hist=False)


hg38_parseTrans, human_geneInfo = load_human_parser()
human_shape = load_human_shape(hg38_parseTrans, 'low_shape.out', rev_hist=False)
human_hl = load_human_degaradation(hg38_parseTrans, rev_hist=False)
human_gini, human_start_coding_gini, human_5UTR_gini, human_3UTR_gini = load_human_gini(human_shape, hg38_parseTrans, rev_hist=False)

