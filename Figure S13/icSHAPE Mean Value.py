
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
import ParseTrans

def calc_ave_SHAPE(shape_dict, trans_list, Parser):
    def my_calc_ave_shape(raw_shape_list):
        valid_shape = [ float(it) for it in raw_shape_list if it != 'NULL' ]
        if len(valid_shape) > 20:
            return sum(valid_shape)/len(valid_shape)
        else:
            return None
    
    import seq
    Shape = {}
    valid_trans = [0]*len(trans_list)
    index = -1
    for trans_id in trans_list:
        index += 1
        ave_shape = my_calc_ave_shape(shape_dict[trans_id])
        if ave_shape:
            try:
                gene_type = Parser.getTransFeature(trans_id, verbose=False)['gene_type']
                gene_type = seq.anno_Methods.gene_type(gene_type)
            except KeyError:
                continue
            if gene_type not in Shape:
                Shape[gene_type] = []
            Shape[gene_type].append( ave_shape )
            valid_trans[index] = 1
    
    return Shape, valid_trans

def draw_icSHAPE_map(ch_shape, np_shape, cy_shape, Parser, will_common=False):
    if will_common:
        common_trans_list = set(ch_shape) & set(np_shape) & set(cy_shape)
        ave_ch, valid_ch = calc_ave_SHAPE(ch_shape, common_trans_list, Parser)
        ave_np, valid_np = calc_ave_SHAPE(np_shape, common_trans_list, Parser)
        ave_cy, valid_cy = calc_ave_SHAPE(cy_shape, common_trans_list, Parser)
        
        valid_trans_list = []
        for i,j,k,tid in zip(valid_ch, valid_np, valid_cy,common_trans_list):
            if 0 not in (i, j, k):
                valid_trans_list.append(tid)
        
        ave_ch, valid_ch = calc_ave_SHAPE(ch_shape, valid_trans_list, Parser)
        ave_np, valid_np = calc_ave_SHAPE(np_shape, valid_trans_list, Parser)
        ave_cy, valid_cy = calc_ave_SHAPE(cy_shape, valid_trans_list, Parser)
    else:
        ave_ch, valid_ch = calc_ave_SHAPE(ch_shape, ch_shape.keys(), Parser)
        ave_np, valid_np = calc_ave_SHAPE(np_shape, np_shape.keys(), Parser)
        ave_cy, valid_cy = calc_ave_SHAPE(cy_shape, cy_shape.keys(), Parser)
    
    RNA_names = ['mRNA', 'snRNA', 'lncRNA', 'pseudogene', 'miRNA', 'snoRNA']
    shape_matrix = tools.init_rect(6, 3, rowNames=RNA_names, colNames=['ch','np','cy'])
    num_matrix = tools.init_rect(6, 3, rowNames=RNA_names, colNames=['ch','np','cy'])
    
    for rna_type in RNA_names:
        ave_ch[rna_type] = ave_ch.get(rna_type, [0.0])
        ave_np[rna_type] = ave_np.get(rna_type, [0.0])
        ave_cy[rna_type] = ave_cy.get(rna_type, [0.0])
    
    for idx in range(6):
        shape_matrix.iloc[idx, 0] = tools.np.mean(ave_ch[ RNA_names[idx] ])
        shape_matrix.iloc[idx, 1] = tools.np.mean(ave_np[ RNA_names[idx] ])
        shape_matrix.iloc[idx, 2] = tools.np.mean(ave_cy[ RNA_names[idx] ])
        
        num_matrix.iloc[idx, 0] = len( ave_ch[ RNA_names[idx] ] )
        num_matrix.iloc[idx, 1] = len( ave_np[ RNA_names[idx] ] )
        num_matrix.iloc[idx, 2] = len( ave_cy[ RNA_names[idx] ] )
    
    import matplotlib
    font = {'family': 'normal', 'weight': 'normal','size': 15}
    matplotlib.rc('font', **font)
    tools.plt.figure(figsize=(7,10))
    tools.plt.subplot(2,1,1)
    tools.sns.heatmap(shape_matrix, annot=True, fmt='.3f', cmap=tools.sns.color_palette("BuGn", 30), vmin=0.1, vmax=0.2)
    tools.plt.yticks(rotation=0)
    tools.plt.subplot(2,1,2)
    tools.sns.heatmap(num_matrix, annot=True, fmt='.1f', cmap=tools.sns.color_palette("Blues", 30))
    tools.plt.yticks(rotation=0)
    tools.plt.tight_layout()

def sub_Shape(raw_shape, trans_list):
    return { k:raw_shape[k] for k in trans_list if k in set(raw_shape.keys()) }


#########################
# 蛋白结合区icSHAPE在三个组分中的差异
#########################


def read_pbr(file_name, SeqDict, pbr_column={'pName': 5, 'score':6, 'transID': 9, 'transStart': 10, 'transEnd': 11}):                                                
    # + 从文件中读取pbr
    # + 默认读取~/lipan/DYNAMIC/modification/HEK293/m6A/choose_good_protein/clipDB.merged.transCoor.bed
    sorted_transList = sorted(SeqDict.keys())
    pbr_ensembl = {}
    IN = open(file_name)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        protein = arr[ pbr_column['pName']-1 ]
        score = float( arr[pbr_column['score']-1] )
        trans_id = arr[ pbr_column['transID']-1 ]
        trans_start = int(arr[ pbr_column['transStart']-1 ]) - 1
        trans_end = int(arr[ pbr_column['transEnd']-1 ])
        if tools.biSearch(trans_id, sorted_transList) and len(SeqDict[trans_id]) > 100:
            if 20 <= trans_start and trans_end <= len( SeqDict[trans_id] ) - 20:
                if protein not in pbr_ensembl:
                    pbr_ensembl[protein] = {}
                if trans_id not in pbr_ensembl[protein]:
                    pbr_ensembl[ protein ][ trans_id ] = []
                pbr_ensembl[ protein ][ trans_id ].append( [trans_start, trans_end, score] )
        line = IN.readline()
    for protein in pbr_ensembl:
        for trans_id in pbr_ensembl[ protein ]:
            pbr_ensembl[ protein ][ trans_id ].sort(key=lambda x: x[0])
    return pbr_ensembl

def filterPBR(raw_pbr_ensembl, minLen=10, maxLen=30):
    """
    过滤掉一些结合区域太短或太长的蛋白结合区域
    2017-5-24
    """
    filt_pbr_ensembl = {}
    for protein in raw_pbr_ensembl:
        filt_pbr_ensembl[protein] = {}
        raw_Count = pbr_count(raw_pbr_ensembl[protein])
        for transID in raw_pbr_ensembl[protein]:
            filt_pbr_ensembl[protein][transID] = [ region for region in raw_pbr_ensembl[protein][transID] if minLen<region[1]-region[0]<maxLen ]
            if len(filt_pbr_ensembl[protein][transID]) == 0: del filt_pbr_ensembl[protein][transID]
        new_Count = pbr_count(filt_pbr_ensembl[protein])
        print "%s is filtered %.2f%%..." % (protein, 100.0*(raw_Count-new_Count)/raw_Count)
    return filt_pbr_ensembl


def calc_RBP_ave_SHAPE_common(RBP, ch_shape_dict, np_shape_dict, cy_shape_dict, Parser):
    def calc_ave_shape(raw_shape_list):
        valid_shape = [ float(it) for it in raw_shape_list if it != 'NULL' ]
        if len(valid_shape) > 10:
            return sum(valid_shape)/len(valid_shape)
        else:
            return None
    
    common_trans_list = set(ch_shape_dict) & set(np_shape_dict) & set(cy_shape_dict) & set(RBP.keys())
    
    import seq
    ch_Shape = {}
    np_Shape = {}
    cy_Shape = {}
    for trans_id in common_trans_list:
        for region in RBP[trans_id]:
            start, end = region[0]+5, region[1]-5
            if start >= end: continue
            try:
                ave_ch_shape = calc_ave_shape(ch_shape_dict[trans_id][start:end])
                ave_np_shape = calc_ave_shape(np_shape_dict[trans_id][start:end])
                ave_cy_shape = calc_ave_shape(cy_shape_dict[trans_id][start:end])
            except KeyError:
                continue
            if ave_ch_shape and ave_np_shape and ave_cy_shape:
                try:
                    gene_type = Parser.getTransFeature(trans_id, verbose=False)['gene_type']
                    gene_type = seq.anno_Methods.gene_type(gene_type)
                except KeyError:
                    continue
                if gene_type not in ch_Shape:
                    ch_Shape[gene_type] = []
                    np_Shape[gene_type] = []
                    cy_Shape[gene_type] = []
                ch_Shape[gene_type].append( ave_ch_shape )
                np_Shape[gene_type].append( ave_np_shape )
                cy_Shape[gene_type].append( ave_cy_shape )
    
    return ch_Shape,np_Shape,cy_Shape

def calc_RBP_ave_SHAPE(RBP, ch_shape_dict, np_shape_dict, cy_shape_dict, Parser):
    def calc_ave_shape(raw_shape_list):
        valid_shape = [ float(it) for it in raw_shape_list if it != 'NULL' ]
        if len(valid_shape) > 10:
            return sum(valid_shape)/len(valid_shape)
        else:
            return None
    
    ch_set = set(ch_shape_dict)
    np_set = set(np_shape_dict)
    cy_set = set(cy_shape_dict)
    
    import seq
    ch_Shape = {}
    np_Shape = {}
    cy_Shape = {}
    for trans_id in RBP:
        
        in_ch = (trans_id in ch_set)
        in_np = (trans_id in np_set)
        in_cy = (trans_id in cy_set)
        
        try:
            gene_type = Parser.getTransFeature(trans_id, verbose=False)['gene_type']
            gene_type = seq.anno_Methods.gene_type(gene_type)
        except KeyError:
            continue
        
        for region in RBP[trans_id]:
            start, end = region[0]+5, region[1]-5
            if start >= end: continue
            
            # for ch
            if in_ch:
                ave_ch_shape = calc_ave_shape(ch_shape_dict[trans_id][start:end])
                if gene_type not in ch_Shape:
                    ch_Shape[gene_type] = []
                if ave_ch_shape: ch_Shape[gene_type].append( ave_ch_shape )
            
            # for np
            if in_np:
                ave_np_shape = calc_ave_shape(np_shape_dict[trans_id][start:end])
                if gene_type not in np_Shape:
                    np_Shape[gene_type] = []
                if ave_np_shape: np_Shape[gene_type].append( ave_np_shape )
            
            # for cy
            if in_cy:
                ave_cy_shape = calc_ave_shape(cy_shape_dict[trans_id][start:end])
                if gene_type not in cy_Shape:
                    cy_Shape[gene_type] = []
                if ave_cy_shape: cy_Shape[gene_type].append( ave_cy_shape )
    
    return ch_Shape, np_Shape, cy_Shape

def draw_RBP_icSHAPE_map(RBP, ch_shape, np_shape, cy_shape, Parser):
    ave_ch, ave_np, ave_cy = calc_RBP_ave_SHAPE(RBP, ch_shape, np_shape, cy_shape, Parser)
    
    RNA_names = ['mRNA', 'snRNA', 'lncRNA', 'pseudogene', 'miRNA', 'snoRNA']
    shape_matrix = tools.init_rect(6, 3, rowNames=RNA_names, colNames=['ch','np','cy'])
    num_matrix = tools.init_rect(6, 3, rowNames=RNA_names, colNames=['ch','np','cy'])
    
    for rna_type in RNA_names:
        ave_ch[rna_type] = ave_ch.get(rna_type, [0.0])
        ave_np[rna_type] = ave_np.get(rna_type, [0.0])
        ave_cy[rna_type] = ave_cy.get(rna_type, [0.0])
    
    for idx in range(6):
        shape_matrix.iloc[idx, 0] = tools.np.mean(ave_ch[ RNA_names[idx] ])
        shape_matrix.iloc[idx, 1] = tools.np.mean(ave_np[ RNA_names[idx] ])
        shape_matrix.iloc[idx, 2] = tools.np.mean(ave_cy[ RNA_names[idx] ])
        
        num_matrix.iloc[idx, 0] = len( ave_ch[ RNA_names[idx] ] )
        num_matrix.iloc[idx, 1] = len( ave_np[ RNA_names[idx] ] )
        num_matrix.iloc[idx, 2] = len( ave_cy[ RNA_names[idx] ] )
    
    import matplotlib
    font = {'family': 'normal', 'weight': 'normal','size': 15}
    matplotlib.rc('font', **font)
    tools.plt.figure(figsize=(7,10))
    tools.plt.subplot(2,1,1)
    tools.sns.heatmap(shape_matrix, annot=True, fmt='.3f', cmap=sns.color_palette("BuGn", 30), vmin=0.1, vmax=0.2)
    tools.plt.yticks(rotation=0)
    tools.plt.subplot(2,1,2)
    tools.sns.heatmap(num_matrix, annot=True, fmt='.1f', cmap=sns.color_palette("Blues", 30))
    tools.plt.yticks(rotation=0)
    tools.plt.tight_layout()


#########################
# 计算三个区域都有值的icSHAPE Heatmap
#########################


def merge_regions(regions_list):
    regions_list = sorted(regions_list, key=lambda x: x[0])
    merged_region = []
    for region in regions_list:
        if len(merged_region) == 0:
            merged_region.append(region)
        elif merged_region[-1][1] >= region[0]:
            merged_region[-1] = (merged_region[-1][0], max(region[1], merged_region[-1][1]) )
        else:
            merged_region.append(region)
    return merged_region


def build_PBR_regions(pbr_ensembl):
    RBP_regions = {}
    RBP_trans = set()
    for protein in pbr_ensembl:
        for tid in pbr_ensembl[protein]:
            if tid not in RBP_trans:
                RBP_trans.add(tid)
                RBP_regions[tid] = []
            RBP_regions[tid] += pbr_ensembl[protein][tid]
    
    # sort
    for tid in RBP_regions:
        RBP_regions[tid].sort(key=lambda x: x[0])
    
    # merge
    for tid in RBP_regions:
        RBP_regions[tid] = merge_regions(RBP_regions[tid])
    
    return RBP_regions

def calc_ave_SHAPE_Sites(icSHAPE, Parser, will_common=True, rem_regions=None, inc_regions=None):
    def classify_shape(tid, ch, np, cy):
        c_ch = []
        c_np = []
        c_cy = []
        idx = 0
        
        if rem_regions:
            rem_tid_set = set(rem_regions)
            tid_in = (tid in rem_tid_set)
        
        if inc_regions:
            inc_tid_set = set(inc_regions)
            tid_in = (tid in inc_tid_set)
        
        for ch_i, np_i, cy_i in zip(ch, np, cy):
            idx += 1
            if will_common and 'NULL' in (ch_i, np_i, cy_i):
                continue
            
            if rem_regions and tid_in:
                in_rem_region = False
                for region in rem_regions[tid]:
                    if region[0] <= idx <= region[1]:
                        in_rem_region = True
                        break
                if in_rem_region:
                    continue
            
            if inc_regions:
                if not tid_in: continue
                in_inc_region = False
                for region in inc_regions[tid]:
                    if region[0] <= idx <= region[1]:
                        in_inc_region = True
                        break
                if not in_inc_region:
                    continue
            
            if ch_i!='NULL': c_ch.append(float(ch_i))
            if np_i!='NULL': c_np.append(float(np_i))
            if cy_i!='NULL': c_cy.append(float(cy_i))
        return c_ch, c_np, c_cy
    
    if rem_regions and inc_regions:
        print "Error: rem_regions and inc_regions are exclusive"
        return
    
    icSHAPE_average = {}
    trans_num = {}
    common_trans_list = set(icSHAPE['ch']) & set(icSHAPE['np']) & set(icSHAPE['cy'])
    for tid in common_trans_list:
        try:
            gene_type = Parser.getTransFeature(tid, verbose=False)['gene_type']
            if gene_type == 'protein_coding': gene_type = 'mRNA'
        except KeyError:
            continue
        if gene_type not in icSHAPE_average:
            icSHAPE_average[gene_type] = {'ch': [], 'np': [], 'cy':[]}
            trans_num[gene_type] = {'ch': 0, 'np': 0, 'cy':0}
        
        c_ch, c_np, c_cy = classify_shape(tid, icSHAPE['ch'][tid], icSHAPE['np'][tid], icSHAPE['cy'][tid])
        if c_ch or c_np or c_cy:
            if c_ch:
                trans_num[gene_type]['ch'] += 1
            if c_np:
                trans_num[gene_type]['np'] += 1
            if c_cy:
                trans_num[gene_type]['cy'] += 1
            
            icSHAPE_average[gene_type]['ch'] += c_ch
            icSHAPE_average[gene_type]['np'] += c_np
            icSHAPE_average[gene_type]['cy'] += c_cy
    
    RNA_names = ['mRNA', 'snRNA', 'lncRNA', 'pseudogene', 'miRNA', 'snoRNA']
    shape_matrix = tools.init_rect(6, 3, rowNames=RNA_names, colNames=['ch','np','cy'])
    nuc_num_matrix = tools.init_rect(6, 3, rowNames=RNA_names, colNames=['ch','np','cy'])
    trans_num_matrix = tools.init_rect(6, 3, rowNames=RNA_names, colNames=['ch','np','cy'])
    
    for rna_type in RNA_names:
        icSHAPE_average[rna_type] = icSHAPE_average.get(rna_type, {'ch':[0.0], 'np':[0.0], 'cy':[0.0]})
        trans_num[rna_type] = trans_num.get(rna_type, {'ch':0, 'np':0, 'cy':0})
    
    for idx in range(6):
        gt = RNA_names[idx]
        
        trans_num_matrix.loc[gt, 'ch'] = trans_num[gt]['ch']
        trans_num_matrix.loc[gt, 'np'] = trans_num[gt]['np']
        trans_num_matrix.loc[gt, 'cy'] = trans_num[gt]['cy']
        
        shape_matrix.loc[gt, 'ch'] = tools.np.mean( icSHAPE_average[gt]['ch'] )
        shape_matrix.loc[gt, 'np'] = tools.np.mean( icSHAPE_average[gt]['np'] )
        shape_matrix.loc[gt, 'cy'] = tools.np.mean( icSHAPE_average[gt]['cy'] )
        
        nuc_num_matrix.loc[gt, 'ch'] = len( icSHAPE_average[gt]['ch'] )
        nuc_num_matrix.loc[gt, 'np'] = len( icSHAPE_average[gt]['np'] )
        nuc_num_matrix.loc[gt, 'cy'] = len( icSHAPE_average[gt]['cy'] )
    
    return shape_matrix, trans_num_matrix, nuc_num_matrix


def plot_matrix(shape_matrix, trans_num_matrix, nuc_num_matrix):
    import matplotlib
    font = {'family': 'normal', 'weight': 'normal','size': 15}
    matplotlib.rc('font', **font)
    tools.plt.figure(figsize=(7,10))
    tools.plt.subplot(3,1,1)
    tools.sns.heatmap(shape_matrix, annot=True, fmt='.3f', cmap=tools.sns.color_palette("BuGn", 30), vmin=0.1, vmax=0.25)
    tools.plt.yticks(rotation=0)
    
    tools.plt.subplot(3,1,2)
    tools.sns.heatmap(trans_num_matrix, annot=True, fmt='.0f', cmap=tools.sns.color_palette("Blues", 30))
    tools.plt.yticks(rotation=0)
    
    tools.plt.subplot(3,1,3)
    tools.sns.heatmap(nuc_num_matrix, annot=True, fmt='.0f', cmap=tools.sns.color_palette("Blues", 30))
    tools.plt.yticks(rotation=0)
    tools.plt.tight_layout()



in_root = "/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/%s_%s_%s/shape.out"

# in vivo

human_ch = tools.loadicSHAPE(in_root % ('hek', 'ch', 'vivo'))
human_np = tools.loadicSHAPE(in_root % ('hek', 'np', 'vivo'))
human_cy = tools.loadicSHAPE(in_root % ('hek', 'cy', 'vivo'))

mouse_ch = tools.loadicSHAPE(in_root % ('mes', 'ch', 'vivo'))
mouse_np = tools.loadicSHAPE(in_root % ('mes', 'np', 'vivo'))
mouse_cy = tools.loadicSHAPE(in_root % ('mes', 'cy', 'vivo'))

mm10_parseTrans = ParseTrans.ParseTransClass(genomeCoorBedFile = "/tmp/mm10.genomeCoor.bed")
hg38_parseTrans = ParseTrans.ParseTransClass(genomeCoorBedFile = "/tmp/hg38.genomeCoor.bed")


human_icSHAPE = {'ch': human_ch, 'np': human_np, 'cy': human_cy}
mouse_icSHAPE = {'ch': mouse_ch, 'np': mouse_np, 'cy': mouse_cy}

### human

shape_matrix, trans_num_matrix, nuc_num_matrix = calc_ave_SHAPE_Sites(human_icSHAPE, hg38_parseTrans, will_common=True)
plot_matrix(shape_matrix, trans_num_matrix, nuc_num_matrix)
tools.plt.savefig('figs/1.pdf')
tools.plt.show()


# remove RBP regions

human_seq = tools.readSeq("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/transcriptome.fa")
clip_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/m6A/choose_good_protein/clipDB.merged.transCoor.bed'
pbr_ensembl = read_pbr(clip_file, human_seq)
merged_RBP_regions = build_PBR_regions(pbr_ensembl)

shape_matrix, trans_num_matrix, nuc_num_matrix = calc_ave_SHAPE_Sites(human_icSHAPE, hg38_parseTrans, will_common=True, rem_regions=merged_RBP_regions)
plot_matrix(shape_matrix, trans_num_matrix, nuc_num_matrix)
tools.plt.savefig('figs/1.pdf')
tools.plt.show()

### mouse

shape_matrix, trans_num_matrix, nuc_num_matrix = calc_ave_SHAPE_Sites(mouse_icSHAPE, mm10_parseTrans, will_common=True)
plot_matrix(shape_matrix, trans_num_matrix, nuc_num_matrix)
tools.plt.savefig('figs/2.pdf')
tools.plt.show()




