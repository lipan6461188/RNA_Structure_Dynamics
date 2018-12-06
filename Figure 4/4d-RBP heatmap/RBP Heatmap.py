
import tools, icSHAPE

def pbr_count(pbr_single):
    return sum([len(v) for v in pbr_single.values()])

def filterPBR(raw_pbr_ensembl, minLen=10, maxLen=30, min_region_num=100):
    filt_pbr_ensembl = {}
    for protein in raw_pbr_ensembl:
        filt_pbr_ensembl[protein] = {}
        raw_Count = pbr_count(raw_pbr_ensembl[protein])
        for transID in raw_pbr_ensembl[protein]:
            filt_pbr_ensembl[protein][transID] = [ region for region in raw_pbr_ensembl[protein][transID] if minLen<region[1]-region[0]<maxLen ]
            if len(filt_pbr_ensembl[protein][transID]) == 0: del filt_pbr_ensembl[protein][transID]
        new_Count = pbr_count(filt_pbr_ensembl[protein])
        if new_Count < min_region_num:
            del filt_pbr_ensembl[protein]
            print "%s is filtered for small dataset: %s" % (protein, new_Count)
        else:
            print "%s is filtered %.2f%%..." % (protein, 100.0*(raw_Count-new_Count)/raw_Count)
    return filt_pbr_ensembl

def generate_PBR_random(pbr_ensembl):
    import random
    
    # get all regions
    dataset = []
    rbp_count = 0
    for protein in pbr_ensembl:
        if protein in ('Random_Sample', 'Genome_Random'): continue
        rbp_count += 1
        for trans_id in pbr_ensembl[protein]:
            for region in pbr_ensembl[protein][trans_id]:
                dataset.append( [trans_id]+region )
    
    # sample from dataset
    number = len(dataset) / rbp_count
    sample_dataset = random.sample(dataset, number)
    
    # generate new set
    sample_dict = {}
    for item in sample_dataset:
        if item[0] not in sample_dict:
            sample_dict[ item[0] ] = []
        sample_dict[ item[0] ].append( item[1:] )
    
    pbr_ensembl['Random_Sample'] = sample_dict

def generate_Wide_random(pbr_ensembl, human_gtf):
    import random
    
    # get all regions
    dataset = []
    rbp_count = 0
    for protein in pbr_ensembl:
        if protein in ('Random_Sample', 'Genome_Random'): continue
        rbp_count += 1
        for trans_id in pbr_ensembl[protein]:
            for region in pbr_ensembl[protein][trans_id]:
                dataset.append( [trans_id]+region )
    
    # randomize
    randomize_dataset = []
    for item in dataset:
        transLen = human_gtf[item[0]]['trans_len']
        PBR_len = item[2]-item[1]
        if transLen > 60 and transLen-PBR_len-30 > 20:
            start = random.randint(20, transLen-PBR_len-30)
            randomize_dataset.append( (item[0], start, start+PBR_len) )
    
    # sample from dataset
    number = len(randomize_dataset) / rbp_count
    sample_dataset = random.sample(randomize_dataset, number)
    
    # generate new set
    sample_dict = {}
    for item in sample_dataset:
        if item[0] not in sample_dict:
            sample_dict[ item[0] ] = []
        sample_dict[ item[0] ].append( item[1:] )
    
    pbr_ensembl['Genome_Random'] = sample_dict

def randomize_all_PBR(pbr_ensembl, human_gtf):
    def random_single_protein(trans_region_list):
        randomize_dataset = {}
        for trans_id in trans_region_list:
            transLen = human_gtf[trans_id]['trans_len']
            
            if transLen <= 60: continue
            
            for region in trans_region_list[trans_id]:
                PBR_len = region[1]-region[0]
                if transLen-PBR_len-30 > 20:
                    start = random.randint(20, transLen-PBR_len-30)
                    new_region = (start, start+PBR_len)
                    if trans_id not in randomize_dataset:
                        randomize_dataset[trans_id] = []
                    randomize_dataset[trans_id].append(new_region)
            randomize_dataset[trans_id].sort(key=lambda x: x[0])
        return randomize_dataset
    
    randomized_PBR = {}
    for protein in pbr_ensembl:
        randomized_PBR[protein] = random_single_protein(pbr_ensembl[protein])
        print "Data Size: %s => %s" % (pbr_count(pbr_ensembl[protein]), pbr_count(randomized_PBR[protein]))
    
    return randomized_PBR

def statistic_single_pbr_shape_mean_2(single_pbr_ensembl, icShape, shrink=5, min_valid=5, will_mean=True):
    import numpy
    def tfl(raw_list):
        # to float list
        valid_list = [ float(it) for it in raw_list if it != 'NULL' ]
        if len(valid_list) >= min_valid:
            return valid_list
        else:
            return None
    
    common_trans_list = set(icShape['ch']) & set(icShape['np']) & set(icShape['cy'])
    
    ch_ShapeList = []; ch_w_ShapeList = []
    np_ShapeList = []; np_w_ShapeList = []
    cy_ShapeList = []; cy_w_ShapeList = []
    
    for transID in single_pbr_ensembl:
        if transID in common_trans_list:
            valid_trans = False
            for region in single_pbr_ensembl[transID]:
                start, end = region[0]+shrink, region[1]-shrink
                if start < end:
                    ch_mean_shape = tfl(icShape['ch'][transID][start:end])
                    np_mean_shape = tfl(icShape['np'][transID][start:end])
                    cy_mean_shape = tfl(icShape['cy'][transID][start:end])
                    if ch_mean_shape and np_mean_shape and cy_mean_shape:
                        valid_trans = True
                        if will_mean:
                            ch_ShapeList.append( numpy.mean( ch_mean_shape ) )
                            np_ShapeList.append( numpy.mean( np_mean_shape ) )
                            cy_ShapeList.append( numpy.mean( cy_mean_shape ) )
                        else:
                            ch_ShapeList.append( (transID, start, end, ch_mean_shape) )
                            np_ShapeList.append( (transID, start, end, np_mean_shape) )
                            cy_ShapeList.append( (transID, start, end, cy_mean_shape) )
            
            if valid_trans:
                ch_shape = tfl(icShape['ch'][transID])
                np_shape = tfl(icShape['np'][transID])
                cy_shape = tfl(icShape['cy'][transID])
                
                ch_w_ShapeList.append( numpy.mean( ch_shape ) )
                np_w_ShapeList.append( numpy.mean( np_shape ) )
                cy_w_ShapeList.append( numpy.mean( cy_shape ) )
    
    return { 'ch': ch_ShapeList, 'np': np_ShapeList, 'cy': cy_ShapeList, 'ch_w': ch_w_ShapeList, 'np_w': np_w_ShapeList, 'cy_w': cy_w_ShapeList }


def statistic_single_pbr_gini(single_pbr_ensembl, icShape, shrink=5, min_valid=5):
    import numpy
    def tfl(raw_list):
        # to float list
        valid_list = [ float(it) for it in raw_list if it != 'NULL' ]
        if len(valid_list) >= min_valid:
            return valid_list
        else:
            return None
    
    common_trans_list = set(icShape['ch']) & set(icShape['np']) & set(icShape['cy'])
    
    ch_ShapeList = []
    np_ShapeList = []
    cy_ShapeList = []
    for transID in single_pbr_ensembl:
        if transID in common_trans_list:
            for region in single_pbr_ensembl[transID]:
                start, end = region[0]+shrink, region[1]-shrink
                if start < end:
                    ch_shape_list = tfl(icShape['ch'][transID][start:end])
                    np_shape_list = tfl(icShape['np'][transID][start:end])
                    cy_shape_list = tfl(icShape['cy'][transID][start:end])
                    if ch_shape_list and np_shape_list and cy_shape_list:
                        ch_gini = tools.calcGINI( ch_shape_list )
                        np_gini = tools.calcGINI( np_shape_list )
                        cy_gini = tools.calcGINI( cy_shape_list )
                        
                        if ch_gini and np_gini and cy_gini:
                            ch_ShapeList.append( ch_gini )
                            np_ShapeList.append( np_gini )
                            cy_ShapeList.append( cy_gini )
    
    return { 'ch': ch_ShapeList, 'np': np_ShapeList, 'cy': cy_ShapeList }

def statistic_all_pbr_shape_mean_2(pbr_ensembl, icSHAPE, shrink=5, min_region=50, min_valid=5, will_mean=True):
    psm = {}; psm_w = {}
    filt_p = set()
    for protein in pbr_ensembl:
        ax = statistic_single_pbr_shape_mean_2(pbr_ensembl[protein], icSHAPE, shrink=shrink, min_valid=min_valid, will_mean=will_mean)
        print protein, len(ax['ch'])
        if len(ax['ch']) > min_region:
            psm[protein] = { 'ch': ax['ch'], 'np': ax['np'], 'cy': ax['cy'] }
            psm_w[protein] = { 'ch': ax['ch_w'], 'np': ax['np_w'], 'cy': ax['cy_w'] }
        else:
            filt_p.add(protein)
    print "...... Filter: ", filt_p
    return psm, psm_w

def statistic_all_pbr_gini(pbr_ensembl, icSHAPE, shrink=5, min_region=50, min_valid=5):
    psm = {}
    filt_p = set()
    for protein in pbr_ensembl:
        ax = statistic_single_pbr_gini(pbr_ensembl[protein], icSHAPE, shrink=shrink, min_valid=min_valid)
        print protein, len(ax['ch'])
        if len(ax['ch']) > min_region:
            psm[protein] = { 'ch': ax['ch'], 'np': ax['np'], 'cy': ax['cy'] }
        else:
            filt_p.add(protein)
    print "...... Filter: ", filt_p
    return psm

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

def pbrShape2DataFrame(pbrShape, balance_num=False, will_mean=True):
    from numpy import mean
    import pandas as pd
    import random
    
    def calc_mean(list_list):
        summ = 0; count = 0;
        for s_list in list_list:
            for it in s_list[3]:
                summ += it
            count += len(s_list[3])
        return 1.0*summ/count
    
    protein_names = pbrShape.keys()
    shape_df = tools.init_rect(rowNum=len(pbrShape), colNum=3, rowNames=protein_names, colNames=['ch','np','cy'])
    count_df = tools.init_rect(rowNum=len(pbrShape), colNum=3, rowNames=protein_names, colNames=['ch','np','cy'])
    
    for protein in protein_names:
        ch = pbrShape[protein]['ch']
        if balance_num: ch = random.sample(ch, balance_num)
        np = pbrShape[protein]['np']
        if balance_num: np = random.sample(np, balance_num)
        cy = pbrShape[protein]['cy']
        if balance_num: cy = random.sample(cy, balance_num)
        
        shape_df.loc[protein, 'ch'] = mean(ch) if will_mean else calc_mean(ch)
        shape_df.loc[protein, 'np'] = mean(np) if will_mean else calc_mean(np)
        shape_df.loc[protein, 'cy'] = mean(cy) if will_mean else calc_mean(cy)
        
        count_df.loc[protein, 'ch'] = len(ch) if will_mean else len(ch)
        count_df.loc[protein, 'np'] = len(np) if will_mean else len(ch)
        count_df.loc[protein, 'cy'] = len(cy) if will_mean else len(ch)
    
    return shape_df, count_df

def pbrShape2Sig(pbrShape, protein_list=[], balance_num=False, will_mean=True):
    import scipy, scipy.stats, random, statsmodels, statsmodels.sandbox.stats.multicomp
    
    def p2n(p):
        if p<1e-5: return 3
        if p<1e-3: return 2
        if p<1e-2: return 1
        return 0
    
    if not protein_list: protein_list = pbrShape.keys()
    
    sigRect = tools.init_rect(len(protein_list), 2, rowNames=protein_list, colNames=['ch_np','np_cy'])
    
    raw_p = []
    for protein in protein_list:
        ch = pbrShape[protein]['ch']
        if not will_mean: ch = [ tools.np.mean(it[3]) for it in ch ]
        if balance_num: ch = random.sample(ch, balance_num)
        
        np = pbrShape[protein]['np']
        if not will_mean: np = [ tools.np.mean(it[3]) for it in np ]
        if balance_num: np = random.sample(np, balance_num)
        
        cy = pbrShape[protein]['cy']
        if not will_mean: cy = [ tools.np.mean(it[3]) for it in cy ]
        if balance_num: cy = random.sample(cy, balance_num)
        
        ch_np = scipy.stats.mannwhitneyu(ch, np)[1]
        np_cy = scipy.stats.mannwhitneyu(np, cy)[1]
        
        raw_p.append(ch_np); raw_p.append(np_cy)
        sigRect.loc[protein, 'ch_np'] = p2n(ch_np)
        sigRect.loc[protein, 'np_cy'] = p2n(np_cy)
    
    adj_p = list(statsmodels.sandbox.stats.multicomp.multipletests(method='bonferroni', pvals=raw_p)[1])
    adj_p.reverse()
    for protein in protein_list:
        sigRect.loc[protein, 'ch_np'] = p2n(adj_p.pop())
        sigRect.loc[protein, 'np_cy'] = p2n(adj_p.pop())
    
    return sigRect



def slice_icSHAPE(icSHAPE, trans_list, wsize=30, wstep=30):
    def tfl(raw_list):
        # to float list
        valid_list = [ float(it) for it in raw_list if it != 'NULL' ]
        if len(valid_list) >= 15:
            return valid_list
        else:
            return None
    ch_slices = []
    np_slices = []
    cy_slices = []
    for tid in trans_list:
        length = len(icSHAPE['ch'][tid])
        start = 0
        while start + wsize <= length:
            ch = tfl(icSHAPE['ch'][tid][start:start+wsize])
            np = tfl(icSHAPE['np'][tid][start:start+wsize])
            cy = tfl(icSHAPE['cy'][tid][start:start+wsize])
            start += wstep
            if ch and np and cy:
                ch_slices += ch
                np_slices += np
                cy_slices += cy
    return ch_slices, np_slices, cy_slices

def get_RBP_slices(pbr_ensembl, icSHAPE, protein_list):
    dataset = {}
    for protein in protein_list:
        dataset[protein] = {'ch':"", 'np':"", 'cy':""}
        process_trans_list = set(pbr_ensembl[protein]) & set(icSHAPE['ch']) & set(icSHAPE['np']) & set(icSHAPE['cy'])
        ch,np,cy = slice_icSHAPE(icSHAPE, process_trans_list, wsize=30, wstep=30)
        dataset[protein]['ch'] = ch
        dataset[protein]['np'] = np
        dataset[protein]['cy'] = cy
    return dataset


###### Load Human Transcriptome

human_seq = tools.readSeq("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/transcriptome.fa")
human_gtf = icSHAPE.loadGTFBed("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Homo_sapiens.transCoordinate.bed")

###### Load Human icSHAPE Data

icSHAPE_ROOT_DIR = '/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/%s/shape.out'
ch_icshape_file = icSHAPE_ROOT_DIR % ("hek_ch_vivo", )
np_icshape_file = icSHAPE_ROOT_DIR % ("hek_np_vivo", )
cy_icshape_file = icSHAPE_ROOT_DIR % ("hek_cy_vivo", )

ch_icshape = tools.loadicSHAPE(ch_icshape_file)
np_icshape = tools.loadicSHAPE(np_icshape_file)
cy_icshape = tools.loadicSHAPE(cy_icshape_file)

Shape = {'ch':ch_icshape, 'np':np_icshape, 'cy':cy_icshape}


###### Load RBP Data

clip_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/m6A/choose_good_protein/clipDB.merged.transCoor.bed'
pbr_ensembl = read_pbr(clip_file, human_seq)

# Filter out 10-50 regions

filt_pbr_ensembl = filterPBR(pbr_ensembl, minLen=10, maxLen=50, min_region_num=100)

# Generate 2 random sets

generate_PBR_random(filt_pbr_ensembl)
generate_Wide_random(filt_pbr_ensembl, human_gtf)

filt_pbr_ensembl['NCBP3'] = filt_pbr_ensembl['C17orf85']; del filt_pbr_ensembl['C17orf85']
rand_pbr_ensembl = randomize_all_PBR(filt_pbr_ensembl, human_gtf)

# Statistic PBR icSHAPE Mean

pbr_shape_mean, pbr_w_shape_mean = statistic_all_pbr_shape_mean_2(filt_pbr_ensembl, Shape, shrink=2, min_region=100, min_valid=15, will_mean=False)
#pbr_gini = statistic_all_pbr_gini(filt_pbr_ensembl, Shape, shrink=2, min_region=100, min_valid=15)

# Covert Shape List to Dataframe

pbrDataFrame, pbrCountFrame = pbrShape2DataFrame(pbr_shape_mean, balance_num=False, will_mean=False)
pbrDataFrame_w, pbrCountFrame_w = pbrShape2DataFrame(pbr_w_shape_mean, balance_num=False, will_mean=True)

# pbrDataFrame = pbrShape2DataFrame(pbr_gini, balance=False)

# Sort Protein Names



protein_loc = """
TAF15   5   0   0
HNRNPD  5   0   0
HNRNPU  5   0   0
HNRNPM  5   0   0
HNRNPA1 5   0   0
HNRNPF  5   0   0
EWSR1   5   0   0
WDR33   5   0   0
CPSF3   5   0   0
CPSF4   5   0   0
RBM15   5   0   0
HNRNPC  5   0   0
CPSF1   5   0   0
CPSF7   5   0   0
RBM15B  5   0   0
ELAVL1  0   5   0
AGO2_miRNA  0   5   5
FXR1    0   5   5
STAU1   0   5   5
CPSF2   5   0   0
MOV10   5   2   0
IGF2BP2 0   5   5
ALKBH5  0   5   5
PTBP1   5   0   0
IGF2BP1 0   5   5
LIN28A  0   5   5
ZC3H7B  0   4   3
FUS 0   5   2
IGF2BP3 0   5   5
LIN28B  0   5   5
YTHDF2  0   5   5
YTHDF1  0   2   2
FXR2    0   5   5
ATXN2   0   3   4
RTCB    0   5   5
CAPRIN1 0   3   5
DGCR8   0   5   3
CSTF2   5   0   0
SRRM4   0   5   2
CSTF2T  5   0   0
FMR1    0   5   5
TARDBP  0   5   3
METTL3  0   5   2
FIP1L1  0   5   2
NUDT21  0   5   1
YTHDC1  0   5   0
CPSF6   5   0   0
NOP58   0   5   5
NCBP3   0   5   2
FBL 0   5   2
"""

def classify_protein():
    plist = protein_loc.strip().split('\n')
    in_nuc = []; both = []; in_cyt = []
    for single_p in plist:
        pName, chr_weight, nuc_weight, cyt_weight = single_p.split()
        if int(chr_weight) > int(cyt_weight) or int(nuc_weight) > int(cyt_weight):
            in_nuc.append(pName)
        elif int(nuc_weight) < int(cyt_weight):
            in_cyt.append(pName)
        else:
            both.append(pName)
    return in_nuc, both, in_cyt

def sorted_protein_list(pbrDataFrame):
    in_nuc, both, in_cyt = classify_protein()
    nuc_pro = []
    both_pro = []
    cyt_pro = []
    for protein in in_nuc:
        ch_shape = pbrDataFrame.loc[protein,'ch']
        np_shape = pbrDataFrame.loc[protein,'np']
        nuc_pro.append( (protein, ch_shape-np_shape) )
    
    for protein in both:
        ch_shape = pbrDataFrame.loc[protein,'ch']
        np_shape = pbrDataFrame.loc[protein,'np']
        both_pro.append( (protein, ch_shape-np_shape) )
    
    for protein in in_cyt:
        ch_shape = pbrDataFrame.loc[protein,'ch']
        np_shape = pbrDataFrame.loc[protein,'np']
        cyt_pro.append( (protein, ch_shape-np_shape) )
    
    nuc_pro = [ it[0] for it in sorted(nuc_pro, key=lambda x: x[1], reverse=True) ]
    both_pro = [ it[0] for it in sorted(both_pro, key=lambda x: x[1], reverse=True) ]
    cyt_pro = [ it[0] for it in sorted(cyt_pro, key=lambda x: x[1], reverse=True) ]
    
    return nuc_pro,both_pro,cyt_pro

nuc_pro,both_pro,cyt_pro = sorted_protein_list(pbrDataFrame)

#pbrDataFrame['ch-np'] = pbrDataFrame['ch']-pbrDataFrame['np']
#pbrDataFrame = pbrDataFrame.sort_values(by="ch-np", inplace=False, ascending=False).iloc[:,0:3]

pbrDataFrame = pbrDataFrame.loc[nuc_pro+both_pro+cyt_pro+['Genome_Random', 'Random_Sample'], :]
pbrCountFrame = pbrCountFrame.loc[pbrDataFrame.index, :]
pbrDataFrame_w = pbrDataFrame_w.loc[pbrDataFrame.index, :]


# Generate Significance Dataframe

sig_dataFrame = pbrShape2Sig(pbr_shape_mean, protein_list=list(pbrDataFrame.index), balance_num=False, will_mean=False)
# sig_dataFrame = pbrShape2Sig(pbr_gini, protein_list=list(pbrDataFrame.index), balance=False)


# Plot icSHAPE heatmap

tools.plt.figure(figsize=(5,15))
tools.sns.heatmap(data=pbrDataFrame, annot=True, cmap=tools.sns.color_palette("RdBu_r", 30), vmin=0.15, vmax=0.30)
tools.plt.yticks(rotation=0)
tools.plt.tight_layout()
tools.plt.savefig("/Share/home/zhangqf8/figs/RBP.pdf")
tools.plt.show()

# Plot Count heatmap

tools.plt.figure(figsize=(5,15))
tools.sns.heatmap(data=pbrCountFrame, annot=True, fmt='.0f', cmap=tools.sns.color_palette("Blues", 30))
tools.plt.yticks(rotation=0)
tools.plt.tight_layout()
tools.plt.savefig("/Share/home/zhangqf8/figs/RBP_count.pdf")
tools.plt.close()

# Plot transcript icSHAPE heatmap

tools.plt.figure(figsize=(5,15))
tools.sns.heatmap(data=pbrDataFrame_w, annot=True, cmap=tools.sns.color_palette("RdBu_r", 30), vmin=0.15, vmax=0.30)
tools.plt.yticks(rotation=0)
tools.plt.tight_layout()
tools.plt.savefig("/Share/home/zhangqf8/figs/RBP_trans.pdf")
tools.plt.close()


# Plot Significance heatmap

anno_frame = sig_dataFrame.replace(3, "***").replace(2, "**").replace(1, "*").replace(0, "")
tools.plt.figure(figsize=(5,15))
tools.sns.heatmap(data=sig_dataFrame, cmap=tools.sns.color_palette("RdBu_r", 30), annot=anno_frame, fmt='s')
tools.plt.yticks(rotation=0)
tools.plt.tight_layout()
tools.plt.savefig("/Share/home/zhangqf8/figs/RBP_sig.pdf")
tools.plt.show()

# print all protein names

for p in pbrDataFrame.index: print p

# plot protein

slices = get_RBP_slices(filt_pbr_ensembl, Shape, pbrDataFrame.index)
shape_df, count_df = pbrShape2DataFrame(slices, balance_num=False, will_mean=True)
shape_df = shape_df.loc[pbrDataFrame.index, :]
count_df = count_df.loc[pbrDataFrame.index, :]

font = {'family': 'normal', 'weight': 'bold','size': 8}
matplotlib.rc('font', **font)

tools.plt.figure(figsize=(5,15))
tools.sns.heatmap(data=shape_df, annot=True, fmt='.3f', cmap=tools.sns.color_palette("RdBu_r", 30), vmin=0.15, vmax=0.30)
tools.plt.yticks(rotation=0)
tools.plt.tight_layout()
tools.plt.savefig("/Share/home/zhangqf8/figs/RBP_slice.pdf")
tools.plt.show()

tools.plt.figure(figsize=(5,15))
tools.sns.heatmap(data=count_df, annot=True, fmt='.3f', cmap=tools.sns.color_palette("RdBu_r", 30))
tools.plt.yticks(rotation=0)
tools.plt.tight_layout()
tools.plt.savefig("/Share/home/zhangqf8/figs/RBP_slice.pdf")
tools.plt.show()











