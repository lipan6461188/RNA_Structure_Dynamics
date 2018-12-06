

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from icSHAPE import *


###################
# Functions
###################

def readDynamics(resultFile):
    dynamicDict = {}
    detected = 0; dynamic = 0;
    IN = open(resultFile)
    line = IN.readline()
    line = IN.readline()
    while line:
        detected += 1
        arr = line.strip().split()
        if arr[-1] == '*':
            dynamic += 1
            if arr[0] not in dynamicDict:
                dynamicDict[arr[0]] = [ int(arr[1]) ]
            else:
                dynamicDict[arr[0]].append( int(arr[1]) )
        line = IN.readline()
    return dynamicDict, detected

def shuffleDynamics(resultFile, dynamicRatio):
    dynamicDict = {}
    IN = open(resultFile)
    line = IN.readline()
    line = IN.readline()
    while line:
        arr = line.strip().split()
        if random.random() < dynamicRatio:       #arr[-1] == '*':
            if arr[0] not in dynamicDict:
                dynamicDict[arr[0]] = [ int(arr[1]) ]
            else:
                dynamicDict[arr[0]].append( int(arr[1]) )
        line = IN.readline()
    return dynamicDict

def get_pbr_sorted_transList(pbr_ensembl):
    pbr_ensembl_transList = {}
    for protein in pbr_ensembl:
        pbr_ensembl_transList[protein] = sorted(pbr_ensembl[protein].keys())
    return pbr_ensembl_transList

def get_pbr_transList(pbr_ensembl):
    trans_list = []
    for protein in pbr_ensembl:
        trans_list += pbr_ensembl[protein].keys()
    trans_list = list(set(trans_list))
    return trans_list

def commonNum(Set1, Set2):
    count = 0
    for item in Set1:
        if item in Set2:
            count += 1
    print count

def selectDynInPBR(dynamic_ensmebl, pbr_ensembl_transList, pbr_ensembl):
    dyn_ensembl_pbr = {}
    for protein in pbr_ensembl:
        pbr_ensembl_transList[protein] = sorted(pbr_ensembl[protein].keys())
    selected_dyn = {}
    for trans_id in dynamic_ensmebl:
        dynamic_Set = copy.copy(dynamic_ensmebl[trans_id])
        for protein in pbr_ensembl_transList:
            if len(dynamic_Set) == 0: break
            if biSearch(trans_id, pbr_ensembl_transList[protein]):
                pbr_Set = pbr_ensembl[protein][trans_id]
                idx = 0
                Length = len(dynamic_Set)
                while idx < Length:
                    find = 0
                    for pbr_region in pbr_Set:
                        if pbr_region[0] < dynamic_Set[idx] < pbr_region[1]:
                            Length -= 1
                            if trans_id not in dyn_ensembl_pbr:
                                dyn_ensembl_pbr[trans_id] = [ dynamic_Set[idx] ]
                            else:
                                dyn_ensembl_pbr[trans_id].append( dynamic_Set[idx] )
                            del dynamic_Set[idx]
                            find = 1
                            break
                    if not find: idx += 1
    return dyn_ensembl_pbr

def readNatureChem_m1A(fileName, human_seq):
    mod_region_dict = {}
    IN = open(fileName)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        trans_id = arr[7]
        start = int(arr[8])
        end = int(arr[9])
        if trans_id not in human_seq:
            line = IN.readline()
            continue
        if trans_id in mod_region_dict:
            mod_region_dict[trans_id].append( (start,end) )
        else:
            mod_region_dict[trans_id] = [ (start,end) ]
        line = IN.readline()
    return mod_region_dict


def selectDynInModification(dynamic_ensmebl, m6A_ensembl, dist=10):
    sorted_List = sorted( m6A_ensembl.keys() )
    selected_dyn = {}
    for trans_id in dynamic_ensmebl:
        dynamic_Set = copy.copy(dynamic_ensmebl[trans_id])
        if biSearch(trans_id, sorted_List):
            m6A_Set = m6A_ensembl[trans_id]
            for dynamic_site in dynamic_Set:
                for m6A_site in m6A_Set:
                    if abs(dynamic_site - m6A_site) < dist:
                        if trans_id not in selected_dyn:
                            selected_dyn[trans_id] = [ dynamic_site ]
                        else:
                            selected_dyn[trans_id].append( dynamic_site )
                        break
    return selected_dyn



def selectDynIn_smallRegion(dynamic_ensmebl, m1A_ensembl):
    sorted_List = sorted( m1A_ensembl.keys() )
    selected_dyn = {}
    for trans_id in dynamic_ensmebl:
        dynamic_Set = copy.copy(dynamic_ensmebl[trans_id])
        if biSearch(trans_id, sorted_List):
            m1A_Set = m1A_ensembl[trans_id]
            for dynamic_site in dynamic_Set:
                for region in m1A_Set:
                    if region[0] <= dynamic_site <= region[1]:
                        if trans_id not in selected_dyn:
                            selected_dyn[trans_id] = [ dynamic_site ]
                        else:
                            selected_dyn[trans_id].append( dynamic_site )
                        break
    return selected_dyn


def intersect_count(ensembl1, ensembl2):
    common = 0
    for trans_id in ensembl1:
        if trans_id in ensembl2:
            for item in ensembl1[trans_id]:
                if item in ensembl2[trans_id]:
                    common += 1
    print '\n=-=-=-=--=-=-=--=-=--='
    print 'Only in Ensembl1:',pbr_count(ensembl1)-common
    print 'In Ensembl1 and Ensembl2', common
    print 'Only in Ensmebl2:',pbr_count(ensembl2)-common
    print '=-=-=-=---=-=-=-=-=-=--\n'


def combineTwoDynEnsembl(ensembl1, ensembl2):
    com_ensembl = copy.deepcopy(ensembl1)
    for trans_id in ensembl2:
        if trans_id in com_ensembl:
            com_ensembl[trans_id] = list(set( com_ensembl[trans_id] + ensembl2[trans_id] ))
        else:
            com_ensembl[trans_id] = copy.deepcopy(ensembl2[trans_id])
    return com_ensembl


def ensemblToList(ensembl, verbose=True):
    List = []
    for trans_id in ensembl:
        for site in ensembl[trans_id]:
            List.append( trans_id+'-'+str(site) )
    if verbose: print 'Data Size:',str(len(List))
    return List

###################
# Load Sequence, CLIP data
###################

human_seq = readSeq(HOME+"/lipan/DYNAMIC/GTF/transcriptome.fa")
clip_file = HOME+'/lipan/DYNAMIC/modification/HEK293/m6A/choose_good_protein/clipDB.merged.transCoor.bed'
pbr_ensembl = read_pbr(clip_file, human_seq)

###################
# Load Dynamics
###################

ch_np_dyn_file = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/HEK293/R_output/ch-np.result'
np_cy_dyn_file = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/HEK293/R_output/np-cy.result'
ch_cy_dyn_file = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/HEK293/R_output/ch-cy.result'

ch_np_dyn, ch_np_detect = readDynamics(ch_np_dyn_file)
np_cy_dyn, np_cy_detect = readDynamics(np_cy_dyn_file)
ch_cy_dyn, ch_cy_detect = readDynamics(ch_cy_dyn_file)

chnp_shuffled = shuffleDynamics(ch_np_dyn_file, 1.0*pbr_count(ch_np_dyn)/ch_np_detect)
npcy_shuffled = shuffleDynamics(np_cy_dyn_file, 1.0*pbr_count(np_cy_dyn)/np_cy_detect)


###################
# RBP & Dynamics
###################

def shuffle_pbr(dyn_FileName, dyn_ratio, pbr_ensembl, times=1000):
    pbr_cover_num_list = []
    pbr_ensembl_transList = get_pbr_sorted_transList(pbr_ensembl)
    idx = 0
    while idx < times:
        if idx%20==0: sys.stdout.write('.'); sys.stdout.flush()
        shuffled_dyn = shuffleDynamics(dyn_FileName, dyn_ratio)
        dyn_ensembl_pbr = selectDynInPBR(shuffled_dyn, pbr_ensembl_transList, pbr_ensembl)
        pbr_List = ensemblToList(dyn_ensembl_pbr, verbose=False)
        pbr_cover_num_list.append( len(pbr_List) )
        idx += 1
    return sorted(pbr_cover_num_list)


pbr_ensembl_transList = get_pbr_sorted_transList(pbr_ensembl)
pbr_trans_list = get_pbr_transList(pbr_ensembl)


# 蛋白覆盖到的转录样本与icSHAPE检测到的转录本之间的交集
commonNum(pbr_trans_list, ch_cy_dyn.keys())
chcy_dyn_ensembl_pbr = selectDynInPBR(ch_cy_dyn, pbr_ensembl_transList, pbr_ensembl)
chcy_pbr_List = ensemblToList(chcy_dyn_ensembl_pbr)

ch_cy_pbr_cover_num_list = shuffle_pbr(ch_cy_dyn_file, 1.0*pbr_count(ch_cy_dyn)/ch_cy_detect, pbr_ensembl, times=1000)

###################
# m6A & Dynamics
###################

def shuffle_m6A(dyn_FileName, dyn_ratio, m6A_ensembl, times=1000):
    m6A_cover_num_list = []
    idx = 0
    while idx < times:
        if idx%20==0: sys.stdout.write('.'); sys.stdout.flush()
        shuffled_dyn = shuffleDynamics(dyn_FileName, dyn_ratio)
        dyn_ensembl_m6A = selectDynInModification(shuffled_dyn, m6A_ensembl, dist=10)
        m6A_List = ensemblToList(dyn_ensembl_m6A, verbose=False)
        m6A_cover_num_list.append( len(m6A_List) )
        idx += 1
    return sorted(m6A_cover_num_list)


m6A_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/m6A/choose_good_protein/m6A_dataset_1.merged.transCoor.bed'
m6A_ensembl = read_m6A(m6A_file, human_seq)

commonNum(ch_cy_dyn, m6A_ensembl)
chcy_dyn_ensmebl_m6A = selectDynInModification(ch_cy_dyn, m6A_ensembl, dist=10)
intersect_count(chcy_dyn_ensembl_pbr, chcy_dyn_ensmebl_m6A)
chcy_m6A_List = ensemblToList(chcy_dyn_ensmebl_m6A)
ch_cy_m6A_cover_num_list = shuffle_m6A(ch_cy_dyn_file, 1.0*pbr_count(ch_cy_dyn)/ch_cy_detect, m6A_ensembl, times=1000)


###################
# pU & Dynamics
###################

def shuffle_pU(dyn_FileName, dyn_ratio, pU_nacent_ensembl, pU_candidate_ensembl, times=1000):
    pU_cover_num_list = []
    idx = 0
    while idx < times:
        if idx%20==0: sys.stdout.write('.'); sys.stdout.flush()
        shuffled_dyn = shuffleDynamics(dyn_FileName, dyn_ratio)
        dyn_ensmebl_nacent = selectDynInModification(shuffled_dyn, pU_nacent_ensembl)
        dyn_ensmebl_candidate = selectDynInModification(shuffled_dyn, pU_candidate_ensembl)
        dyn_ensmebl_pU = combineTwoDynEnsembl( dyn_ensmebl_nacent, dyn_ensmebl_candidate)
        pU_List = ensemblToList(dyn_ensmebl_pU, verbose=False)
        pU_cover_num_list.append( len(pU_List) )
        idx += 1
    return sorted(pU_cover_num_list)


pU_nacent_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/pU/newData/nascent_pU.transCoor.bed'
pU_candidate_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/pU/newData/pU_candidate.transCoor.bed'

pU_nacent_ensembl = read_m6A(pU_nacent_file, human_seq)
pU_candidate_ensembl = read_m6A(pU_candidate_file, human_seq)

chcy_dyn_ensmebl_old_pU = selectDynInModification(ch_cy_dyn, pU_old_ensembl)

chcy_dyn_ensmebl_nacent_pU = selectDynInModification(ch_cy_dyn, pU_nacent_ensembl)
chcy_dyn_ensmebl_candidate_pU = selectDynInModification(ch_cy_dyn, pU_candidate_ensembl)

intersect_count(chcy_dyn_ensembl_pbr, chcy_dyn_ensmebl_nacent_pU)
intersect_count(chcy_dyn_ensembl_pbr, chcy_dyn_ensmebl_candidate_pU)

chcy_dyn_ensmebl_pU = combineTwoDynEnsembl( chcy_dyn_ensmebl_nacent_pU, chcy_dyn_ensmebl_candidate_pU)
chcy_pU_List = ensemblToList(chcy_dyn_ensmebl_pU)

ch_cy_pU_cover_num_list = shuffle_pU(ch_cy_dyn_file, 1.0*pbr_count(ch_cy_dyn)/ch_cy_detect, pU_nacent_ensembl, pU_candidate_ensembl, times=1000)


###################
# m1A & Dynamics
###################


def shuffle_m1A(dyn_FileName, dyn_ratio, m1A_nature_ensembl, times=1000):
    m1A_cover_num_list = []
    idx = 0
    while idx < times:
        if idx%20==0: sys.stdout.write('.'); sys.stdout.flush()
        shuffled_dyn = shuffleDynamics(dyn_FileName, dyn_ratio)
        dyn_ensmebl_nature = selectDynInModification(shuffled_dyn, m1A_nature_ensembl)
        m1A_List = ensemblToList(dyn_ensmebl_nature, verbose=False)
        m1A_cover_num_list.append( len(m1A_List) )
        idx += 1
    return sorted(m1A_cover_num_list)



m1A_nature_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/m1A/m1A.GSE70485.valid.transCoor.bed'
m1A_nature_ensembl = read_m6A(m1A_nature_file, human_seq)

chcy_dyn_ensmebl_nature_m1A = selectDynInModification(ch_cy_dyn, m1A_nature_ensembl, dist=10)

intersect_count(chcy_dyn_ensembl_pbr, chcy_dyn_ensmebl_nature_m1A)
chcy_dyn_ensmebl_m1A = chcy_dyn_ensmebl_nature_m1A
chcy_m1A_List = ensemblToList(chcy_dyn_ensmebl_m1A)

ch_cy_m1A_cover_num_list = shuffle_m1A(ch_cy_dyn_file, 1.0*pbr_count(ch_cy_dyn)/ch_cy_detect, m1A_nature_ensembl, times=1000)


"""
###################
# HNRNPC & Dynamics
###################

RBP_list = ['HNRNPC', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'STAU1', 'FXR1', 'HNRNPU', 'TAF15', 'PTBP1', 'CPSF1', 'YTHDF1', 'EWSR1']
RBP_list.sort()

shuffle_list = []
true_list = []
for RBP in RBP_list:
    print RBP
    HNRNPC = { RBP: pbr_ensembl[RBP] }
    pbr_ensembl_transList = get_pbr_sorted_transList(HNRNPC)
    pbr_trans_list = get_pbr_transList(HNRNPC)
    len(set(pbr_trans_list)&set(ch_cy_dyn))
    chcy_dyn_ensembl_pbr = selectDynInPBR(ch_cy_dyn, pbr_ensembl_transList, HNRNPC)
    chcy_pbr_List = ensemblToList(chcy_dyn_ensembl_pbr)
    ch_cy_pbr_cover_num_list = shuffle_pbr(ch_cy_dyn_file, 1.0*pbr_count(ch_cy_dyn)/ch_cy_detect, HNRNPC, times=1000)
    
    shuffle_list.append(ch_cy_pbr_cover_num_list)
    true_list.append(len(chcy_pbr_List))


plt.figure(figsize=(20,3))
for idx, RBP in enumerate(RBP_list):
    plt.subplot(1, len(RBP_list), idx+1)
    tools.shuffle_plot(shuffle_list[idx], true_list[idx])
    plt.title(RBP)

plt.tight_layout()
plt.savefig("figs/rbp.pdf")
plt.show()
"""

