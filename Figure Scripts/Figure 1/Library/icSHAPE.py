#-*- coding:utf-8 -*-

import os,commands,sys,re,random
import seaborn as sns
import pandas as pd
import numpy
import numpy as np
import scipy
from multiprocessing import Pool
from functools import partial
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
import copy
import statsmodels
import statsmodels.sandbox.stats.multicomp
import random
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import commands

# Home Path
HOME = os.environ.get('HOME')

def readSeq(seq_name, removeVersion=True):
    # + 读取转录组序列
    transSeq = {}
    IN = open(seq_name)
    line = IN.readline()
    cur_trans = ''
    while line:
        if line[0] == '>':
            cur_trans = line[1:].split()[0]
            if removeVersion: cur_trans = cur_trans.split(".")[0]
            transSeq[ cur_trans ] = ''
        else:
            transSeq[ cur_trans ] += line.strip()
        line = IN.readline()
    return transSeq

def loadGeneName(file_name):
    "加载基因的名字"
    gene_name = {}
    IN = open(file_name)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        gene_name[ arr[0] ] = arr[1]
        line = IN.readline()
    return gene_name

def loadGTFBed(annoBedFile):
    # + 加载文件为 oganism.transCoordinate.bed
    # + Return Dict: trans_id ==> gene_id trans_len utr ....
    def parseUTR(utr_str, trans_len):
        # + 解析UTR字符串
        # + 没有5'UTR的转录本 utr_5_start=utr_5_end=0
        # + 没有3'UTR的转录本 utr_3_start=utr_3_end=trans_len+1
        utr_5_start = utr_5_end = 0 
        utr_3_start = utr_3_end = trans_len + 1
        cds_start = cds_end = 0
        if utr_str:
            utr_bak = utr_str.split(',')
            utr = []
            for utr_item in utr_bak:
                utr.append([int(i) for i in utr_item.split('-')])
            if utr[0][0] == 1:
                # + 有5'UTR
                utr_5_start = 1
                for i in range(len(utr)-1):
                    if utr[i][1] != utr[i+1][0] - 1:
                        utr_5_end = utr[i][1]
                        utr_3_start = utr[i+1][0]
                        utr_3_end = utr[-1][1]
                        break
                if utr_5_end == 0: 
                    utr_5_end = utr[-1][1]    
            else:
                # + 没有5'UTR，只有3'UTR
                utr_3_start = utr[0][0]
                utr_3_end = utr[-1][1]
        cds_start = utr_5_end + 1
        cds_end = utr_3_start - 1
        elem_coor = {'utr_5_start':utr_5_start, 'utr_5_end':utr_5_end, 'utr_3_start':utr_3_start, 'utr_3_end':utr_3_end, 'cds_start':cds_start, 'cds_end':cds_end}
        return elem_coor
    # + 解析函数
    H = open(annoBedFile)
    annoBed = dict()
    line = H.readline()
    while line:
        arr = line.strip().split()
        trans_id = arr[0]
        gene_id = arr[1]
        gene_type = arr[2]
        trans_len = int(arr[3])
        utr = ''
        if len(arr) == 6:
            utr = arr[5]
        utr_coor = parseUTR(utr, trans_len)
        annoBed[trans_id] = {}
        annoBed[trans_id]['gene_id'] = gene_id
        annoBed[trans_id]['gene_type'] = gene_type
        annoBed[trans_id]['trans_len'] = trans_len
        for it in utr_coor:
            annoBed[trans_id][it] = utr_coor[it]
        line = H.readline()
    print 'Features For Transcript: '
    print '  *'+'\n  *'.join( annoBed[annoBed.keys()[0]].keys() )
    return annoBed

def makeGene2TransDict(annoBed):
    # + 制作一个从Gene->Trans的字典
    # + 输入参数是loadGTFBed()函数返回的值
    Gene2TransDict = {}
    for trans_id in annoBed:
        gene_id = annoBed[trans_id]['gene_id']
        if gene_id not in Gene2TransDict:
            Gene2TransDict[gene_id] = []
        Gene2TransDict[gene_id].append(trans_id)
    return Gene2TransDict

def getRPKMLifePairs(trans_rpkm, gene_half_life, gtf):
    # + 获得RPKM和基因Half-Life之间的对
    # + 一个基因有多个转录本，取RPKM最大的转录本
    gene2transDict = makeGene2TransDict(gtf)
    RPKMLifePairs = []
    for gene_id in gene_half_life:
        if gene_id in gene2transDict:
            max_rpkm_trans_id = ''
            max_rpkm = 0
            for trans_id in gene2transDict[gene_id]:
                if trans_id in trans_rpkm:
                    if trans_rpkm[trans_id] > max_rpkm:
                        max_rpkm_trans_id = trans_id
                        max_rpkm = trans_rpkm[trans_id]
            if max_rpkm != 0:
                RPKMLifePairs.append( [gene_id, max_rpkm, gene_half_life[gene_id]] )
    return pd.DataFrame(RPKMLifePairs, columns=['gene_id', 'RPKM', 'life'])

def read_m6A(file_name, SeqDict):
    # + 读取m6A数据，并用转录组序列校正
    m6A_ensembl = {}
    IN = open(file_name)
    line = IN.readline()
    total_num = 0
    drop_num = 0
    while line:
        total_num += 1
        arr = line.strip().split()
        trans_id = arr[7]
        trans_loc = int(arr[8]) - 1
        if trans_id in SeqDict and len(SeqDict[trans_id]) > 100:
            if 20 <= trans_loc <= len( SeqDict[trans_id] ) - 20:
                if trans_id not in m6A_ensembl:
                    m6A_ensembl[ trans_id ] = []
                m6A_ensembl[ trans_id ].append( trans_loc )
            else:
                drop_num += 1
        else:
            drop_num += 1
        line = IN.readline()
    for trans_id in m6A_ensembl:
        m6A_ensembl[trans_id].sort()
    print "Total:", total_num, "Drop:",drop_num
    return m6A_ensembl


def read_m6A_withRatio(file_name, SeqDict, motif=False, A_Site=3):
    " 读取m6A数据，并用转录组序列校正 "
    " 返回m6A的位点和位点处的confidence "
    m6A_ensembl = {}
    IN = open(file_name)
    noMotif = True
    if motif:
        noMotif = False
        motif_len = len(motif[0])
        post_len = motif_len - A_Site
    line = IN.readline()
    total_num = 0
    drop_num = 0
    while line:
        total_num += 1
        arr = line.strip().split()
        trans_id = arr[7]
        trans_loc = int(arr[8]) - 1
        conf = float(arr[4])
        if trans_id in SeqDict and len(SeqDict[trans_id]) > 100:
            if 20 <= trans_loc <= len( SeqDict[trans_id] ) - 20:
                if noMotif or SeqDict[trans_id][(trans_loc-A_Site+1):(trans_loc+post_len+1)] in motif:
                    if trans_id not in m6A_ensembl:
                        m6A_ensembl[ trans_id ] = []
                    m6A_ensembl[ trans_id ].append( [trans_loc, conf] )
                elif not noMotif:
                    drop_num += 1
            else:
                drop_num += 1
        else:
            drop_num += 1
        line = IN.readline()
    for trans_id in m6A_ensembl:
        m6A_ensembl[trans_id].sort(key=lambda x: x[0])
    print "Total:", total_num, "Drop:",drop_num
    return m6A_ensembl


def read_m6A_withMotif(file_name, SeqDict, motif=['GGAC','GAAC'], A_Site=3):
    # + 读取m6A数据，并用转录组序列校正
    # + 只保留修饰位点处的序列motif是GGAC或GAAC的修饰位点
    # + motif中A的位置必须一致，且长度相同
    m6A_ensembl = {}
    IN = open(file_name)
    line = IN.readline()
    motif_len = len(motif[0])
    pre_len = A_Site - 1
    post_len = motif_len - A_Site
    total_num = 0
    drop_num = 0
    while line:
        total_num += 1
        arr = line.strip().split()
        trans_id = arr[7]
        trans_loc = int(arr[8]) - 1
        if trans_id in SeqDict and len(SeqDict[trans_id]) > 100:
            if 20 <= trans_loc <= len( SeqDict[trans_id] ) - 20:
                if SeqDict[trans_id][(trans_loc-A_Site+1):(trans_loc+post_len+1)] not in motif:
                    line = IN.readline()
                    drop_num += 1
                    continue
                if trans_id not in m6A_ensembl:
                    m6A_ensembl[ trans_id ] = []
                m6A_ensembl[ trans_id ].append( trans_loc )
            else:
                drop_num += 1
        else:
            drop_num += 1
        line = IN.readline()
    for trans_id in m6A_ensembl:
        m6A_ensembl[trans_id].sort()
    print "Total:", total_num, "Drop:",drop_num
    return m6A_ensembl

def load_m6A_region(file_name, SeqDict, motif='GGACT'):
    # + 载入mES的m6A数据，这些数据都是一个个的区域，需要到区域中去扫描motif
    m6A_Sites = {}
    IN = open(file_name)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        trans_id = arr[7]
        trans_start = int(arr[8]) - 1
        trans_end = int(arr[9])
        if trans_id in SeqDict:
            seq = SeqDict[trans_id][trans_start:trans_end]
            for i in range(len(seq)-len(motif)+1):
                if seq[i:(i+len(motif))] == motif:
                    if trans_id not in m6A_Sites:
                        m6A_Sites[trans_id] = []
                    m6A_Sites[trans_id].append( trans_start+i+2 )
        line = IN.readline()
    for trans_id in m6A_Sites:
        m6A_Sites[trans_id].sort()
    return m6A_Sites

def read_biased_m6A(file_name, SeqDict, min_data_size=2500, Ratio_start=0, Ratio_end=1):
    # + 读取m6A数据，并用转录组序列校正
    # + 选出其中排序从低到高[Ratio_start, Ratio_start]部分的修饰位点
    # + 并保证数据量在min_data_size之上
    m6A_raw_ensembl = []
    IN = open(file_name)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        score = float(arr[4])
        trans_id = arr[7]
        trans_loc = int(arr[8]) - 1
        if trans_id in SeqDict and len(SeqDict[trans_id]) > 100:
            if 20 <= trans_loc <= len( SeqDict[trans_id] ) - 20:
                m6A_raw_ensembl.append( [score, trans_id, trans_loc] )
        line = IN.readline()
    m6A_raw_ensembl.sort(key=lambda x: x[0])
    DataSize = len(m6A_raw_ensembl)
    sys.stdout.write("Total Data Size: "+str(DataSize)+'\n')
    if min_data_size > DataSize:
        sys.stdout.write("min_data_size shouldn't larger than DataSize\n")
        return {}
    start_loc = int(DataSize * Ratio_start)
    end_loc = int(DataSize * Ratio_end)
    sample_DataSize = end_loc - start_loc
    sys.stdout.write("Raw Sample Data Size: "+str(sample_DataSize)+'\n')
    if sample_DataSize < min_data_size:
        amount_left = min_data_size - sample_DataSize
        if start_loc < amount_left / 2:
            end_loc += amount_left - start_loc
            start_loc = 0
        elif DataSize - end_loc < amount_left / 2:
            start_loc -= amount_left - (DataSize - end_loc)
            end_loc = DataSize
        else:
            start_loc -= amount_left / 2
            end_loc += amount_left / 2
    sys.stdout.write("start_loc: "+str(start_loc)+'\n')
    sys.stdout.write("end_loc: "+str(end_loc)+'\n')
    sampled_m6A = m6A_raw_ensembl[ start_loc:end_loc ]
    m6A_ensembl = {}
    for tp in sampled_m6A:
        trans_id = tp[1]
        trans_loc = tp[2]
        if trans_id not in m6A_ensembl:
            m6A_ensembl[ trans_id ] = []
        m6A_ensembl[ trans_id ].append( trans_loc )
    for trans_id in m6A_ensembl:
        m6A_ensembl[trans_id].sort()
    return m6A_ensembl


def loadicSHAPE(file_name, SeqDict=None, removeVersion=True):
    # + 读取icSHAPE数据
    # + 并用长度信息校正
    icSHAPE = {}
    IN = open(file_name)
    line = IN.readline()
    while line:
        arr = line.strip().split()
	trans_id = arr[0]
        if removeVersion: trans_id = arr[0].split(".")[0]
        shape = arr[3:]
        if SeqDict:
            if trans_id in SeqDict:
                if len(SeqDict[trans_id]) == len(shape):
                    icSHAPE[ trans_id ] = shape
        else:
            icSHAPE[ trans_id ] = shape
        line = IN.readline()
    return icSHAPE


def readRPKM(file_name, icshape):
    # + 从icshape文件或者.rt文件中读取转录本的表达量
    # + 只保留存在于icshape字典中的转录本
    RPKM = {}
    IN = open(file_name)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        trans_id = arr[0]
        rpkm = float(arr[2])
        if trans_id in icshape:
            RPKM[trans_id] = rpkm
        line = IN.readline()
    return RPKM

def readRPKM_DMSO(file_name_1, file_name_2, icshape):
    # + 从两个.rt文件中读取转录本的表达量，取均值
    # + 只保留存在于icshape字典中的转录本
    RPKM = {}
    IN = open(file_name_1)
    line = IN.readline()
    while line:
        if line.startswith('#'):
            line = IN.readline()
            continue
        arr = line.strip().split()
        trans_id = arr[0]
        rpkm = float(arr[2])
        if trans_id in icshape:
            RPKM[trans_id] = rpkm
        line = IN.readline()
        if line:
            line = IN.readline()
    IN = open(file_name_2)
    line = IN.readline()
    while line:
        if line.startswith('#'):
            line = IN.readline()
            continue
        arr = line.strip().split()
        trans_id = arr[0]
        rpkm = float(arr[2])
        if trans_id in icshape:
            if trans_id in RPKM:
                RPKM[trans_id] = (rpkm + RPKM[trans_id])/2
            else:
                RPKM[trans_id] = rpkm
        line = IN.readline()
        if line:
            line = IN.readline()
    return RPKM

def computeGINI(my_list, valid_cutoff=10):
    # + 计算一段列表值的GINI系数，列表中可以允许有NULL
    # + 要求有效值的个数大于等于valid_cutoff
    def GINI(list_of_values):
        # + 输入值全都是非NULL值
        length = len(list_of_values)
        total = sum(list_of_values)
        if total == 0: 
            return -1
        Sorted_Array = sorted(list_of_values)
        accum, giniB = 0, 0
        for i in Sorted_Array:
            accum += i
            giniB += accum - i / 2.0
        fair_area = accum * length / 2.0
        return (fair_area - giniB) / fair_area
    floatArr = []
    for item in my_list:
        if item != 'NULL':
            floatArr.append(float(item))
    gini_index = -1
    if len(floatArr) > valid_cutoff:
        gini_index = round(GINI(floatArr),3)
    return gini_index

def computeGINI_splice(my_list, valid_cutoff=10, splice=20):
    # 把序列切成一个个片段来计算
    GINI = []
    idx = 0
    while len(my_list)-idx > splice/2:
        gini = computeGINI(my_list[idx:idx+splice], valid_cutoff)
        if gini != -1:
            GINI.append(gini)
        idx += splice
    return GINI


def computeRNAGINI_sortbyRPKM(icshape, gtf, rpkm):
    # + 对不同类别的转录本进行分类
    # + 每类转录本单独记录rpkm和gini index
    # + icshape: icshape字典
    # + gtf: loadGTFBed()读入的注释信息
    # + rpkm: readRPKM()函数读入的RPKM信息
    def mean(List):
        vList = []
        for it in List:
            if it != 'NULL':
                vList.append( float(it) )
        return numpy.mean(vList)
    transGINIDict = {}
    for trans_id in rpkm:
        if trans_id not in gtf: continue
        gene_type = gtf[trans_id]['gene_type']
        if gene_type not in transGINIDict:
            transGINIDict[gene_type]= []
        shape_v = icshape[trans_id]
        gini = computeGINI(shape_v)
        if gini != -1:
            transGINIDict[gene_type].append( [ gini, rpkm[trans_id], mean(shape_v) ] )
    for gene_type in transGINIDict:
        transGINIDict[gene_type].sort(key=lambda x: x[1])
    return transGINIDict

def generate_randA(m6A_ensembl, SeqDict, base='A'):
    # + 产生随机地A位点
    # + 调整base参数也可以生成随机地其他位点
    def allASites(seq):
        # + 返回所有的位置是A的位点坐标
        # + 0-based
        A_index = []
        for i in range(len(seq)):
            if seq[i] == base:
                A_index.append( i )
        return A_index
    rand_A_ensembl = {}
    for trans_id in m6A_ensembl:
        rand_A_ensembl[ trans_id ] = []
        seq = SeqDict[trans_id]
        ASites = allASites(seq)
        m6A_nums = len(m6A_ensembl[trans_id])
        i = 0
        while i < m6A_nums:
            randASite = random.choice(ASites)
            if 20 <= randASite <= len( SeqDict[trans_id] ) - 20:
                rand_A_ensembl[ trans_id ].append( randASite )
                i += 1
        rand_A_ensembl[ trans_id ].sort()
    return rand_A_ensembl


def generate_randA_global(m6A_ensembl, SeqDict, base='A'):
    # + 从所有的存在于m6A_ensembl的转录本中产生随机地A位点
    # + 从所有的转录本中抽样这些A碱基
    # + 调整base参数也可以生成随机地其他位点
    def allASites(seq):
        # + 返回所有的位置是A的位点坐标
        # + 0-based
        A_index = []
        for i in range(len(seq)):
            if seq[i] == base:
                A_index.append( i )
        return A_index
    all_ensembl = []
    for trans_id in SeqDict:
        seq = SeqDict[trans_id]
        ASites = allASites(seq)
        all_ensembl += [ [trans_id, i] for i in ASites if  20 <= i <= len( SeqDict[trans_id] ) - 20]
    rand_A_ensembl = {}
    sampled_ensembl = random.sample(all_ensembl, pbr_count(m6A_ensembl))
    for itPair in sampled_ensembl:
        trans_id = itPair[0]
        site = itPair[1]
        if trans_id not in rand_A_ensembl:
            rand_A_ensembl[trans_id] = []
        rand_A_ensembl[trans_id].append( site )
    return rand_A_ensembl

def generate_randA_withMotif(m6A_ensembl, SeqDict, motif=['GGAC','GAAC'], A_Site=3):
    # + 在带有m6A修饰的转录本内部产生随机地带有motif的A位点
    # + GGAC或GAAC
    # + motif中A的位置必须一致，且长度相同
    # + 避免Pos和Neg没有交叉
    def allASites_withMotif(seq):
        # + 返回所有的位置是A的位点坐标
        # + 0-based
        A_index = []
        motif_len = len(motif[0])
        for i in range(len(seq)-motif_len+1):
            if seq[i:(i+motif_len)] in motif:
                A_index.append( i + A_Site - 1 )
        return A_index
    all_randA = []
    sample_num = sum([ len(i) for i in m6A_ensembl.values() ])
    for trans_id in m6A_ensembl:
        seq = SeqDict[trans_id]
        ASites = allASites_withMotif(seq)
        for site in ASites:
            if 20 <= site <= len( SeqDict[trans_id] ) - 20 and site not in m6A_ensembl[trans_id]:
                all_randA.append( (trans_id, site) )
    sample_List = random.sample(all_randA, sample_num)
    rand_A_ensembl = {}
    for sl_item in sample_List:
        if sl_item[0] not in rand_A_ensembl:
            rand_A_ensembl[sl_item[0]] = []
        rand_A_ensembl[sl_item[0]].append( sl_item[1] )
    for trans_id in rand_A_ensembl:
        rand_A_ensembl[ trans_id ].sort()
    return rand_A_ensembl

def generate_randA_withMotif_global(m6A_ensembl, SeqDict, motif=['GGAC','GAAC'], A_Site=3):
    # + 在整个转录组上产生随机地带有motif的A位点
    # + GGAC或GAAC
    # + motif中A的位置必须一致，且长度相同
    # + 避免Pos和Neg没有交叉
    def allASites_withMotif(seq):
        # + 返回所有的位置是A的位点坐标
        # + 0-based
        A_index = []
        motif_len = len(motif[0])
        for i in range(len(seq)-motif_len+1):
            if seq[i:(i+motif_len)] in motif:
                A_index.append( i + A_Site - 1 )
        return A_index
    all_randA = []
    m6A_transList = sorted(m6A_ensembl.keys())
    sample_num = sum([ len(i) for i in m6A_ensembl.values() ])
    for trans_id in SeqDict:
        seq = SeqDict[trans_id]
        ASites = allASites_withMotif(seq)
        for site in ASites:
            if 20 <= site <= len( SeqDict[trans_id] ) - 20:
                if biSearch(trans_id, m6A_transList):
                    if site not in m6A_ensembl[trans_id]:
                        all_randA.append( (trans_id, site) )
                else:
                    all_randA.append( (trans_id, site) )
    sample_List = random.sample(all_randA, sample_num*3)
    rand_A_ensembl = {}
    for sl_item in sample_List:
        if sl_item[0] not in rand_A_ensembl:
            rand_A_ensembl[sl_item[0]] = []
        rand_A_ensembl[sl_item[0]].append( sl_item[1] )
    for trans_id in rand_A_ensembl:
        rand_A_ensembl[ trans_id ].sort()
    return rand_A_ensembl


def generate_randN(m6A_ensembl, SeqDict):
    # + 产生随机地N位点
    rand_N_ensembl = {}
    for trans_id in m6A_ensembl:
        rand_N_ensembl[ trans_id ] = []
        seq = SeqDict[trans_id]
        m6A_nums = len(m6A_ensembl[trans_id])
        for i in range(m6A_nums*2):
            randNSite = random.randint(20,len(seq)-20)
            rand_N_ensembl[ trans_id ].append( randNSite )
        rand_N_ensembl[ trans_id ].sort()
    print 'RandN Nums:', pbr_count(rand_N_ensembl)
    return rand_N_ensembl

def get_line_plot_matrix(ensembl_pos, ensembl_neg, pos_label, neg_label, length=21):
    # + Default Length 21
    # + 获得画折线图所需的DataFrame
    xx = []
    sig_points = []
    for idx_1 in range(length):
        for idx_2 in range(len(ensembl_pos[idx_1])):
            xx.append( [idx_1-int(length/2), pos_label, ensembl_pos[idx_1][idx_2] ] )
        for idx_2 in range(len(ensembl_neg[idx_1])):
            xx.append( [idx_1-int(length/2), neg_label, ensembl_neg[idx_1][idx_2] ] )
        p = scipy.stats.ttest_ind(ensembl_pos[idx_1], ensembl_neg[idx_1])[1]
        if p < 0.01:
            my_max = max( numpy.mean(ensembl_pos[idx_1]), numpy.mean(ensembl_neg[idx_1]) )
            sig_points.append( [idx_1, my_max, idx_1 ] )
    return pd.DataFrame(xx, columns=['site', 'type', 'icshape']), sig_points




def diff_scatter_points(ensembl, icshape, start=-1, end=2, valid_cutoff=3):
    # + 取出修饰位点对应的两个组分下的icSHAPE值
    # + 以便后面画散点图
    # + start包括，end不包括
    # + valid_cutoff表示当非NULL值在这个数以上时才保留
    def countNoNULL(List):
        num_of_noNULL = 0
        for it in List:
            if it != 'NULL':
                num_of_noNULL += 1
        return num_of_noNULL
    def mean(List):
        vList = []
        for it in List:
            if it != 'NULL':
                vList.append( float(it) )
        return numpy.mean(vList)
    icshape_point = {'chnp':[], 'npcy':[]}
    for trans_id in ensembl:
        for m6ASite in ensembl[trans_id]:
            ch_shape = []
            np_shape = []
            cy_shape = []
            if trans_id in icshape['ch']:
                ch_shape = icshape['ch'][trans_id][(m6ASite+start):(m6ASite+end)]
            if trans_id in icshape['np']:
                np_shape = icshape['np'][trans_id][(m6ASite+start):(m6ASite+end)]
            if trans_id in icshape['cy']:
                cy_shape = icshape['cy'][trans_id][(m6ASite+start):(m6ASite+end)]
            if countNoNULL(ch_shape) >= valid_cutoff and countNoNULL(np_shape) >= valid_cutoff:
                icshape_point['chnp'].append( (mean(ch_shape),mean(np_shape)) )
            if countNoNULL(np_shape) >= valid_cutoff and countNoNULL(cy_shape) >= valid_cutoff:
                icshape_point['npcy'].append( (mean(np_shape),mean(cy_shape)) )
    return icshape_point

def fetch_icshape_single(ensembl, single_icshape, start=0, end=5, valid_cutoff=5):
    # + 单个组分
    # + 从修饰位点的周围取出icshape值
    # + start包括，end不包括
    # + valid_cutoff表示当非NULL值在这个数以上时才保留
    def countNoNULL(List):
        num_of_noNULL = 0
        for it in List:
            if it != 'NULL':
                num_of_noNULL += 1
        return num_of_noNULL
    icshape_ensembl = []
    for trans_id in ensembl:
        for m6ASite in ensembl[trans_id]:
            shape = []
            if trans_id in single_icshape:
                shape = single_icshape[trans_id][(m6ASite+start):(m6ASite+end)]
            if countNoNULL(shape) >= valid_cutoff:
                icshape_ensembl.append(shape)
    print 'Total Nums:', len(icshape_ensembl)
    return icshape_ensembl

def fetch_icshape(ensembl, icshape, start=0, end=5, valid_cutoff=5):
    # + 从修饰位点的周围取出icshape值
    # + start包括，end不包括
    # + valid_cutoff表示当非NULL值在这个数以上时才保留
    icshape_ensembl = {'ch':[], 'np':[], 'cy':[], 'ch_vitro':[], 'np_vitro':[], 'cy_vitro':[] }
    icshape_ensembl['ch'] = fetch_icshape_single(ensembl, icshape['ch'], start, end, valid_cutoff)
    icshape_ensembl['np'] = fetch_icshape_single(ensembl, icshape['np'], start, end, valid_cutoff)
    icshape_ensembl['cy'] = fetch_icshape_single(ensembl, icshape['cy'], start, end, valid_cutoff)
    icshape_ensembl['ch_vitro'] = fetch_icshape_single(ensembl, icshape['ch_vitro'], start, end, valid_cutoff)
    icshape_ensembl['np_vitro'] = fetch_icshape_single(ensembl, icshape['np_vitro'], start, end, valid_cutoff)
    icshape_ensembl['cy_vitro'] = fetch_icshape_single(ensembl, icshape['cy_vitro'], start, end, valid_cutoff)
    return icshape_ensembl

def permutation(icshape_ensembl, sample_size=50, sample_times=200):
    # + permutate icshape from icshap_emsembl
    # + sample_size表示每一次抽取的样本数
    # + sample_times表示抽取的次数
    def computeMean(sampled_List):
        ensembl = []
        for List in sampled_List:
            for it in List:
                if it != 'NULL':
                    ensembl.append( float(it) )
        return numpy.mean(ensembl)
    if len(icshape_ensembl) < sample_size*2:
        print 'icshape_ensembl size is too small'
        return [0]*sample_times
    permutated_values = []
    for i in range(sample_times):
        sampled_values = random.sample(icshape_ensembl, sample_size)
        permutated_values.append( computeMean(sampled_values) )
    return permutated_values


def fetch_diff_icshape(ensembl, icshape, start=0, end=5, valid_cutoff=5):
    # + 从修饰位点的周围取出两个组分下icshape的差值
    # + start包括，end不包括
    # + valid_cutoff表示当非NULL值在这个数以上时才保留
    def countNoNULL(List):
        return len(List) - pd.Series(List).value_counts().get('NULL', 0)
    def diff(List1, List2):
        diff_values = []
        for tt in zip(List1, List2):
            if tt[0] != 'NULL' and tt[1] != 'NULL':
                diff_values.append( abs(float(tt[0])-float(tt[1])) )
            else:
                diff_values.append('NULL')
        return diff_values
    icshape_ensembl = {'chnp':[], 'npcy':[]}
    for trans_id in ensembl:
        for m6ASite in ensembl[trans_id]:
            ch_shape = []
            np_shape = []
            cy_shape = []
            if trans_id in icshape['ch']:
                ch_shape = icshape['ch'][trans_id][(m6ASite+start):(m6ASite+end)]
            if trans_id in icshape['np']:
                np_shape = icshape['np'][trans_id][(m6ASite+start):(m6ASite+end)]
            if trans_id in icshape['cy']:
                cy_shape = icshape['cy'][trans_id][(m6ASite+start):(m6ASite+end)]
            if countNoNULL(ch_shape) >= valid_cutoff and countNoNULL(np_shape) >= valid_cutoff:
                diff_vs = diff(ch_shape, np_shape)
                if countNoNULL(diff_vs) > valid_cutoff:
                    icshape_ensembl['chnp'].append( diff_vs )
            if countNoNULL(np_shape) >= valid_cutoff and countNoNULL(cy_shape) >= valid_cutoff:
                diff_vs = diff(np_shape, cy_shape)
                if countNoNULL(diff_vs) > valid_cutoff:
                    icshape_ensembl['npcy'].append( diff_vs )
                icshape_ensembl['npcy'].append( diff(np_shape, cy_shape) )
    return icshape_ensembl

def fetch_sites_ensembl_single(icshape_ensembl_single, length):
    # + 单个组分
    # + 得到icSHAP_ensembl各个位点上的值的集合
    site_ensembl = []
    for idx in range(length):
        site_ensembl.append( [] )
    if length != len(icshape_ensembl_single[0]):
        print "The length is not equal"
        return site_ensembl
    for shape_v in icshape_ensembl_single:
        for idx in range(len(shape_v)):
            if shape_v[idx] != 'NULL':
                site_ensembl[idx].append( float(shape_v[idx]) )
    return site_ensembl



def fetch_sites_ensembl(icshape_ensembl, length):
    # + 得到icSHAP_ensembl各个位点上的值的集合
    site_ensembl = {'ch':[], 'np':[], 'cy':[], 'ch_vitro':[], 'np_vitro':[], 'cy_vitro':[]}
    site_ensembl['ch'] = fetch_sites_ensembl_single(icshape_ensembl['ch'], length=length)
    site_ensembl['np'] = fetch_sites_ensembl_single(icshape_ensembl['np'], length=length)
    site_ensembl['cy'] = fetch_sites_ensembl_single(icshape_ensembl['cy'], length=length)
    site_ensembl['ch_vitro'] = fetch_sites_ensembl_single(icshape_ensembl['ch_vitro'], length=length)
    site_ensembl['np_vitro'] = fetch_sites_ensembl_single(icshape_ensembl['np_vitro'], length=length)
    site_ensembl['cy_vitro'] = fetch_sites_ensembl_single(icshape_ensembl['cy_vitro'], length=length)
    return site_ensembl


def fetch_diff_DataFrame_and_sig(pos_sites_ensembl, neg_sites_ensembl, pos_label, neg_label):
    # + pos_sites_ensembl是正样本（或其中一个样本的各个位点的Ensembl）
    # + neg_sites_ensembl是负样本（或另外一个样本的各个位点的Ensembl）
    # + pos_label/neg_label是正负样本标签
    xx = []
    sig_points = []
    if len(pos_sites_ensembl) != len(neg_sites_ensembl):
        sys.stdout.write("Pos != Neg\n")
        return [], []
    Seq_Len = len(pos_sites_ensembl)
    for idx_1 in range(Seq_Len):
        for idx_2 in range(len( pos_sites_ensembl[idx_1] )):
            xx.append( [idx_1-Seq_Len/2, pos_label, pos_sites_ensembl[idx_1][idx_2]] )
        for idx_2 in range(len( neg_sites_ensembl[idx_1] )):
            xx.append( [idx_1-Seq_Len/2, neg_label, neg_sites_ensembl[idx_1][idx_2]] )
        p = scipy.stats.ttest_ind(pos_sites_ensembl[idx_1], neg_sites_ensembl[idx_1])[1]
        if p < 0.01:
            my_max = max( numpy.mean(pos_sites_ensembl[idx_1]), numpy.mean(neg_sites_ensembl[idx_1]) )
            sig_points.append( [idx_1, my_max+0.04] )
    return pd.DataFrame(xx, columns=['site', 'type', 'icshape']), sig_points


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
        if biSearch(trans_id, sorted_transList) and len(SeqDict[trans_id]) > 100:
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


def generate_randPBR_single(single_pbr_ensembl, SeqDict):
    # + 产生随机的PBR(单个蛋白)
    rand_pbr_ensembl = {}
    for trans_id in single_pbr_ensembl:
        rand_pbr_ensembl[ trans_id ] = []
        seq = SeqDict[trans_id]
        pbr_nums = len(single_pbr_ensembl[trans_id])
        for i in range(pbr_nums):
            region_len = single_pbr_ensembl[trans_id][i][1] - single_pbr_ensembl[trans_id][i][0]
            randPBRStart = random.randint(20,len(seq)-20-region_len)
            rand_pbr_ensembl[ trans_id ].append( [randPBRStart, randPBRStart+region_len] )
        rand_pbr_ensembl[ trans_id ].sort(key=lambda x: x[0])
    return rand_pbr_ensembl

def generate_randPBR(pbr_ensembl, SeqDict):
    # + 产生随机的PBR(所有蛋白)
    randPBR_ensembl = {}
    for protein in pbr_ensembl:
        randPBR_ensembl[protein] = generate_randPBR_single( pbr_ensembl[protein], SeqDict )
    return randPBR_ensembl


def biSearch(item, Set):
    # + 二分查找
    # + Set必需事先从小到大排序
    start = 0
    end = len(Set) - 1
    while start <= end:
        middle = (start + end) / 2
        if Set[middle] < item:
            start = middle + 1
        elif Set[middle] > item:
            end = middle - 1
        else:
            return True
    return False

def fetch_icshape_pbr_single(single_pbr_ensembl, protein_name, icshape, hang=0, hangRatio=0, valid_cutoff=10):
    # + 从蛋白结合区取出icSHAPE值（单个蛋白）
    # + valid_cutoff表示当非NULL值在这个数以上时才保留
    # + hang表示蛋白结合区两端的N个碱基值不算在里面，缩小PBR
    # + hangRatio表示蛋白结合区两端的hangRatio这个比例的碱基值不算在里面
    # + 两个都有时，先hang，然后再hangRatio
    def countNoNULL(List):
        num_of_noNULL = 0
        for it in List:
            if it != 'NULL':
                num_of_noNULL += 1
        return num_of_noNULL
    icshape_ensembl = {'ch':[], 'np':[], 'cy':[], 'ch_vitro':[], 'np_vitro':[], 'cy_vitro':[]}
    totol_num = len(single_pbr_ensembl)
    i = 0
    for trans_id in single_pbr_ensembl:
        if i % 1000 == 0:
            print protein_name+": %.2f%%" % (1.0*i/totol_num*100, )
        i += 1
        for pbr_s in single_pbr_ensembl[trans_id]:
            ch_shape = []
            np_shape = []
            cy_shape = []
            ch_vitro_shape = []
            np_vitro_shape = []
            cy_vitro_shape = []
            if pbr_s[1] - pbr_s[0] <= 2*hang:
                continue
            start = pbr_s[0] + hang
            end = pbr_s[1] - hang
            start = int(start + (end-start)*hangRatio)
            end = int(end - (end-start)*hangRatio)
            "vivo"
            if biSearch(trans_id, icshape['trans_ch']):# trans_id in icshape['ch']:
                ch_shape = icshape['ch'][trans_id][start:end]
            if biSearch(trans_id, icshape['trans_np']):# trans_id in icshape['np']:
                np_shape = icshape['np'][trans_id][start:end]
            if biSearch(trans_id, icshape['trans_cy']):# trans_id in icshape['cy']:
                cy_shape = icshape['cy'][trans_id][start:end]
            "vitro"
            if biSearch(trans_id, icshape['trans_ch_vitro']):# trans_id in icshape['cy']:
                ch_vitro_shape = icshape['ch_vitro'][trans_id][start:end]
            if biSearch(trans_id, icshape['trans_np_vitro']):# trans_id in icshape['cy']:
                np_vitro_shape = icshape['np_vitro'][trans_id][start:end]
            if biSearch(trans_id, icshape['trans_cy_vitro']):# trans_id in icshape['cy']:
                cy_vitro_shape = icshape['cy_vitro'][trans_id][start:end]
            "vivo"
            if countNoNULL(ch_shape) >= valid_cutoff:
                icshape_ensembl['ch'].append(ch_shape)
            if countNoNULL(np_shape) >= valid_cutoff:
                icshape_ensembl['np'].append(np_shape)
            if countNoNULL(cy_shape) >= valid_cutoff:
                icshape_ensembl['cy'].append(cy_shape)
            "vitro"
            if countNoNULL(ch_vitro_shape) >= valid_cutoff:
                icshape_ensembl['ch_vitro'].append(ch_vitro_shape)
            if countNoNULL(np_vitro_shape) >= valid_cutoff:
                icshape_ensembl['np_vitro'].append(np_vitro_shape)
            if countNoNULL(cy_vitro_shape) >= valid_cutoff:
                icshape_ensembl['cy_vitro'].append(cy_vitro_shape)
    print protein_name+" finished!"
    return icshape_ensembl

# hnRNPC = fetch_icshape_pbr_single(pbr_ensembl['HNRNPC'], 'HNRNPC', icshape)

def get_pbr_cover_m6A(pbr_ensembl, m6A_ensembl, extend=0):
    # + 获得所有与m6A有交集的PBR
    # + extend表示把m6A位点左右延长成一个区域
    m6A_trans_Set = sorted(m6A_ensembl.keys())
    pbr_ensembl_cover_m6A = {}
    for protein in pbr_ensembl:
        pbr_ensembl_cover_m6A[protein] = {}
        for trans_id in pbr_ensembl[protein]:
            if not biSearch(trans_id, m6A_trans_Set):
                # + 这个转录本不是m6A修饰的其中一个转录本
                continue
            for pbr in pbr_ensembl[protein][trans_id]:
                pbr_start = pbr[0]
                pbr_end = pbr[1]
                for mod_site in m6A_ensembl[trans_id]:
                    mod_start = mod_site - extend
                    mod_end = mod_site + extend
                    if mod_start <= pbr_end and pbr_start <= mod_end:
                        # + 覆盖修饰位点
                        if trans_id not in pbr_ensembl_cover_m6A[protein]:
                            pbr_ensembl_cover_m6A[protein][trans_id] = []
                        pbr_ensembl_cover_m6A[protein][trans_id].append(pbr)
                        break
    return pbr_ensembl_cover_m6A

def buildMotifIndex(SeqDict, motif=['GGAC','GAAC','AGAC','AAAC'], Plus_Site=3):
    # + 从所有的转录本中获得motif的位置
    # + Plus_Site表示返回的坐标从motif的起点加上这个偏移
    # + motif的长度必需一致
    def allMotifSites(seq):
        # + 返回所有的位置是motif的位点坐标加上偏移量
        # + 0-based
        motif_index = []
        motif_len = len(motif[0])
        for i in range(len(seq)-motif_len+1):
            if seq[i:(i+motif_len)] in motif:
                motif_index.append( i + Plus_Site - 1 )
        return motif_index
    MotifIndex = {}
    for trans_id in SeqDict:
        motif_idx = allMotifSites(SeqDict[trans_id])
        if len(motif_idx) >= 1:
            MotifIndex[trans_id] = motif_idx
    return MotifIndex

def get_allMotif_fromTransSet(TransList, SeqDict, motifIndex, sample_num=False):
    # + 从一套转录本的集合中得到某个motif的集合
    # + sample_num表示从所有的结果中采样的数量
    motifDict = {}
    for trans_id in TransList:
        if trans_id in motifIndex:
            sites = motifIndex[trans_id]
            if len(sites) > 0:
                motifDict[trans_id] = sites
    total_count = sum([len(i) for i in motifDict.values()])
    print 'Total Size: ' + str(total_count)
    if sample_num:
        if sample_num >= total_count:
            print 'Warning: sample size is equal to total size'
            return motifDict
        sample_prob = 1.0 * sample_num / total_count
        sampledDict = {}
        for trans_id in motifDict:
            for site in motifDict[trans_id]:
                if random.random() <= sample_prob:
                    if trans_id not in sampledDict:
                        sampledDict[trans_id] = []
                    sampledDict[trans_id].append( site )
        sampled_count = sum([len(i) for i in sampledDict.values()])
        print 'Sampled Size: '+str(sampled_count)
        return sampledDict
    else:
        return motifDict

def get_pbr_MotifSet(pbr_MinusTransList, SeqDict, motifIndex):
    # + 获得所有蛋白结合转录本中的motif位置
    psudo_m6A_ensembl = {}
    for protein in pbr_MinusTransList:
        psudo_m6A_ensembl[protein] = get_allMotif_fromTransSet(pbr_MinusTransList[protein], SeqDict, motifIndex=motifIndex, sample_num=False)
    return psudo_m6A_ensembl

def get_pbr_cover_m6A_sep_m6A_single(pbr, m6A_ensembl, extend=0):
    # + 获得所有与m6A有交集的PBR
    # + extend表示把m6A位点左右延长成一个区域
    m6A_trans_Set = sorted(m6A_ensembl.keys())
    pbr_cover_m6A = {}
    for trans_id in pbr:
        if not biSearch(trans_id, m6A_trans_Set):
            # + 这个转录本不是m6A修饰的其中一个转录本
            continue
        for region in pbr[trans_id]:
            pbr_start = region[0]
            pbr_end = region[1]
            for mod_site in m6A_ensembl[trans_id]:
                mod_start = mod_site - extend
                mod_end = mod_site + extend
                if mod_start <= pbr_end and pbr_start <= mod_end:
                    # + 覆盖修饰位点
                    if trans_id not in pbr_cover_m6A:
                        pbr_cover_m6A[trans_id] = []
                    pbr_cover_m6A[trans_id].append(region)
                    break
    return pbr_cover_m6A


def get_pbr_cover_m6A_sep_m6A(pbr_ensembl, psudo_m6A_ensembls, extend=0):
    # + psudo_m6A_ensembls是每一个蛋白都有一个单独的List
    # + 获得所有与m6A有交集的PBR
    # + extend表示把m6A位点左右延长成一个区域
    pbr_ensembl_cover_m6A = {}
    for protein in pbr_ensembl:
        pbr_ensembl_cover_m6A[protein] = get_pbr_cover_m6A_sep_m6A_single(pbr_ensembl[protein], psudo_m6A_ensembls[protein], extend=extend)
    return pbr_ensembl_cover_m6A


def getout_m6A_overlap_with_pbr(pbr_ensembl_single, m6A_ensembl, conf_cutoff=2.8 ,extend=60):
    # + 从数据集中取出与置信度高于某个阈值的pbr有交集的m6A位点
    pbr_overlap_m6A_ensembl = {}
    pbr_tranSet = sorted(pbr_ensembl_single.keys())
    for trans_id in m6A_ensembl:
        if biSearch(trans_id, pbr_tranSet):
            for site in m6A_ensembl[trans_id]:
                m6A_start = site - extend
                m6A_end = site + extend
                confSet = []
                for region in pbr_ensembl_single[trans_id]:
                    pbr_start = region[0]
                    pbr_end = region[1]
                    pbr_conf = region[2]
                    if pbr_conf < conf_cutoff:
                        continue
                    if m6A_start <= pbr_end and pbr_start <= m6A_end:
                        if trans_id not in pbr_overlap_m6A_ensembl:
                            pbr_overlap_m6A_ensembl[trans_id] = []
                        confSet.append( pbr_conf )
                if len(confSet) >= 1:
                    pbr_overlap_m6A_ensembl[trans_id].append([site, max(confSet) ])
    return pbr_overlap_m6A_ensembl

def get_icshape_df_m6A_Site(m6A_ensembl, ShapeDict, extend=5):
    # + 从一套m6A修饰位点中取出两端的icSHAPE值
    def noNULL(List):
        for it in List:
            if it == 'NULL':
                return False
        return True
    
    icshape_v = []
    Confs = []
    for trans_id in m6A_ensembl:
        sites = m6A_ensembl[trans_id]
        for siteInfo in sites:
            site = siteInfo[0]
            pbrConf = siteInfo[1]
            shape = ShapeDict[trans_id][(site-extend):(site+extend+1)]
            if noNULL(shape) and len(shape) == 2*extend+1:
                icshape_v.append([float(i) for i in shape])
                Confs.append(pbrConf)
    
    xx = pd.DataFrame(icshape_v, columns=[str(i-extend) for i in range(2*extend+1)])
    return xx, Confs

def ColorMap(floatList):
    # + 给出一段数字，返回一串颜色
    fList = pd.Series(floatList)
    fList = fList/sorted(floatList)[int(0.9*len(floatList))]
    fList[fList>1] = 1
    colors = []
    for it in list(fList):
        colors.append([1,1-it,1-it])
    return pd.Series(colors)


def pbr_count(pbr_single):
    # + 得到某一个pbr的总结合区域个数
    return sum([len(v) for v in pbr_single.values()])

def get_OnlyPBR_trans_single(pbr, TransList):
    # + 获取单个只有PBR结合而m6A没有结合的转录本
    sort_TransList = sorted(TransList)
    pbr_MinusTransList = []
    for trans_id in pbr:
        if not biSearch(trans_id, sort_TransList):
            pbr_MinusTransList.append(trans_id)
    return pbr_MinusTransList

def get_OnlyPBR_trans(pbr_ensembl, TransList):
    # + 获取只有PBR结合而m6A没有结合的转录本
    pbr_MinusTransList = {}
    for protein in pbr_ensembl:
        pbr_MinusTransList[protein] = get_OnlyPBR_trans_single(pbr_ensembl[protein], TransList)
    return pbr_MinusTransList

def balance_two_dataset(ref_pbr_ensembl, psudo_pbr_ensembl):
    # + 使用ref_pbr_ensembl来平衡数据
    # + ref_pbr_ensembl是参考的pbr数据，psudo_pbr_ensembl是我们要调整数据量的数据
    balanced_pbr_ensembl = {}
    for protein in ref_pbr_ensembl:
        pbr_num = sum([ len(i) for i in ref_pbr_ensembl[protein].values() ])
        psudo_pbr_num = sum([ len(i) for i in psudo_pbr_ensembl[protein].values() ])
        if psudo_pbr_num > pbr_num*1.5:
            sample_prob = 1.0 * pbr_num / psudo_pbr_num
            balanced_pbr_ensembl[protein] = {}
            for trans_id in psudo_pbr_ensembl[protein]:
                for region in psudo_pbr_ensembl[protein][trans_id]:
                    if random.random() <= sample_prob:
                        if trans_id not in balanced_pbr_ensembl[protein]:
                            balanced_pbr_ensembl[protein][trans_id] = []
                        balanced_pbr_ensembl[protein][trans_id].append( region )
        else:
            balanced_pbr_ensembl[protein] = psudo_pbr_ensembl[protein]
    return balanced_pbr_ensembl


def fetch_conf_from_PBR_Ensembl_single(pbr, func):
    conf = []
    for trans_id in pbr:
        conf += [ func(region[2]) for region in pbr[trans_id] ]
    conf.sort()
    return conf

def fetch_conf_from_PBR_Ensembl(pbr_ensembl, func=numpy.log):
    # + 从PBR中取出所有的置信度
    # + func表示对数据做转换
    conf = {}
    for protein in pbr_ensembl:
        conf[protein] = fetch_conf_from_PBR_Ensembl_single(pbr_ensembl[protein], func=func )
    return conf

def fetch_icshape_pbr(pbr_ensembl, icshape, hang=0, hangRatio=0, valid_cutoff=10):
    # + 从蛋白结合区取出icSHAPE值（所有蛋白）
    # + 非常耗时
    icshape_ensembl = {}
    for protein in pbr_ensembl:
        print protein
        icshape_ensembl[protein] = fetch_icshape_pbr_single(pbr_ensembl[protein], protein, icshape, hang=hang, hangRatio=hangRatio, valid_cutoff=valid_cutoff)
    return icshape_ensembl

def fetch_icshape_pbr_paralle(pbr_ensembl, icshape, hang=0, hangRatio=0, valid_cutoff=10):
    # + 从蛋白结合区取出icSHAPE值（所有蛋白并行版）
    # + 非常耗时
    num_cores = 3       # ? number of cores on your machine
    pbr_ensembl_split = []
    for protein in pbr_ensembl:
        pbr_ensembl_split.append( {protein: pbr_ensembl[protein]} )
    fetch_icshape_pbr_model = partial(fetch_icshape_pbr, icshape=icshape, valid_cutoff=valid_cutoff, hang=hang, hangRatio=hangRatio)
    pool = Pool(num_cores)
    print 'start...'
    icshape_ensembl_list = pool.map(fetch_icshape_pbr_model, pbr_ensembl_split, 1)
    print 'finished...'
    pool.close()
    pool.join()
    icshape_ensembl = {}
    for it in icshape_ensembl_list:
        if len(it) != 1:
            print 'Error'
            return 0
        protein = it.keys()[0]
        icshape_ensembl[ protein ] = it[protein]
    return icshape_ensembl

def save_pbr_icshape(pbr_icshape, file_name):
    # + 把上面fetch_icshape_pbr得到的结果存储在文件中
    OUT = open("/Share/home/zhangqf/lipan/tmp/"+file_name, 'w')
    for protein in pbr_icshape:
        print >>OUT, '>>>'+protein
        for compartment in pbr_icshape[protein]:
            print >>OUT, '$$$'+compartment
            for icshape_v in pbr_icshape[protein][compartment]:
                print >>OUT, '\t'.join(icshape_v)
    OUT.close()


def read_pbr_icshape(file_name):
    # + 从save_pbr_icshape存储的文件中读取数据
    pbr_icshape = {}
    IN = open("/Share/home/zhangqf/lipan/tmp/"+file_name)
    line = IN.readline()
    cur_protein = ''
    cur_compartment = ''
    while line:
        if line.startswith('>>>'):
            cur_protein = line.strip('\n>')
            pbr_icshape[cur_protein] = {}
        elif line.startswith('$$$'):
            cur_compartment = line.strip('\n$')
            pbr_icshape[cur_protein][cur_compartment] = []
        else:
            pbr_icshape[cur_protein][cur_compartment].append( line.strip().split() )
        line = IN.readline()
    return pbr_icshape


def permutation_pbr(pbr_icshape_ensembl, sample_size=50, sample_times=200):
    # + 这个针对蛋白的permutation
    # + 与m6A不同
    # ? 未测试
    pbr_permutated_values = {}
    for protein in pbr_icshape_ensembl:
        sys.stdout.write(protein+': ')
        pbr_permutated_values[protein] = {}
        for compartment in pbr_icshape_ensembl[protein]:
            sys.stdout.write(compartment+'\t')
            pbr_permutated_values[protein][compartment] = permutation(pbr_icshape_ensembl[protein][compartment], sample_size, sample_times)
        sys.stdout.write('\n')
        sys.stdout.flush()
    return pbr_permutated_values


def permutation_global_pbr(pbr_icshape_ensembl, sample_size=50, sample_times=200):
    # + 这个函数从所有的蛋白结合位点中抽取出一部分结合位点做permutation
    # + 计算所有蛋白结合位点的单双链偏好性
    sampled_pbr = {}
    # + 首先统计pbr总量
    total_data = {'ch':[], 'np':[], 'cy':[], 'ch_vitro':[], 'np_vitro':[], 'cy_vitro':[]}
    for protein in pbr_icshape_ensembl:
        for compartment in pbr_icshape_ensembl[protein]:
            total_data[compartment] += pbr_icshape_ensembl[protein][compartment]
    # + 抽样
    protein_num = len(pbr_icshape_ensembl)
    "vivo"
    sampled_pbr['ch'] = random.sample(total_data['ch'], len(total_data['ch'])/protein_num)
    sampled_pbr['np'] = random.sample(total_data['np'], len(total_data['np'])/protein_num)
    sampled_pbr['cy'] = random.sample(total_data['cy'], len(total_data['cy'])/protein_num)
    "vitro"
    sampled_pbr['ch_vitro'] = random.sample(total_data['ch_vitro'], len(total_data['ch_vitro'])/protein_num)
    sampled_pbr['np_vitro'] = random.sample(total_data['np_vitro'], len(total_data['np_vitro'])/protein_num)
    sampled_pbr['cy_vitro'] = random.sample(total_data['cy_vitro'], len(total_data['cy_vitro'])/protein_num)
    pbr_permutated_values = {}
    "vivo"
    pbr_permutated_values['ch'] = permutation(sampled_pbr['ch'], sample_size, sample_times)
    pbr_permutated_values['np'] = permutation(sampled_pbr['np'], sample_size, sample_times)
    pbr_permutated_values['cy'] = permutation(sampled_pbr['cy'], sample_size, sample_times)
    "vitro"
    pbr_permutated_values['ch_vitro'] = permutation(sampled_pbr['ch_vitro'], sample_size, sample_times)
    pbr_permutated_values['np_vitro'] = permutation(sampled_pbr['np_vitro'], sample_size, sample_times)
    pbr_permutated_values['cy_vitro'] = permutation(sampled_pbr['cy_vitro'], sample_size, sample_times)
    return pbr_permutated_values



def fetch_diff_icshape_pbr_single(single_pbr_ensembl, protein_name, icshape, hang=0, hangRatio=0, valid_cutoff=10):
    # + 把蛋白结合区域的不同组分下的icSHAPE差值取出
    # + valid_cutoff表示当非NULL值在这个数以上时才保留
    # + hang表示蛋白结合区两端的N个碱基值不算在里面，缩小PBR
    # + hangRatio表示蛋白结合区两端的hangRatio这个比例的碱基值不算在里面
    # + 两个都有时，先hang，然后再hangRatio
    def countNoNULL(List):
        num_of_noNULL = 0
        for it in List:
            if it != 'NULL':
                num_of_noNULL += 1
        return num_of_noNULL
    def diff(List1, List2):
        diff_values = []
        for tt in zip(List1, List2):
            if tt[0] != 'NULL' and tt[1] != 'NULL':
                diff_values.append( abs(float(tt[0])-float(tt[1])) )
            else:
                diff_values.append('NULL')
        return diff_values
    def biSearch(item, Set):
        start = 0
        end = len(Set) - 1
        while start <= end:
            middle = (start + end) / 2
            if Set[middle] < item:
                start = middle + 1
            elif Set[middle] > item:
                end = middle - 1
            else:
                return True
        return False
    icshape_ensembl = {'chnp':[], 'npcy':[]}
    totol_num = len(single_pbr_ensembl)
    i = 0
    for trans_id in single_pbr_ensembl:
        if i % 1000 == 0:
            print protein_name+": %.2f%%" % (1.0*i/totol_num*100, )
        i += 1
        for pbr_s in single_pbr_ensembl[trans_id]:
            ch_shape = []
            np_shape = []
            cy_shape = []
            if pbr_s[1] - pbr_s[0] <= 2*hang:
                continue
            start = pbr_s[0] + hang
            end = pbr_s[1] - hang
            start = int(start + (end-start)*hangRatio)
            end = int(end - (end-start)*hangRatio)
            if biSearch(trans_id, icshape['trans_ch']):
                ch_shape = icshape['ch'][trans_id][start:end]
            if biSearch(trans_id, icshape['trans_np']):
                np_shape = icshape['np'][trans_id][start:end]
            if biSearch(trans_id, icshape['trans_cy']):
                cy_shape = icshape['cy'][trans_id][start:end]
            if countNoNULL(ch_shape) >= valid_cutoff and countNoNULL(np_shape) >= valid_cutoff:
                diff_vs = diff(ch_shape, np_shape)
                if countNoNULL(diff_vs) > valid_cutoff:
                    icshape_ensembl['chnp'].append( diff_vs )
            if countNoNULL(np_shape) >= valid_cutoff and countNoNULL(cy_shape) >= valid_cutoff:
                diff_vs = diff(np_shape, cy_shape)
                if countNoNULL(diff_vs) > valid_cutoff:
                    icshape_ensembl['npcy'].append( diff_vs )
                icshape_ensembl['npcy'].append( diff(np_shape, cy_shape) )
    return icshape_ensembl

# hnRNPC = fetch_diff_icshape_pbr_single(pbr_ensembl['HNRNPC'],'HNRNPC' , icshape, hangRatio=0.2, valid_cutoff=8)

def fetch_diff_icshape_pbr(pbr_ensembl, icshape, hang=0, hangRatio=0, valid_cutoff=10):
    # + 从蛋白结合区取出两个组分下icSHAPE差值（所有蛋白）
    # + 非常耗时
    diff_icshape_ensembl = {}
    for protein in pbr_ensembl:
        print protein
        diff_icshape_ensembl[protein] = fetch_diff_icshape_pbr_single(pbr_ensembl[protein], protein, icshape, hang=hang, hangRatio=hangRatio, valid_cutoff=valid_cutoff)
    return diff_icshape_ensembl


def pickUpSig(List, min_len=100, max_portion=0.4, step_len=10, mode="mann"):
    # + 从一串有序的数组中计算出两端差异最显著的组合即显著值
    p_v_method = scipy.stats.mannwhitneyu if mode == 'mann' else scipy.stats.ttest_ind
    list_len = len(List)
    if list_len * max_portion <= min_len:
        return None, None
    cur_len = min_len
    min_pValue = 1
    len_with_mPV = cur_len
    while list_len * max_portion > cur_len:
        preList = List[:cur_len]
        postList = List[-cur_len:]
        pV = p_v_method(preList, postList)[1]
        #print pV, cur_len
        if pV < min_pValue:
            min_pValue = pV
            len_with_mPV = cur_len
        cur_len += step_len
    return min_pValue, len_with_mPV

def percentile(List, percent):
    "相当于numpy.percentile"
    sort_list = sorted(List)
    list_len = len(sort_list)
    sample_site = int(percent * list_len)
    return sort_list[sample_site]

def normalize_pbr_conf(pbr_ensembl, up_ratio=0.95, down_ratio=0.05, mode='normal'):
    "把所有蛋白的结合强度进行归一化成 0-1 "
    for protein in pbr_ensembl:
        print 'Process '+protein+'...'
        conf_Set = []
        for trans_id in pbr_ensembl[protein]:
            for region in pbr_ensembl[protein][trans_id]:
                conf_Set.append( region[2] )
        conf_Set.sort()
        if mode == 'normal':
            upper = percentile(conf_Set, up_ratio)
            lower = percentile(conf_Set, down_ratio)
            for trans_id in pbr_ensembl[protein]:
                for idx in range(len(pbr_ensembl[protein][trans_id])):
                    conf = pbr_ensembl[protein][trans_id][idx][2]
                    if conf > upper:
                        pbr_ensembl[protein][trans_id][idx][2] = upper
                    elif conf < lower:
                        pbr_ensembl[protein][trans_id][idx][2] = lower
                    pbr_ensembl[protein][trans_id][idx][2] = 1.0 * (pbr_ensembl[protein][trans_id][idx][2] - lower) / (upper - lower + 0.001)
        elif mode == 'scalling':
            scaling_factor = numpy.mean( conf_Set[int(18.0/20*len(conf_Set)):int(19.0/20*len(conf_Set))] )
            upper = percentile( [float(i)/scaling_factor for i in conf_Set] , up_ratio)
            lower = percentile( [float(i)/scaling_factor for i in conf_Set] , down_ratio)
            for trans_id in pbr_ensembl[protein]:
                for idx in range(len(pbr_ensembl[protein][trans_id])):
                    conf = 1.0 * pbr_ensembl[protein][trans_id][idx][2] / scaling_factor
                    if conf > upper:
                        pbr_ensembl[protein][trans_id][idx][2] = upper
                    elif conf < lower:
                        pbr_ensembl[protein][trans_id][idx][2] = lower
                    pbr_ensembl[protein][trans_id][idx][2] = 1.0 * (conf - lower) / (upper - lower + 0.001)


def sub_pbr_ensembl(raw_pbr_ensembl, SeqDict):
    "取pbr集合中的转录本子集, 即存在于SeqDict中的集合"
    sorted_Seq_transList = sorted(SeqDict.keys())
    small_pbr_ensembl = {}
    for protein in raw_pbr_ensembl:
        small_pbr_ensembl[protein] = {}
        for trans_id in raw_pbr_ensembl[protein]:
            if biSearch(trans_id, sorted_Seq_transList):
                small_pbr_ensembl[protein][trans_id] = raw_pbr_ensembl[protein][trans_id]
    return small_pbr_ensembl

def vennTwoSet(List1, List2):
    only_1 = []
    common = []
    List2_sorted = sorted(List2)
    only_2 = copy.copy(List2)
    for item in List1:
        if biSearch(item, List2_sorted):
            common.append(item)
            only_2.remove(item)
        else:
            only_1.append(item)
    return (only_1, common, only_2)

# (only_1, common, only_2) = vennTwoSet([1,2,3,4,5], [4,5,6,7,8])

def cmp(Item1, Item2):
    if Item1[0] > Item2[0]:
        return 1
    elif Item1[0] == Item2[0]:
        if Item1[1] > Item2[1]:
            return 1
        else:
            return -1
    else:
        return -1

"""
merged_region_list = [ [1,5, 0.2], [1,5, 0.4], [2, 6, 0.6], [9,20, 12], [7,9,0.4] ]
merged_region_list.sort(cmp=cmp)
"""

def merge_region_list(regionList1, regionList2):
    merged_region_list = regionList1 + regionList2
    merged_region_list.sort(cmp=cmp)
    idx = 0
    Length = len(merged_region_list)
    while idx < Length - 1:
        if merged_region_list[idx][1] >= merged_region_list[idx+1][0]:
            merged_region_list[idx][1] = merged_region_list[idx+1][1]
            if type(merged_region_list[idx][2]) == list:
                merged_region_list[idx][2].append( merged_region_list[idx+1][2] )
            else:
                merged_region_list[idx][2] = [ merged_region_list[idx][2], merged_region_list[idx+1][2] ]
            del merged_region_list[idx+1]
            Length -= 1
        else:
            idx += 1
    for item in merged_region_list:
        item[2] = round( numpy.mean(item[2]), 3)
    return merged_region_list

# merge_region_list( [ [1, 10, 0.1], [14, 20, 0.3] ], [ [7, 12, 0.3], [18, 25, 0.1] ] )

def combineTwoPBREnsembl(pbr_ensembl1, pbr_ensembl2, merge=True):
    combine_ensembl = {}
    for protein in pbr_ensembl1:
        if protein in pbr_ensembl2:
            print 'Processing ', protein, '...'
            (only_1, common, only_2) = vennTwoSet(pbr_ensembl1[protein].keys(), pbr_ensembl2[protein].keys())
            print 'Point 1'
            combine_ensembl[protein] = {}
            for trans_id in only_1:
                combine_ensembl[protein][trans_id] = pbr_ensembl1[protein][trans_id]
            for trans_id in only_2:
                combine_ensembl[protein][trans_id] = pbr_ensembl2[protein][trans_id]
            for trans_id in common:
                if merge:
                    combine_ensembl[protein][trans_id] = merge_region_list(pbr_ensembl1[protein][trans_id], pbr_ensembl2[protein][trans_id])
                else:
                    combine_ensembl[protein][trans_id] = pbr_ensembl1[protein][trans_id] + pbr_ensembl2[protein][trans_id]
        else:
            combine_ensembl[protein] = pbr_ensembl1[protein]
    for protein in pbr_ensembl2:
        if protein not in pbr_ensembl1:
            combine_ensembl[protein] = pbr_ensembl2[protein]
    return combine_ensembl

# pbr_ensembl = combineTwoPBREnsembl(Piranha_pbr_ensembl, PARalyzer_pbr_ensembl)


def permutate_choice(chooseSet, chooseNum, chooseTime=1000, method=numpy.mean):
    "从数据集chooseSet中每次抽取chooseNum个用函数method得到一个统计量，重复chooseTime次"
    sample_items = []
    for i in range(chooseTime):
        sample_items.append( method(random.sample(chooseSet, chooseNum)) )
    return sample_items


def fetch_random_icshape(single_icshape, length, capacity):
    "从icshape字典中随机地抽取capacity条长度为length的"
    "序列片段，当做List返回"
    valid_icshape_List = []
    for List in single_icshape.values():
        valid_icshape_List += [float(i) for i in List if i != 'NULL']
    print 'Total icSHAPE Values:', len(valid_icshape_List)
    sampled_site = random.sample(range(len(valid_icshape_List)-length), capacity)
    sampled_icshape_List = [ valid_icshape_List[start:start+length] for start in sampled_site ]
    return sampled_icshape_List

def zipEnsembl(ensembl):
    "把定义的ensembl变成列表"
    parsedEnsembl = []
    for trans_id in ensembl:
        for site in ensembl[trans_id]:
            parsedEnsembl.append( [trans_id, site] )
    return parsedEnsembl

def read_refSeqID_2_ensemblID(file_name):
    "读取refSeq和Ensembl ID的对应文件"
    RefSeq2Ensembl = {}
    IN = open(file_name)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        ensembl_ID = arr[0].split('.')[0]
        refSeq_ID = arr[1].split('.')[0]
        RefSeq2Ensembl[refSeq_ID] = ensembl_ID
        line = IN.readline()
    return RefSeq2Ensembl


def permutationList(DataList, sample_size=50, sample_times=200):
    if len(DataList) < sample_size*2:
        print 'Warning: DataList size is too small!!'
        sample_size = len(DataList) / 2
    permutated_values = []
    for i in range(sample_times):
        sampled_values = random.sample(DataList, sample_size)
        permutated_values.append( numpy.mean(sampled_values) )
    return permutated_values


def FDR(x):
    """
    Assumes a list or numpy array x which contains p-values for multiple tests
    Copied from p.adjust function from R  
    """
    o = [i[0] for i in sorted(enumerate(x), key=lambda v:v[1],reverse=True)]
    ro = [i[0] for i in sorted(enumerate(o), key=lambda v:v[1])]
    q = sum([1.0/i for i in xrange(1,len(x)+1)])
    l = [q*len(x)/i*x[j] for i,j in zip(reversed(xrange(1,len(x)+1)),o)]
    l = [l[k] if l[k] < 1.0 else 1.0 for k in ro]
    return l

def calc_benjamini_hochberg_corrections(p_values, num_total_tests):
    """
    Calculates the Benjamini-Hochberg correction for multiple hypothesis
    testing from a list of p-values *sorted in ascending order*.

    See
    http://en.wikipedia.org/wiki/False_discovery_rate#Independent_tests
    for more detail on the theory behind the correction.

    **NOTE:** This is a generator, not a function. It will yield values
    until all calculations have completed.

    :Parameters:
    - `p_values`: a list or iterable of p-values sorted in ascending
      order
    - `num_total_tests`: the total number of tests (p-values)

    """
    prev_bh_value = 0
    for i, p_value in enumerate(p_values):
        bh_value = p_value * num_total_tests / (i + 1)
        # Sometimes this correction can give values greater than 1,
        # so we set those values at 1
        bh_value = min(bh_value, 1)

        # To preserve monotonicity in the values, we take the
        # maximum of the previous value or this one, so that we
        # don't yield a value less than the previous.
        bh_value = max(bh_value, prev_bh_value)
        prev_bh_value = bh_value
        yield bh_value



def predictStructure_RNAfold(Seq, Shape=""):
    """ Predict RNA Structure using RNAfold
    Seq: sequence
    Shape: shape
    The length of each seq must be equal to shape
    """
    if Shape and len(Shape) != len(Seq):
        print >>sys.stderr, "Shape and Seq must have same length"
        raise Exception("Different Length")
    randID = random.randint(10000,99999)
    tmp_fa_file = "/tmp/tmp_%s.fa" % (randID, )
    tmp_shape_file = "/tmp/tmp_%s.shape" % (randID, )
    # prepare tmp.fa
    FA = open(tmp_fa_file, 'w')
    print >>FA, Seq
    FA.close()
    if Shape:
        # prepare tmp.shape
        SHAPE = open(tmp_shape_file, 'w')
        for idx in range(len(Seq)):
            if Shape[idx] != "NULL":
                print >>SHAPE, str(idx+1)+"\t"+str(float(Shape[idx])*2)
                #print str(idx+1)+"\t"+str(float(Shape[idx])*2)
            else:
                print >>SHAPE, str(idx+1)+"\t-1"
                #print str(idx+1)+"\t-1"
        SHAPE.close()
        # Predict Secondary Structure
        RNAFold_CMD = 'RNAfold --noPS --shape=%s < %s' % (tmp_shape_file, tmp_fa_file)
        return_code, predicted_structure = commands.getstatusoutput(RNAFold_CMD)
        os.remove(tmp_shape_file)
    else:
        RNAFold_CMD = 'RNAfold --noPS < %s' % (tmp_fa_file, )
        return_code, predicted_structure = commands.getstatusoutput(RNAFold_CMD)
    os.remove(tmp_fa_file)
    predicted_structure = predicted_structure.strip().split()[1]
    return predicted_structure



def plotRNAStructure(Seq, Shape, Structure, title, fileName):
    VARNA_CMD = """java -cp /Share/home/zhangqf8/usr/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd \
        -sequenceDBN "%s" \
        -structureDBN "%s" \
        -colorMap "%s" \
        -resolution "10" \
        -title "%s" \
        -o %s \
        -algorithm naview \
        -titleSize 10
    """
    color_map = ";".join([shape if shape != "NULL" else "0" for shape in Shape])
    print >>sys.stdout, VARNA_CMD % (Seq, Structure, color_map, title, fileName)
    os.system(VARNA_CMD % (Seq, Structure, color_map, title, fileName))

def bi_fold(seq_1, seq_2, permit_local_pairing=False):
    OUT = open("/tmp/seq_1.fa", "w")
    print >>OUT, ">seq_1\n%s" % (seq_1, )
    OUT.close()
    OUT = open("/tmp/seq_2.fa", "w")
    print >>OUT, ">seq_2\n%s" % (seq_2, )
    OUT.close()
    if not permit_local_pairing:
        CMD = "bifold -i /tmp/seq_1.fa /tmp/seq_2.fa /tmp/bifold.ct"
    else:
        CMD = "bifold /tmp/seq_1.fa /tmp/seq_2.fa /tmp/bifold.ct"
    #print CMD
    os.system(CMD)
    return_code, return_string = commands.getstatusoutput( "grep ENERGY /tmp/bifold.ct | wc -l" )
    structure_number = int( return_string.strip() )
    structure_list = []
    ct2dot_cmd = "ct2dot /tmp/bifold.ct %d /dev/stdout"
    for idx in range(structure_number):
        return_code, return_string = commands.getstatusoutput( ct2dot_cmd % (idx+1, ) )
        lines = return_string.split('\n')
        energy = float(lines[0].strip().split()[-2])
        structure = return_string.split()[8]
        structure_list.append( (energy, lines[2]) )
    cur_seq = seq_1 + "III" + seq_2
    return cur_seq, structure_list


def prepare_VARNA_colormap_file(shape_list, file_name, base=0):
    OUT = open(file_name, "w")
    for init_idx in range(base):
        print >>OUT, "%d\t%s" % (init_idx+1, "0.0")
    for idx,shape in enumerate(shape_list):
        if shape != "NULL":
            print >>OUT, "%d\t%s" % (idx+1+base, shape)
        else:
            print >>OUT, "%d\t%s" % (idx+1+base, 0.0)
    OUT.close()


def prepare_structure_prediction(shape_list, file_name):
    OUT = open(file_name, "w")
    for idx,shape in enumerate(shape_list):
        if shape != "NULL":
            print >>OUT, "%d\t%s" % (idx+1, shape)
        else:
            print >>OUT, "%d\t%s" % (idx+1, -999)
    OUT.close()

