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

def readSeq(seq_name, removeVersion=False):
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

def writeSeq(seqDict, outFile):
    from seq import cutSeq
    OUT = open(outFile, 'w')
    for trans_id in seqDict:
        print >>OUT, '>%s\n%s' % (trans_id, cutSeq(seqDict[trans_id]))
    OUT.close()


def readDot(inFile):
    structure = {}
    IN = open(inFile)
    line = IN.readline()
    cur_trans_id = ""
    while line:
        if line[0] == '>':
            cur_trans_id = line[1:].strip().split()[0]
            structure[cur_trans_id] = [ IN.readline().strip(), IN.readline().strip() ]
        line = IN.readline()
    IN.close()
    return structure

def loadicSHAPE(file_name, removeVersion=False):
    icSHAPE = {}
    IN = open(file_name)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        trans_id = arr[0]
        if removeVersion: trans_id = arr[0].split(".")[0]
        shape = arr[3:]
        icSHAPE[ trans_id ] = shape
        line = IN.readline()
    return icSHAPE

def init_rect(rowNum, colNum, rowNames=[], colNames=[], init_value=None):
    import pandas as pd
    import numpy as np
    
    if colNames:
        assert(len(colNames)==colNum)
    else:
        colNames = np.arange(colNum)
    
    if rowNames:
        assert(len(rowNames)==rowNum)
    else:
        rowNames = np.arange(rowNum)
    
    df = pd.DataFrame(np.zeros((rowNum, colNum)), index=rowNames, columns=colNames)
    if init_value == None:
        return df
    else:
        df.iloc[:,:] = init_value
        return df


def loadCodon(codon_file):
    """
        Each line: Ser=UCU|UCC|UCA|UCG|AGU|AGC
    """
    codon = {}
    IN = open(codon_file)
    line = IN.readline()
    while line:
        data = line.strip().split('=')
        for cur_codon in data[1].split('|'):
            codon[cur_codon] = data[0]
            codon[cur_codon.replace('U','T')] = data[0]
        line = IN.readline()
    IN.close()
    return codon

def sample_file(inFile, outFile, ratio, skip=''):
    import random
    
    assert( 0 <= ratio <= 1 )
    
    IN = open(inFile)
    OUT = open(outFile, 'w')
    line = IN.readline()
    while line:
        if skip and line.startswith(skip):
            pass
        else:
            if random.random() < ratio:
                OUT.writelines(line)
        line = IN.readline()
    IN.close()
    OUT.close()

def permutate_pValue(shuffle_list, true_value, mode='low'):
    assert( mode in ('low', 'high') )
    if mode == 'high':
        shuffle_list = sorted(shuffle_list, reverse=True)
    elif mode == 'low':
        shuffle_list = sorted(shuffle_list, reverse=False)
    sites = len(shuffle_list)
    for idx in range(len(shuffle_list)):
        if mode == 'high':
            if true_value > shuffle_list[idx]:
                sites = idx
                break
        if mode == 'low':
            if true_value < shuffle_list[idx]:
                sites = idx
                break
    return 1.0*sites/len(shuffle_list)


def find_all_match(pattern, string):
    import re
    matches = []
    for item in re.finditer(pattern, string):
        s, e = item.start(), item.end()
        matches.append( ((s+1, e), string[s:e]) )
    return matches

def uniq_list(raw_list):
    uniqued_list = []
    for item in raw_list:
        if item not in uniqued_list:
            uniqued_list.append(item)
    return uniqued_list

def biSearch(item, Set):
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

def cluster_dataFrame(raw_dataFrame, row_cluster=True, col_cluster=True):
    import matplotlib.pyplot as plt
    
    data2d = sns.clustermap(raw_dataFrame, row_cluster=row_cluster, col_cluster=col_cluster).data2d
    plt.close()
    
    return data2d

def shuffle_plot(shuffle_list, true_value, ylim=()):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    
    sns.set_style('white')
    ax = sns.stripplot(y=shuffle_list, jitter=True, alpha=0.15)
    sns.boxplot(y=shuffle_list)
    plt.plot(0, true_value, '^', markersize=15, color="#e72d27")
    
    if ylim:
        plt.ylim(ylim[0], ylim[1])

def calcGINI(my_list, valid_cutoff=10):
    def GINI(list_of_values):
        length = len(list_of_values)
        total = sum(list_of_values)
        if total == 0: 
            return None
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
    
    if len(floatArr) > valid_cutoff:
        return GINI(floatArr)
    return None

def scatter_plot(dict_x, dict_y, color="", xlabel="", ylabel="", title="", filt_list=None, xlim=None, ylim=None, Parser=None, xlog=False, ylog=False):
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
    ax = sns.jointplot(data=df, x='x', y='y', kind='reg', color=color)
    
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    if title:
        plt.title(title)
    
    return ax

def readDot(dot_name):
    transDot = {}
    IN = open(dot_name)
    line = IN.readline()
    while line:
        if line[0] == '>':
            tid = line.strip()[1:]
            seq = IN.readline().strip()
            ss = IN.readline().strip()
            transDot[tid] = (seq, ss)
        elif line.strip() == '':
            pass
        else:
            print "Error: \n\tundifined line: %s" % (line, )
            return
        line = IN.readline()
    IN.close()
    return transDot



# <======== Calculate Structure-Shape consistency =========>

def Shape_positive_rate(ss_code, shape_list, cutoff):
    Pos_Num = 0
    True_Pos = 0
    False_Pos = 0
    Neg_Num = 0
    for idx, code in enumerate(list(ss_code)):
        if shape_list[idx] != 'NULL':
            if code != ".":
                Pos_Num += 1
                if float(shape_list[idx]) <= cutoff:
                    True_Pos += 1
                else:
                    pass
            else:
                Neg_Num += 1
                if float(shape_list[idx]) <= cutoff:
                    False_Pos += 1
                else:
                    pass
    return 1.0*True_Pos/Pos_Num, 1.0*False_Pos/Neg_Num

def calc_shape_ROC(ss_code, shape_list, step=0.01):
    assert(len(ss_code)==len(shape_list))
    ROC = []
    cutoff = -step
    while cutoff < 1.0 + step:
        TPR, FPR = Shape_positive_rate(ss_code, shape_list, cutoff)
        ROC.append( (FPR, TPR) )
        cutoff += step
    return ROC

def calc_AUC(ROC):
    area = 0.5*ROC[0][0]*ROC[0][1]
    for idx in range(1, len(ROC)-1):
        area += (ROC[idx+1][0]-ROC[idx][0]) * ROC[idx+1][1] + 0.5 * (ROC[idx+1][0]-ROC[idx][0]) * (ROC[idx+1][1]-ROC[idx][1])
    return area

def plot_ROC(ROC):
    x = [ i[0] for i in ROC ]
    y = [ i[1] for i in ROC ]
    plt.plot(x, y, '-')
    plt.xlim(0,1)
    plt.ylim(0,1)


#<=====================================>


