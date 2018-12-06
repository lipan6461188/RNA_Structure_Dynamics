import re, sys, os, getopt, time
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy


def gene_type(raw_type):
    "Convert Raw Gene Type to Our Defined Gene Type"
    valid_gene_type = ('pseudogene', 'snoRNA', 'snRNA', 'miRNA', 'misc_RNA', 'rRNA')
    #lncRNA_class = ('antisense','lincRNA','processed_transcript','sense_intronic','TEC','sense_overlapping')
    lncRNA_class = ('3prime_overlapping_ncrna','antisense','lincRNA','non_coding','sense_intronic','sense_overlapping','processed_transcript')
    if raw_type in valid_gene_type: return raw_type;
    if re.match('.*pseudogene',raw_type): return 'pseudogene';
    if raw_type == 'protein_coding': return 'mRNA';
    if raw_type in lncRNA_class: return 'lncRNA';
    return 'other'

def loadGeneType(file_name):
    "加载基因的类型"
    typeDict = {}
    IN = open(file_name)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        typeDict[ arr[5] ] = gene_type(arr[6])
        line = IN.readline()
    return typeDict

def statisticDataSize(file_name, geneTypeDict, transList=None):
    "统计每种类型转录本的数据量"
    def transType(trans_id):
        if trans_id in geneTypeDict:
            return geneTypeDict[trans_id]
        else:
            return 'unKnown'
    table = pd.read_csv(file_name,sep="\t")
    if transList:
        sorted_transList = sorted(transList)
        in_List = (table['trans_id'] + "_" + table['site_1_based'].apply(lambda x:str(x))).apply(lambda x: biSearch(x, sorted_transList))
        table = table.loc[in_List,:]
    table['type'] = table['trans_id'].apply(transType)
    typeSiteDict = dict( table['trans_id'].apply(transType).value_counts() )
    typeTransNumDict = dict( pd.Series(table['trans_id'].unique()).apply(transType).value_counts() )
    typeDynSiteDict = dict( table[table['label']=="*"]['trans_id'].apply(transType).value_counts() )
    #print typeDynSiteDict
    typeDynRatio = {}
    for mtype in typeSiteDict:
        if mtype == 'unKnown':
            continue
        try:
            typeDynRatio[mtype] = 1.0*typeDynSiteDict[mtype]/typeSiteDict[mtype]
        except KeyError:
            typeDynRatio[mtype] = 0.0
    return typeSiteDict, typeTransNumDict, typeDynRatio


def transList(file_name):
    "获得文件中的转录本列表"
    table = pd.read_csv(file_name,sep="\t")
    return list(table['trans_id'] + "_" + table['site_1_based'].apply(lambda x:str(x)))

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

def intersect(List1, List2):
    "取交集"
    sorted_List2 = sorted(List2)
    common_Set = []
    for item in List1:
        if biSearch(item, sorted_List2):
            common_Set.append(item)
    return common_Set


########## 统计数量


def no_overlap_statistics(fileDict):
    typeSiteDict = {}
    typeTransNumDict = {}
    typeDynRatio = {}
    for Component in fileDict:
        typeSiteDict[Component], typeTransNumDict[Component], typeDynRatio[Component] = statisticDataSize(fileDict[Component], transTypeDict)
    return typeSiteDict, typeTransNumDict, typeDynRatio

def overlap_statistics(fileDict):
    typeSiteDict = {}
    typeTransNumDict = {}
    typeDynRatio = {}
    common_transList = []
    for Component in fileDict:
        print 'intersect '+Component+'...'
        if Component == fileDict.keys()[0]:
            common_transList = transList(fileDict[Component])
        else:
            common_transList = intersect(transList(fileDict[Component]), common_transList)
    for Component in fileDict:
        typeSiteDict[Component], typeTransNumDict[Component], typeDynRatio[Component] = statisticDataSize(fileDict[Component], transTypeDict, transList=common_transList)
    return typeSiteDict, typeTransNumDict, typeDynRatio

def plot(dataSet, title, vv=False):
    typeSiteDict, typeTransNumDict, typeDynRatio = no_overlap_statistics(dataSet)
    ""
    RNA_type_List = ['mRNA', 'snRNA', 'lncRNA', 'pseudogene', 'miRNA', 'snoRNA']
    #mes_RNA_type_List = ['snoRNA', 'mRNA','pseudogene','lncRNA', 'snRNA', 'miRNA']
    typeSite_df = pd.DataFrame(typeSiteDict).loc[RNA_type_List, :]
    typeTransNum_df = pd.DataFrame(typeTransNumDict).loc[RNA_type_List, :]
    typeDynRatio_df = pd.DataFrame(typeDynRatio).loc[RNA_type_List, :]
    print typeDynRatio_df
    ""
    anno = []
    for row in typeDynRatio_df.values:
        anno.append( [ "%.2f%%" % (i*100,) for i in row ] )
    ""
    anno = pd.DataFrame(anno, columns=typeDynRatio_df.columns)
    anno.index = typeDynRatio_df.index
    ""
    row_index = RNA_type_List #typeDynRatio_df.sum(axis=1).sort_values(ascending=False).index
    col_index = ['chcy', 'chnp', 'npcy'] #typeDynRatio_df.sum(axis=0).sort_values(ascending=False).index
    #row_index = typeDynRatio_df.sum(axis=1).sort_values(ascending=False).index
    #col_index = typeDynRatio_df.sum(axis=0).sort_values(ascending=False).index
    if vv:
        col_index = ['ch', 'np', 'cy']
    "画图"
    fig = plt.figure(figsize=(15,10))
    ""
    log_typeDynRatio_df = numpy.log10(typeDynRatio_df*100)
    plt.subplot(2,2,1)
    ax = sns.heatmap(data=typeDynRatio_df.loc[row_index, col_index], annot=anno.loc[row_index, col_index], fmt='5s', cmap="YlGnBu", annot_kws={'size':'large', 'weight':'bold'})
    plt.setp(ax.yaxis.get_majorticklabels(), rotation=0)
    plt.title("Dynamic Ratio")
    ""
    plt.subplot(2,2,2)
    ax = sns.heatmap(data=typeDynRatio_df.loc[row_index, col_index], annot=anno.loc[row_index, col_index], fmt='5s', cmap="YlGnBu", vmin=0.00, vmax=0.35, center=0.12, annot_kws={'size':'large', 'weight':'bold'})
    #ax = sns.heatmap(data=log_typeDynRatio_df.loc[row_index, col_index], annot=anno.loc[row_index, col_index], fmt='5s', cmap="YlGnBu", vmin=0.00, vmax=1.30, center=0.8, annot_kws={'size':'large', 'weight':'bold'})
    plt.setp(ax.yaxis.get_majorticklabels(), rotation=0)
    plt.title("Dynamic Ratio")
    ""
    plt.subplot(2,2,3)
    ax = sns.heatmap(data=typeSite_df.loc[row_index, col_index], annot=True, fmt='g', cmap="YlGnBu", annot_kws={'size':'large', 'weight':'bold'})
    plt.setp(ax.yaxis.get_majorticklabels(), rotation=0)
    plt.title("Data Size(Sites)")
    ""
    plt.subplot(2,2,4)
    ax = sns.heatmap(data=typeTransNum_df.loc[row_index, col_index], annot=True, fmt='g', cmap="YlGnBu", annot_kws={'size':'large', 'weight':'bold'})
    plt.setp(ax.yaxis.get_majorticklabels(), rotation=0)
    plt.title("Data Size(Transcripts)")
    ""
    plt.savefig(os.environ.get( "HOME" ) +"/figs/%s.pdf" % (title, ))
    plt.close()


##### HEK293

HEK293_refbed = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Homo_sapiens.trans.bed'
HEK293_vivo = {'chnp':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/HEK293/R_output/ch-np.result',
'npcy':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/HEK293/R_output/np-cy.result',
'chcy':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/HEK293/R_output/ch-cy.result'}

HEK293_vitro = {'chnp':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/HEK293/R_output/chT-npT.result',
'npcy':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/HEK293/R_output/npT-cyT.result',
'chcy':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/HEK293/R_output/chT-cyT.result'}

HEK293_vivo_vitro = {'ch':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/HEK293/R_output/ch-chT.result',
'np':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/HEK293/R_output/np-npT.result',
'cy':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/HEK293/R_output/cy-cyT.result'}

transTypeDict = loadGeneType(HEK293_refbed)

dataSet = HEK293_vivo
title = 'HEK293_vivo'
plot(dataSet, title)

dataSet = HEK293_vitro
title = 'HEK293_vitro'
plot(dataSet, title)

dataSet = HEK293_vivo_vitro
title = 'HEK293_vivo_vitro'
plot(dataSet, title, True)

###### mES

mES_refbed = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Mus_musculus.trans.bed'
mES_vivo = {'chnp':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/mES/R_output/ch-np.result',
'npcy':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/mES/R_output/np-cy.result',
'chcy':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/mES/R_output/ch-cy.result'}

mES_vitro = {'chnp':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/mES/R_output/chT-npT.result',
'npcy':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/mES/R_output/npT-cyT.result',
'chcy':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/mES/R_output/chT-cyT.result'}

mES_vivo_vitro = {'ch':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/mES/R_output/ch-chT.result',
'np':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/mES/R_output/np-npT.result',
'cy':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/mES/R_output/cy-cyT.result'}

transTypeDict = loadGeneType(mES_refbed)

dataSet = mES_vivo
title = 'mES_vivo'
plot(dataSet, title)

dataSet = mES_vitro
title = 'mES_vitro'
plot(dataSet, title)

dataSet = mES_vivo_vitro
title = 'mES_vivo_vitro'
plot(dataSet, title, True)

