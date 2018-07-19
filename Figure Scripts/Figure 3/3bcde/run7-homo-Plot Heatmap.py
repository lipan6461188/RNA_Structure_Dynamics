
##### homo


HEK293_refbed = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Homo_sapiens.trans.bed'
mES_refbed = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Mus_musculus.trans.bed'

vivo = {'ch':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/homo/R_output/hg38_mm10_ch.result',
'np':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/homo/R_output/hg38_mm10_np.result',
'cy':'/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/homo/R_output/hg38_mm10_cy.result'}

hek_geneTypeDict = loadGeneType(HEK293_refbed)
mes_geneTypeDict = loadGeneType(mES_refbed)




def no_overlap_statistics(fileDict):
    typeSiteDict = {}
    typeTransNumDict = {}
    typeDynRatio = {}
    for Component in fileDict:
        typeSiteDict[Component], typeTransNumDict[Component], typeDynRatio[Component] = statisticDataSize(fileDict[Component], hek_geneTypeDict, mes_geneTypeDict)
    return typeSiteDict, typeTransNumDict, typeDynRatio


def statisticDataSize(file_name, hek_geneTypeDict, mes_geneTypeDict, transList=None):
    "统计每种类型转录本的数据量"
    def transType(trans_id, geneTypeDict):
        if trans_id in geneTypeDict:
            return geneTypeDict[trans_id]
        else:
            return 'unKnown'
    table = pd.read_csv(file_name,sep="\t")
    if transList:
        sorted_transList = sorted(transList)
        in_List = (table['hekID'] + "_" + table['hek_1_based'].apply(lambda x:str(x))).apply(lambda x: biSearch(x, sorted_transList))
        table = table.loc[in_List,:]
    " 过滤掉两个物种中RNA类型的比较 "
    hek_type = table['hekID'].apply(transType, args=(hek_geneTypeDict,))
    mes_type = table['mesID'].apply(transType, args=(mes_geneTypeDict,))
    table = table.loc[hek_type == mes_type,:]; table.index = range(len(table))
    typeSiteDict = dict( table['hekID'].apply(transType, args=(hek_geneTypeDict,)).value_counts() )
    typeTransNumDict = dict( pd.Series(table['hekID'].unique()).apply(transType, args=(hek_geneTypeDict,)).value_counts() )
    typeDynSiteDict = dict( table[table['label']=="*"]['hekID'].apply(transType, args=(hek_geneTypeDict,)).value_counts() )
    typeDynRatio = {}
    for mtype in typeSiteDict:
        if mtype == 'unKnown':
            continue
        typeDynRatio[mtype] = 1.0*typeDynSiteDict[mtype]/typeSiteDict[mtype]
    return typeSiteDict, typeTransNumDict, typeDynRatio

typeSiteDict, typeTransNumDict, typeDynRatio = no_overlap_statistics(vivo)
typeSite_df = pd.DataFrame(typeSiteDict).loc[['snoRNA', 'mRNA','lncRNA', 'snRNA', 'miRNA'], :]
typeTransNum_df = pd.DataFrame(typeTransNumDict).loc[['snoRNA', 'mRNA','lncRNA', 'snRNA', 'miRNA'], :]
typeDynRatio_df = pd.DataFrame(typeDynRatio).loc[['snoRNA', 'mRNA','lncRNA', 'snRNA', 'miRNA'], :]

anno = []
for row in typeDynRatio_df.values:
    anno.append( [ "%.2f%%" % (i*100,) for i in row ] )

anno = pd.DataFrame(anno, columns=typeDynRatio_df.columns)
anno.index = typeDynRatio_df.index

row_index = ['mRNA', 'snRNA', 'lncRNA', 'miRNA', 'snoRNA'] #typeDynRatio_df.mean(axis=1).sort_values(ascending=False).index
col_index = ['ch', 'np', 'cy']#typeDynRatio_df.mean(axis=0).sort_values(ascending=False).index


#######################################
#######################################
################ 画图 #################
#######################################
#######################################

plt.close()
fig = plt.figure(figsize=(15,10))

plt.subplot(2,2,1)
ax = sns.heatmap(data=typeDynRatio_df.loc[row_index, col_index], annot=anno.loc[row_index, col_index], fmt='5s', cmap="YlGnBu", annot_kws={'size':'large', 'weight':'bold'})
plt.setp(ax.yaxis.get_majorticklabels(), rotation=0)
plt.title("Dynamic Ratio")

plt.subplot(2,2,2)
ax = sns.heatmap(data=typeDynRatio_df.loc[row_index, col_index], annot=anno.loc[row_index, col_index], fmt='5s', cmap="YlGnBu", vmin=0, vmax=0.35, center=0.12, annot_kws={'size':'large', 'weight':'bold'})
plt.setp(ax.yaxis.get_majorticklabels(), rotation=0)
plt.title("Dynamic Ratio")

plt.subplot(2,2,3)
ax = sns.heatmap(data=typeSite_df.loc[row_index, col_index], annot=True, fmt='g', cmap="YlGnBu", annot_kws={'size':'large', 'weight':'bold'})
plt.setp(ax.yaxis.get_majorticklabels(), rotation=0)
plt.title("Data Size(Sites)")

plt.subplot(2,2,4)
ax = sns.heatmap(data=typeTransNum_df.loc[row_index, col_index], annot=True, fmt='g', cmap="YlGnBu", annot_kws={'size':'large', 'weight':'bold'})
plt.setp(ax.yaxis.get_majorticklabels(), rotation=0)
plt.title("Data Size(Transcripts)")

plt.savefig(os.environ.get( "HOME" ) +"/figs/hek293_mes_homology.pdf")

