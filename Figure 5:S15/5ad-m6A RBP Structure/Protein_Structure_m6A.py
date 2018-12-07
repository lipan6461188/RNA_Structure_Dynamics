
from icSHAPE import *

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

def random_color():
    return random.choice(['#4c72b0', '#55a868', '#c44e52', '#8172b2', '#ccb974', '#64b5cd'])

human_seq = readSeq("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/transcriptome.fa")
m6A_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/m6A/choose_good_protein/m6A_dataset_1.merged.transCoor.bed'

icSHAPE_ROOT_DIR = "/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/%s/shape.out"
np_icshape_file = icSHAPE_ROOT_DIR % ("hek_np_vivo", )
np_icshape = loadicSHAPE(np_icshape_file, human_seq)

m6A_ensembl = read_m6A(m6A_file, human_seq)
new_pbr_columns = {'pName': 5, 'score':6, 'transID': 10, 'transStart': 11, 'transEnd': 12}

# 1. Piranha

clip_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/m6A/choose_good_protein/test/allProteins/Piranha_0.01/HEK293.Piranha.merge.transCoor.txt'
Piranha_pbr_ensembl = read_pbr(clip_file, human_seq, new_pbr_columns)
normalize_pbr_conf(Piranha_pbr_ensembl, mode = 'scalling')

# 2. PARalyzer

clip_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/m6A/choose_good_protein/test/allProteins/PARalyzer/HEK293.PARalyzer.merge.transCoor.txt'
PARalyzer_pbr_ensembl = read_pbr(clip_file, human_seq, new_pbr_columns)
normalize_pbr_conf(PARalyzer_pbr_ensembl, mode = 'scalling')

# 3. PARalyzer + Piranha YTHD + Piranha HNRNPC

add_protein = {'HNRNPC': Piranha_pbr_ensembl['HNRNPC']}
for Piranha_protein in Piranha_pbr_ensembl:
    if 'YTHD' in Piranha_protein and Piranha_protein not in PARalyzer_pbr_ensembl:
        add_protein[Piranha_protein] = Piranha_pbr_ensembl[Piranha_protein]
    if 'METTL' in Piranha_protein and Piranha_protein not in PARalyzer_pbr_ensembl:
        add_protein[Piranha_protein] = Piranha_pbr_ensembl[Piranha_protein]

add_protein_Sub = sub_pbr_ensembl(add_protein, sub_human_seq)
PARalyzerSub_pbr_ensembl = sub_pbr_ensembl(PARalyzer_pbr_ensembl, sub_human_seq)

pbr_ensembl = combineTwoPBREnsembl(add_protein_Sub, PARalyzerSub_pbr_ensembl)
dataSizeName = 'PARalyzer_YTHD'


#####################
###    Figure 5a
#####################

motif_index = buildMotifIndex(human_seq, motif=['GGACT'], Plus_Site=3)
motif_m6A_ensembl = read_m6A_withMotif(m6A_file, human_seq, motif=['GGACT'], A_Site=3)

sigSet_step1 = []

#@@@@@@@@@@@
EXTEND = 0
#@@@@@@@@@@@

DataSize = ''
dataSet = []
for protein in pbr_ensembl:
    # + 获得没有m6A修饰的转录本
    hnRNPC_trans_without_m6A = get_OnlyPBR_trans_single(pbr_ensembl[protein], m6A_ensembl.keys())
    # + 获得没有m6A修饰的转录本中的motif
    hnRNPC_rand_m6A_ensembls_A = get_allMotif_fromTransSet(hnRNPC_trans_without_m6A, sub_human_seq, motifIndex=motif_index, sample_num=False)
    # + 与没有修饰的转录本中的motif有交集的蛋白结合区域
    hnRNPC_rand_pbr_ensembl_cover_m6A = get_pbr_cover_m6A_sep_m6A_single(pbr_ensembl[protein], hnRNPC_rand_m6A_ensembls_A, extend=EXTEND)
    hnRNPC_m6A_pbr_ensembl_cover_m6A = get_pbr_cover_m6A_sep_m6A_single(pbr_ensembl[protein], motif_m6A_ensembl, extend=EXTEND)
    rand_m6A_conf = fetch_conf_from_PBR_Ensembl_single(hnRNPC_rand_pbr_ensembl_cover_m6A, func=float)
    m6A_conf = fetch_conf_from_PBR_Ensembl_single(hnRNPC_m6A_pbr_ensembl_cover_m6A, func=float)
    
    DataSize += "%s\t%d\t%d\n" % (protein, len(rand_m6A_conf), len(m6A_conf) )
    if len(rand_m6A_conf) > 30 and len(m6A_conf) > 30:
        p_value = scipy.stats.mannwhitneyu(rand_m6A_conf, m6A_conf)[1]
        diff = numpy.mean(m6A_conf) - numpy.mean(rand_m6A_conf)
        dataSet.append([p_value, diff, protein])


print DataSize
print >>open("datasize.txt", 'w'), DataSize
dataSet.sort(key=lambda x:x[1])
dataSet_xx = pd.DataFrame(dataSet, columns=('p_value', 'diff', 'protein'))
dataSet_xx['corrected_p_value'] = statsmodels.sandbox.stats.multicomp.multipletests(method='fdr_bh', pvals=dataSet_xx['p_value'])[1]

dataSet_xx.to_csv("/Share/home/zhangqf8/lipan/tmp/Step_1.txt", sep="\t")

dataSet_xx['log_p_value'] = -numpy.log10(dataSet_xx['p_value']) # 直接对P-value求log
dataSet_xx['color'] = '#9f9fa0'
dataSet_xx.loc[dataSet_xx.log_p_value >= -numpy.log10(0.05),'color'] = '#4C72B0'

sigSet_step1 = list(dataSet_xx.loc[dataSet_xx.corrected_p_value < 0.05,'protein'])
print sigSet_step1


# 画成散点图
sns.set_style("ticks")
sns.despine(offset=10, trim=True);
ax=sns.lmplot('log_p_value', 'diff',
           data=dataSet_xx,
           fit_reg=False,
           hue="protein",
           palette=sns.color_palette(dataSet_xx.color),
           legend=False,
           scatter_kws={"marker": "D", "s": 70}
           )
# 蛋白名字
for point in list(dataSet_xx.values):
    x = point[3]
    y = point[1]
    text = point[2]
    if x < numpy.log10(0.05):
        plt.text(x, y, text, horizontalalignment='right',verticalalignment='top', fontsize=4, color=random_color())
    else:
        plt.text(x, y, text, horizontalalignment='left',verticalalignment='top', fontsize=3, color=random_color())

plt.axvline(x=numpy.log10(0.05), linestyle='dashed', color="#C44E52")
plt.axhline(y=0.0, linestyle='dashed', color="black")
plt.title('Step 1 - Protein Binding Preference of m6A\n'+dataSizeName)
plt.xlabel('lg(p-value)')
plt.ylabel('Diff of Normalized Protein Bind Strength')
plt.xlim(0, 15) # plt.xlim(0, 60)
plt.ylim(-0.35, 0.45)
plt.tight_layout()
plt.savefig(HOME+"/figs/Step1-%s.pdf" % (dataSizeName, ))
plt.show()




#####################
###    Figure 5d
#####################


sub_human_seq = {}
for trans_id in human_seq:
    if trans_id in np_icshape:
        sub_human_seq[trans_id] = human_seq[trans_id]

m6A_ensembl = read_m6A(m6A_file, sub_human_seq)
motif_index = buildMotifIndex(sub_human_seq, motif=['GGACT'], Plus_Site=3)

DataSize = ''
dataSet = []

plt.close()
index = 1
for my_protein in ['HNRNPC', 'IGF2BP3']:
    hnRNPC_trans_without_m6A = get_OnlyPBR_trans_single(pbr_ensembl[my_protein], m6A_ensembl.keys())
    hnRNPC_rand_m6A_ensembls_A = get_allMotif_fromTransSet(hnRNPC_trans_without_m6A, sub_human_seq, motifIndex=motif_index, sample_num=False)
    hnRNPC_rand_pbr_ensembl_cover_m6A = get_pbr_cover_m6A_sep_m6A_single(pbr_ensembl[my_protein], hnRNPC_rand_m6A_ensembls_A, extend=0)
    hnRNPC_rand_m6A_overConfPBR = getout_m6A_overlap_with_pbr(hnRNPC_rand_pbr_ensembl_cover_m6A, hnRNPC_rand_m6A_ensembls_A, conf_cutoff=0, extend=0)
    hnRNPC_rand_xx, hnRNPC_rand_conf = get_icshape_df_m6A_Site(hnRNPC_rand_m6A_overConfPBR, np_icshape)
    DataSize += '%s\t%d\n' % (my_protein, len(hnRNPC_rand_conf))
    
    #index_mean = hnRNPC_rand_xx.loc[:,['-3','-2','-1','0','1','2','3']].apply(numpy.mean, axis=1).sort_values().index
    index_mean = hnRNPC_rand_xx.loc[:,['-3','-2','-1','0','1','2','3']].apply(numpy.mean, axis=1).sort_values().index
    conf_sort_by_meanShape = list(pd.Series(hnRNPC_rand_conf)[index_mean])
    log_conf = [float(i) for i in conf_sort_by_meanShape]
        
    Len = int( len(log_conf) * 0.15 )
    preList = log_conf[:Len]
    postList = log_conf[-Len:]
    #pV = scipy.stats.ttest_ind(preList, postList)[1]
    pV = scipy.stats.mannwhitneyu(preList, postList)[1]
    plt.subplot(1,2,index); index+=1
    plt.violinplot([preList, postList], showmedians=True)
    plt.title(my_protein)
    plt.ylim(0.0, 1.6)
    
    diff = numpy.median(log_conf[-Len:]) - numpy.median(log_conf[:Len])
    dataSet.append([ pV, diff, my_protein])

print DataSize
plt.savefig(HOME+"/rbp.pdf")
plt.show()


dataSet = pd.DataFrame(dataSet, columns=['p_value', 'diff', 'protein'])
print dataSet





