
import icSHAPE, tools

###### Load Human Transcriptome

human_seq = tools.readSeq("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/transcriptome.fa")
mouse_seq = tools.readSeq("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/mouse_transcriptome.fa")


###### Load Human icSHAPE Data

icSHAPE_ROOT_DIR = '/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/%s/shape.out'
ch_icshape_file = icSHAPE_ROOT_DIR % ("hek_ch_vivo", )
np_icshape_file = icSHAPE_ROOT_DIR % ("hek_np_vivo", )
cy_icshape_file = icSHAPE_ROOT_DIR % ("hek_cy_vivo", )

ch_icshape = tools.loadicSHAPE(ch_icshape_file)
np_icshape = tools.loadicSHAPE(np_icshape_file)
cy_icshape = tools.loadicSHAPE(cy_icshape_file)

human_shape = {'ch':ch_icshape, 'np':np_icshape, 'cy':cy_icshape}

###### Load Mouse icSHAPE Data

icSHAPE_ROOT_DIR = '/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/%s/shape.out'
ch_icshape_file = icSHAPE_ROOT_DIR % ("mes_ch_vivo", )
np_icshape_file = icSHAPE_ROOT_DIR % ("mes_np_vivo", )
cy_icshape_file = icSHAPE_ROOT_DIR % ("mes_cy_vivo", )

ch_icshape = tools.loadicSHAPE(ch_icshape_file)
np_icshape = tools.loadicSHAPE(np_icshape_file)
cy_icshape = tools.loadicSHAPE(cy_icshape_file)

mouse_shape = {'ch':ch_icshape, 'np':np_icshape, 'cy':cy_icshape}

###### Compare GINI

def get_common_gini_list(shape_comp_1, shape_comp_2, valid_cutoff=15):
    gini_list = []
    
    for tid in set(shape_comp_1) & set(shape_comp_2):
        valid_shape_1 = []
        valid_shape_2 = []
        for s_1, s_2 in zip(shape_comp_1[tid], shape_comp_2[tid]):
            if 'NULL' not in (s_1,s_2):
                valid_shape_1.append( float(s_1) )
                valid_shape_2.append( float(s_2) )
        gini_1 = tools.calcGINI(valid_shape_1, valid_cutoff=valid_cutoff)
        gini_2 = tools.calcGINI(valid_shape_2, valid_cutoff=valid_cutoff)
        if gini_1 != -1 and gini_2 != -1:
            gini_list.append((tid, gini_1, gini_2))
    
    gini_df = tools.pd.DataFrame(gini_list, columns=['tid', 'gini_1', 'gini_2'])
    print gini_df.shape
    return gini_df


gini_ch_np = get_common_gini_list(mouse_shape['ch'], mouse_shape['np'], valid_cutoff=15)
gini_ch_cy = get_common_gini_list(mouse_shape['ch'], mouse_shape['cy'], valid_cutoff=15)
gini_np_cy = get_common_gini_list(mouse_shape['np'], mouse_shape['cy'], valid_cutoff=15)

gini_ch_np.columns = ['tid', 'gini_ch', 'gini_np']
gini_ch_cy.columns = ['tid', 'gini_ch', 'gini_cy']
gini_np_cy.columns = ['tid', 'gini_np', 'gini_cy']

ax = tools.sns.jointplot(data=gini_ch_np, x='gini_ch', y='gini_np', xlim=(0.2, 1), ylim=(0.2, 1), kind='reg' )
tools.plt.savefig("figs/ch_np.pdf")
tools.plt.close()

ax = tools.sns.jointplot(data=gini_ch_cy, x='gini_ch', y='gini_cy', xlim=(0.2, 1), ylim=(0.2, 1), kind='reg' )
tools.plt.savefig("figs/ch_cy.pdf")
tools.plt.close()

ax = tools.sns.jointplot(data=gini_np_cy, x='gini_np', y='gini_cy', xlim=(0.2, 1), ylim=(0.2, 1), kind='reg' )
tools.plt.savefig("figs/np_cy.pdf")
tools.plt.close()






