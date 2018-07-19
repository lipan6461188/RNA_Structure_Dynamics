
from icSHAPE import *
import tools


##########################
##    Some Functions
##########################

def get_lineplot_matrix(pos_sites_ensembl, neg_sites_ensembl, pos_sites_label, neg_sites_label):
    xx = []
    sig_points = []
    for idx_1 in range(21):
        for idx_2 in range(len(pos_sites_ensembl[idx_1])):
            xx.append( [idx_1-10, pos_sites_label, pos_sites_ensembl[idx_1][idx_2] ] )
        for idx_2 in range(len(neg_sites_ensembl[idx_1])):
            xx.append( [idx_1-10, neg_sites_label, neg_sites_ensembl[idx_1][idx_2] ] )
        p = scipy.stats.mannwhitneyu(pos_sites_ensembl[idx_1], neg_sites_ensembl[idx_1])[1]
        if p < 0.01:
            my_max = max( numpy.mean(pos_sites_ensembl[idx_1]), numpy.mean(neg_sites_ensembl[idx_1]) )
            sig_points.append( [idx_1, my_max+0.08 ] )
    xx = pd.DataFrame(xx, columns=['site', 'type', 'icshape'])
    return xx, sig_points

def plot_line(xx, sig_points, title, ylim=(0.10, 0.50)):
    ax = sns.pointplot(x="site", y="icshape", hue="type", data=xx)
    plt.title(title)
    for pair in sig_points:
        ax.text(pair[0], pair[1], s="*", ha ='center', fontsize = 25, color=(202/255.0, 75/255.0, 78/255.0,))
    if ylim: plt.ylim(*ylim)

def fetch_common_icSHAPE_of_modification(m6A, icSHAPE, start=-10, end=11, M=2):
    
    ch_trans_set = set(icSHAPE['ch'].keys())
    np_trans_set = set(icSHAPE['np'].keys())
    cy_trans_set = set(icSHAPE['cy'].keys())
    
    ch = []
    np_1 = []
    np_2 = []
    cy = []
    
    ch_np_count = 0
    np_cy_count = 0
    
    for trans_id in m6A:
        in_ch = (trans_id in ch_trans_set)
        in_np = (trans_id in np_trans_set)
        in_cy = (trans_id in cy_trans_set)
        
        for pos in m6A[trans_id]:
            s = pos+start
            e = pos+end
            if s<0: continue
            
            # for ch/np
            if in_ch and in_np:
                if e < len(icSHAPE['ch'][trans_id]):
                    ch_shape = icSHAPE['ch'][trans_id][s:e]
                    np_shape = icSHAPE['np'][trans_id][s:e]
                    if ch_shape.count('NULL')<=M and np_shape.count('NULL')<=M:
                        ch_np_count += 1
                        for shape,idx in zip( ch_shape, range(start, end) ):
                            if shape != 'NULL':
                                ch.append( (float(shape), idx, 'ch') )
                        for shape,idx in zip( np_shape, range(start, end) ):
                            if shape != 'NULL':
                                np_1.append( (float(shape), idx, 'np') )
            
            # for np/cy
            if in_np and in_cy:
                if e < len(icSHAPE['np'][trans_id]):
                    np_shape = icSHAPE['np'][trans_id][s:e]
                    cy_shape = icSHAPE['cy'][trans_id][s:e]
                    if np_shape.count('NULL')<=M and cy_shape.count('NULL')<=M:
                        np_cy_count += 1
                        for shape,idx in zip( np_shape, range(start, end) ):
                            if shape != 'NULL':
                                np_2.append( (float(shape), idx, 'np') )
                        for shape,idx in zip( cy_shape, range(start, end) ):
                            if shape != 'NULL':
                                cy.append( (float(shape), idx, 'cy') )
    
    # calc significance
    sig_ch_np = calcSig(ch, np_1, start=start, end=end)
    sig_np_cy = calcSig(np_2, cy, start=start, end=end)
    
    ch_np_df = pd.DataFrame(ch+np_1, columns=['icshape', 'site', 'type'])
    np_cy_df = pd.DataFrame(np_2+cy, columns=['icshape', 'site', 'type'])
    
    print "ch_np data size: ", ch_np_count
    print "np_cy data size: ", np_cy_count
    return ch_np_df, np_cy_df, sig_ch_np, sig_np_cy

def calcSig(compart_1, compart_2, start, end):
    import numpy
    
    SigPoints = []
    for idx in range(start, end):
        shape_list_1 = [ it[0] for it in compart_1 if it[1] == idx ]
        shape_list_2 = [ it[0] for it in compart_2 if it[1] == idx ]
        p = scipy.stats.mannwhitneyu(shape_list_1, shape_list_2)[1]
        if p < 0.01:
            my_max = max( numpy.mean(shape_list_1), numpy.mean(shape_list_2) )
            SigPoints.append( [idx-start, my_max+0.08, idx] )
    
    return SigPoints

def fetch_icSHAPE_of_modification(m6A, icSHAPE, start=-10, end=11, M=2):
    
    ch_trans_set = set(icSHAPE['ch'].keys())
    np_trans_set = set(icSHAPE['np'].keys())
    cy_trans_set = set(icSHAPE['cy'].keys())
    
    ch = []
    np = []
    cy = []
    
    ch_count = 0
    np_count = 0
    cy_count = 0
    
    for trans_id in m6A:
        in_ch = (trans_id in ch_trans_set)
        in_np = (trans_id in np_trans_set)
        in_cy = (trans_id in cy_trans_set)
        
        for pos in m6A[trans_id]:
            s = pos+start
            e = pos+end
            if s<0: continue
            
            # for ch
            if in_ch:
                if e < len(icSHAPE['ch'][trans_id]):
                    ch_shape = icSHAPE['ch'][trans_id][s:e]
                    if ch_shape.count('NULL')<=M:
                        ch_count += 1
                    for shape,idx in zip( ch_shape, range(start, end) ):
                            if shape != 'NULL':
                                ch.append( (float(shape), idx, 'ch') )
            
            # for np
            if in_np:
                if e < len(icSHAPE['np'][trans_id]):
                    np_shape = icSHAPE['np'][trans_id][s:e]
                    if np_shape.count('NULL')<=M:
                        np_count += 1
                    for shape,idx in zip( np_shape, range(start, end) ):
                            if shape != 'NULL':
                                np.append( (float(shape), idx, 'np') )
            
            # for cy
            if in_cy:
                if e < len(icSHAPE['cy'][trans_id]):
                    cy_shape = icSHAPE['cy'][trans_id][s:e]
                    if cy_shape.count('NULL')<=M:
                        cy_count += 1
                    for shape,idx in zip( cy_shape, range(start, end) ):
                            if shape != 'NULL':
                                cy.append( (float(shape), idx, 'cy') )
    
    # calc significance
    sig_ch_np = calcSig(ch, np, start=start, end=end)
    sig_np_cy = calcSig(np, cy, start=start, end=end)
    
    ch_np_df = pd.DataFrame(ch+np, columns=['icshape', 'site', 'type'])
    np_cy_df = pd.DataFrame(np+cy, columns=['icshape', 'site', 'type'])
    
    print "ch data size: ", ch_count
    print "np data size: ", np_count
    print "cy data size: ", cy_count
    
    return ch_np_df, np_cy_df, sig_ch_np, sig_np_cy


def fetch_icSHAPE_random_of_modification(modification, random_mod, icSHAPE, start=-10, end=11, M=2, Sequence={}):
    
    ch_trans_set = set(icSHAPE['ch'].keys())
    np_trans_set = set(icSHAPE['np'].keys())
    cy_trans_set = set(icSHAPE['cy'].keys())
    
    ch = []
    ch_random = []
    np = []
    np_random = []
    cy = []
    cy_random = []
    
    ch_count = 0; ch_rand_count = 0
    np_count = 0; np_rand_count = 0
    cy_count = 0; cy_rand_count = 0
    
    for trans_id in modification:
        in_ch = (trans_id in ch_trans_set)
        in_np = (trans_id in np_trans_set)
        in_cy = (trans_id in cy_trans_set)
        
        for pos in modification[trans_id]:
            s = pos+start
            e = pos+end
            if s<0: continue
            
            # Check Seuquence
            if Sequence:
                sub_seq = Sequence[trans_id][pos-2:pos+3]
                assert sub_seq == 'GGACT'
            
            # for ch
            if in_ch:
                if e < len(icSHAPE['ch'][trans_id]):
                    ch_shape = icSHAPE['ch'][trans_id][s:e]
                    if ch_shape.count('NULL')<=M:
                        ch_count += 1
                        for shape,idx in zip( ch_shape, range(start, end) ):
                            if shape != 'NULL':
                                ch.append( (float(shape), idx, 'ch') )
            
            # for np
            if in_np:
                if e < len(icSHAPE['np'][trans_id]):
                    np_shape = icSHAPE['np'][trans_id][s:e]
                    if np_shape.count('NULL')<=M:
                        np_count += 1
                        for shape,idx in zip( np_shape, range(start, end) ):
                            if shape != 'NULL':
                                np.append( (float(shape), idx, 'np') )
            
            # for cy
            if in_cy:
                if e < len(icSHAPE['cy'][trans_id]):
                    cy_shape = icSHAPE['cy'][trans_id][s:e]
                    if cy_shape.count('NULL')<=M:
                        cy_count += 1
                        for shape,idx in zip( cy_shape, range(start, end) ):
                            if shape != 'NULL':
                                cy.append( (float(shape), idx, 'cy') )
    
    for trans_id in random_mod:
        in_ch = (trans_id in ch_trans_set)
        in_np = (trans_id in np_trans_set)
        in_cy = (trans_id in cy_trans_set)
        
        for pos in random_mod[trans_id]:
            s = pos+start
            e = pos+end
            if s<0: continue
            
            # Check Sequence
            if Sequence:
                sub_seq = Sequence[trans_id][pos-2:pos+3]
                assert sub_seq == 'GGACT'
            
            # for ch
            if in_ch:
                if e < len(icSHAPE['ch'][trans_id]):
                    ch_shape = icSHAPE['ch'][trans_id][s:e]
                    if ch_shape.count('NULL')<=M:
                        ch_rand_count += 1
                        for shape,idx in zip( ch_shape, range(start, end) ):
                            if shape != 'NULL':
                                ch_random.append( (float(shape), idx, 'rand_ch') )
            
            # for np
            if in_np:
                if e < len(icSHAPE['np'][trans_id]):
                    np_shape = icSHAPE['np'][trans_id][s:e]
                    if np_shape.count('NULL')<=M:
                        np_rand_count += 1
                        for shape,idx in zip( np_shape, range(start, end) ):
                            if shape != 'NULL':
                                np_random.append( (float(shape), idx, 'rand_np') )
            
            # for cy
            if in_cy:
                if e < len(icSHAPE['cy'][trans_id]):
                    cy_shape = icSHAPE['cy'][trans_id][s:e]
                    if cy_shape.count('NULL')<=M:
                        cy_rand_count += 1
                        for shape,idx in zip( cy_shape, range(start, end) ):
                            if shape != 'NULL':
                                cy_random.append( (float(shape), idx, 'rand_cy') )
    
    # calc significance
    sig_ch = calcSig(ch, ch_random, start=start, end=end)
    sig_np = calcSig(np, np_random, start=start, end=end)
    sig_cy = calcSig(cy, cy_random, start=start, end=end)
    
    ch_df = pd.DataFrame(ch+ch_random, columns=['icshape', 'site', 'type'])
    np_df = pd.DataFrame(np+np_random, columns=['icshape', 'site', 'type'])
    cy_df = pd.DataFrame(cy+cy_random, columns=['icshape', 'site', 'type'])
    
    print "ch/random data size: ", ch_count, ch_rand_count
    print "np/random data size: ", np_count, np_rand_count
    print "cy/random data size: ", cy_count, cy_rand_count
    
    return ch_df, np_df, cy_df, sig_ch, sig_np, sig_cy


def filter_m6A(m6A, Sequence):
    new_m6A = {}
    for tid in m6A:
        for pos in m6A[tid]:
            sub_seq = Sequence[tid][pos-2:pos+3]
            if sub_seq == 'GGACT':
                if tid not in new_m6A:
                    new_m6A[tid] = []
                new_m6A[tid].append(pos)
    return new_m6A

def generate_randA_withMotif(m6A_ensembl, SeqDict, motif='GGACT', A_Site=3, sample_times=2):
    import tools, random
    
    def allASites_withMotif(seq):
        A_index = []
        all_m6A = tools.find_all_match(motif, seq)
        for it in all_m6A:
            A_index.append( it[0][0]+A_Site-2 )
        return A_index    
    
    all_randA = []
    sample_num = sum([ len(i) for i in m6A_ensembl.values() ])
    for trans_id in m6A_ensembl:
        seq = SeqDict[trans_id]
        ASites = allASites_withMotif(seq)
        for site in ASites:
            if 20 <= site <= len( SeqDict[trans_id] ) - 20 and site not in m6A_ensembl[trans_id]:
                all_randA.append( (trans_id, site) )
    
    if sample_times*sample_num > len(all_randA):
        sample_List = all_randA
    else:
        sample_List = random.sample(all_randA, sample_times*sample_num)
    
    rand_A_ensembl = {}
    for sl_item in sample_List:
        if sl_item[0] not in rand_A_ensembl:
            rand_A_ensembl[sl_item[0]] = []
        rand_A_ensembl[sl_item[0]].append( sl_item[1] )
    
    for trans_id in rand_A_ensembl:
        rand_A_ensembl[ trans_id ].sort()
    
    return rand_A_ensembl

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
        while i < 2*m6A_nums:
            randASite = random.choice(ASites)
            if 20 <= randASite <= len( SeqDict[trans_id] ) - 20:
                rand_A_ensembl[ trans_id ].append( randASite )
                i += 1
        rand_A_ensembl[ trans_id ].sort()
    return rand_A_ensembl


def show_modiBase(modification, Fasta):
    for tid in modification:
        if tid not in Fasta: continue
        for pos in modification[tid]:
            print tid, pos, Fasta[tid][pos]


##########################
##    Load Human icSHAPE
##########################

human_seq = readSeq("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/transcriptome.fa")

icSHAPE_ROOT_DIR = '/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/%s/shape.out'
ch_icshape_file = icSHAPE_ROOT_DIR % ("hek_ch_vivo", )
np_icshape_file = icSHAPE_ROOT_DIR % ("hek_np_vivo", )
cy_icshape_file = icSHAPE_ROOT_DIR % ("hek_cy_vivo", )
ch_vitro_icshape_file = icSHAPE_ROOT_DIR % ("hek_ch_vitro", )
np_vitro_icshape_file = icSHAPE_ROOT_DIR % ("hek_np_vitro", )
cy_vitro_icshape_file = icSHAPE_ROOT_DIR % ("hek_cy_vitro", )

ch_icshape = tools.loadicSHAPE(ch_icshape_file)
np_icshape = tools.loadicSHAPE(np_icshape_file)
cy_icshape = tools.loadicSHAPE(cy_icshape_file)
ch_vitro_icshape = tools.loadicSHAPE(ch_vitro_icshape_file)
np_vitro_icshape = tools.loadicSHAPE(np_vitro_icshape_file)
cy_vitro_icshape = tools.loadicSHAPE(cy_vitro_icshape_file)

human_icshape = {'ch':ch_icshape, 'np':np_icshape , 'cy':cy_icshape, 
        'trans_ch': sorted(ch_icshape.keys()), 'trans_np': sorted(np_icshape.keys()),'trans_cy': sorted(cy_icshape.keys()),
        'ch_vitro': ch_vitro_icshape, 'np_vitro': np_vitro_icshape, 'cy_vitro': cy_vitro_icshape, 
        'trans_ch_vitro': sorted(ch_vitro_icshape.keys()), 'trans_np_vitro': sorted(np_vitro_icshape.keys()),'trans_cy_vitro': sorted(cy_vitro_icshape.keys()) }


##########################
##    Load Human pU Dataset
##########################

pU_nacent_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/pU/newData/nascent_pU.transCoor.bed'
pU_candidate_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/pU/newData/pU_candidate.transCoor.bed'

pU_nacent_ensembl = read_m6A(pU_nacent_file, human_seq)
pU_candidate_ensembl = read_m6A(pU_candidate_file, human_seq)




##########################
##    Load Human m6A Dataset
##########################

m6A_human_file = "/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/m6A/choose_good_protein/m6A_dataset_1.merged.transCoor.bed"

m6A_ensembl_raw = read_m6A(m6A_human_file, human_seq); print pbr_count(m6A_ensembl_raw)
m6A_ensembl = filter_m6A(m6A_ensembl_raw, human_seq); print pbr_count(m6A_ensembl)

rand_m6A = generate_randA_withMotif(m6A_ensembl, human_seq, motif='GGACT', A_Site=3); print pbr_count(rand_m6A)


##########################
## Generate Random U
##########################

randT_ensembl = generate_randA(pU_candidate_ensembl, human_seq, base='T')

pU_icshape = fetch_icshape(pU_candidate_ensembl, human_icshape, start=-10, end=11, valid_cutoff=10)
randT_icshape = fetch_icshape(randT_ensembl, human_icshape, start=-10, end=11, valid_cutoff=10)

pU_sites_ensembl = fetch_sites_ensembl(pU_icshape, length=21)
randT_sites_ensembl = fetch_sites_ensembl(randT_icshape, length=21)


##########################
## Plot pU
##########################

##### Plot in vivo pU profile

ch_xx, ch_sig_points = get_lineplot_matrix(pU_sites_ensembl['ch'], randT_sites_ensembl['ch'], pos_sites_label='pU', neg_sites_label='randU')
np_xx, np_sig_points = get_lineplot_matrix(pU_sites_ensembl['np'], randT_sites_ensembl['np'], pos_sites_label='pU', neg_sites_label='randU')
cy_xx, cy_sig_points = get_lineplot_matrix(pU_sites_ensembl['cy'], randT_sites_ensembl['cy'], pos_sites_label='pU', neg_sites_label='randU')

## Plot

fig = plt.figure(figsize=(30,30))

plt.subplot(3,1,1)
plot_line( ch_xx, ch_sig_points, "CH pU %d %d" % (len(pU_sites_ensembl['ch'][0]), len(randT_sites_ensembl['ch'][1])) )
plt.subplot(3,1,2)
plot_line( np_xx, np_sig_points, "NP pU %d %d" % (len(pU_sites_ensembl['np'][0]), len(randT_sites_ensembl['np'][1])) )
plt.subplot(3,1,3)
plot_line( cy_xx, cy_sig_points, "CY pU %d %d" % (len(pU_sites_ensembl['cy'][0]), len(randT_sites_ensembl['cy'][1])) )

plt.savefig(HOME+"/figs/pU-candidate-vivo.pdf")
plt.show()


##########################
## Plot human m6A with Random (AC)
##########################

m6A_human_file = "/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/m6A/choose_good_protein/m6A_dataset_1.merged.transCoor.bed"
m6A_ensembl = read_m6A(m6A_human_file, human_seq); print pbr_count(m6A_ensembl_raw)
rand_m6A = generate_randA_withMotif(m6A_ensembl, human_seq, motif='AC', A_Site=1); print pbr_count(rand_m6A)

ch_df, np_df, cy_df, sig_ch, sig_np, sig_cy = fetch_icSHAPE_random_of_modification(m6A_ensembl, rand_m6A, human_icshape, start=-10, end=11, M=2, Sequence={})

plt.figure(figsize=(10,10))
plt.subplot(3,1,1)
plot_line(ch_df, sig_ch, title="ch", ylim=(0,0.5))
plt.subplot(3,1,2)
plot_line(np_df, sig_np, title="np", ylim=(0,0.5))
plt.subplot(3,1,3)
plot_line(cy_df, sig_cy, title="cy", ylim=(0,0.5))
plt.savefig("figs/m6A.pdf")
plt.show()


##########################
## Plot human pU with Random (U)
##########################

pU_nacent_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/pU/newData/nascent_pU.transCoor.bed'
pU_candidate_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/pU/newData/pU_candidate.transCoor.bed'

pU_nacent_ensembl = read_m6A(pU_nacent_file, human_seq)
pU_candidate_ensembl = read_m6A(pU_candidate_file, human_seq)

rand_pU = generate_randA(pU_candidate_ensembl, human_seq, base='T'); print pbr_count(rand_pU)

ch_df, np_df, cy_df, sig_ch, sig_np, sig_cy = fetch_icSHAPE_random_of_modification(pU_candidate_ensembl, rand_pU, human_icshape, start=-10, end=11, M=2, Sequence={})

plt.figure(figsize=(10,10))
plt.subplot(3,1,1)
plot_line(ch_df, sig_ch, title="ch", ylim=(0,0.5))
plt.subplot(3,1,2)
plot_line(np_df, sig_np, title="np", ylim=(0,0.5))
plt.subplot(3,1,3)
plot_line(cy_df, sig_cy, title="cy", ylim=(0,0.5))
plt.savefig("figs/pU.pdf")
plt.show()


##########################
## Plot Human HNRNPC with Random (UUUU)
##########################

def parse_site_from_region(RBP_Regions, pattern, Sequence, findBest=False):
    def find_best_m6A(raw_m6A_locus, length):
        middle_site = length/2
        sorted_m6A_locus = sorted( raw_m6A_locus, key=lambda x: abs(x[0][0]-middle_site) )
        return [ sorted_m6A_locus[0] ]
    
    import tools
    bind_sites = {}
    search_count = 0
    for tid in RBP_Regions:
        if tid not in Sequence: continue
        for region in RBP_Regions[tid]:
            start, end = region[0], region[1]
            subSeq = Sequence[tid][start:end]
            locus = tools.find_all_match(pattern, subSeq)
            
            if len(locus) == 0: continue
            search_count += 1
            if findBest: locus = find_best_m6A(locus, len(subSeq))
            
            if tid not in bind_sites:
                bind_sites[tid] = []
            for cur_locus in locus:
                s = start + cur_locus[0][0] - 1
                e = start + cur_locus[0][1] - 1
                bind_sites[tid].append((s, e))
    
    print "Pattern Search Count: %s, Ratio: %.2f%%" % (search_count, 100.0*search_count/pbr_count(RBP_Regions))
    return bind_sites



clip_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/m6A/choose_good_protein/clipDB.merged.transCoor.bed'
pbr_ensembl = read_pbr(clip_file, human_seq)

##### HNRNPC

hnrnpc = parse_site_from_region(pbr_ensembl['HNRNPC'], 'TTTT', human_seq)
rand_hnrnpc = generate_randA_withMotif(hnrnpc, human_seq, motif='TTTT', A_Site=1); print pbr_count(rand_hnrnpc)

hnrnpc = {k:[it[0] for it in hnrnpc[k]] for k in hnrnpc}
ch_df, np_df, cy_df, sig_ch, sig_np, sig_cy = fetch_icSHAPE_random_of_modification(hnrnpc, rand_hnrnpc, human_icshape, start=-10, end=11, M=4, Sequence={})

plt.figure(figsize=(10,10))
plt.subplot(3,1,1)
plot_line(ch_df, sig_ch, title="ch", ylim=(0,0.6))
plt.subplot(3,1,2)
plot_line(np_df, sig_np, title="np", ylim=(0,0.6))
plt.subplot(3,1,3)
plot_line(cy_df, sig_cy, title="cy", ylim=(0,0.6))
plt.savefig("figs/hnrnpc.pdf")
plt.show()


