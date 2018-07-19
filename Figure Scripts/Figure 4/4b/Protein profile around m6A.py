

import tools, icSHAPE, ParseTrans


def pbr_count(pbr_single):
    return sum([len(v) for v in pbr_single.values()])

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


###### Load huamn Seq

human_seq = tools.readSeq("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/transcriptome.fa")

###### Load Human RBP Data

clip_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/m6A/choose_good_protein/clipDB.merged.transCoor.bed'
pbr_ensembl = read_pbr(clip_file, human_seq)

###### Load human Parser

hg38_parseTrans = ParseTrans.ParseTransClass(genomeCoorBedFile = "/tmp/hg38.genomeCoor.bed")

###### Load m6A

m6A_file = '/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/m6A/choose_good_protein/m6A_dataset_1.merged.transCoor.bed'
m6A_ensembl = icSHAPE.read_m6A(m6A_file, human_seq)




def count_rbp(m6A, rbp, shrink=5):
    rbp_cov = [0]*200
    common_set = set(m6A.keys()) & set(rbp.keys())
    for tid in common_set:
        for m6A_pos in m6A[tid]:
            for rbp_region in rbp[tid]:
                start = rbp_region[0] - m6A_pos + shrink
                end = rbp_region[1] - m6A_pos - shrink
                if start >= end: continue
                for pos in range(start, end):
                    if -100 <= pos < 100:
                        rbp_cov[100+pos] += 1
    
    return rbp_cov

def shuffle_rbp(rbp, Parser, times=2):
    shuffled_rbp = {}
    for tid in rbp:
        try:
            length = Parser.getTransFeature(tid, verbose=False)['trans_len']
            if length < 100:
                continue
        except KeyError:
            continue
        for region in rbp[tid]:
            s,e = region[:2]
            if s >= e: continue
            span = e-s
            if span > 80: continue
            for i in range(times):
                rand_s = tools.random.randint(0, length-span)
                rand_e = rand_s + span
                if tid not in shuffled_rbp:
                    shuffled_rbp[tid] = []
                shuffled_rbp[tid].append( [rand_s, rand_e] )
    return shuffled_rbp

def generate_shuffled_pbr_ensembl(pbr_ensembl, Parser, times=2):
    shuffled_pbr_ensembl = {}
    for protein in pbr_ensembl:
        print 'Process %s ....' % (protein, )
        shuffled_pbr_ensembl[protein] = shuffle_rbp(pbr_ensembl[protein], Parser, times=times)
    return shuffled_pbr_ensembl

def generate_shuffled_m6A_ensembl(m6A_ensembl, Sequence, times=2):
    shuffled_m6A_ensembl = {}
    no_found = 0
    for tid in m6A_ensembl:
        try: 
            seq = Sequence[tid]
        except KeyError:
            continue
        for m6A_pos in m6A_ensembl[tid]:
            sub_seq = seq[m6A_pos-1:m6A_pos+2]
            locus = tools.find_all_match(sub_seq, seq)
            my_possible = [ it[0][0]+1-1 for it in locus ]
            
            if m6A_pos in my_possible: 
                my_possible.remove(m6A_pos)
            else:
                no_found += 1
            
            if len(my_possible) == 0: continue
            sampled_m6A = random.sample( my_possible, min(times, len(my_possible)) )
            if tid not in shuffled_m6A_ensembl:
                shuffled_m6A_ensembl[tid] = []
            shuffled_m6A_ensembl[tid] += sampled_m6A
    
    print "not found num: %s" % (no_found, )
    return shuffled_m6A_ensembl

shuffled_pbr_ensembl = generate_shuffled_pbr_ensembl(pbr_ensembl, hg38_parseTrans, times=1)
#shuffled_pbr_ensembl_2 = generate_shuffled_pbr_ensembl(pbr_ensembl, hg38_parseTrans, times=2)

shuffled_m6A_ensembl = generate_shuffled_m6A_ensembl(m6A_ensembl, human_seq, times=1)



###########################
## 画出所有的蛋白结合profile
###########################

# 按照数据量排序所有的蛋白

def get_protein_count(m6A, rbp):
    pbr_count_dict = {}
    for i, protein in enumerate(rbp.keys()):
        pbr_count = count_rbp(m6A, rbp[protein], shrink=8)
        pbr_count_dict[protein] = pbr_count
    return pbr_count_dict

pbr_count_dict = get_protein_count(m6A_ensembl, pbr_ensembl)
pbr_count_list = [ [pn, sum(pbr_count_dict[pn])] for pn in pbr_count_dict ]
pbr_count_list.sort(key=lambda x: x[1], reverse=True)
trans_list = [ it[0] for it in pbr_count_list ]

# 画图

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def draw_profile(trans_list, m6A_ensembl, pbr_ensembl, shuffled_pbr_ensembl=None, shuffled_m6A_ensembl=None, shrink=5):
    tools.plt.close()
    tools.plt.figure(figsize=(30,15))
    index = 1
    for i, protein in enumerate(trans_list[:64]):
        if protein not in ('YTHDF2', 'IGF2BP3'): continue
        pbr_count = count_rbp(m6A_ensembl, pbr_ensembl[protein], shrink=shrink)
        
        if shuffled_pbr_ensembl:
            shuffled_pbr_count = count_rbp(m6A_ensembl, shuffled_pbr_ensembl[protein], shrink=shrink)
        
        if shuffled_m6A_ensembl:
            shuffled_m6A_count = count_rbp(shuffled_m6A_ensembl, pbr_ensembl[protein], shrink=shrink)
        
        tools.plt.subplot(8,8,index)
        index += 1
        tools.plt.plot(range(-100,100), pbr_count, color='#55A868')
        
        if shuffled_pbr_ensembl:
            tools.plt.plot(range(-100,100), shuffled_pbr_count, color='#C44E52')
        
        if shuffled_m6A_ensembl:
            tools.plt.plot(range(-100,100), shuffled_m6A_count, color='#8172B2')
        
        tools.plt.axvline(x=0, ymin=0, ymax = 1000, linewidth=0.5, color='k')
        tools.plt.title(protein)
    
    tools.plt.tight_layout()
    tools.plt.savefig(icSHAPE.HOME+"/figs/rbp_profile.pdf")
    tools.plt.show()


draw_profile(trans_list, m6A_ensembl, pbr_ensembl, shuffled_pbr_ensembl=shuffled_pbr_ensembl, shuffled_m6A_ensembl=shuffled_m6A_ensembl, shrink=1)

