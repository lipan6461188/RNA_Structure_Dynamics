

from icSHAPE import *
import tools
from ParseTrans import *

font = {'family': 'normal', 'weight': 'bold','size': 15}
matplotlib.rc('font', **font)

def readPeak(peak_file, max_p_cutoff=1):
    peaks = {}
    IN = open(peak_file)
    line = IN.readline()
    while line:
        chr_id, start, end, X, nums, strand, p_value = line.strip().split()
        if float(p_value) <= max_p_cutoff:
            locus = (chr_id, strand)
            if locus not in peaks:
                peaks[ locus ] = []
            peaks[ locus ].append( (int(start), int(end), float(nums), float(p_value)) )
        line = IN.readline()
    IN.close()
    return peaks

def peak_nums(peak):
    nums = 0
    for locus in peak:
        nums += len(peak[locus])
    return nums

def peak_overlap(peaks_1, peaks_2, locus_list=[], verbose=True):
    uniq_1 = {}
    uniq_2 = {}
    common = {}
    label_2 = {}
    
    ## update peaks_1 and peaks_2
    if locus_list:
        peaks_1 = { k:peaks_1[k] for k in peaks_1 if k in locus_list }
        peaks_2 = { k:peaks_2[k] for k in peaks_2 if k in locus_list }
    
    for locus in peaks_2:
        label_2[locus] = [0]*len(peaks_2[locus])
    for locus in peaks_1:
        if locus not in peaks_2:
            uniq_1[locus] = peaks_1[locus]
            continue
        index = 0
        for peak1 in peaks_1[locus]:
            s,e = peak1[:2]
            find = False
            while index < len(peaks_2[locus]):
                s_2,e_2 = peaks_2[locus][index][:2]
                if e<s_2:
                    index = max(0, index-1)
                    break
                if s_2 < e and s < e_2:
                    if locus not in common:
                        common[locus] = []
                    common[locus].append( (peak1, peaks_2[locus][index]) )
                    label_2[locus][index] = 1
                    find = True
                index += 1
            if not find:
                if locus not in uniq_1:
                    uniq_1[locus] = []
                uniq_1[locus].append(peak1)
    for locus in peaks_2:
        if locus not in peaks_1:
            uniq_2[locus] = peaks_2[locus]
            continue
        for index in range(len(peaks_2[locus])):
            if label_2[locus][index] == 0:
                if locus not in uniq_2:
                    uniq_2[locus] = []
                uniq_2[locus].append(peaks_2[locus][index])
    if verbose:
        print peak_nums(uniq_1), peak_nums(common), peak_nums(uniq_2)
    return uniq_1, common, uniq_2

def uniq_list(raw_list):
    my_uniq_list = []
    for item in raw_list:
        if item not in my_uniq_list:
            my_uniq_list.append(item)
    return my_uniq_list

def parse_common_peaks(common_peaks):
    common_1 = {}
    common_2 = {}
    for locus in common_peaks:
        if locus not in common_1:
            common_1[locus] = []
        if locus not in common_2:
            common_2[locus] = []
        for item in common_peaks[locus]:
            common_1[locus].append( item[0] )
            common_2[locus].append( item[1] )
    for locus in common_1: common_1[locus] = uniq_list(common_1[locus])
    for locus in common_2: common_2[locus] = uniq_list(common_2[locus])
    return common_1, common_2


def overlap_3_peaks(peak_1, peak_2, peak_3):
    uniq_1,common_12,uniq_2 = peak_overlap(peak_1, peak_2, verbose=False)
    print peak_nums(uniq_1), peak_nums(uniq_2), peak_nums(common_12)
    
    uniq_1,common_13,uniq_2 = peak_overlap(peak_1, peak_3, verbose=False)
    print peak_nums(uniq_1), peak_nums(uniq_2), peak_nums(common_13)
    
    uniq_1,common_23,uniq_2 = peak_overlap(peak_2, peak_3, verbose=False)
    print peak_nums(uniq_1), peak_nums(uniq_2), peak_nums(common_23)

def common_peaks_pair(common_peaks):
    common_list_1 = []
    common_list_2 = []
    for locus in common_peaks:
        common_list_1 += [item[0] for item in common_peaks[locus]]
        common_list_2 += [item[1] for item in common_peaks[locus]]
    return common_list_1, common_list_2

def read_m6A(inFile, expand=0):
    m6A = {}
    for line in open(inFile):
        chr_id,start,end,dot,fold_change,strand = line.strip().split()
        start,end,fold_change = int(start),int(end), fold_change
        locus = ( chr_id, strand )
        if locus not in m6A:
            m6A[locus] = []
        m6A[locus].append( (start-expand, end+expand, fold_change) )
    for locus in m6A:
        m6A[locus].sort(key=lambda x: x[0])
    return m6A

def filter_m6A(m6A, max_cutoff):
    new_m6A = {}
    for locus in m6A:
        new_m6A[locus] = [ item for item in m6A[locus] if float(item[2]) <= max_cutoff ]
    return new_m6A

def parse_site_from_region(m6A, m6A_pattern, SeqHandle, findBest=False):
    def find_best_m6A(raw_m6A_locus, length):
        middle_site = length/2
        sorted_m6A_locus = sorted( raw_m6A_locus, key=lambda x: abs(x[0][0]-middle_site) )
        return [ sorted_m6A_locus[0] ]
    
    import tools
    sub_locus = {}
    for locus in m6A:
        for m6A_site in m6A[locus]:
            start, end = int(m6A_site[0]), int(m6A_site[1])
            sequence = SeqHandle.fetch(locus[0], Start=start-1, End=end, Strand=locus[1])
            m6A_locus = tools.find_all_match(m6A_pattern, sequence)
            if len(m6A_locus) == 0: continue
            if findBest: m6A_locus = find_best_m6A(m6A_locus, len(sequence))
            for cur_m6A_locus in m6A_locus:
                if locus not in sub_locus:
                    sub_locus[locus] = []
                #s,e = cur_m6A_locus[0][0] + start - 1, cur_m6A_locus[0][1] + start - 1
                if locus[1] == '+':
                    s = start + cur_m6A_locus[0][0] - 1
                    e = start + cur_m6A_locus[0][1] - 1
                else:
                    length = end - start + 1
                    s = start + length - cur_m6A_locus[0][1]
                    e = start + length - cur_m6A_locus[0][0]
                sub_locus[locus].append((s, e, cur_m6A_locus[1]))
    return sub_locus

def shuffle_m6A(m6A, max_dist):
    shuffled_m6A = {}
    for locus in m6A:
        shuffled_m6A[locus] = []
        for site in m6A[locus]:
            randV = random.randint(-max_dist, max_dist)
            new_site = ( max(0, site[0]+randV), max(0, site[1]+randV) )
            shuffled_m6A[locus].append( new_site )
        shuffled_m6A[locus].sort(key=lambda x: x[0])
    return shuffled_m6A

def prepare_shuffle_regions(m6A, WT, KO):
    shuffled_regions = {}
    shuffled_regions.update(m6A)
    shuffled_regions.update(WT)
    shuffled_regions.update(KO)
    
    # intergrate_regions
    intergrate_region_list = []
    for locus in shuffled_regions:
        for site in shuffled_regions[locus]:
            intergrate_region_list.append( (locus, site) )
    
    return intergrate_region_list

def shuffle_m6A_2(shuffle_num, intergrate_region_list, max_dist):
    import random
    
    # Shuffle m6A Sites
    shuffled_m6A = {}
    
    sampled_m6A = random.sample(intergrate_region_list, shuffle_num)
    for s_site in sampled_m6A:
        locus, site = s_site
        if locus not in shuffled_m6A:
            shuffled_m6A[locus] = []
        randV = random.randint(-max_dist, max_dist)
        new_site = ( max(0, site[0]+randV), max(0, site[1]+randV) )
        shuffled_m6A[locus].append( new_site )
    
    # Clean my Shuffled dataset
    for locus in shuffled_m6A:
        shuffled_m6A[locus].sort(key=lambda x: x[0])
    
    return shuffled_m6A

def fetch_peak_around_m6A(peaks, m6A, tolerate=50):
    new_peaks = {}
    for locus in peaks:
        if locus not in m6A:
            continue
        for peak in peaks[locus]:
            for m6A_site in m6A[locus]:
                if peak[0]<=m6A_site[1]+tolerate and m6A_site[0]<=peak[1]+tolerate:
                    if locus not in new_peaks:
                        new_peaks[locus] = []
                    new_peaks[locus].append(peak)
                    break
    return new_peaks

def choose_longest_transcript(transCoor_list, Parser):
    gene2trans = {}
    for item in transCoor_list:
        chr_id, ch_start, ch_end, trans_id, trans_start, trans_end = item
        trans_info = Parser.getTransFeature(transID=trans_id)
        geneName = trans_info['gene_name']
        transLen = trans_info['trans_len']
        if geneName not in gene2trans:
            gene2trans[geneName] = ( trans_id, transLen )
        else:
            if transLen > gene2trans[geneName][1]:
                gene2trans[geneName] = ( trans_id, transLen )
    preserve_trans_list = [ gene2trans[k][0] for k in gene2trans ]
    return [ item for item in transCoor_list if item[3] in preserve_trans_list ]

def locate_transcript(peaks, Parser, min_range=20):
    Trans_Peak = {}
    for locus in peaks:
        chr_id, strand = locus
        for region in peaks[locus]:
            start, end = region[:2]
            left = region[2:]
            try:
                #print chr_id, start, end, strand
                trans_locus = Parser.genomeCoor2transCoor(chr_id, start+1, end, strand)
                trans_locus = choose_longest_transcript(trans_locus, Parser)
            except KeyError:
                continue
            for tl in trans_locus:
                trans_id, start, end = tl[3:6]
                if end-start<min_range:
                    continue
                if trans_id not in Trans_Peak:
                    Trans_Peak[trans_id] = []
                Trans_Peak[trans_id].append((start, end, left))
    return Trans_Peak

def prepare_m6A_shuffle_db(m6A_trans, Parser, mode='randomA'):
    import random, tools
    assert mode in ('randomA', 'randomAC', 'randomMotif')
    
    randomDB = {}
    
    for trans_id in m6A_trans:
        try:
            transLen = Parser.getTransFeature(transID=trans_id)['trans_len']
        except KeyError:
            continue
        for m6A_pos in m6A_trans[trans_id]:
            start, end = m6A_pos[0], m6A_pos[1]
            length = end - start
            
            if 10 >= transLen-10-length: continue
            if trans_id not in randomDB:
                randomDB[trans_id] = []
            
            sequence = Parser.getTransFeatureSeq(transID=trans_id, feature='whole', source='GENCODE', verbose=False)
            if mode == 'randomA':
                AList = tools.find_all_match('A', sequence[10:transLen-10-length])
            elif mode == 'randomAC':
                AList = tools.find_all_match('AC', sequence[10:transLen-10-length])
            elif mode == 'randomMotif':
                AList = tools.find_all_match('[TAG][GAT][GA]AC[TAC][GCT]', sequence[10:transLen-10-length])
            randomDB[trans_id].extend( [ it[0][0]+10 for it in AList ] )
    
    return randomDB


def shuffle_m6A_trans(m6A_trans, Parser, mode='random', Shuffle_db=None):    
    import random, tools
    assert mode in ('random', 'randomDB')
    
    if mode == 'randomDB' and not Shuffle_db:
        print "Error: You should specify Shuffle_db if you use randomDB mode"
    
    if mode == 'random' and Shuffle_db:
        print "Warning: Shuffle_db is unusable in randomDB mode"
    
    shuffled_m6A_trans = {}
    
    for trans_id in m6A_trans:
        try:
            transLen = Parser.getTransFeature(transID=trans_id)['trans_len']
        except KeyError:
            print trans_id, 'not found'
            continue
        
        if mode == 'random':
            for m6A_pos in m6A_trans[trans_id]:
                start, end = m6A_pos[0], m6A_pos[1]
                length = end - start
                
                if 10 >= transLen-10-length: continue
                
                if trans_id not in shuffled_m6A_trans:
                    shuffled_m6A_trans[trans_id] = []
                
                if mode == 'random':
                    randPos = random.randint(10, transLen-10-length)
                    shuffled_m6A_trans[trans_id].append((randPos, randPos+length))
        elif mode == 'randomDB':
            num = len(m6A_trans[trans_id])
            if len(Shuffle_db[trans_id]) < num: num = len(Shuffle_db[trans_id])
            len_list = [ it[1]-it[0] for it in m6A_trans[trans_id] ]
            try:
                randPos = random.sample(Shuffle_db[trans_id], num)
            except KeyError:
                print 'Warning: %s not in Shuffle_db' % (trans_id, )
            shuffled_m6A_trans[trans_id] = [ ( randPos[i], randPos[i]+len_list[i] ) for i in range(num) ]
    
        shuffled_m6A_trans[trans_id].sort(key=lambda x: x[0])
    
    return shuffled_m6A_trans

def out_common_reads_num(common_peaks, outFile):
    OUT = open(outFile, 'w')
    for locus in common_peaks:
        for site in common_peaks[locus]:
            print >>OUT, "%s\t%s" % ( site[0][2], site[1][2] )
    OUT.close()

def shuffle_diff_num(KO, WT, m6A):
    
    ## Shuffle m6A in m6A region and LIN28A regions(in Genome)
    
    m6A_KO_peaks = fetch_peak_around_m6A(KO, m6A, tolerate=5);
    m6A_WT_peaks = fetch_peak_around_m6A(WT, m6A, tolerate=5);
    uniq_KO,common,uniq_WT = peak_overlap(m6A_KO_peaks, m6A_WT_peaks, verbose=False)
    true_common = peak_nums(common)
    true_uniq_sum = peak_nums(uniq_KO)+peak_nums(uniq_WT)
    true_diff = peak_nums(uniq_KO)-peak_nums(uniq_WT)
    true_ratio = 1.0*peak_nums(uniq_KO)/(peak_nums(uniq_KO)+peak_nums(uniq_WT))
    diff_list = []; KO_ratio = []
    
    intergrate_region_list = prepare_shuffle_regions(m6A, WT, KO)
    
    while len(KO_ratio) < 100:
        s_m6A = shuffle_m6A_2(peak_nums(m6A), intergrate_region_list, max_dist=50000)
        m6A_KO_peaks = fetch_peak_around_m6A(KO, s_m6A, tolerate=5)
        m6A_WT_peaks = fetch_peak_around_m6A(WT, s_m6A, tolerate=5)
        uniq_KO,common,uniq_WT = peak_overlap(m6A_KO_peaks, m6A_WT_peaks, verbose=False)
        KO_num = peak_nums(uniq_KO); WT_num = peak_nums(uniq_WT)
        if KO_num+WT_num >= true_uniq_sum:
            KO_ratio.append( 1.0*KO_num/(KO_num+WT_num) )
            print len(KO_ratio), KO_ratio[-1]
    
    print tools.permutate_pValue(KO_ratio, true_ratio, mode='high')
    return true_ratio, sorted(KO_ratio, reverse=True)

def sample_peaks(raw_peaks, sample_num):
    # Got a random subset from peaks
    if peak_nums(raw_peaks) <= sample_num:
        print "Error: FATAL Error"
        return None
    peak_list = []
    for trans_id in raw_peaks:
        for peak in raw_peaks[trans_id]:
            peak_list.append((trans_id, peak))
    sampled_peak_list = random.sample(peak_list, sample_num)
    sampled_peaks = {}
    for trans_id, peak in sampled_peak_list:
        if trans_id not in sampled_peaks:
            sampled_peaks[ trans_id ] = []
        sampled_peaks[ trans_id ].append(peak)
    return sampled_peaks

def calc_peaks_m6A_sig(Peaks, m6A, Parser, randomDB, tolerate=5, times=1000):
    # Overlap m6A and Peaks to calc the significance between them
    random_num = []
    for i in range(times):
        if i%100==0: sys.stdout.writelines('.'); sys.stdout.flush()
        random_m6A_trans = shuffle_m6A_trans(m6A, Parser, mode='randomDB', Shuffle_db=randomDB)
        cur_peaks = fetch_peak_around_m6A(Peaks, random_m6A_trans, tolerate=tolerate)
        random_num.append( peak_nums(cur_peaks) )
    
    true_peaks = fetch_peak_around_m6A(Peaks, m6A, tolerate=tolerate)
    print peak_nums(true_peaks), tools.permutate_pValue(sorted(random_num), peak_nums(true_peaks), mode='low')
    return peak_nums(true_peaks), random_num

def show_peaks_info(Peaks, Parser):
    index = 1
    geneParser = Parser.getGeneParser()
    for locus in Peaks:
        for position in Peaks[locus]:
            gene_list = Parser.genomeCoor2geneCoor(locus[0], position[0], position[1], locus[1])
            geneNameList = [ geneParser[it[3]]['gene_name'] for it in gene_list ]
            if len(geneNameList) == 0: continue
            geneNameStr = ";".join(geneNameList)
            print "%s %s:%s-%s %s %s %s %s" % (index, locus[0], position[0], position[1], locus[1], position[2], position[3], geneNameStr)
            index += 1

def covert_m6A_2_bedGraph(m6A, outFile, SeqHandle, m6A_pattern):
    OUT = open(outFile, 'w')
    for locus in m6A:
        for item in m6A[locus]:
            start, end = item[:2]
            sequence = SeqHandle.fetch(locus[0], Start=start-1, End=end, Strand=locus[1])
            m6A_sites = tools.find_all_match(m6A_pattern, sequence)
            print >>OUT, "%s\t%s\t%s\t%s" % (locus[0], start-1, end, 0.1)
            for region,subseq in m6A_sites:
                if locus[1] == '+':
                    s = start + region[0] - 1
                    e = start + region[1] - 1
                else:
                    length = end - start + 1
                    s = start + length - region[1]
                    e = start + length - region[0]
                    #s = start + region[0] + 1
                    #e = start + region[1] + 1
                print >>OUT, "%s\t%s\t%s\t%s" % (locus[0], s-1, e, 1.0)
    OUT.close()

def covert_Peak_2_bedGraph(peaks, outFile):
    OUT = open(outFile, 'w')
    for locus in peaks:
        for peak in peaks[locus]:
            print >>OUT, "%s\t%s\t%s\t%s" % (locus[0], peak[0]-1, peak[1], 1)
    OUT.close()



######################################
#  Plot KO_Peaks/WT_Peaks around m6A Sites
######################################

def generate_m6A_Peaks(m6A, window_size=50, extend=100):
    m6A_Peaks = {}
    for locus in m6A:
        chr_id, strand = locus
        m6A_Peaks[locus] = []
        for pos in m6A[locus]:
            pos = pos[0]
            start = pos-extend
            end = start+window_size
            while end<=pos+extend:
                m6A_Peaks[locus].append( (start, end) )
                start += window_size
                end = start+window_size
    return m6A_Peaks


def load_bam_bed(inFile):
    Reads = {}
    for line in open(inFile):
        data = line.strip().split()
        locus = (data[0], data[5])
        start = int(data[1])
        end = int(data[2])
        if locus not in Reads: Reads[locus] = []
        Reads[locus].append( (start, end) )
    for locus in Reads: Reads[locus].sort( key=lambda x: x[0] )
    return Reads

def count_reads_num(BED, chr_id, strand, start, end, min_overlap=10):
    count = 0
    locus = (chr_id, strand)
    
    for region in BED[locus]:
        if start<region[1] and region[0]<end:
            if min(region[1], end) - max(start, region[0]) >= min_overlap:
                count += 1
    
    return count

def get_peak_reads_num_ratio(BED_1, BED_2, Peaks, min_load=20, min_overlap=10):
    read_list_1 = []
    read_list_2 = []
    index = 0
    print len(Peaks)
    for locus in Peaks:
        sys.stdout.writelines('.'); sys.stdout.flush()
        chr_id, strand = locus
        if locus not in BED_1 or locus not in BED_2: continue
        for peak in Peaks[locus]:
            start, end = peak[:2]
            read_num_1 = count_reads_num(BED_1, chr_id, strand, start, end, min_overlap=min_overlap)
            read_num_2 = count_reads_num(BED_2, chr_id, strand, start, end, min_overlap=min_overlap)
            if read_num_1+read_num_2>=min_load:
                read_list_1.append(read_num_1)
                read_list_2.append(read_num_2)
    sys.stdout.writelines('\n'); sys.stdout.flush()
    return read_list_1, read_list_2

def plot_KO_WT_dist(KO_reads, WT_reads, bins=30, min_load=0):
    log_Ratio = []; raw_Ratio = []
    for r_KO, r_WT in zip(KO_reads, WT_reads):
        if r_KO+r_WT >= min_load:
            if r_WT != 0:
                ratio = 1.0*r_KO/r_WT
                log_Ratio.append(np.log(ratio))
                raw_Ratio.append(ratio)
            else:
                log_Ratio.append(np.log(10.0))
                raw_Ratio.append(10.0)
    
    print 'x<=-1: ', len([it for it in log_Ratio if it<=-1])
    print '-1<x<0: ', len([it for it in log_Ratio if -1<it<0])
    print '0<x<1: ', len([it for it in log_Ratio if 0<it<1])
    print '1<=x: ', len([it for it in log_Ratio if 1<=it])
    sns.distplot(log_Ratio, bins=bins, hist=False)
    plt.axvline(x=0, ymin=-1, ymax=1, linewidth=2, color='k')
    plt.axvline(x=1, ymin=-1, ymax=1, linewidth=2, color='k')
    plt.axvline(x=-1, ymin=-1, ymax=1, linewidth=2, color='k')
    return raw_Ratio

def combine_peaks(Peaks_1, Peaks_2):
    locus_set = set(Peaks_1) | set(Peaks_2)
    
    all_Peaks = {}
    for locus in locus_set:
        if (locus in Peaks_1) and (locus in Peaks_2):
            all_Peaks[locus] = []
            peaks_1_label = [0]*len(Peaks_2[locus])
            for region_1 in Peaks_1[locus]:
                find = False
                for idx, region_2 in enumerate(Peaks_2[locus]):
                    if region_1[0] <= region_2[1] and region_2[0] <= region_1[1]:
                        peaks_1_label[idx] = 1
                        start, end = min(region_1[0], region_2[0]), max(region_1[1], region_2[1])
                        all_Peaks[locus].append((start, end))
                        find = True
                        break
                if not find:
                    all_Peaks[locus].append(region_1)
            for idx,region_2 in enumerate(Peaks_2[locus]):
                if peaks_1_label[idx] == 0:
                    all_Peaks[locus].append(region_2)
        elif locus in Peaks_1:
            all_Peaks[locus] = Peaks_1[locus]
        elif locus in Peaks_2:
            all_Peaks[locus] = Peaks_2[locus]
        all_Peaks[locus].sort(key=lambda x: x[0])
    
    return all_Peaks



#####################
# Load Parser File
#####################

mm10_parser = ParseTransClass(  genomeCoorBedFile = "/150T/zhangqf/GenomeAnnotation/Gencode/mm10.genomeCoor.bed",
                                seqFileName="/150T/zhangqf/GenomeAnnotation/Gencode/mm10_transcriptome.fa")

#####################
# Load Peak Files
#####################

ROOT_KO = "/Share/home/zhangqf8/lipan/DYNAMIC/LIN28A/KO/"
ROOT_WT = "/Share/home/zhangqf8/lipan/DYNAMIC/LIN28A/WT/"

KO = readPeak(ROOT_KO+"KO_all_pir.bed", max_p_cutoff=1e-3); print peak_nums(KO)
WT = readPeak(ROOT_WT+"WT_all_pir.bed", max_p_cutoff=1e-3); print peak_nums(WT)

KO_trans = locate_transcript(KO, mm10_parser, min_range=20); print peak_nums(KO_trans)
WT_trans = locate_transcript(WT, mm10_parser, min_range=20); print peak_nums(WT_trans)




#####################
# Read m6A Data
#####################

m6A_Pattern = "[TAG][GAT][GA]AC[TAC][GCT]"
import seq
mm10_seq = seq.seqClass("/150T/zhangqf/GenomeAnnotation/genome/GRCm38.p5.genome.fa")

m6A = read_m6A("/Share/home/zhangqf8/lipan/DYNAMIC/LIN28A/m6A/mm10_m6A.bed", expand=0); print peak_nums(m6A)
m6A = filter_m6A(m6A, max_cutoff=0.8); print peak_nums(m6A)

m6A = parse_site_from_region(m6A, m6A_Pattern, mm10_seq, findBest=False); print peak_nums(m6A)
m6A_trans = locate_transcript(m6A, mm10_parser, min_range=1); print peak_nums(m6A_trans)



####### Read Raw Reads Coverage

bamBedFile_KO = "/Share/home/zhangqf8/lipan/DYNAMIC/LIN28A/KO/KO_all_sample_2.bed"
bamBedFile_WT = "/Share/home/zhangqf8/lipan/DYNAMIC/LIN28A/WT/WT_all.bed"

KO_bed = load_bam_bed(bamBedFile_KO)
WT_bed = load_bam_bed(bamBedFile_WT)


####### Get LIN28A Peaks around m6A

m6A_KO_peaks = fetch_peak_around_m6A(KO, m6A, tolerate=50); print peak_nums(m6A_KO_peaks)
m6A_WT_peaks = fetch_peak_around_m6A(WT, m6A, tolerate=50); print peak_nums(m6A_WT_peaks)
uniq_KO,common,uniq_WT = peak_overlap(m6A_KO_peaks, m6A_WT_peaks)

####### Distribution

cPeaks = combine_peaks(m6A_KO_peaks, m6A_WT_peaks); print peak_nums(cPeaks)

KO_reads, WT_reads = get_peak_reads_num_ratio(KO_bed, WT_bed, cPeaks, min_load=10, min_overlap=10)
Ratio = plot_KO_WT_dist(KO_reads, WT_reads, bins=30)
plt.xlim(-2,2)
plt.savefig("figs/dist.pdf")
plt.show()

scipy.stats.ttest_1samp(np.log(Ratio), 0.0)


#####################
# RNA-Seq
#####################


RNA_Seq_KO = load_bam_bed("/Share/home/zhangqf8/lipan/DYNAMIC/LIN28A/RNA-Seq/KO_sample.bed")
RNA_Seq_WT = load_bam_bed("/Share/home/zhangqf8/lipan/DYNAMIC/LIN28A/RNA-Seq/WT_sample.bed")

KO_bed = load_bam_bed("/Share/home/zhangqf8/lipan/DYNAMIC/LIN28A/KO/KO_all_sample_2.bed")
WT_bed = load_bam_bed("/Share/home/zhangqf8/lipan/DYNAMIC/LIN28A/WT/WT_all.bed")

cPeaks = combine_peaks(m6A_KO_peaks, m6A_WT_peaks); print peak_nums(cPeaks)

m6A_KO_peaks = fetch_peak_around_m6A(KO, m6A, tolerate=50); print peak_nums(m6A_KO_peaks)
m6A_WT_peaks = fetch_peak_around_m6A(WT, m6A, tolerate=50); print peak_nums(m6A_WT_peaks)
uniq_KO,common,uniq_WT = peak_overlap(m6A_KO_peaks, m6A_WT_peaks)

cPeaks = combine_peaks(m6A_KO_peaks, m6A_WT_peaks); print peak_nums(cPeaks)

KO_reads_RNASeq, WT_reads_RNASeq = get_peak_reads_num_ratio(RNA_Seq_KO, RNA_Seq_WT, cPeaks, min_load=0, min_overlap=10)
KO_reads_CLIP, WT_reads_CLIP = get_peak_reads_num_ratio(KO_bed, WT_bed, cPeaks, min_load=0, min_overlap=10)


rSeq_ratio = []; Clip_ratio = []
for r_KO_rSeq, r_WT_rSeq, r_KO_Clip, r_WT_Clip in zip(KO_reads_RNASeq, WT_reads_RNASeq, KO_reads_CLIP, WT_reads_CLIP):
    if r_KO_rSeq+r_WT_rSeq>=10 and r_KO_Clip+r_WT_Clip>=10:
        if r_WT_rSeq != 0: ratio = 1.0*r_KO_rSeq/r_WT_rSeq
        else: ratio = 10.0
        rSeq_ratio.append( np.log2(ratio) )
        if r_WT_Clip != 0: ratio = 1.0*r_KO_Clip/r_WT_Clip
        else: ratio = 10.0
        Clip_ratio.append( np.log2(ratio) )


scipy.stats.pearsonr(rSeq_ratio, Clip_ratio)
plt.scatter(rSeq_ratio, Clip_ratio)
plt.xlabel("log2(RNA_seq_KO/RNA_seq_WT)")
plt.ylabel("log2(CLIP_seq_KO/CLIP_seq_WT)")
plt.tight_layout()
plt.show()






