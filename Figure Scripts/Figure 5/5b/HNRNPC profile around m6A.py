
###### Input Data: /Users/lee/Seafile/icSHAPE_Paper/Nature\ Figure/Data/HNRNPC_m6A.hg19.txt

####################
# 1. No shuffle
####################

def statistic_HNRNPC_m6A(inFile, seqHandle):
    count = [0] * 200
    line_count = 0
    for line in open(inFile):
        line_count += 1
        data = line.strip().split()
        chr_id, strand = data[0][:-1], data[0][-1]
        hnrnpc_s, hnrnpc_e = int(data[1]), int(data[2])+1
        if strand == '+':
            m6A_site = int(data[3])+2
            AC = seqHandle.fetch(chr_id, m6A_site, m6A_site+2, strand)
        elif strand == '-':
            m6A_site = int(data[3])+2
            AC = seqHandle.fetch(chr_id, m6A_site-1, m6A_site+1, strand)
        
        ### 检测
        if AC != 'AC':
            print 'A Error'
            print line.strip(), line_count
            print AC
            return
        if seqHandle.fetch(chr_id, hnrnpc_s, hnrnpc_e, strand) != 'T'*(hnrnpc_e-hnrnpc_s):
            print 'T Error'
            return line.strip(), line_count
            print seqHandle.fetch(chr_id, hnrnpc_s, hnrnpc_e, strand)
            return
        
        #print m6A_site-hnrnpc_s
        ### 计数
        bias_1 = abs(hnrnpc_s - m6A_site)
        bias_2 = abs(hnrnpc_e - m6A_site)
        if hnrnpc_s < hnrnpc_e < m6A_site:
            if strand == '+':
                s = 100 - (m6A_site - hnrnpc_s)
                e = 100 - (m6A_site - hnrnpc_e)
            if strand == '-':
                s = 100 + (m6A_site - hnrnpc_e)
                e = 100 + (m6A_site - hnrnpc_s)
        elif m6A_site < hnrnpc_s < hnrnpc_e:
            if strand == '+':
                s = 100 + (hnrnpc_s - m6A_site)
                e = 100 + (hnrnpc_e - m6A_site)
            if strand == '-':
                s = 100 - (hnrnpc_e - m6A_site)
                e = 100 - (hnrnpc_s - m6A_site)
        for i in range(s, e+1):
            count[i] += 1
    
    return count

import tools, seq

handle = seq.seqClass("/150T/zhangqf/GenomeAnnotation/genome/GRCh37.p13.genome.fa")
count = statistic_HNRNPC_m6A("/tmp/HNRNPC_m6A.hg19.txt", handle)


tools.plt.plot(range(-100,100), count)
tools.plt.savefig("figs/hnrnpc_m6A.pdf")
tools.plt.show()



####################
# 2. Shuffle m6A or PBR
####################



def read_HNRNPC_m6A(inFile):
    HNRNPC = {}
    m6A = {}
    
    for line in open(inFile):
        chr_strand, start, end, m6A_start = line.strip().split()
        locus = (chr_strand[:-1], chr_strand[-1])
        start, end, m6A_start = int(start), int(end), int(m6A_start)
        if locus not in HNRNPC:
            HNRNPC[locus] = []
            m6A[locus] = []
        HNRNPC[locus].append( (start, end) )
        
        m6A_pos = m6A_start+2
        m6A[locus].append( m6A_pos )
    
    for locus in HNRNPC:
        HNRNPC[locus].sort( key=lambda x: x[0] )
        m6A[locus].sort()
    
    return HNRNPC, m6A


def generate_random_m6A(m6A, seqHandle, max_dist=1000):
    rand_m6A = {}
    
    for locus in m6A:
        rand_m6A[locus] = []
        for pos in m6A[locus]:
            motif = seqHandle.fetch(locus[0], pos-1, pos+2, locus[1])
            flank_seq = seqHandle.fetch(locus[0], pos-max_dist, pos+max_dist, locus[1])
            all_motif = tools.find_all_match(motif, flank_seq)
            ss, se = random.sample(all_motif, 1)[0][0]
            #print ss, se
            if locus[1] == '+':
                abs_s, abs_e = pos-max_dist+ss-1, pos-max_dist+se-1
            elif locus[1] == '-':
                abs_s, abs_e = pos+max_dist-se, pos+max_dist-ss
            
            #print (locus, abs_s, abs_e)
            #sample_motif = seqHandle.fetch(locus[0], abs_s, abs_e+1, locus[1])
            #print locus, motif, sample_motif
            rand_m6A[locus].append( abs_s+1 )
    
    return rand_m6A


def generate_random_HNRNPC(HNRNPC, seqHandle, max_dist=1000):
    rand_HNRNPC = {}
    
    for locus in HNRNPC:
        rand_HNRNPC[locus] = []
        for s,e in HNRNPC[locus]:
            motif = seqHandle.fetch(locus[0], s, e+1, locus[1])
            
            flank_seq = seqHandle.fetch(locus[0], s-max_dist, e+max_dist, locus[1])
            all_motif = tools.find_all_match(motif, flank_seq)
            ss, se = random.sample(all_motif, 1)[0][0]
            #print ss, se
            if locus[1] == '+':
                abs_s, abs_e = s-max_dist+ss-1, s-max_dist+se-1
            elif locus[1] == '-':
                abs_s, abs_e = e+max_dist-se, e+max_dist-ss
            
            #print (locus, abs_s, abs_e)
            #sample_motif = seqHandle.fetch(locus[0], abs_s, abs_e+1, locus[1])
            #if s == abs_s:
            #    print locus, motif, sample_motif, s, abs_s, 11111111
            #else:
            #    print locus, motif, sample_motif, s, abs_s
            rand_HNRNPC[locus].append( (abs_s, abs_e) )
    
    return rand_HNRNPC


def count_cov(m6A, HNRNPC):
    
    count = [0] * 200
    
    for locus in set(m6A) & set(HNRNPC):
        strand = locus[1]
        for m6A_site in m6A[locus]:
            for hnrnpc_s, hnrnpc_e in HNRNPC[locus]:
                ### 计数
                bias_1 = abs(hnrnpc_s - m6A_site)
                bias_2 = abs(hnrnpc_e - m6A_site)
                if hnrnpc_s < hnrnpc_e < m6A_site:
                    if strand == '+':
                        s = 100 - (m6A_site - hnrnpc_s)
                        e = 100 - (m6A_site - hnrnpc_e)
                    if strand == '-':
                        s = 100 + (m6A_site - hnrnpc_e)
                        e = 100 + (m6A_site - hnrnpc_s)
                elif m6A_site < hnrnpc_s < hnrnpc_e:
                    if strand == '+':
                        s = 100 + (hnrnpc_s - m6A_site)
                        e = 100 + (hnrnpc_e - m6A_site)
                    if strand == '-':
                        s = 100 - (hnrnpc_e - m6A_site)
                        e = 100 - (hnrnpc_s - m6A_site)
                for i in range(s, e+1):
                    if 0 <= i < 200:
                        count[i] += 1
    
    return count


import tools, seq

handle = seq.seqClass("/150T/zhangqf/GenomeAnnotation/genome/GRCh37.p13.genome.fa")
HNRNPC, m6A = read_HNRNPC_m6A("/tmp/HNRNPC_m6A.hg19.txt")

rand_m6A = generate_random_m6A(m6A, handle, max_dist=2000)
rand_HNRNPC = generate_random_HNRNPC(HNRNPC, handle, max_dist=2000)


true_count = count_cov(m6A, HNRNPC)
shuffle_m6A_count = count_cov(rand_m6A, HNRNPC)
shuffle_HNRNPC_count = count_cov(m6A, rand_HNRNPC)

tools.plt.plot(range(-100,100), true_count, color='#55A868')
tools.plt.plot(range(-100,100), shuffle_HNRNPC_count, color='#C44E52')
tools.plt.plot(range(-100,100), shuffle_m6A_count, color='#8172B2')
tools.plt.savefig("figs/HNRNPC.pdf")
tools.plt.show()












