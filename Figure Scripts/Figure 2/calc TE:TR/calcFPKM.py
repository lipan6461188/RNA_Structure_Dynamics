#!/Share/home/zhangqf/usr/anaconda/bin/python
#-*- coding:utf-8 -*-

""" Compute FPKM From SAM Files

Limitation:
    1. You Must Sort the SAM File Before Calculate FPKM according to Reads ID so that Same Reads are together
        You Can Sort Your SAM File Like This:
            samtools view -hb test.sam | samtools sort -m 1G -n - | samtools view -h > test.sorted.sam
    
    2. The SAM File Must Include @SQ Tag to Give the information of Length of Transcripts, below is a example:
        @SQ     SN:ENST00000311027_ENSG00000110514      LN:5990
    
    3. The Script Can not be used to Statistic Gene RPKM in SAM file mapped to genome.

Version:
    2017-6-16
"""

import re, sys, os, getopt, time, datetime

Usage = """
## --------------------------------------
Calculate Gene FPKM From SAM File

Command:
%s -i map.sam [-o FPKM_file] [-m 1] [-g 1] [-s 1] [-r 20] [--ribo_occup transCDSfile] [--offset=12] [--rtFile rtFile] 1> file.fpkm

# what it is:
 -i         mapped sam file(must be sorted by read names)
 -o         output fpkm to file instead of standard output ( default standard output )
 -m 0/1     if to consider multi-map, each multi-mapped reads will be equally distributed to each genes it mapped
            (default: 0)
 -g 0/1     if to consider gapped reads(with N in Cigar)
            (default: 0)
 -s 0/1     if strand-specific, won't consider reverse complementary reads if on
            (default: 1)
 -r int     only output genes with more than [int] reads mapped to
            (default: 20)
 --ribo_occup <file>    A File include Transcript CDS Region, Format like this.
                        transID[\t]transLen[\t]1_based_cds_start[\t]1_based_cds_end
                        (default: null, calculate ribosome occupancy if file provided)
 --offset <int>         5' ends of reads will be offset <int> nucleotided to match P-sites location of ribosome
                        (default: 12)
 --rtFile <file>        A File include base 5' pos counts, specific 'stdout' will output into stdout, it's mutex 
                        with --ribo_occup option
                        (default: null)
 --bdFile <file>        A File include base density counts, specific 'stdout' will output into stdout, it's mutex 
                        with --ribo_occup option
                        (default: null)

Limitation:
    1. You Must Sort the SAM File Before Calculate FPKM according to Reads ID so that Same Reads are together
        You Can Sort Your SAM File Like This:
            samtools view -hb test.sam | samtools sort -m 1G -n - | samtools view -h > test.sorted.sam
    
    2. The SAM File Must Include @SQ Tag to Give the information of Length of Transcripts, below is a example:
        @SQ     SN:ENST00000311027_ENSG00000110514      LN:5990
    
    3. The Script Can not be used to Statistic Gene RPKM in SAM file mapped to genome.

Version:
    2017-6-16
    2017-7-4 Add BaseDensity Count Function

""" % (sys.argv[0], )



def init():
    params = {'sam_file': '', 'fpkm_file': '', 'consider_multimap': False, 'consider_gapped_reads': False, 'strand_specific': True, 'minReads': 20, 'ribo_occup': '', 'offset': 12, 'rtFile': None, 'bdFile':None}
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:m:g:s:r:', ['ribo_occup=', 'offset=', 'rtFile=', 'bdFile='])
    for op, value in opts:
        if op == '-h':
            print >>sys.stderr, Usage;
            sys.exit(-1)
        elif op == '-i':
            params['sam_file'] = os.path.abspath(value)
        elif op == '-o':
            params['fpkm_file'] = os.path.abspath(value)
        elif op == '-m':
            params['consider_multimap'] = True if int(value) else False
        elif op == '-g':
            params['consider_gapped_reads'] = True if int(value) else False
        elif op == '-s':
            params['strand_specific'] = True if int(value) else False
        elif op == '-r':
            params['minReads'] = int(value)
        elif op == '--ribo_occup':
            params['ribo_occup'] = os.path.abspath(value)
        elif op == '--offset':
            params['offset'] = int(value)
        elif op == '--rtFile':
            params['rtFile'] = value
        elif op == '--bdFile':
            params['bdFile'] = value
    if not params['sam_file']:
        print >>sys.stderr, Usage;
        sys.exit(-1)
    if params['ribo_occup'] and params['rtFile'] and params['rtFile'] != 'stdout':
        print >>sys.stderr, 'Error: --ribo_occup is mutex with --rtFile option'
        print >>sys.stderr, Usage
        sys.exit(-1)
    return params

def main():
    params = init()
    if params['ribo_occup']:
        transCDSInfo = readTransCDSRegion(params['ribo_occup'])
    else:
        transCDSInfo = None
    GeneCount = calcReadsCount_in_SAM(params['sam_file'], \
        considerMultipleMapReads=params['consider_multimap'], \
        considerStranded=params['strand_specific'], \
        considerGappedReads=params['consider_gapped_reads'], \
        transCDSInfo=transCDSInfo, offset=params['offset'], \
        rtFile=params['rtFile'], bdFile=params['bdFile'])
    if transCDSInfo: calcOccupancyForEachGene(GeneCount, transCDSInfo)
    else: calcFPKMForEachGene(GeneCount)
    outPutFPKM(GeneCount, fileName=params['fpkm_file'], minReads=params['minReads'])
    if params['rtFile']:
        outPutGeneRT(GeneCount, fileName=params['rtFile'], minReads=params['minReads'])
    if params['bdFile']:
        outPutGeneBD(GeneCount, fileName=params['bdFile'], minReads=params['minReads'])

"""

sys.argv = ['bin', '-i', '/Share/home/zhangqf5/tanglei/paris/result/huh7_766_72_1/Aligned.toTranscriptome.out.sam', '--bdFile', 'bd.txt', '-o', './rpkm.txt', '-g', '1', '-m', '1']
params = init()
if params['ribo_occup']:
    transCDSInfo = readTransCDSRegion(params['ribo_occup'])
else:
    transCDSInfo = None
GeneCount = calcReadsCount_in_SAM(params['sam_file'], considerMultipleMapReads=params['consider_multimap'], considerStranded=params['strand_specific'], considerGappedReads=params['consider_gapped_reads'], transCDSInfo=transCDSInfo, offset=params['offset'], rtFile=params['rtFile'], bdFile=params['bdFile'])
if transCDSInfo: calcOccupancyForEachGene(GeneCount, transCDSInfo)
else: calcFPKMForEachGene(GeneCount)

outPutFPKM(GeneCount, fileName=params['fpkm_file'], minReads=params['minReads'])
if params['rtFile']:
    outPutGeneRT(GeneCount, fileName=params['rtFile'], minReads=params['minReads'])

if params['bdFile']:
    outPutGeneBD(GeneCount, fileName=params['bdFile'], minReads=params['minReads'])

"""

def readTransCDSRegion(transCDSFileName):
    """ Read a File include Transcript CDS Region
    The File Should Format like this:
        transID[\t]transLen[\t]1_based_cds_start[\t]1_based_cds_end
    """
    transCDSInfo = {}
    IN = open(transCDSFileName)
    line = IN.readline()
    while line:
        if line.startswith('#') or not line.strip():
            line = IN.readline()
            continue
        (transID, transLen, cds_start, cds_end) = line.strip().split()
        transCDSInfo[transID] = [int(transLen), int(cds_start), int(cds_end)] #{'length': transLen, 'cds_start': cds_start, 'cds_end': cds_end}
        line = IN.readline()
    return transCDSInfo

def FPKM(mapped_reads_in_the_gene, gene_len, total_mapped_reads):
    """ Calculate the FPKM/RPKM
        mapped_reads_in_the_gene / ( (gene_len/pow(10, 3)) * (total_mapped_reads/pow(10, 6)) )
    """
    return 1.0 * mapped_reads_in_the_gene / ( (1.0*gene_len/pow(10, 3)) * (1.0*total_mapped_reads/pow(10, 6)) + 1 )

def addGeneCount_fpkm(GeneCount, Reads, considerMultipleMapReads, considerStranded, considerGappedReads):
    """ Add Read Count To GeneCount Dict
    Reads => [ [flag, geneID, cigar, pos] ]
    """
    # First Filter Some Reads
    idx = 0
    ReadsNums = len(Reads)
    while idx < ReadsNums:
        (flag, geneID, cigar, pos) = Reads[idx]
        # If unmapped, remove it
        if flag & 4:
            del Reads[idx]
            ReadsNums -= 1
        # If stranded, remove reverse complementary reads
        elif considerStranded and flag & 16:
            del Reads[idx]
            ReadsNums -= 1
        # If Gapped-map and we don't consider gapped-reads, skip this read
        elif not considerGappedReads and 'N' in cigar:
            del Reads[idx]
            ReadsNums -= 1
        else:
            idx += 1
    # If Still multi-map and we don't consider multi-map, skip this read
    if len(Reads) > 1 and not considerMultipleMapReads:
        return
    if len(Reads) == 0:
        return
    elif len(Reads) == 1:
        (flag, geneID, cigar, pos) = Reads[0]
        GeneCount[geneID]['uniq'] += 1
    else:
        for (flag, geneID, cigar, pos) in Reads:
            GeneCount[geneID]['multi'] += [len(Reads)]
    if len(Reads) == 1 or ( len(Reads) > 1 and considerMultipleMapReads):
        GeneCount['MAPPED_READS'] += 1

def addGeneCount_occupancy(GeneCount, Reads, geneCDSInfo, offset, considerMultipleMapReads, considerStranded, considerGappedReads, noCDSInfoTrans):
    """ Add Read Count(Only Mapped to CDS Region) To GeneCount Dict
    Reads => [ [flag, geneID, cigar, pos] ]
    """
    def inGeneCDS(read, GeneCount, geneCDSInfo, noCDSInfoTrans):
        (flag, geneID, cigar, pos) = read
        try:
            (gene_length, gene_cds_start, gene_cds_end) = geneCDSInfo[geneID]
            if GeneCount[geneID]['length'] != gene_length:
                if geneID not in noCDSInfoTrans: 
                    print >>sys.stderr, "Warning: %s has different length in SAM(%d) and geneCDSInfo(%d). Skip it!" % (geneID, GeneCount[geneID]['length'], gene_length)
                    noCDSInfoTrans.append(geneID)
                return 0 # error
            true_pos = pos + offset
            if gene_cds_start < true_pos < gene_cds_end:
                return 1 # plus 1
            return 2 # no need to plus, no error
        except KeyError:
            if geneID not in noCDSInfoTrans: 
                print >>sys.stderr, "Warning: %s not in geneCDSInfo. Skip it!" % (geneID, )
                noCDSInfoTrans.append(geneID)
            return 0 # error
    # First Filter Some Reads
    idx = 0
    ReadsNums = len(Reads)
    while idx < ReadsNums:
        (flag, geneID, cigar, pos) = Reads[idx]
        # If unmapped, remove it
        if flag & 4:
            del Reads[idx]
            ReadsNums -= 1
        # If stranded, remove reverse complementary reads
        elif considerStranded and flag & 16:
            del Reads[idx]
            ReadsNums -= 1
        # If Gapped-map and we don't consider gapped-reads, skip this read
        elif not considerGappedReads and 'N' in cigar:
            del Reads[idx]
            ReadsNums -= 1
        else:
            idx += 1
    # If Still multi-map and we don't consider multi-map, skip this read
    if len(Reads) > 1 and not considerMultipleMapReads:
        return
    if len(Reads) == 0:
        return
    elif len(Reads) == 1:
        in_gene_state = inGeneCDS(Reads[0], GeneCount, geneCDSInfo, noCDSInfoTrans)
        if in_gene_state == 1:
            GeneCount[geneID]['uniq'] += 1
        elif in_gene_state == 0: # error
            return
        elif in_gene_state == 2: # no error, in utr, mapped reads should plus 1
            return
    else:
        for read in Reads:
            in_gene_state = inGeneCDS(read, GeneCount, geneCDSInfo, noCDSInfoTrans)
            if in_gene_state == 1:
                GeneCount[geneID]['multi'] += [len(Reads)]
            elif in_gene_state == 0: # error
                return
            elif in_gene_state == 2: # no error, in utr, mapped reads should plus 1
                return
    if len(Reads) == 1 or ( len(Reads) > 1 and considerMultipleMapReads):
        GeneCount['MAPPED_READS'] += 1

def addGeneCount_5pRT(GeneCount, Reads, considerMultipleMapReads, considerStranded, considerGappedReads):
    """ Add 5 prime RT Count(Only Mapped to CDS Region) To GeneCount Dict
    Reads => [ [flag, geneID, cigar, pos] ]
    """
    # First Filter Some Reads
    idx = 0
    ReadsNums = len(Reads)
    while idx < ReadsNums:
        (flag, geneID, cigar, pos) = Reads[idx]
        # If unmapped, remove it
        if flag & 4:
            del Reads[idx]
            ReadsNums -= 1
        # If stranded, remove reverse complementary reads
        elif considerStranded and flag & 16:
            del Reads[idx]
            ReadsNums -= 1
        # If Gapped-map and we don't consider gapped-reads, skip this read
        elif not considerGappedReads and 'N' in cigar:
            del Reads[idx]
            ReadsNums -= 1
        else:
            idx += 1
    # If Still multi-map and we don't consider multi-map, skip this read
    if len(Reads) > 1 and not considerMultipleMapReads:
        return
    if len(Reads) == 0:
        return
    elif len(Reads) == 1:
        (flag, geneID, cigar, pos) = Reads[0]
        GeneCount[geneID]['RT'][pos-1] += 1
    else:
        for (flag, geneID, cigar, pos) in Reads:
            GeneCount[geneID]['RT'][pos-1] += 1
    #if len(Reads) == 1 or ( len(Reads) > 1 and considerMultipleMapReads):
    #    GeneCount['MAPPED_READS'] += 1



def getMatchPos(Cigar, pos):
    """ Get Match region from cigar and pos
    Example:
        getMatchPos(Cigar='2S2M12N3M2M2S', pos=8211398)
        getMatchPos(Cigar='2S2M12N1N3M2M2S', pos=8211398)
    """
    def remove(List, elem):
        return filter(lambda a: a != elem, List)
    matchPos = []
    numVec = [ int(a) for a in remove(re.split('[MIDNSHP=X]', Cigar), '') ]
    symVec = remove(re.split('\d', Cigar), '')
    start = pos
    end = start
    has_value = False
    # remove head
    idx = 0
    while symVec[idx] not in ('M', 'D', '=', 'X'):
        idx += 1
    while idx < len(numVec):
        sym = symVec[idx]
        if sym in ('M', 'D', '=', 'X'):
            end += numVec[idx]
        else:
            if start != end:
                matchPos.append([start, end-1])
            if sym not in ('I', 'P', 'S', 'H'):
                start = end = end + numVec[idx]
            else:
                start = end
        idx += 1
    if start != end:
        matchPos.append([start, end-1])
    return matchPos


def addGeneCount_baseDensity(GeneCount, Reads, considerMultipleMapReads, considerStranded, considerGappedReads):
    """ Add Base Density(Only Mapped to CDS Region) To GeneCount Dict
    Reads => [ [flag, geneID, cigar, pos] ]
    """
    # First Filter Some Reads
    idx = 0
    ReadsNums = len(Reads)
    while idx < ReadsNums:
        (flag, geneID, cigar, pos) = Reads[idx]
        # If unmapped, remove it
        if flag & 4:
            del Reads[idx]
            ReadsNums -= 1
        # If stranded, remove reverse complementary reads
        elif considerStranded and flag & 16:
            del Reads[idx]
            ReadsNums -= 1
        # If Gapped-map and we don't consider gapped-reads, skip this read
        elif not considerGappedReads and 'N' in cigar:
            del Reads[idx]
            ReadsNums -= 1
        else:
            idx += 1
    # If Still multi-map and we don't consider multi-map, skip this read
    if len(Reads) > 1 and not considerMultipleMapReads:
        return
    if len(Reads) == 0:
        return
    elif len(Reads) == 1:
        (flag, geneID, cigar, pos) = Reads[0]
        match_region = getMatchPos(cigar, pos)
        for region in match_region:
            for base_pos in range(region[0], region[1]+1):
                try:
                    GeneCount[geneID]['BD'][base_pos-1] += 1
                except IndexError:
                    print Reads, base_pos
    else:
        for (flag, geneID, cigar, pos) in Reads:
            match_region = getMatchPos(cigar, pos)
            for region in match_region:
                for base_pos in range(region[0], region[1]+1):
                    try:
                        GeneCount[geneID]['BD'][base_pos-1] += 1
                    except IndexError:
                        print Reads, base_pos


def calcReadsCount_in_SAM(samFileName, considerMultipleMapReads=False, considerStranded=True, considerGappedReads=False, transCDSInfo=None, offset=12, rtFile=None, bdFile=None):
    """ Statistic Reads Count in Each Gene
    """
    # Some Neccessary Variable
    GeneCount = {}
    GeneCount['MAPPED_READS'] = 0
    SAM = open(samFileName)
    line = SAM.readline()
    last_reads_id = ''
    last_reads_list = []
    start_geneCount = False
    noCDSInfoTrans = []
    while line:
        if line.startswith('@'):
            if line.startswith('@SQ'):
                (tag, geneID, length) = line.strip().split()
                geneID = geneID.split(':')[1]
                length = int(length.split(':')[1])
                GeneCount[geneID] = {'uniq': 0, 'multi': [], 'length': length}
                if rtFile:
                    GeneCount[geneID]['RT'] = [0]*length
                if bdFile:
                    GeneCount[geneID]['BD'] = [0]*length
            line = SAM.readline()
            continue
        if not start_geneCount:
            print >>sys.stderr, '\tStart to Statistic Gene Count...'
            start_geneCount = True
        SAMLineItems = line.strip().split('\t')
        (ReadName, Flag, Chr, Pos, MapQ, Cigar, Rnext, Pnext, Tlen, Seq, Quality) = SAMLineItems[:11]
        if ReadName != last_reads_id:
            if last_reads_id == '':
                last_reads_id = ReadName
                last_reads_list = [ [int(Flag), Chr, Cigar, int(Pos)] ]
            else:
                if not transCDSInfo:
                    addGeneCount_fpkm(GeneCount, last_reads_list, considerMultipleMapReads, considerStranded, considerGappedReads)
                    if rtFile:
                        addGeneCount_5pRT(GeneCount, last_reads_list, considerMultipleMapReads, considerStranded, considerGappedReads)
                    if bdFile:
                        addGeneCount_baseDensity(GeneCount, last_reads_list, considerMultipleMapReads, considerStranded, considerGappedReads)
                else:
                    addGeneCount_occupancy(GeneCount, last_reads_list, transCDSInfo, offset, considerMultipleMapReads, considerStranded, considerGappedReads, noCDSInfoTrans)
                last_reads_id = ReadName
                last_reads_list = [ [int(Flag), Chr, Cigar, int(Pos)] ]
        else:
            last_reads_list.append([int(Flag), Chr, Cigar, int(Pos)])
        line = SAM.readline()
    if not transCDSInfo:
        addGeneCount_fpkm(GeneCount, last_reads_list, considerMultipleMapReads, considerStranded, considerGappedReads)
        if rtFile:
            addGeneCount_5pRT(GeneCount, last_reads_list, considerMultipleMapReads, considerStranded, considerGappedReads)
        if bdFile:
            addGeneCount_baseDensity(GeneCount, last_reads_list, considerMultipleMapReads, considerStranded, considerGappedReads)
    else:
        addGeneCount_occupancy(GeneCount, last_reads_list, transCDSInfo, offset, considerMultipleMapReads, considerStranded, considerGappedReads, noCDSInfoTrans)
    return GeneCount

def calcFPKMForEachGene(GeneCount):
    for GeneID in GeneCount:
        if GeneID != 'MAPPED_READS':
            mapped_reads = GeneCount[GeneID]['uniq'] + sum([ 1.0/mapped_locs for mapped_locs in GeneCount[GeneID]['multi'] ])
            if GeneCount[GeneID]['length'] != 0:
                GeneCount[GeneID]['FPKM'] = FPKM(mapped_reads_in_the_gene=mapped_reads, gene_len=GeneCount[GeneID]['length'], total_mapped_reads=GeneCount['MAPPED_READS'])

def calcOccupancyForEachGene(GeneCount, geneCDSInfo):
    for GeneID in GeneCount:
        if GeneID != 'MAPPED_READS':
            mapped_reads = GeneCount[GeneID]['uniq'] + sum([ 1.0/mapped_locs for mapped_locs in GeneCount[GeneID]['multi'] ])
            try:
                (gene_length, gene_cds_start, gene_cds_end) = geneCDSInfo[GeneID]
                cds_length = gene_cds_end - gene_cds_start + 1
                if GeneCount[GeneID]['length'] != 0:
                    GeneCount[GeneID]['FPKM'] = FPKM(mapped_reads_in_the_gene=mapped_reads, gene_len=cds_length, total_mapped_reads=GeneCount['MAPPED_READS'])
            except KeyError:
                continue

def outPutFPKM(GeneCount, fileName='', minReads=20):
    if fileName:
        OUT = open(fileName, 'w')
    else:
        OUT = sys.stdout
    print >>OUT, "# Total Mapped Reads: %d" % ( GeneCount["MAPPED_READS"], )
    print >>OUT, "#%s\t%s\t%s\t%s\t%s" % ( 'id_name', 'length', 'uniq_map', 'multi_map', 'fpkm' )
    for GeneID in sorted(GeneCount.keys()):
        if GeneID != 'MAPPED_READS':
            if GeneCount[GeneID]['uniq']+len(GeneCount[GeneID]['multi']) >= minReads and 'FPKM' in GeneCount[GeneID]:
                gene_rpkm_info = "%s\t%d\t%d\t%d\t%.3f" % (GeneID, GeneCount[GeneID]['length'], GeneCount[GeneID]['uniq'], len(GeneCount[GeneID]['multi']), GeneCount[GeneID]['FPKM'])
                print >>OUT, gene_rpkm_info
    OUT.close()

def outPutGeneRT(GeneCount, fileName='stdout', minReads=20):
    if fileName != 'stdout':
        OUT = open(fileName, 'w')
    else:
        OUT = sys.stdout
    print >>OUT, "# Total Mapped Reads: %d" % ( GeneCount["MAPPED_READS"], )
    print >>OUT, "#%s\t%s\t%s\t%s\t%s" % ( 'id_name', 'length', 'uniq_map', 'multi_map', 'RT...' )
    for GeneID in sorted(GeneCount.keys()):
        if GeneID != 'MAPPED_READS':
            if sum(GeneCount[GeneID]['RT']) >= minReads:
                print >>OUT, "%s\t%d\t%d\t%d\t%s" % (GeneID, GeneCount[GeneID]['length'], GeneCount[GeneID]['uniq'], len(GeneCount[GeneID]['multi']), '\t'.join([ str(rt) for rt in GeneCount[GeneID]['RT'] ]))
    OUT.close()

def outPutGeneBD(GeneCount, fileName='stdout', minReads=20):
    if fileName != 'stdout':
        OUT = open(fileName, 'w')
    else:
        OUT = sys.stdout
    print >>OUT, "# Total Mapped Reads: %d" % ( GeneCount["MAPPED_READS"], )
    print >>OUT, "#%s\t%s\t%s\t%s\t%s" % ( 'id_name', 'length', 'uniq_map', 'multi_map', 'BD...' )
    for GeneID in sorted(GeneCount.keys()):
        if GeneID != 'MAPPED_READS':
            if GeneCount[GeneID]['uniq']+len(GeneCount[GeneID]['multi']) >= minReads:
                print >>OUT, "%s\t%d\t%d\t%d\t%s" % (GeneID, GeneCount[GeneID]['length'], GeneCount[GeneID]['uniq'], len(GeneCount[GeneID]['multi']), '\t'.join([ str(rt) for rt in GeneCount[GeneID]['BD'] ]))
    OUT.close()





main();





