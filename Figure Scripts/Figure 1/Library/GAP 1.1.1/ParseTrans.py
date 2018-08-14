#-*- coding:utf-8 -*-

""" Get transcript Feature from genomeCoorBedFile

Help:

    from ParseTrans import *

    Parser = ParseTransClass(genomeCoorBedFile = "hg19.genomeCoor.bed")
    or:
    Parser = ParseTransClass(genomeCoorBedFile = "hg19.genomeCoor.bed", seqFileName="hg19_transcriptome.fa")

    Parser.getTransFeature('ENST00000456328.2_1')
    Parser.getTransFeatureSeq(transID='ENST00000456328.2_1', feature='start_codon', source='GENCODE')

    Parser.addSeq(seqFileName = "hg19_transcriptome.fa")

    GeneParser = Parser.getGeneParser()
    
    geneIntron = Parser.getGeneIntron('ENSG00000223972.5_2')['ENST00000456328.2_1']
    geneExon = Parser.getGeneExon('ENSG00000223972.5_2')['ENST00000456328.2_1']

    Parser.getTransFeatureSeq(transID='XM_005244749', feature='stop_codon', source='NCBI')

Coordinate Covert:

    Gene => Genome/Trans
        Parser.geneCoor2genomeCoor("ENSG00000227232.5_2", 10)
        Parser.geneCoor2transCoor("ENSG00000227232.5_2", 1, 11)

    Trans => Genome/Gene
        Parser.transCoor2genomeCoor("ENST00000456328.2_1", 1)
        Parser.transCoor2geneCoor("ENST00000456328.2_1", 1, 11)

    Genome => Trans/Gene
        Parser.genomeCoor2geneCoor("chr1", 11869, 11869+10, '+')
        Parser.genomeCoor2transCoor("chr1", 11869, 11869+10, '+')

feature：must be one of
        *   whole
        *   5utr(utr5, 5_utr, utr_5)
        *   3utr(utr3, 3_utr, utr_3)
        *   cds(CDS)
        *   start_codon(startcodon, startCodon)
        *   stop_codon(stopcodon, stopCodon)

Attention:
    NCBI/Gencode Stop Codon is different:
        NCBI Stop Codon is the last three nucleotides in CDS
        Gencode Stop Codon is the first three nucleotides in 3_utr

"""

from NCBI_Genome import *
import CoorFunc
import sys

___ParseTrans = {
    'first_biuldGeneParser': 1
}

def showBeautifulExample(exmaple, title):
    print >>sys.stderr, '# =-=-=-=-=-=-=-=-=-=-=-=-=-%s=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=' % (title, )
    Pattern = "\t%15s:\t%-25s"
    for attr in sorted(exmaple.keys()):
        print >>sys.stderr, Pattern % (attr, exmaple[attr])
    print >>sys.stderr, '# =-=-=-=-=-=-=-=-=-=-=-=-=-%s=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=' % (len(title)/2 * "=-", )

def format_Exon_UTR_str(raw_str, strand):
    regions_list = [ (int(region.split('-')[0]), int(region.split('-')[1])) for region in raw_str.split(',') ]
    if strand == '+':
        regions_list.sort()
    else:
        regions_list.sort(reverse=True)
    return ",".join([str(region[0])+"-"+str(region[1]) for region in regions_list])

def getTransIParser(genomeCoorBedFile, showAttr=True, verbose=True):
    """ 获得转录本区段解析句柄
    """
    def parseUTR(utr_str, trans_len):
        # + 解析UTR字符串
        # + 没有5'UTR的转录本 utr_5_start=utr_5_end=0
        # + 没有3'UTR的转录本 utr_3_start=utr_3_end=trans_len+1
        utr_5_start = utr_5_end = 0 
        utr_3_start = utr_3_end = trans_len + 1
        cds_start = cds_end = 0
        if utr_str:
            utr_bak = utr_str.split(',')
            utr = []
            for utr_item in utr_bak:
                utr.append([int(i) for i in utr_item.split('-')])
            if utr[0][0] == 1:
                # + 有5'UTR
                utr_5_start = 1
                for i in range(len(utr)-1):
                    if utr[i][1] != utr[i+1][0] - 1:
                        utr_5_end = utr[i][1]
                        utr_3_start = utr[i+1][0]
                        utr_3_end = utr[-1][1]
                        break
                if utr_5_end == 0: 
                    utr_5_end = utr[-1][1]    
            else:
                # + 没有5'UTR，只有3'UTR
                utr_3_start = utr[0][0]
                utr_3_end = utr[-1][1]
        cds_start = utr_5_end + 1
        cds_end = utr_3_start - 1
        elem_coor = {'utr_5_start':utr_5_start, 'utr_5_end':utr_5_end, 'utr_3_start':utr_3_start, 'utr_3_end':utr_3_end, 'cds_start':cds_start, 'cds_end':cds_end}
        return elem_coor
    # + 解析函数
    IN = open(genomeCoorBedFile)
    transIParser = dict();
    line = IN.readline()
    while line:
        arr = line.strip().split('\t')
        (Chr, Start, End, Strand, GeneName_GeneID, TransID, GeneType, Exon) = arr[:8]
        Exon = format_Exon_UTR_str(Exon, Strand)
        try:
            (GeneName, GeneID) = (GeneName_GeneID.split('=')[0], GeneName_GeneID.split('=')[1])
        except:
            print >>sys.stderr, GeneName_GeneID
            return
        exonList = NCBI_Genome_Class.norm_exons(Exon)[0]
        TransLen = int(exonList[-1][1])
        utrString = ''
        if len(arr) == 9:
            arr[8] = format_Exon_UTR_str(arr[8], Strand)
            utrList = NCBI_Genome_Class.norm_utr(Exon, arr[8], Strand)[1]
            utrString = ','.join([str(item[0])+'-'+str(item[1]) for item in utrList])
        utr_coor = parseUTR(utrString, TransLen)
        transIParser[TransID] = {}
        transIParser[TransID]['gene_id'] = GeneID
        transIParser[TransID]['gene_name'] = GeneName
        transIParser[TransID]['gene_type'] = GeneType
        transIParser[TransID]['trans_len'] = TransLen
        transIParser[TransID]['chr'] = Chr
        transIParser[TransID]['start'] = int(Start)
        transIParser[TransID]['end'] = int(End)
        transIParser[TransID]['strand'] = Strand
        transIParser[TransID]['exon_str'] = Exon
        for it in utr_coor:
            transIParser[TransID][it] = utr_coor[it]
        line = IN.readline()
    # show attribuets
    if showAttr and len(transIParser) != 0:
        exmaple = transIParser[transIParser.keys()[0]]
        showBeautifulExample(exmaple, title='transIParser Object Example')
    return transIParser



def readSeq(seqFileName):
    """ 读取Fasta文件
    测试：
        1. Gencode
            hg38_Seq = readSeq(os.environ.get('HOME') + '/lipan/DYNAMIC/GTF/Gencode/gencode/hg38_transcriptome.fa')
        2. NCBI
            hg38_Seq = readSeq(os.environ.get('HOME') + '/lipan/DYNAMIC/GTF/ncbi/ncbi/hg38_transcriptome.fa')
    """
    Seq = {}
    IN = open(seqFileName)
    line = IN.readline()
    cur_trans = ''
    while line:
        if line[0] == '>':
            cur_trans = line[1:].split()[0]
            Seq[ cur_trans ] = ''
        else:
            Seq[ cur_trans ] += line.strip()
        line = IN.readline()
    return Seq

def transSeq(Seq, transIParser, transID, feature, source):
    from seq import biSearch
    """
    feature is one of whole/5utr(utr5, 5_utr, utr_5)/3utr(utr3, 3_utr, utr_3)/cds(CDS)/start_codon(startcodon, startCodon)/stop_codon(stopcodon, stopCodon)
    测试：
        transSeq(Seq=hg38_Seq, transIParser=hg38_transIParser, transID=transID, feature='stop_codon')
    """
    try:
        trans_seq = Seq[transID]
    except KeyError:
        print >>sys.stderr, 'Warning: no this transID: %s in Seq' % (transID, )
        return ''
    try:
        parser = transIParser[transID]
    except KeyError:
        print >>sys.stderr, 'Warning: no this transID: %s in transIParser' % (transID, )
        return ''
    if feature in ('utr5', '5utr', '5_utr', 'utr_5'):
        start = parser['utr_5_start'] - 1
        end = parser['utr_5_end']
    elif feature in ('utr3', '3utr', '3_utr', 'utr_3'):
        start = parser['utr_3_start'] + 2
        end = parser['utr_3_end']
    elif feature in ('cds', 'CDS'):
        start = parser['cds_start'] - 1
        end = parser['cds_end']
    elif feature in ('start_codon', 'startcodon', 'startCodon'):
        start = parser['cds_start'] - 1
        end = parser['cds_start'] + 2
    elif feature in ('stop_codon', 'stopcodon', 'stopCodon'):
        if source in ('NCBI', 'ncbi'):
            start = parser['cds_end'] - 3
            end = parser['cds_end']
        elif source in ('Gencode', 'GENCODE', 'Ensembl', 'ENSEMBL'):
            start = parser['cds_end']
            end = parser['cds_end'] + 3
        else:
            print >>sys.stderr, 'source must be one of NCBI/ncbi/Gencode/GENCODE/Emsembl'
    elif feature == 'whole':
        start = 0
        end = parser['utr_3_end']
    else:
        print >>sys.stderr, 'Warning: feature is one of whole/5utr(utr5, 5_utr, utr_5)/3utr(utr3, 3_utr, utr_3)/cds(CDS)/start_codon(startcodon, startCodon)/stop_codon(stopcodon, stopCodon)'
        return ''
    if parser['gene_type'] not in ('mRNA', 'protein_coding'):
        print >>sys.stderr, 'Warining: %s is a %s RNA' % (transID, parser['gene_type'])
    return trans_seq[start:end]

def biuldGeneParser(transIParser, showAttr=True):
    """ 构建一个基因ID => 转录本ID的字典
    测试：
        geneInfo = biuldGeneParser(hg38_transIParser)
    注意：One Gene may have multiple gene types
    """
    if ___ParseTrans['first_biuldGeneParser']:
        print >>sys.stderr, 'Warning: Gene Coordinate System: [1_based_start, 1_based_end]'
        ___ParseTrans['first_biuldGeneParser'] = 0
    geneInfo = {}
    for transID in transIParser:
        geneID = transIParser[transID]['gene_id']
        Chr = transIParser[transID]['chr']
        Start = transIParser[transID]['start']
        End = transIParser[transID]['end']
        Strand = transIParser[transID]['strand']
        geneName = transIParser[transID]['gene_name']
        geneType = transIParser[transID]['gene_type']
        try:
            geneInfo[geneID]['transcript'] += [transID] # = geneInfo[geneID].get(geneID, []) + [transID]
            if geneType not in geneInfo[geneID]['gene_type']: geneInfo[geneID]['gene_type'].append(geneType)
            if Start < geneInfo[geneID]['start']: 
                geneInfo[geneID]['start'] = Start
                geneInfo[geneID]['length'] = geneInfo[geneID]['end'] - geneInfo[geneID]['start'] + 1
            if End > geneInfo[geneID]['end']: 
                geneInfo[geneID]['end'] = End
                geneInfo[geneID]['length'] = geneInfo[geneID]['end'] - geneInfo[geneID]['start'] + 1
        except KeyError:
            geneInfo[geneID] = {}
            geneInfo[geneID]['transcript'] = [ transID ]
            geneInfo[geneID]['chr'] = Chr
            geneInfo[geneID]['start'] = Start
            geneInfo[geneID]['end'] = End
            geneInfo[geneID]['length'] = End - Start + 1
            geneInfo[geneID]['strand'] = Strand
            geneInfo[geneID]['gene_name'] = geneName
            geneInfo[geneID]['gene_type'] = [ geneType ]
    # show attribuets
    if showAttr and len(geneInfo) != 0:
        exmaple = geneInfo[geneInfo.keys()[0]]
        showBeautifulExample(exmaple, title='geneInfo Object Example')
    return geneInfo

def GeneIntron(transIParser, geneParser, geneID):
    class NoExonError(Exception):
        pass
    def Intron_From_Exom(exon_str, start, end, strand):
        """ Get Intron Area From Exon String
        Example:
            exon_str = "489361-489710,485040-485208,476887-476945,476738-476882"
            start = 476738
            end = 489710
            Intron_From_Exom(exon_str, start, end, strand="-")
        """
        import copy
        genomeCoorIntron = []
        exonList = [ [int(exon_region.split('-')[0]), int(exon_region.split('-')[1])] for exon_region in exon_str.split(',') ]
        if len(exonList) == 0: raise NoExonError
        exonList.sort(key=lambda x: x[0], reverse=False)
        if start < exonList[0][0]:
            genomeCoorIntron.append( [start, exonList[0][0]-1] )
        for idx in range(len(exonList)-1):
            exon_region = exonList[idx]
            genomeCoorIntron.append( [exonList[idx][1]+1, exonList[idx+1][0]-1] )
        if end > exonList[-1][1]:
            genomeCoorIntron.append( [exonList[-1][1]+1, end] )
        GeneCoorIntron = copy.deepcopy(genomeCoorIntron)
        if strand == '+':
            for idx in range(len(GeneCoorIntron)):
                GeneCoorIntron[idx][0] = GeneCoorIntron[idx][0] - start + 1
                GeneCoorIntron[idx][1] = GeneCoorIntron[idx][1] - start + 1
        elif strand == '-':
            for idx in range(len(GeneCoorIntron)):
                GeneCoorIntron[idx][0] = end - GeneCoorIntron[idx][0] + 1
                GeneCoorIntron[idx][1] = end - GeneCoorIntron[idx][1] + 1
                GeneCoorIntron[idx][0], GeneCoorIntron[idx][1] = GeneCoorIntron[idx][1], GeneCoorIntron[idx][0]
        genomeCoorIntron.sort(key=lambda x: x[0], reverse=True if strand == '-' else False)
        GeneCoorIntron.sort(key=lambda x: x[0])
        return genomeCoorIntron, GeneCoorIntron
    gene_start = geneParser[geneID]['start']
    gene_end = geneParser[geneID]['end']
    strand = geneParser[geneID]['strand']
    IsoformIntrons = {}
    for transID in geneParser[geneID]['transcript']:
        exon_str = transIParser[transID]['exon_str']
        trans_start = transIParser[transID]['start']
        trans_end = transIParser[transID]['end']
        try:
            #IsoformIntrons[transID] = Intron_From_Exom(exon_str, gene_start, gene_end, strand)
            IsoformIntrons[transID] = Intron_From_Exom(exon_str, trans_start, trans_end, strand)
        except NoExonError:
            print >>sys.stderr, "Warning: %s have no Exon_str, Skip it" % (transID, )
    return IsoformIntrons


def GeneExon(transIParser, geneParser, geneID):
    class NoExonError(Exception):
        pass
    def Exom_from_str(exon_str, start, end, strand):
        """ Get Intron Area From Exon String
        Example:
            exon_str = "489361-489710,485040-485208,476887-476945,476738-476882"
            start = 476738
            end = 489710
            Intron_From_Exom(exon_str, start, end, strand="-")
        """
        import copy
        genomeCoorIntron = []
        exonList = [ [int(exon_region.split('-')[0]), int(exon_region.split('-')[1])] for exon_region in exon_str.split(',') ]
        if len(exonList) == 0: raise NoExonError
        exonList.sort(key=lambda x: x[0], reverse=False)
        for idx in range(len(exonList)):
            if strand == '+':
                gene_start, gene_end = exonList[idx][0]-start+1, exonList[idx][1]-start+1
            elif strand == '-':
                gene_end, gene_start = end-exonList[idx][0]+1, end-exonList[idx][1]+1
            genomeCoorIntron.append([gene_start, gene_end])
        genomeCoorIntron.sort(key=lambda x:x[0])
        return exonList, genomeCoorIntron
    gene_start = geneParser[geneID]['start']
    gene_end = geneParser[geneID]['end']
    strand = geneParser[geneID]['strand']
    IsoformExons = {}
    for transID in geneParser[geneID]['transcript']:
        exon_str = transIParser[transID]['exon_str']
        try:
            IsoformExons[transID] = Exom_from_str(exon_str, gene_start, gene_end, strand)
        except NoExonError:
            print >>sys.stderr,"Warning: %s have no Exon_str, Skip it" % (transID, )
    return IsoformExons


class ParseTransClass(object):
    """ 解析转录本元件 """
    def __init__(self, genomeCoorBedFile, seqFileName='', showAttr=True):
        """ 
        genomeCoorBedFile：hg38.genomeCoor.bed / mm10.genomeCoor.bed
        seqFileName[option]: hg38_transcriptome.fa / mm10_transcriptome.fa
        """
        self.seqFileName = seqFileName
        self.genomeCoorBedFile = genomeCoorBedFile
        self.TransIParser = getTransIParser(self.genomeCoorBedFile, showAttr=showAttr)
        if self.seqFileName != '':
            self.Seq = readSeq(self.seqFileName)
        print >>sys.stderr, 'Warning: Trans/Gene Coorninate System: [1_based_start, 1_based_end]'
    def getTransFeatureSeq(self, transID, feature, source):
        """ get transcription seq by this function """
        if self.seqFileName == "":
            print >>sys.stderr, 'Not define Seq! Using addSeq() to add Seq'
            return ''
        return transSeq(self.Seq, self.TransIParser, transID, feature, source)
    def addSeq(self, seqFileName):
        """ If not define seqFileName in construct, add it using this function """
        self.seqFileName = seqFileName
        self.Seq = readSeq(self.seqFileName)
    def getTransFeature(self, transID, showAttr=False, verbose=True):
        """ return the feature of a transcript. Mind that 1-based coordinate system """
        try:
            parser = self.TransIParser[transID]
        except KeyError:
            if verbose: print >>sys.stderr, 'Warning: no this transID: %s in transIParser' % (transID, )
            raise KeyError
        if showAttr:
            showBeautifulExample(parser, title='transID Info Show')
        return parser
    def getGeneParser(self, showAttr=True):
        if hasattr(self, 'GeneParser'):
            pass
        else:
            self.GeneParser = biuldGeneParser(self.TransIParser, showAttr=showAttr)
        return self.GeneParser
    def getGeneInfo(self, showAttr=True):
        return self.getGeneParser(showAttr)
    def getGeneIntron(self, geneID, showAttr=False):
        if hasattr(self, 'GeneParser'):
            pass
        else:
            self.GeneParser = biuldGeneParser(self.TransIParser, showAttr=showAttr)
        return GeneIntron(self.TransIParser, self.GeneParser, geneID)
    def getGeneExon(self, geneID, showAttr=False):
        if hasattr(self, 'GeneParser'):
            pass
        else:
            self.GeneParser = biuldGeneParser(self.TransIParser, showAttr=showAttr)
        return GeneExon(self.TransIParser, self.GeneParser, geneID)
    def getGeneBin(self):
        if hasattr(self, 'GeneBin'):
            pass
        else:
            self.GeneBin = CoorFunc.binize_Gene(self, bw=100000)
        return self.GeneBin
    def getTransBin(self):
        if hasattr(self, 'TransBin'):
            pass
        else:
            self.TransBin = CoorFunc.binize_Trans(self, bw=100000)
        return self.TransBin
    def geneCoor2genomeCoor(self, GeneID, Pos):
        return CoorFunc.geneCoor2genomeCoor(self, GeneID, Pos)
    def genomeCoor2geneCoor(self, Chr, Start, End, Strand):
        GeneBin = self.getGeneBin()
        return CoorFunc.genomeCoor2geneCoor(GeneBin, self, Chr, Start, End, Strand)
    def genomeCoor2transCoor(self, Chr, Start, End, Strand):
        TransBin = self.getTransBin()
        return CoorFunc.genomeCoor2transCoor(TransBin, self, Chr, Start, End, Strand)
    def transCoor2genomeCoor(self, TransID, Pos):
        return CoorFunc.transCoor2genomeCoor(self, TransID, Pos)
    def geneCoor2transCoor(self, GeneID, Start, End):
        TransBin = self.getTransBin()
        return CoorFunc.geneCoor2transCoor(TransBin, self, GeneID, Start, End)
    def transCoor2geneCoor(self, TransID, Start, End):
        GeneBin = self.getGeneBin()
        return CoorFunc.transCoor2geneCoor(GeneBin, self, TransID, Start, End)













