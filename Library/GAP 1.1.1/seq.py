#-*- coding:utf-8 -*-

"""
这个模块主要是一些与序列操作有关的函数

"""

import re

def biSearch(item, Set):
    """ 二分查找函数，Set必须是已经从小到大排序的
    使用方法：
        biSearch( 3, [-1,0,1,3,4] ) # True
    """
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


def reverse_comp(seq):
        """逆序互补序列"""
        reverseDict = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N',
                   'a':'t', 'c':'g', 't':'a', 'g':'c'}
        return ''.join(map(lambda x: reverseDict[x], list(seq[::-1])))


def cutSeq(rawSeq, lineLen=60):
    """把完整连续的序列切割成分割的片段
    测试：
        print cutSeq('ATCG', lineLen=60)
        print cutSeq('ATCG'*60, lineLen=60)
        print cutSeq(('ATCG'*60)[:-1], lineLen=60)
        print cutSeq('ATCG', lineLen=1)
    """
    cut_seq = ''
    idx = 0
    while idx < len(rawSeq):
        cut_seq += rawSeq[idx:idx+lineLen]+'\n'
        idx += lineLen
    return cut_seq[:-1]

class anno_Methods(object):
    """
    与基因组注释相关的方法
    """
    @staticmethod
    def gene_type(raw_type):
        """ 常见的基因类型分类方法:
            Convert Raw Gene Type to Our Defined Gene Type
        """
        valid_gene_type = ('pseudogene', 'snoRNA', 'snRNA', 'miRNA', 'misc_RNA', 'rRNA')
        #lncRNA_class = ('antisense','lincRNA','processed_transcript','sense_intronic','TEC','sense_overlapping')
        lncRNA_class = ('3prime_overlapping_ncrna','antisense','lincRNA','non_coding','sense_intronic','sense_overlapping','processed_transcript')
        if raw_type in valid_gene_type: return raw_type;
        if re.match('.*pseudogene',raw_type): return 'pseudogene';
        if raw_type == 'protein_coding': return 'mRNA';
        if raw_type in lncRNA_class: return 'lncRNA';
        return 'other'
    @staticmethod
    def loadGTFBed(annoBedFile):
        """ 加载一个基因组处理过的注释文件: oganism.transCoordinate.bed
        使用方法：
            hg38Anno = anno_Methods.loadGTFBed("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Homo_sapiens.transCoordinate.bed")
            mm10Anno = anno_Methods.loadGTFBed("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Mus_musculus.transCoordinate.bed")

        返回：
             dict: trans_id ==> gene_id trans_len utr ....
        """
        def parseUTR(utr_str, trans_len):
            """
                解析UTR字符串
                没有5'UTR的转录本 utr_5_start=utr_5_end=0
                没有3'UTR的转录本 utr_3_start=utr_3_end=trans_len+1
            """
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
        H = open(annoBedFile)
        annoBed = dict()
        line = H.readline()
        while line:
            arr = line.strip().split()
            trans_id = arr[0]
            gene_id = arr[1]
            gene_type = arr[2]
            trans_len = int(arr[3])
            utr = ''
            if len(arr) == 6:
                utr = arr[5]
            utr_coor = parseUTR(utr, trans_len)
            annoBed[trans_id] = {}
            annoBed[trans_id]['gene_id'] = gene_id
            annoBed[trans_id]['gene_type'] = gene_type
            annoBed[trans_id]['trans_len'] = trans_len
            for it in utr_coor:
                annoBed[trans_id][it] = utr_coor[it]
            line = H.readline()
        print 'Features For Transcript: '
        print '  *'+'\n  *'.join( annoBed[annoBed.keys()[0]].keys() )
        return annoBed


import pysam
class seqClass(object):
    """获取大基因组的序列
    使用方法：
        hg38_seq = seq.seqClass(hg38_genome_file)
        hg38_seq.fetch('chr1', 56554, 56554+2370, '+') 
        hg38_seq.fetch('chr1', 56554, 56554+2370, '-') complementary reverse sequence
    常用的基因组文件：
        /Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.fa
        /Share/home/zhangqf8/lipan/DYNAMIC/GTF/mm10.fa
        /Share/home/zhangqf8/lipan/DYNAMIC/GTF/transcriptome.fa
        /Share/home/zhangqf8/lipan/DYNAMIC/GTF/mouse_transcriptome.fa   
    """
    def __init__(self, fa_fileName):
        self.fileName = fa_fileName
        self.genome = pysam.Fastafile(fa_fileName)
        print "seqClass: 输入坐标前闭后开，0-based"
    def fetch(self, Chr, Start, End, Strand="+"):
        """获取序列"""
        if Strand == '+':
            return self.genome.fetch(Chr, Start, End)
        if Strand == '-':
            return reverse_comp(self.genome.fetch(Chr, Start, End))
    def has(self, Chr):
        return Chr in self.genome.references


from pyliftover import LiftOver
class liftOverClass(LiftOver):
    """
    基因组版本转换：
        hg19_2_hg38 = liftOverClass(from_db='hg19', to_db='hg38')
        hg19_2_hg38.convert_coor(Chr='chr1', Pos=109303388, Strand='+')
        # [('chr1', 108760766, '+', 20851231461)]
    """
    def __init__(self, from_db, to_db):
        LiftOver.__init__(self, from_db=from_db, to_db=to_db)
    def convert_coor(self, Chr, Pos, Strand='+'):
        return self.convert_coordinate(chromosome=Chr, position=Pos, strand=Strand)


enviroment = {
    "maxChr": 999999999999, 
    "matchEnergy": {
        'AA': -999, 'AC': -999, 'AG': -999, 'AT': 1,
        'CA': -999, 'CC': -999, 'CG': 1, 'CT': -999,
        'GA': -999, 'GC': 1, 'GG': -999, 'GT': 1,
        'TA': 1, 'TC': -999, 'TG': 1, 'TT': -999},
    "gapPenalty": -1
}

def localAlignment(seq1, seq2):
    """
    local Alignment two sequence
    test:
        seq1 = 'CAAATACTAGGGAAAGACTAGGAGGATGAGCCAGGGTTGCTACTA'
        seq2 = 'TTTGTTGCTGCTTCCCT'#[::-1]
        ( maxScore, alignment, seq1_s, seq1_e, seq2_s, seq2_e ) = localAlignment(seq1, seq2[::-1]);
        print seq1[seq1_s-1:seq1_e]
        print seq2[seq2_s-1:seq2_e]
        print alignment
    """
    # initialzation
    seq1_len = len(seq1)
    seq2_len = len(seq2)
    matrix = []
    for i in range(seq2_len):
        matrix.append( [] )
        for j in range(seq1_len):
            matrix[i].append( {} )
            matrix[i][j]['score'] = 0
            matrix[i][j]['pointer'] = 'none'
    # fill
    maxScore = 0; maxI = 0; maxJ = 0;
    for i in range(1, seq2_len):
        for j in range(1, seq1_len):
            # calculate match score
            pair = seq1[j-1] + seq2[i-1]
            diagonalScore = matrix[i-1][j-1]['score'] + enviroment['matchEnergy'][pair]
            # calculate gap score
            upScore = matrix[i-1][j]['score'] + enviroment['gapPenalty']
            leftScore = matrix[i][j-1]['score'] + enviroment['gapPenalty']
            if diagonalScore <= 0 and upScore <= 0 and leftScore <= 0:
                matrix[i][j]['score'] = 0
                matrix[i][j]['pointer'] = 'none'
                continue
            # choose best score
            if diagonalScore >= upScore:
                if diagonalScore >= leftScore:
                    matrix[i][j]['score'] = diagonalScore
                    matrix[i][j]['pointer'] = 'diagonal'
                else:
                    matrix[i][j]['score'] = leftScore
                    matrix[i][j]['pointer'] = 'left'
            else:
                if upScore >= leftScore:
                    matrix[i][j]['score'] = upScore
                    matrix[i][j]['pointer'] = 'up'
                else:
                    matrix[i][j]['score'] = leftScore
                    matrix[i][j]['pointer'] = 'left'
            if matrix[i][j]['score'] > maxScore:
                maxI = i
                maxJ = j
                maxScore = matrix[i][j]['score']
    # trace-back
    align1 = align2 = ''
    i = maxI; j = maxJ;
    while 1:
        if matrix[i][j]['pointer'] == 'none': break
        if matrix[i][j]['pointer'] == 'diagonal':
            align1 += seq1[j-1]
            align2 += seq2[i-1]
            i -= 1; j -= 1;
        elif matrix[i][j]['pointer'] == 'left':
            align1 += seq1[j-1]
            align2 += '-'
            j -= 1
        elif matrix[i][j]['pointer'] == 'up':
            align1 += '-'
            align2 += seq2[i-1]
            i -= 1
    align = align1[::-1] + ':' + align2
    return maxScore, align, j+1, maxJ, i+1, maxI


