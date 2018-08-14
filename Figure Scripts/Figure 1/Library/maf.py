#-*- coding:utf-8 -*-

"""
这个模块主要与 .maf 文件的相关操作有关，下面是目前实现的最重要的几个方法：

        # maf 转成 tMaf
        mafFile = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/hg38.mm10.synNet.maf'
        hg38_GTF_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Homo_sapiens.gtf'
        mm10_GTF_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Mus_musculus.GRCm38.81.gtf'
        maf = mafClass(mafFile)
        maf.gmaf2tMaf(hg38_GTF_file, mm10_GTF_file, outFileName='/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/chrY.tmaf', Chr_to_processed='chrY')
        或者分别同时处理多个染色体：
        for chrSym in ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']:
            maf.gmaf2tMaf(hg38_GTF_file, mm10_GTF_file, outFileName='/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/chrs/'+chrSym+'.tmaf', Chr_to_processed=chrSym)

        hg38_genome_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.fa'
        hg38_transcriptome_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/transcriptome.fa'
        mm10_genome_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/mm10.fa'
        mm10_transcriptome_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/mouse_transcriptome.fa'

        # 对某一个 tmaf 文件的操作
        tMaf_FileName = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/chrs/chr1.tmaf'
        ctMaf_FileName = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/chrs/chr1.ctmaf'
        mafClass.check_tMaf_format( tMaf_FileName, hg38_genome_file, hg38_transcriptome_file, mm10_genome_file, mm10_transcriptome_file ) # 检查tMaf的正确性
        mafClass.tMaf2ctMaf( tMaf_FileName, ctMaf_FileName ) # tMaf 转成 ctMaf 文件

        # 针对目录下的所有 tmaf 文件进行操作
        tMaf_Dir = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/chrs/'
        ctMaf_Dir = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/chrs/tcMaf/'
        mafClass.check_tMaf_format_inDir( tMaf_Dir, hg38_genome_file, hg38_transcriptome_file, mm10_genome_file, mm10_transcriptome_file )
        mafClass.tMaf2ctMaf_inDir(tMaf_Dir, ctMaf_Dir)

"""

"""

测试专用：

# 文件
hg38_GTF_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Homo_sapiens.gtf'
mm10_GTF_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Mus_musculus.GRCm38.81.gtf'

hg38_genome_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.fa'
mm10_genome_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/mm10.fa'

hg38_transcriptome_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/transcriptome.fa'
mm10_transcriptome_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/mouse_transcriptome.fa'


# 加载Handle
hg38_transcriptome = seq.seqClass(hg38_transcriptome_file)
mm10_transcriptome = seq.seqClass(mm10_transcriptome_file)

hg38_genome = seq.seqClass(hg38_genome_file)
mm10_genome = seq.seqClass(mm10_genome_file)

hg38_gt = genome2Trans.genome2TransClass(hg38_GTF_file)
mm10_gt = genome2Trans.genome2TransClass(mm10_GTF_file)

"""



import genome2Trans
import seq
import threading
import os
import re


class gMaf_2_tMaf_Methods(object):
    "把 gmaf 转成 tmaf 相关的函数"
    @staticmethod
    def sparse_mafS_Line(Line):
        """
        解析maf文件中的s行

        比如： sparse_S_Line（"s hg38.chr1 61943 29 + 248956422 ctccgcctcccggggtcaagctattctcc")

        """
        (label, species_chr, chrStart, length, strand, chrLength, Seq) = Line.strip().split()
        (species, Chr) = species_chr.split('.')
        chrStart = int(chrStart); length = int(length); chrLength = int(chrLength)
        if strand == '-':
            chrEnd = chrLength-chrStart
            chrStart = chrLength-chrStart-length
        else:
            chrEnd = chrStart+length
        return Chr, species, chrStart, chrEnd, strand, Seq
    @staticmethod
    def index_Base_Coor(seq, base_index):
        """从seq序列中找出坐标为base_index的全局坐标
        比如从 ATC--ATG-GTA 中找到碱基坐标为6的全局坐标
        base_index为1-based，返回0-based

        比如： index_Base_Coor('ATC--AGCTA--AGT', 10)

        """
        global_index = -1
        alphaCount = 0
        for idx in range(len(seq)):
            global_index += 1
            if seq[idx] != '-':
                alphaCount += 1 
                if alphaCount == base_index:
                    return global_index
        print 'Error: ',seq,'doesn\'t have a base_index', base_index
        return False
    @staticmethod
    def find_Base_Coor(seq, global_index, isStart):
        """从seq序列中找出全局坐标为global_index的碱基坐标
        比如从 ATC--ATG-GTA 中找到全局坐标为6的碱基坐标为5
        global_index为0-based，返回1-based
        isStart表示当global_index遇到的字符为-时，应该把base_index加一还是不变

        比如： find_Base_Coor('ATC--AGCTA--AGT', 3, isStart=1)

        """
        base_index = 0
        for idx in range(len(seq[:global_index+1])):
            if seq[idx] != '-':
                base_index += 1
        if seq[global_index] == '-':
            if isStart:
                return base_index+1
        return base_index
    @staticmethod
    def SeqSeg(seq1, seq2, alphaPosStart1, alphaPosEnd1, alphaPosStart2, alphaPosEnd2):
        """ [1_based_start, 1_based_end]
            seq1: ATC--AGCTA--AGT
            seq2: T-GGTTC-ATGTT-A
            获取seq1中[start1, end1]和seq2中[start2, end2]刚好重叠的那一部分并计算碱基偏移值

            比如：SeqSeg('ATC--AGCTA--AGT', 'T-GGTTC-ATGTT-A', 3, 10, 7, 10)

        """
        seqDict = {}
        seqPosStart1 = gMaf_2_tMaf_Methods.index_Base_Coor(seq1, alphaPosStart1)
        seqPosEnd1 = gMaf_2_tMaf_Methods.index_Base_Coor(seq1, alphaPosEnd1)
        seqPosStart2 = gMaf_2_tMaf_Methods.index_Base_Coor(seq2, alphaPosStart2)
        seqPosEnd2 = gMaf_2_tMaf_Methods.index_Base_Coor(seq2, alphaPosEnd2)
        if (not seqPosStart1) or (not seqPosStart2) or (not seqPosEnd1) or (not seqPosEnd2):
            return False
        if seqPosStart1 <= seqPosEnd2 and seqPosStart2 <= seqPosEnd1:
            Start = max(seqPosStart1, seqPosStart2)
            End = min(seqPosEnd1, seqPosEnd2)
            seqDict['seq1'] = seq1[Start:End+1]
            seqDict['seq2'] = seq2[Start:End+1]
            seqDict['start_bias1'] = gMaf_2_tMaf_Methods.find_Base_Coor(seq1, Start, isStart=1) - alphaPosStart1
            seqDict['start_bias2'] = gMaf_2_tMaf_Methods.find_Base_Coor(seq2, Start, isStart=1) - alphaPosStart2
            seqDict['end_bias1'] = gMaf_2_tMaf_Methods.find_Base_Coor(seq1, End, isStart=0) - alphaPosEnd1
            seqDict['end_bias2'] = gMaf_2_tMaf_Methods.find_Base_Coor(seq2, End, isStart=0) - alphaPosEnd2
            return seqDict
        return False
    @staticmethod
    def MAF2TransMAF(gMaf_FileName, tMaf_FileName, species1_gt, species2_gt, Chr_to_processed=None):
        """
        把基因组坐标的MAF转成转录组坐标的MAF
        Chr_to_processed是指定要处理的染色体，None表示全都处理
        """
        IN = open(gMaf_FileName)
        OUT = open(tMaf_FileName, 'w')
        line = IN.readline()
        groupCount = 0
        last_chr = ''
        while line:
            if line[0] != 's':
                line = IN.readline()
                continue
            groupCount += 1
            (Chr1, species1, chrStart1, chrEnd1, strand1, Seq1) = gMaf_2_tMaf_Methods.sparse_mafS_Line(line)
            if Chr_to_processed:
                if Chr1 != Chr_to_processed:
                    line = IN.readline()
                    line = IN.readline()
                    continue
            if Chr1 != last_chr:
                print 'Now Processing '+Chr1+'...'
                last_chr = Chr1
            line = IN.readline()
            (Chr2, species2, chrStart2, chrEnd2, strand2, Seq2) = gMaf_2_tMaf_Methods.sparse_mafS_Line(line)
            try:
                sense_species1_trans = species1_gt.genomeCoor2TransCoor(Chr1, chrStart1, chrEnd1, strand1)
                sense_species2_trans = species2_gt.genomeCoor2TransCoor(Chr2, chrStart2, chrEnd2, strand2)
                antisense_species1_trans = species1_gt.genomeCoor2TransCoor(Chr1, chrStart1, chrEnd1, seq.reverse_strand(strand1))
                antisense_species2_trans = species2_gt.genomeCoor2TransCoor(Chr2, chrStart2, chrEnd2, seq.reverse_strand(strand2))
            except KeyError:
                line = IN.readline()
                continue
            if len(sense_species1_trans)>0 and len(sense_species2_trans)>0:
                for item1 in sense_species1_trans:
                    for item2 in sense_species2_trans:
                        if strand1 == '-':
                            seqStart1 = chrEnd1 - item1[2] + 1
                            seqEnd1 = chrEnd1 - item1[1] + 1
                        else:
                            seqStart1 = item1[1] - chrStart1
                            seqEnd1 = item1[2] - chrStart1
                        if strand2 == '-':
                            seqStart2 = chrEnd2 - item2[2] + 1
                            seqEnd2 = chrEnd2 - item2[1] + 1
                        else:
                            seqStart2 = item2[1] - chrStart2
                            seqEnd2 = item2[2] - chrStart2
                        seqDict = gMaf_2_tMaf_Methods.SeqSeg(Seq1, Seq2, seqStart1, seqEnd1, seqStart2, seqEnd2)
                        if not seqDict: continue
                        if strand1 == '+':
                            print >>OUT, "%s.%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s" % (species1, Chr1, item1[1]-1+seqDict['start_bias1'], item1[2]+seqDict['end_bias1'], strand1, item1[3], item1[4]+seqDict['start_bias1'], item1[5]+seqDict['end_bias1'], seqDict['seq1'])
                        else:
                            print >>OUT, "%s.%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s" % (species1, Chr1, item1[1]-1-seqDict['end_bias1'], item1[2]-seqDict['start_bias1'], strand1, item1[3], item1[4]+seqDict['start_bias1'], item1[5]+seqDict['end_bias1'], seqDict['seq1'])
                        if strand2 == '+':
                            print >>OUT, "%s.%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s" % (species2, Chr2, item2[1]-1+seqDict['start_bias2'], item2[2]+seqDict['end_bias2'], strand2, item2[3], item2[4]+seqDict['start_bias2'], item2[5]+seqDict['end_bias2'], seqDict['seq2'])
                        else:
                            print >>OUT, "%s.%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s" % (species2, Chr2, item2[1]-1-seqDict['end_bias2'], item2[2]-seqDict['start_bias2'], strand2, item2[3], item2[4]+seqDict['start_bias2'], item2[5]+seqDict['end_bias2'], seqDict['seq2'])
                        print >>OUT, ""
            if len(antisense_species1_trans)>0 and len(antisense_species2_trans)>0:
                for item1 in antisense_species1_trans:
                    for item2 in antisense_species2_trans:
                        if seq.reverse_strand(strand1) == '-':
                            seqStart1 = chrEnd1 - item1[2] + 1
                            seqEnd1 = chrEnd1 - item1[1] + 1
                        else:
                            seqStart1 = item1[1] - chrStart1
                            seqEnd1 = item1[2] - chrStart1
                        if seq.reverse_strand(strand2) == '-':
                            seqStart2 = chrEnd2 - item2[2] + 1
                            seqEnd2 = chrEnd2 - item2[1] + 1
                        else:
                            seqStart2 = item2[1] - chrStart2
                            seqEnd2 = item2[2] - chrStart2
                        seqDict = gMaf_2_tMaf_Methods.SeqSeg(seq.reverse_comp(Seq1), seq.reverse_comp(Seq2), seqStart1, seqEnd1, seqStart2, seqEnd2)
                        if not seqDict: continue
                        if seq.reverse_strand(strand1) == '-':
                            print >>OUT, "%s.%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s" % (species1, Chr1, item1[1]-1-seqDict['end_bias1'], item1[2]-seqDict['start_bias1'], seq.reverse_strand(strand1), item1[3], item1[4]+seqDict['start_bias1'], item1[5]+seqDict['end_bias1'], seqDict['seq1'])
                        else:
                            print >>OUT, "%s.%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s" % (species1, Chr1, item1[1]-1+seqDict['start_bias1'], item1[2]+seqDict['end_bias1'], seq.reverse_strand(strand1), item1[3], item1[4]+seqDict['start_bias1'], item1[5]+seqDict['end_bias1'], seqDict['seq1'])
                        if seq.reverse_strand(strand2) == '-':
                            print >>OUT, "%s.%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s" % (species2, Chr2, item2[1]-1-seqDict['end_bias2'], item2[2]-seqDict['start_bias2'], seq.reverse_strand(strand2), item2[3], item2[4]+seqDict['start_bias2'], item2[5]+seqDict['end_bias2'], seqDict['seq2'])
                        else:
                            print >>OUT, "%s.%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s" % (species2, Chr2, item2[1]-1+seqDict['start_bias2'], item2[2]+seqDict['end_bias2'], seq.reverse_strand(strand2), item2[3], item2[4]+seqDict['start_bias2'], item2[5]+seqDict['end_bias2'], seqDict['seq2'])
                        print >>OUT, ""
                        break
                    break
            line = IN.readline()
        print 'Finished '+last_chr+' Processed'
        OUT.close()


class check_tMaf_Methods(object):
    """
    一套检查输出的tMaf结果的函数集合
    """
    @staticmethod
    def sparse_tMaf_Line(Line):
        """
        解析tMaf文件中的非空行
        """
        (species_chr, genome_s, genome_e, strand, transID, trans_s, trans_e, seq) = Line.strip().split()
        Chr = species_chr.split('.')[1]
        return Chr, int(genome_s), int(genome_e), strand, transID, int(trans_s), int(trans_e), seq
    @staticmethod
    def base_len(seq):
        """
        纯碱基长度：

        比如：base_len('A-GTCGAGTCG-TAGT--CGATC')

        """
        bLen = 0
        for idx in range(len(seq)):
            if seq[idx] != '-':
                bLen += 1
        return bLen
    @staticmethod
    def pure_seq(rawAlignSeq):
        """
        去除原始比对序列中的 - 号，得到纯序列
        
        比如：pure_seq("GCAAC--AGC-----GGCG--ACAGGC")
        
        """
        pureSeq = ''
        for idx in range(len(rawAlignSeq)):
            if rawAlignSeq[idx] != '-':
                pureSeq += rawAlignSeq[idx]
        return pureSeq
    @staticmethod
    def check_transcriptome_MAF(tMaf_FileName, sp1_genome_FileName, sp1_transcriptome_FileName, sp2_genome_FileName, sp2_transcriptome_FileName, verbose=False):
        sp1_genomeSeq = seq.seqClass(sp1_genome_FileName)
        sp1_transcriptomeSeq = seq.seqClass(sp1_transcriptome_FileName)
        sp2_genomeSeq = seq.seqClass(sp2_genome_FileName)
        sp2_transcriptomeSeq = seq.seqClass(sp2_transcriptome_FileName)
        IN = open(tMaf_FileName)
        human = IN.readline()
        mouse = IN.readline()
        line = IN.readline()
        count = 1
        last_chr = ''
        while line:
            (Chr1, genome_s1, genome_e1, strand1, transID1, trans_s1, trans_e1, seq1) = check_tMaf_Methods.sparse_tMaf_Line(human)
            (Chr2, genome_s2, genome_e2, strand2, transID2, trans_s2, trans_e2, seq2) = check_tMaf_Methods.sparse_tMaf_Line(mouse)
            if Chr1 != last_chr:
                print 'Now Checking '+Chr1
                last_chr = Chr1
            try:
                if len(seq1) != len(seq2):
                    print 'Error 1', count
                    print human
                    print mouse
                    return
                if check_tMaf_Methods.base_len(seq1) != trans_e1 - trans_s1 + 1:
                    print 'Error 2', count
                    print human
                    print mouse
                    return
                if check_tMaf_Methods.base_len(seq2) != trans_e2 - trans_s2 + 1:
                    print 'Error 3', count
                    print human
                    print mouse
                    return
                if check_tMaf_Methods.base_len(seq1) != genome_e1 - genome_s1:
                    print 'Error 4', count
                    print human
                    print mouse
                    return
                if check_tMaf_Methods.base_len(seq2) != genome_e2 - genome_s2:
                    print 'Error 5', count
                    print human
                    print mouse
                    return
                if sp1_genomeSeq.fetch(Chr1, genome_s1, genome_e1, strand1).upper() != sp1_transcriptomeSeq.fetch(transID1, trans_s1-1, trans_e1).upper() or sp1_transcriptomeSeq.fetch(transID1, trans_s1-1, trans_e1).upper() != check_tMaf_Methods.pure_seq(seq1).upper():
                    print 'Error 6', count
                    print human
                    print mouse
                    return
                if sp2_genomeSeq.fetch(Chr2, genome_s2, genome_e2, strand2).upper() != sp2_transcriptomeSeq.fetch(transID2, trans_s2-1, trans_e2).upper() or sp2_transcriptomeSeq.fetch(transID2, trans_s2-1, trans_e2).upper() != check_tMaf_Methods.pure_seq(seq2).upper():
                    print 'Error 7', count
                    print human
                    print mouse
                    return
            except KeyError:
                if verbose:
                    print 'Check tMaf: Key Error -> '
                    print '\t\tsp_1: ', (Chr1, genome_s1, genome_e1, strand1, transID1, trans_s1, trans_e1)
                    print '\t\tsp_2: ', (Chr2, genome_s2, genome_e2, strand2, transID2, trans_s2, trans_e2)
                else:
                    pass
            human = IN.readline()
            mouse = IN.readline()
            line = IN.readline()
            count += 1
        del sp1_genomeSeq
        del sp2_genomeSeq
        del sp1_transcriptomeSeq
        del sp2_transcriptomeSeq
        print 'All Right'





class tMaf_ctMaf_Methods(object):
    """
    tMaf 转成 ctMaf 的相关函数

    使用方法：
        tMaf_FileName = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/chrs/chr1.tmaf'
        ctMaf_FileName = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/chrs/chr1.ctmaf'
        tMafContainer = tMaf_ctMaf_Methods.read_tMaf(tMaf_FileName)
        combined_tMafContainer = tMaf_ctMaf_Methods.combine_tMaf(tMafContainer)
        tMaf_ctMaf_Methods.write_combined_tMaf( combined_tMafContainer, ctMaf_FileName )

    """
    @staticmethod
    def read_tMaf(tMaf_FileName):
        "读取 tMaf 格式的函数"
        tMafContainer = {}
        IN = open(tMaf_FileName)
        line = IN.readline()
        while line:
            (speciels_chr1, chr_start1, chr_end1, strand1, transID1, trans_start1, trans_end1, seq1) = line.strip().split()
            line = IN.readline()
            (speciels_chr2, chr_start2, chr_end2, strand2, transID2, trans_start2, trans_end2, seq2) = line.strip().split()
            line = IN.readline()
            line = IN.readline()
            trans_key = transID1+'-'+transID2
            seqPair = { 'seq1':[ int(trans_start1), int(trans_end1), seq1 ], 'seq2':[ int(trans_start2), int(trans_end2), seq2 ] }
            if trans_key in tMafContainer:
                tMafContainer[ trans_key ].append( seqPair )
            else:
                tMafContainer[ trans_key ] = [ seqPair ]
        for trans_key in tMafContainer:
            tMafContainer[trans_key].sort(key=lambda x: x['seq1'][0])
        return tMafContainer
    @staticmethod
    def _uniq_ListPair(trans1_ListPair, trans2_ListPair):
        " Unique List Pair "
        length = len(trans1_ListPair)
        idx = 0
        while idx < length - 1:
            if trans1_ListPair[idx] == trans1_ListPair[idx+1] and trans2_ListPair[idx] == trans2_ListPair[idx+1]:
                del trans1_ListPair[idx]
                del trans2_ListPair[idx]
                length -= 1
            else:
                idx += 1
    @staticmethod
    def _combine_ListPair(trans1_ListPair, trans2_ListPair):
        " 把几个ListPair合并在一起 "
        length = len(trans1_ListPair)
        idx = 0
        while idx < length - 1:
            if trans1_ListPair[idx][1] == trans1_ListPair[idx+1][0]-1 and trans2_ListPair[idx][1] == trans2_ListPair[idx+1][0]-1:
                trans1_ListPair[idx][1] = trans1_ListPair[idx+1][1]
                trans2_ListPair[idx][1] = trans2_ListPair[idx+1][1]
                trans1_ListPair[idx][2] += trans1_ListPair[idx+1][2]
                trans2_ListPair[idx][2] += trans2_ListPair[idx+1][2]
                del trans1_ListPair[idx+1]
                del trans2_ListPair[idx+1]
                length -= 1
            else:
                idx += 1
    @staticmethod
    def combine_tMaf(tMafContainer):
        " 合并 tMaf 文件 "
        filtered_combined_tMaf = {}
        for trans_key in tMafContainer:
            trans1_List = [ item['seq1'] for item in tMafContainer[trans_key] ]
            trans2_List = [ item['seq2'] for item in tMafContainer[trans_key] ]
            tMaf_ctMaf_Methods._uniq_ListPair(trans1_List, trans2_List)
            tMaf_ctMaf_Methods._combine_ListPair(trans1_List, trans2_List)
            (transID1, transID2) = trans_key.split('-')
            filtered_combined_tMaf[ trans_key ] = [ trans1_List, trans2_List ]
        return filtered_combined_tMaf
    @staticmethod
    def write_combined_tMaf( combined_tMafContainer, ctMaf_FileName ):
        " 把 tMaf 写到一个 ctMaf 文件中 "
        OUT = open(ctMaf_FileName, 'w')
        for trans_key in combined_tMafContainer:
            (transID1, transID2) = trans_key.split('-')
            print >>OUT, '>' + transID1 + '\t' + transID2
            for idx in range(len(combined_tMafContainer[trans_key][0])):
                segment1 = combined_tMafContainer[trans_key][0][idx]
                segment2 = combined_tMafContainer[trans_key][1][idx]
                print >>OUT, str(segment1[0])+'-'+str(segment1[1])+'\t'+segment1[2]+'\t'
                print >>OUT, str(segment2[0])+'-'+str(segment2[1])+'\t'+segment2[2]+'\t'
                print >>OUT, ''
        OUT.close()


class mafClass(object):
    """
    maf文件相关的操作：
    1. 读取maf文件（未完成）
    2. 把基因组坐标的maf文件转成转录组坐标的maf文件: gmaf2tMaf（2017-4-25）
    3. 检查tMaf文件的正确性：check_tMaf_format（2017-5-4）
    4. 合并 tMaf 文件为 ctMaf 文件：tMaf2ctMaf（2017-5-4）
    5. 把基因组坐标的maf文件转成基因坐标的maf文件（未完成）
    6. 计算两个物种间的相似程度（未完成）
    
    使用方法：

        # maf 转成 tMaf
        mafFile = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/hg38.mm10.synNet.maf'
        hg38_GTF_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Homo_sapiens.gtf'
        mm10_GTF_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Mus_musculus.GRCm38.81.gtf'
        maf = mafClass(mafFile)
        maf.gmaf2tMaf(hg38_GTF_file, mm10_GTF_file, outFileName='/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/chrY.tmaf', Chr_to_processed='chrY')
        或者分别同时处理多个染色体：
        for chrSym in ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']:
            maf.gmaf2tMaf(hg38_GTF_file, mm10_GTF_file, outFileName='/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/chrs/'+chrSym+'.tmaf', Chr_to_processed=chrSym)

        hg38_genome_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.fa'
        hg38_transcriptome_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/transcriptome.fa'
        mm10_genome_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/mm10.fa'
        mm10_transcriptome_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/mouse_transcriptome.fa'

        # 对某一个 tmaf 文件的操作
        tMaf_FileName = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/chrs/chr1.tmaf'
        ctMaf_FileName = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/chrs/chr1.ctmaf'
        mafClass.check_tMaf_format( tMaf_FileName, hg38_genome_file, hg38_transcriptome_file, mm10_genome_file, mm10_transcriptome_file ) # 检查tMaf的正确性
        mafClass.tMaf2ctMaf( tMaf_FileName, ctMaf_FileName ) # tMaf 转成 ctMaf 文件

        # 针对目录下的所有 tmaf 文件进行操作
        tMaf_Dir = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/chrs/'
        ctMaf_Dir = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/chrs/tcMaf/'
        mafClass.check_tMaf_format_inDir( tMaf_Dir, hg38_genome_file, hg38_transcriptome_file, mm10_genome_file, mm10_transcriptome_file )
        mafClass.tMaf2ctMaf_inDir(tMaf_Dir, ctMaf_Dir)

    """
    def __init__(self, gmap_FileName=None):
        self.gmaf_fileName = gmap_FileName
        self.species1_gtfFile = ''
        self.species2_gtfFile = ''
        print 'mafClass: maf文件是0-based左闭右开'
    def gmaf2tMaf(self, species1_gtfFile, species2_gtfFile, outFileName, Chr_to_processed=None):
        """
        把genome coordinated maf文件转成transcriptome coordinated maf文件
        """
        "读文件"
        if not self.gmaf_fileName:
            print "Error: 没有定义 gMaf 文件"
            return
        if self.species1_gtfFile != species1_gtfFile:
            self.species1_gtfFile = species1_gtfFile
            self.species1_gt = genome2Trans.genome2TransClass(species1_gtfFile)
        if self.species2_gtfFile != species2_gtfFile:
            self.species2_gtfFile = species2_gtfFile
            self.species2_gt = genome2Trans.genome2TransClass(species2_gtfFile)
        "转坐标"
        trd = threading.Thread(target=gMaf_2_tMaf_Methods.MAF2TransMAF,args=(self.gmaf_fileName, outFileName, self.species1_gt, self.species2_gt, Chr_to_processed))
        trd.setDaemon(True)
        trd.start()
    @staticmethod
    def check_tMaf_format(tMaf_FileName, sp1_genome_FileName, sp1_transcriptome_FileName, sp2_genome_FileName, sp2_transcriptome_FileName):
        "检查 tMaf 文件的正确性"
        check_tMaf_Methods.check_transcriptome_MAF(tMaf_FileName, sp1_genome_FileName, sp1_transcriptome_FileName, sp2_genome_FileName, sp2_transcriptome_FileName)
    @staticmethod
    def check_tMaf_format_inDir(tMaf_Dir, sp1_genome_FileName, sp1_transcriptome_FileName, sp2_genome_FileName, sp2_transcriptome_FileName):
        "检查目录下所有.tmaf文件的正确性"
        tMaf_Dir = os.path.abspath(tMaf_Dir)
        tMaf_Files = os.listdir(tMaf_Dir)
        for tMaf in tMaf_Files:
            if re.search("\.tmaf$", tMaf):
                mafClass.check_tMaf_format(tMaf_Dir+'/'+tMaf, sp1_genome_FileName, sp1_transcriptome_FileName, sp2_genome_FileName, sp2_transcriptome_FileName)
    @staticmethod
    def tMaf2ctMaf(tMaf_FileName, ctMaf_FileName):
        " tMaf 转成 combined tMaf "
        tMafContainer = tMaf_ctMaf_Methods.read_tMaf(tMaf_FileName)
        combined_tMafContainer = tMaf_ctMaf_Methods.combine_tMaf(tMafContainer)
        tMaf_ctMaf_Methods.write_combined_tMaf( combined_tMafContainer, ctMaf_FileName )
    @staticmethod
    def tMaf2ctMaf_inDir(tMaf_Dir, ctMaf_Dir):
        " 把目录下的所有 tMaf 文件转成 ctMaf 文件 "
        tMaf_Dir = os.path.abspath(tMaf_Dir)
        ctMaf_Dir = os.path.abspath(ctMaf_Dir)
        tMaf_Files = os.listdir(tMaf_Dir)
        for tMaf in tMaf_Files:
            if re.search("\.tmaf$", tMaf):
                print 'Now tMaf To ctMaf '+tMaf+'...'
                ctMaf = '.'.join(tMaf.split('.')[:-1]) + '.ctmaf'
                mafClass.tMaf2ctMaf( tMaf_Dir+'/'+tMaf, ctMaf_Dir+'/'+ctMaf )


class ctMaf_Homology_Methods(object):
    """
    通过 ctmaf 文件获取一个物种的转录本同源与另一个物种某转录本的信息。

    使用方法：
        ctMaf_FileName = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/hg38_mm10.ctmaf'
        ctMaf_Container = ctMaf_Homology_Methods.read_ctMaf_master_slave(ctMaf_FileName)
        ctMaf_Homology_Methods.get_slave_base_to_master_base(ctMaf_Container, 'ENST00000427245', 2, 7)

    """
    @staticmethod
    def read_ctMaf_master_slave(ctMaf_FileName):
        " 以 master-slave 方式读取 ctmaf "
        ctMaf_Container = {}
        IN = open(ctMaf_FileName)
        line = IN.readline()
        while line:
            if line[0] == '>':
                (transID1, transID2) = line[1:].strip().split()
                if transID1 in ctMaf_Container:
                    ctMaf_Container[ transID1 ][ transID2 ] = []
                else:
                    ctMaf_Container[ transID1 ] = {}
                    ctMaf_Container[ transID1 ][ transID2 ] = []
                line = IN.readline()
                while line and line[0] != '>':
                    (start_end1, seq1) = line.strip().split()
                    line = IN.readline()
                    (start_end2, seq2) = line.strip().split()
                    (start1, end1) = start_end1.split('-')
                    (start2, end2) = start_end2.split('-')
                    line = IN.readline() # 空行
                    line = IN.readline() # >行 或者 新行
                    ctMaf_Container[ transID1 ][ transID2 ].append( [ [int(start1), int(end1), seq1], [int(start2), int(end2), seq2] ] )
            else:
                print "Unexpected Error"
        return ctMaf_Container
    @staticmethod
    def get_slave_base_to_master_base(ctMaf_Container, master_TransID, query_Master_Start, query_Master_End, return_slave_gap=False):
        """ 
            获取 Slave 物种对应 Master 物种的同源区段序列 
            ===> 请提供 1-based 转录本起点和 1-based终点，左闭右闭 <===
            当 return_slave_gap 为 True 时，如果slave_base是一个gap，那就不返回
        """
        slave_master = {}
        if master_TransID in ctMaf_Container:
            for slave_TransID in ctMaf_Container[master_TransID]:
                for pair in ctMaf_Container[master_TransID][slave_TransID]:
                    master_start = pair[0][0]
                    master_end = pair[0][1]
                    master_seq = pair[0][2]
                    slave_start = pair[1][0]
                    slave_end = pair[1][1]
                    slave_seq = pair[1][2]
                    Start = max( master_start, query_Master_Start )
                    End = min( master_end, query_Master_End )
                    if Start <= End:
                        if slave_TransID not in slave_master:
                            slave_master[slave_TransID] = []
                        base_coor_start = Start - master_start + 1
                        base_coor_end = End - master_start + 1
                        for master_base_seq_index in range( base_coor_start, base_coor_end+1 ):
                            global_index = gMaf_2_tMaf_Methods.index_Base_Coor(master_seq, master_base_seq_index)
                            master_base = master_seq[ global_index ].upper()
                            slave_base = slave_seq[ global_index ].upper()
                            if not return_slave_gap and slave_base == '-': continue
                            slave_base_seq_index = gMaf_2_tMaf_Methods.find_Base_Coor( slave_seq, global_index, isStart=1 )
                            slave_master[slave_TransID].append( [ master_base_seq_index+master_start-1, slave_base_seq_index+slave_start-1, master_base+'-'+slave_base ] )
        return slave_master
    @staticmethod
    def generate_base_correlate_file(ctMaf_Container, species1_anno, species2_anno, outFileName):
	import seq as SeqModule
        """
        把ctmaf文件转成下面的格式：
            ENST00000341423 30      ENSMUST00000085546      1       A-A
            ENST00000341423 31      ENSMUST00000085546      2       A-A
            ENST00000341423 32      ENSMUST00000085546      3       T-T
            ENST00000341423 33      ENSMUST00000085546      4       G-G
            ENST00000341423 34      ENSMUST00000085546      5       T-T
            ENST00000341423 35      ENSMUST00000085546      6       T-T
            ENST00000341423 36      ENSMUST00000085546      7       A-A

        使用方法：
	    from seq import *

            species1_refbed_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Homo_sapiens.transCoordinate.bed'
            species2_refbed_file = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Mus_musculus.transCoordinate.bed'

            species1_anno = anno_Methods.loadGTFBed(species1_refbed_file)
            species2_anno = anno_Methods.loadGTFBed(species2_refbed_file)

            ctMaf_FileName = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/hg38_mm10.ctmaf'
            ctMaf_Container = ctMaf_Homology_Methods.read_ctMaf_master_slave(ctMaf_FileName)

            outFileName = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/hg38_mm10.basecoor'
            ctMaf_Homology_Methods.generate_base_correlate_file(ctMaf_Container, species1_anno, species2_anno, outFileName)
        """
        import sys
        species1_transList = sorted( species1_anno.keys() )
        species2_transList = sorted( species2_anno.keys() )
        OUT = open(outFileName, 'w')
        total = len(ctMaf_Container)
        count = 0
        for sp1_transID in ctMaf_Container:
            if count % 10 == 0: sys.stdout.write('.'); sys.stdout.flush()
            if count % 1000 == 0: sys.stdout.write("%.3f%%" % ( 100.0*count/total )); sys.stdout.flush()
            count += 1
            if not SeqModule.biSearch(sp1_transID, species1_transList): continue
            transLen = species1_anno[sp1_transID]['trans_len']
            homo = ctMaf_Homology_Methods.get_slave_base_to_master_base(ctMaf_Container, sp1_transID, 1, transLen)
            PairTransSelectList = []
            for sp2_transID in homo:
                if not SeqModule.biSearch(sp2_transID, species2_transList): continue
                sp1_gene_type = SeqModule.anno_Methods.gene_type(species1_anno[sp1_transID]['gene_type'])
                sp2_gene_type = SeqModule.anno_Methods.gene_type(species2_anno[sp2_transID]['gene_type'])
                if sp1_gene_type == sp2_gene_type:
                    PairTransSelectList.append( [ sp2_transID, len(homo[sp2_transID])] )
            if len(PairTransSelectList) == 0: continue
            sp2_transID = sorted(PairTransSelectList, key=lambda x: x[1])[-1][0]
            for basePair in homo[sp2_transID]:
                hek_coor = basePair[0] - 1
                mes_coor = basePair[1] - 1
                item = [ sp1_transID, str(hek_coor+1), sp2_transID, str(mes_coor+1), basePair[2] ]
                print >>OUT, '\t'.join(item)
        OUT.close()









