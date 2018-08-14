#-*- coding:utf-8 -*-

"""

有待完成的功能：
    1. 读取随机顺序的GFF3文件
    2. 有序输出GFF3文件
    3. NCBI/Gencode之间基因名字互转功能

Version:
    2017-5-29
    2017-5-30 add __getRNATransID/__build_transID/__getRNAGeneType functions
"""

import sys, os

class NCBI_Genome_Class(object):
    """这是一个操作NCBI注释相关的类
    主要包括下面的功能：
        1. 读取NCBI的GFF3注释文件： 
            read_ncbi_gff3(gff3FileName)
        2. 构建 NC 标志染色体和经典的 chr 标志染色体之间的转换：
            build_chr_To_NC_dict(pureNCID=False)
            build_NC_To_chr_dict(pureNCID=False)
        3. 把原始的GFF3转换成基因组坐标的Bed文件
            write_genomeCoor_bed(genomeCoorFileName, onlyChr=False, pureTransID=True)
        4. 把基因组坐标的Bed文件转成转录组坐标的Bed文件
            static genomeCoorBed_To_transCoorBed(genomeCoorFileName, transCoorFileName)
        5. 根据GFF3注释从基因组中解析出转录组
            writeTranscriptome(genomeFileName, transcriptomeFileName, genomeChrSym='chr', pureTransID=False, verbose=True)
    
    注意：
        NCBI官方下载的 ==转录组== 文件与这里用函数writeTranscriptome取出的不对应，因为官方的文件在mRNA的末尾加上了polyA

    使用方法：
        from NCBI_Genome import *

        # read gff3 annotation and build a NCBIClass Object
        NCBI_hg38 = NCBI_Genome_Class('/Share/home/zhangqf8/lipan/DYNAMIC/GTF/ncbi/GRCh38.p7.gff3')
        NCBI_hg19 = NCBI_Genome_Class('/Share/home/zhangqf8/lipan/DYNAMIC/GTF/ncbi/GRCh37.p13.gff3')
        NCBI_mm10 = NCBI_Genome_Class('/Share/home/zhangqf8/lipan/DYNAMIC/GTF/ncbi/GRCm38.p4.gff3')

        # convert gff3 transcripts infomation to genome-coordinated bed file
        NCBI_hg38.write_genomeCoor_bed(genomeCoorFileName="~/hg38.genomeCoor.bed", onlyChr=True, pureTransID=True)
        NCBI_hg19.write_genomeCoor_bed(genomeCoorFileName="~/hg19.genomeCoor.bed", onlyChr=True, pureTransID=True)
        NCBI_mm10.write_genomeCoor_bed(genomeCoorFileName="~/mm10.genomeCoor.bed", onlyChr=True, pureTransID=True)

        # covert genome-coordinated bed file to transcript-coordinated bed file
        NCBI_Genome_Class.genomeCoorBed_To_transCoorBed(os.environ.get("HOME")+"/hg38.genomeCoor.bed", os.environ.get("HOME")+"/hg38.transCoor.bed")
        NCBI_Genome_Class.genomeCoorBed_To_transCoorBed(os.environ.get("HOME")+"/hg19.genomeCoor.bed", os.environ.get("HOME")+"/hg19.transCoor.bed")
        NCBI_Genome_Class.genomeCoorBed_To_transCoorBed(os.environ.get("HOME")+"/mm10.genomeCoor.bed", os.environ.get("HOME")+"/mm10.transCoor.bed")

        # build chr => NC dict
        chrToNC_hg38 = NCBI_hg38.build_chr_To_NC_dict(pureNCID=False)
        chrToNC_hg19 = NCBI_hg19.build_chr_To_NC_dict(pureNCID=True)
        chrToNC_mm10 = NCBI_mm10.build_chr_To_NC_dict(pureNCID=True)

        # build NC => chr dict
        NCToChr_hg38 = NCBI_hg38.build_NC_To_chr_dict(pureNCID=True)
        NCToChr_hg19 = NCBI_hg19.build_NC_To_chr_dict(pureNCID=False)
        NCToChr_mm10 = NCBI_mm10.build_NC_To_chr_dict(pureNCID=True)

        # parse transcriptome from genome fasta file
        # Gencode下载的基因组
        NCBI_hg38.writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/hg38.fa", os.environ.get("HOME")+"/hg38_transcriptome.fa", genomeChrSym='chr', pureTransID=True, verbose=True)
        NCBI_hg19.writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/hg19.fa", os.environ.get("HOME")+"/hg19_transcriptome.fa", genomeChrSym='chr', pureTransID=True, verbose=True)
        NCBI_mm10.writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/mm10.fa", os.environ.get("HOME")+"/mm10_transcriptome.fa", genomeChrSym='chr', pureTransID=True, verbose=True)
        # NCBI 下载的基因组
        NCBI_hg38.writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/ncbi/GRCh38.p7.fa", os.environ.get("HOME")+"/hg38_ncbi.fa", genomeChrSym='chr', pureTransID=True, verbose=True)
        NCBI_mm10.writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/ncbi/GRCm38.p4.fa", os.environ.get("HOME")+"/mm10_ncbi.fa", genomeChrSym='chr', pureTransID=True, verbose=True)
    """
    def __init__(self, gff3FileName):
        self.gff3FileName = gff3FileName
        self.gff3_container = self.read_ncbi_gff3(gff3FileName)
        print '\tbuild transID Set...'
        self.__build_transID()
    def __getRNATransID(self, rna_ID):
        """获得一个RNA的转录本ID
        测试
            __getRNATransID('rna192')
        """
        RNA = self.gff3_container['RNA'][rna_ID]
        if 'transcript_id' in RNA:
            return RNA['transcript_id']
        else:
            return rna_ID
        """
        if 'miRBase' in RNA:
            return RNA['gene'] + '_' + RNA['miRBase']
        if RNA['gbkey'] in ('tRNA', 'rRNA'):
            return RNA['gene']+'_'+RNA['gbkey']
        """
        return ""
    def __build_transID(self):
        """生成所有RNA的ID
        测试：
            __build_transID()
        """
        tmp_TransIDDict = {}
        RNAIDDict = {}
        for rna_ID in self.gff3_container['RNA']:
            TransID = self.__getRNATransID(rna_ID)
            """
            try:
                tmp_TransIDDict[TransID]['inIt'] = 0
                TransID = self.gff3_container['RNA'][rna_ID]['GeneID'] + '_' + TransID
            except KeyError:
                tmp_TransIDDict[TransID] = {}
            """
            RNAIDDict[rna_ID] = TransID
        self.rnaid_2_transID = RNAIDDict
    def __getRNAGeneType(self, rna_ID):
        """获得一个RNA的基因类型
        测试
            __getRNAGeneType('rna192')
        """
        RNA = self.gff3_container['RNA'][rna_ID]
        try:
            #GeneType = self.gff3_container['gene'][RNA['GeneID']]['gene_biotype']
            GeneType = RNA['ncrna_class']
        except KeyError:
            GeneType = RNA['gbkey']
        return GeneType
    def __parse_attributes_gff3(self, attributesString):
        """ 解析attributes列
        测试：
            attri1 = 'ID=rna0;Parent=gene0;Dbxref=GeneID:497097,Genbank:XM_006495550.3,MGI:MGI:3528744;Name=XM_006495550.3;gbkey=mRNA;gene=Xkr4;model_evidence=Supporting evidence includes similarity to: 3 mRNAs%2C 32 ESTs%2C 4 Proteins%2C 110 long SRA reads%2C and 99%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 7 samples with support for all annotated introns;product=X-linked Kx blood group related 4%2C transcript variant X1;transcript_id=XM_006495550.3'
            __parse_attributes_gff3(attri1)
        """
        attributes = {}
        for Tuple in attributesString.split(';'):
            (Key, Value) = Tuple.split('=')
            if Key == 'Dbxref':
                for miniTuple in Value.split(','):
                    miniKey_miniValue = miniTuple.split(':')
                    if miniKey_miniValue[0] in attributes:
                         attributes[ miniKey_miniValue[0]+'_2' ] = miniKey_miniValue[-1]
                    else:
                        attributes[ miniKey_miniValue[0] ] = miniKey_miniValue[-1]
            else:
                attributes[ Key ] = Value
        return attributes
    def read_ncbi_gff3(self, gff3FileName):
        """ 读取 ncbi gff3 文件
        测试：
            hg19_gff3_container = read_ncbi_gff3('/Share/home/zhangqf8/lipan/DYNAMIC/GTF/ncbi/GRCh37.p13.gff3')
            hg38_gff3_container = read_ncbi_gff3('/Share/home/zhangqf8/lipan/DYNAMIC/GTF/ncbi/GRCh38.p7.gff3')
            mm10_gff3_container = read_ncbi_gff3('/Share/home/zhangqf8/lipan/DYNAMIC/GTF/ncbi/GRCm38.p4.gff3')
        """
        print '\nRead ncbi GFF3 Annotation: ', gff3FileName
        gff3_container = {'region': {}, 'gene':{}, 'transcript':{}, 'mRNA':{}, 'ncRNA':{}, 'rRNA':{}, 'tRNA':{}, 'CDS':{}, 'exon':{}, 'RNA':{}}
        IN = open(gff3FileName)
        count = 0
        line = IN.readline()
        while line:
            count += 1
            if count % 100000 == 0: print '\tRead %d Lines...' % (count,)
            if line[0] == '#':
                line = IN.readline()
                continue
            (Chr, Source, Type, Chr_Start, Chr_End, Score, Strand, Phase, Attributes) = line.strip().split('\t')
            attributes = self.__parse_attributes_gff3(Attributes)
            attributes['chr'] = Chr; attributes['start'] = Chr_Start; attributes['end'] = Chr_End; attributes['strand'] = Strand
            if Type == 'region':
                if attributes['start'] == '1':
                    gff3_container[ 'region' ][ Chr ] = attributes
            if Type == 'gene':
                if 'gene' not in attributes:
                    attributes['gene'] = 'NULL'
                gff3_container[ 'gene' ][ attributes['GeneID'] ] = attributes
            if Type in ('transcript', 'mRNA', 'ncRNA', 'rRNA', 'primary_transcript'):
                if 'gene' not in attributes:
                    attributes['gene'] = 'NULL'
                if 'GeneID' not in attributes:
                    attributes['GeneID'] = 'NULL'
                gff3_container[ 'RNA' ][ attributes['ID'] ] = attributes
            if Type in ('lnc_RNA', 'miRNA', 'pseudogene', 'snoRNA'):
                if 'gene' not in attributes:
                    attributes['gene'] = 'NULL'
                if 'GeneID' not in attributes:
                    attributes['GeneID'] = 'NULL'
                attributes['ncrna_class'] = Type
                gff3_container[ 'RNA' ][ attributes['ID'] ] = attributes
            if Type == 'tRNA':
                attributes['GeneID'] = 'NULL'
                attributes['gene'] = attributes['product']
                gff3_container[ 'RNA' ][ attributes['ID'] ] = attributes
            if Type == 'exon':
                gff3_container[ 'exon' ][ attributes['Parent'] ] = gff3_container[ 'exon' ].get( attributes['Parent'], []) + [attributes]
            if Type == 'CDS':
                gff3_container[ 'CDS' ][ attributes['Parent'] ] = gff3_container[ 'CDS' ].get( attributes['Parent'], []) + [attributes]
            line = IN.readline()
        return gff3_container
    def build_chr_To_NC_dict(self, pureNCID=False):
        """ 构建 chrxxx => NCxxx 的字典
        测试：
            chrToNC_hg19 = build_chr_To_NC_dict(hg19_gff3_container, pureNCID=True)
            chrToNC_hg38 = build_chr_To_NC_dict(hg38_gff3_container, pureNCID=True)
            chrToNC_mm10 = build_chr_To_NC_dict(mm10_gff3_container, pureNCID=True)
        """
        chrToNC = {}
        regions = self.gff3_container['region']
        for ChrID in regions:
            if not ChrID.startswith('NC'): continue
            pureChrID = ChrID.split('.')[0] if pureNCID else ChrID
            if 'chromosome' in regions[ChrID]:
                chrToNC['chr'+regions[ChrID]['chromosome']] = pureChrID
            else:
                chrToNC['chrM'] = pureChrID
        return chrToNC
    def build_NC_To_chr_dict(self, pureNCID=False):
        """构建 NCxxx => chrxxx 的字典
        测试：
            NCToChr_hg19 = build_chr_To_NC_dict(hg19_gff3_container, pureNCID=False)
            NCToChr_hg38 = build_chr_To_NC_dict(hg38_gff3_container, pureNCID=True)
            NCToChr_mm10 = build_chr_To_NC_dict(mm10_gff3_container, pureNCID=True)
        """
        NCToChr = {}
        regions = self.gff3_container['region']
        for ChrID in regions:
            if not ChrID.startswith('NC'): continue
            pureChrID = ChrID.split('.')[0] if pureNCID else ChrID
            if 'chromosome' in regions[ChrID]:
                NCToChr[pureChrID] =  'chr'+regions[ChrID]['chromosome']
            else:
                NCToChr[pureChrID] = 'chrM'
        return NCToChr
    def __utr_String(self, exonString, cdsString, strand):
        """ 从 exonString 中解析出 cdsString
        测试：
            "plus strand"
            cdsString = "9710456-9710596,9715541-9715769,9715849-9716078,9716440-9716619,9716959-9717108,9717537-9717626,9718694-9718915,9719921-9720017,9720112-9720242,9720434-9720661,9720742-9720909,9721127-9721248,9721444-9721587,9721761-9721860,9721975-9722153,9722244-9722356,9722528-9722606,9723125-9723292,9723969-9724092,9724276-9724421,9724804-9724936,9726909-9727046"
            exonString = "9629889-9629915,9630138-9630224,9691467-9691571,9710424-9710596,9715541-9715769,9715849-9716078,9716440-9716619,9716959-9717108,9717537-9717626,9718694-9718915,9719921-9720017,9720112-9720242,9720434-9720661,9720742-9720909,9721127-9721248,9721444-9721587,9721761-9721860,9721975-9722153,9722244-9722356,9722528-9722606,9723125-9723292,9723969-9724092,9724276-9724421,9724804-9724936,9726909-9728920"
            strand = '+'
            __utr_String(exonString, cdsString, strand)
            "minus strand"
            cdsString = "9069504-9069536,9058152-9058250,9057448-9057608,9047610-9047734,9041785-9041937,9040064-9040189,9039800-9039987,9039552-9039662,9038828-9038929,9038431-9038506,9037897-9038024,9037586-9037789"
            exonString = "9073076-9073140,9069504-9069594,9058152-9058250,9057448-9057608,9047610-9047734,9041785-9041937,9040064-9040189,9039800-9039987,9039552-9039662,9038828-9038929,9038431-9038506,9037897-9038024,9036941-9037789"
            strand = '-'
            __utr_String(exonString, cdsString, strand)
        """
        def correctList(List, strand):
            """校正有重叠的序列
            测试：
                plus_error_list = [[1,29],[28,90],[102, 203]]
                correctList(plus_error_list, strand='+')
                plus_error_list

                minus_error_list = [[102, 203],[28,90],[1,29]]
                correctList(minus_error_list, strand='-')
                minus_error_list
            """
            if strand == '+':
                idx = 0
                while idx < len(List) - 1:
                    if List[idx][1] >= List[idx+1][0]:
                        List[idx] = [ List[idx][0], List[idx][1] ]
                        del List[idx+1]
                    else:
                        idx += 1
            if strand == '-':
                idx = 0
                while idx < len(List) - 1:
                    if List[idx][0] <= List[idx+1][1]:
                        List[idx] = [ List[idx+1][0], List[idx][1] ]
                        del List[idx+1]
                    else:
                        idx += 1
        exonList = exonString.split(',')
        exonList = [ [int(exon.split('-')[0]), int(exon.split('-')[1])] for exon in exonList ]
        correctList(exonList, strand)
        cdsList = cdsString.split(',')
        cdsList = [ [int(cds.split('-')[0]), int(cds.split('-')[1])] for cds in cdsList ]
        correctList(cdsList, strand)
        utr_5 = []; utr_3 = []
        if strand == '+':
            for (exon_start, exon_end) in exonList:
                if exon_end < cdsList[0][0]:
                    utr_5.append( (exon_start, exon_end) )
                elif (exon_start < cdsList[0][0] and exon_end == cdsList[0][1]) or (exon_start < cdsList[0][0] and len(cdsList) == 1):
                    utr_5.append( (exon_start, cdsList[0][0]-1) )
                    break
                elif exon_start == cdsList[0][0]:
                    break
                else:
                    print 'Impossible Event 1'
                    print 'exonString:',exonString
                    print 'cdsString:',cdsString
                    raise Exception("Impossible Event 1")
                    break
            for (exon_start, exon_end) in exonList[::-1]:
                if exon_start > cdsList[-1][1]:
                    utr_3.append( (exon_start, exon_end) )
                elif (cdsList[-1][1] < exon_end and exon_start == cdsList[-1][0]) or (cdsList[-1][1] < exon_end and len(cdsList) == 1):
                    utr_3.append( (cdsList[-1][1]+1, exon_end) )
                    break
                elif exon_end == cdsList[-1][1]:
                    break
                else:
                    print 'Impossible Event 2'
                    print 'exonString:',exonString
                    print 'cdsString:',cdsString
                    raise Exception("Impossible Event 1")
                    break
            utr_3.reverse()
        else:
            for (exon_start, exon_end) in exonList:
                if exon_start > cdsList[0][1]:
                    utr_5.append( (exon_start, exon_end) )
                elif (cdsList[0][1] < exon_end and exon_start == cdsList[0][0]) or (cdsList[0][1] < exon_end and len(cdsList) == 1):
                    utr_5.append( (cdsList[0][1]+1, exon_end) )
                    break
                elif exon_end == cdsList[0][1]:
                    break
                else:
                    print 'Impossible Event 3'
                    print 'exonString:',exonString
                    print 'cdsString:',cdsString
                    raise Exception("Impossible Event 1")
                    break
            for (exon_start, exon_end) in exonList[::-1]:
                if exon_end < cdsList[-1][0]:
                    utr_3.append( (exon_start, exon_end) )
                elif (exon_start < cdsList[-1][0] and exon_end == cdsList[-1][1]) or (exon_start < cdsList[-1][0] and len(cdsList) == 1):
                    utr_3.append( (exon_start, cdsList[-1][0]-1) )
                    break
                elif exon_start == cdsList[-1][0]:
                    break
                else:
                    print 'Impossible Event 4'
                    print 'exonString:',exonString
                    print 'cdsString:',cdsString
                    raise Exception("Impossible Event 1")
                    break
            utr_3.reverse()
        utr_5_string = ','.join([str(item[0])+'-'+str(item[1]) for item in utr_5])
        utr_3_string = ','.join([str(item[0])+'-'+str(item[1]) for item in utr_3])
        return utr_5_string, utr_3_string
    def __formatRNALine(self, rna_ID, pureTransID=False):
        """格式化一个转录本的信息
        测试
            print __formatRNALine('rna192', True) # plus strand mRMA
            print __formatRNALine('rna139', True) # plus strand mRMA
        """
        RNA = self.gff3_container['RNA'][rna_ID]
        "chromosome"
        chrID = RNA['chr']
        Chr = chrID
        if chrID.startswith('NC'):
            if 'chromosome' in self.gff3_container['region'][chrID]:
                Chr = 'chr' + self.gff3_container['region'][chrID]['chromosome']
            else:
                Chr = 'chrM'
        "Start, End, Strand, GeneID, TransID, GeneType"
        """
        try:
            GeneType = self.gff3_container['gene'][RNA['GeneID']]['gene_biotype']
        except KeyError:
            if 'ncrna_class' in RNA:
                GeneType = RNA['ncrna_class']
            else:
                GeneType = RNA['gbkey']
        """
        GeneType = self.__getRNAGeneType(rna_ID)
        """
        try:
            TransID = RNA['transcript_id']
        except KeyError:
            if RNA['gbkey'] in ('tRNA', 'rRNA'):
                TransID = RNA['gene']+'_'+RNA['gbkey']
            elif 'miRBase' in RNA:
                TransID = RNA['miRBase']
            else:
                return ""
        """
        TransID = self.rnaid_2_transID[rna_ID]
        if TransID == '': return ""
        (Start, End, Strand, GeneID, TransID, GeneType) = (RNA['start'], RNA['end'], RNA['strand'], RNA['gene']+'='+RNA['GeneID'], TransID, GeneType)
        if pureTransID: TransID = TransID.split('.')[0]
        "exon"
        exons = self.gff3_container['exon'][rna_ID]
        exonsList = [ exon['start']+'-'+exon['end'] for exon in exons ]
        "CDS"
        utrString = ''
        if rna_ID in self.gff3_container['CDS']:
            CDSs = self.gff3_container['CDS'][rna_ID]
            cdsList = [ cds['start']+'-'+cds['end'] for cds in CDSs ]
            cdsString = ','.join(cdsList)
            exonString = ','.join(exonsList)
            try:
                (utr_5_string, utr_3_string) = self.__utr_String(exonString, cdsString, Strand)
            except:
                print 'Skip this transcript: ',rna_ID
                print 'exonString:',exonString
                print 'cdsString:',cdsString
                print 'strand:',Strand
                print '\n'
                return ''
            utrString = (utr_5_string+','+utr_3_string).strip(',')
        formatString = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (Chr, Start, End, Strand, GeneID, TransID, GeneType, ','.join(exonsList))
        if utrString: formatString += '\t'+utrString
        return formatString
    def write_genomeCoor_bed(self, genomeCoorFileName, onlyChr=False, pureTransID=True):
        """ 把GFF3格式的注释转成基因组坐标的bed文件
            chr1    975205  981029  -       ENSG00000187642 ENST00000433179 protein_coding  978881-981029,976499-976624,975205-976269       975205-976171
        测试：
            write_genomeCoor_bed(hg19_gff3_container, genomeCoorFileName="~/hg19.bed", onlyChr=True, pureTransID=True)
            write_genomeCoor_bed(hg38_gff3_container, genomeCoorFileName="~/hg38.bed", onlyChr=True, pureTransID=True)
            write_genomeCoor_bed(mm10_gff3_container, genomeCoorFileName="~/mm10.bed", onlyChr=True, pureTransID=True)
        参数：
            onlyChr: 只输出染色体以chr开头的转录本
            pureTransID: 转录本的ID去除 .version
        """
        import random, os, commands
        tmpFileName = '/tmp/genomeCoor'+str(random.randint(1,1000))+'.bed'
        TMP = open(tmpFileName, 'w')
        for rna_ID in self.gff3_container['RNA']:
            TranscriptLine = self.__formatRNALine(rna_ID, pureTransID=pureTransID)
            if TranscriptLine == '': continue
            if onlyChr:
                if TranscriptLine.startswith('chr'):
                    print >>TMP, TranscriptLine
            else:
                print >>TMP, TranscriptLine
        TMP.close()
        return_code, output = commands.getstatusoutput( "cat %s | sort -k 1,1 -k2n,3n > %s" % (tmpFileName, genomeCoorFileName) )
        if return_code != 0:
            print 'Error: cat and sort step is broken. Please check the file: ', tmpFileName
        else:
            os.remove(tmpFileName)
            print 'Success!'
    def writeTranscriptome(self, genomeFileName, transcriptomeFileName, genomeChrSym='chr', pureTransID=False, verbose=True):
        """从注释文件个基因组文件中解析出转录组文件
        测试：
            # Gencode下载的基因组
            writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/hg38.fa", os.environ.get("HOME")+"/hg38.fa", hg38_gff3_container, genomeChrSym='chr', pureTransID=True)
            writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/mm10.fa", os.environ.get("HOME")+"/mm10.fa", mm10_gff3_container, genomeChrSym='chr', pureTransID=True)
            # NCBI 下载的基因组
            writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/ncbi/GRCh38.p7.fa", os.environ.get("HOME")+"/hg38_ncbi.fa", hg38_gff3_container, genomeChrSym='chr', pureTransID=True)
            writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/ncbi/GRCm38.p4.fa", os.environ.get("HOME")+"/mm10_ncbi.fa", mm10_gff3_container, genomeChrSym='chr', pureTransID=True)
        注意：这个文件可能和NCBI上直接下载的文件不一样，NCBI上的mRNA喜欢3'UTR加polyA
        """
        import seq as SeqFunc
        if genomeChrSym == 'chr':
            NCToChr = self.build_NC_To_chr_dict(pureNCID=False)
        OUT = open(transcriptomeFileName, 'w')
        Seq = SeqFunc.seqClass(genomeFileName)
        count = 0
        no_chr_list = []
        for rna_ID in self.gff3_container['RNA']:
            count += 1
            if count % 1000 == 0: print '\tProcessed %.2f%% ...' % (100.0*count/len(self.gff3_container['RNA']))
            "trancript information"
            GeneID = self.gff3_container['RNA'][rna_ID]['GeneID']
            GeneName = self.gff3_container['RNA'][rna_ID]['gene']
            GeneType = self.__getRNAGeneType(rna_ID)     #self.gff3_container['RNA'][rna_ID]['trans_type']
            RNA_info = GeneID+'|'+GeneName+'|'+GeneType
            "strand, chromosome, "
            strand = self.gff3_container['RNA'][rna_ID]['strand']
            Chr = self.gff3_container['RNA'][rna_ID]['chr']
            #if not Chr.startswith('NC'): continue
            if genomeChrSym == 'chr': 
                try:
                    Chr = NCToChr[Chr]
                except KeyError:
                    pass
            if not Seq.has(Chr):
                if Chr not in no_chr_list:
                    no_chr_list.append(Chr)
                    if verbose:
                        print >>sys.stderr, 'Warning: %s not in %s' % (Chr, genomeFileName)
                continue
            "Transcript ID"
            """
            try:
                TransID = self.gff3_container['RNA'][rna_ID]['transcript_id']
            except KeyError:
                if self.gff3_container['RNA'][rna_ID]['gbkey'] in ('tRNA', 'rRNA'):
                    TransID = self.gff3_container['RNA'][rna_ID]['gene']+'_'+self.gff3_container['RNA'][rna_ID]['gbkey']
                else:
                    continue
            """
            TransID = self.rnaid_2_transID[rna_ID]
            if pureTransID: TransID = TransID.split('.')[0]
            RNA_seq = ''
            exons = self.gff3_container['exon'][rna_ID]
            exonsList = [ [int(exon['start']), int(exon['end'])] for exon in exons ]
            for exon in exonsList:
                try:
                    RNA_seq += Seq.fetch(Chr, exon[0]-1, exon[1], strand)
                except KeyError:
                    if verbose:
                        print >>sys.stderr, 'Warning: KeyError -> %s\t%d\t%d\t%s' % (Chr, exon[0]-1, exon[1], strand)
                    continue
            print >>OUT, '>%s\t%s\n%s' % (TransID, RNA_info, SeqFunc.cutSeq(RNA_seq))
        OUT.close()
    """
    下面是static类型的方法
    """
    @staticmethod
    def norm_exons(string):
        """把 ExonString 基因组坐标归一化成转录组坐标
        测试：
            norm_exons( "173753-174414,172557-172688,169049-169264,168100-168165,165884-165942,164263-164791,158593-158674,155767-155831,146300-146339" )
        """
        tuples = list(); raw_tuple = list()
        arr = string.split(',')
        for indx in range(len(arr)):
            aarr = arr[indx].split('-')
            if indx == 0:  tuples.append( (1, abs(int(aarr[1])-int(aarr[0]))+1 ) )
            else: tuples.append( (tuples[-1][1]+1, abs(int(aarr[1])-int(aarr[0]))+tuples[-1][1]+1 ) )
            raw_tuple.append( (int(aarr[0]), int(aarr[1])) )
        return tuples, raw_tuple
    @staticmethod
    def norm_utr(exon_str, utr_str, strand):
        """把 UTRString 基因组坐标归一化成转录组坐标
        测试：
            norm_utr("859993-860328,861302-861393,865535-865716,866419-866469,871152-871276,874420-874509,874655-874840,876524-876686,877516-877631,877790-877868,877939-878438,878633-878757,879078-879954", "859993-860328,861302-861321,879534-879954", '+')
        """
        norm_tuple, raw_tuple = NCBI_Genome_Class.norm_exons(exon_str)
        utr_tuple = list()
        arr = utr_str.split(',')
        for indx in range(len(arr)):
            aarr = [ int(i) for i in arr[indx].split('-') ]
            flag = False
            for i in range(len(raw_tuple)):
                if strand == '-':
                    if raw_tuple[i][0] == aarr[0] <= aarr[1] == raw_tuple[i][1]:
                        utr_tuple.append( ( norm_tuple[i][0], norm_tuple[i][1] ) ); flag = True; break
                    elif raw_tuple[i][0] == aarr[0] <= aarr[1] < raw_tuple[i][1]:
                        utr_tuple.append( ( raw_tuple[i][1]-aarr[1]+norm_tuple[i][0], norm_tuple[i][1] ) ); flag = True; break
                    elif raw_tuple[i][0] < aarr[0] <= aarr[1] == raw_tuple[i][1]:
                        utr_tuple.append( ( norm_tuple[i][0], raw_tuple[i][1]-aarr[0]+norm_tuple[i][0] ) ); flag = True; break
                if strand == '+':
                    if raw_tuple[i][0] == aarr[0] <= aarr[1] == raw_tuple[i][1]:
                        utr_tuple.append( ( norm_tuple[i][0], norm_tuple[i][1] ) ); flag = True; break
                    elif raw_tuple[i][0] == aarr[0] <= aarr[1] < raw_tuple[i][1]:
                        utr_tuple.append( ( norm_tuple[i][0], norm_tuple[i][0]+aarr[1]-aarr[0] ) ); flag = True; break
                    elif raw_tuple[i][0] < aarr[0] <= aarr[1] == raw_tuple[i][1]:
                        utr_tuple.append( ( norm_tuple[i][1]-(aarr[1]-aarr[0]), norm_tuple[i][1] ) ); flag = True; break
            if not flag: print 'Unexpected Result'
        return norm_tuple, utr_tuple
    @staticmethod
    def genomeCoorBed_To_transCoorBed(genomeCoorFileName, transCoorFileName):
        """把基因组坐标的 bed 文件转成转录组坐标的 bed 文件
            ENST00000496112 ENSG00000106636 protein_coding  1270    1-272,273-355,356-456,457-561,562-627,628-1270  1-168,664-1270
        """
        IN = open(genomeCoorFileName)
        OUT = open(transCoorFileName, 'w')
        line = IN.readline()
        while line:
            arr = line.strip().split('\t')
            exonList = NCBI_Genome_Class.norm_exons(arr[7])[0]
            Len = exonList[-1][1]
            exonString = ','.join([str(item[0])+'-'+str(item[1]) for item in exonList])
            transCoorString = arr[5]+'\t'+arr[4]+'\t'+arr[6]+'\t'+str(Len)+'\t'+exonString
            if len(arr) == 9:
                utrList = NCBI_Genome_Class.norm_utr(arr[7], arr[8], arr[3])[1]
                utrString = ','.join([str(item[0])+'-'+str(item[1]) for item in utrList])
                transCoorString += '\t'+utrString
            print >>OUT, transCoorString
            line = IN.readline()
        OUT.close()

