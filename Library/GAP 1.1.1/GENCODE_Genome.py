#-*- coding:utf-8 -*-


"""

有待完成的功能：
    1. 读取随机顺序的GTF文件
    2. 有序输出GTF文件
    3. NCBI/Gencode之间基因名字互转功能

Version:
    2017-5-29

"""

import sys, os

class No_GeneName_Error(Exception):
    pass

class GENCODE_Genome_Class(object):
    """这是一个操作GENCODE注释相关的类
    主要包括下面的功能：
        1. 读取Gencode的GFF3注释文件： 
            read_gencode_gtf(gtfFileName)
        2. 把原始的GTF转换成基因组坐标的Bed文件
            write_genomeCoor_bed(genomeCoorFileName, onlyChr=False, pureTransID=False, pureGeneID=False)
        3. 把基因组坐标的Bed文件转成转录组坐标的Bed文件
            import NCBI_Genome
            NCBI_Genome.NCBI_Genome_Class.genomeCoorBed_To_transCoorBed(genomeCoorFileName, transCoorFileName)
        4. 根据GTF注释从基因组中解析出转录组
            writeTranscriptome(genomeFileName, transcriptomeFileName, pureTransID=False, onlyChr=True, verbose=True)
    
    注意：
        Gencode官方下载的 ==转录组== 文件与这里用函数writeTranscriptome取出的转录本对应度非常好，可代替使用

    使用方法：
        from GENCODE_Genome import *

        # read gff3 annotation and build a GENCODEClass Object
        GENCODE_hg38 = GENCODE_Genome_Class('/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Gencode/GRCh38.gtf')
        GENCODE_hg19 = GENCODE_Genome_Class('/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Gencode/GRCh37.gtf')
        GENCODE_mm10 = GENCODE_Genome_Class('/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Gencode/GRCm38.gtf')
        GENCODE_mm9 = GENCODE_Genome_Class('/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Gencode/GRCm37.gtf')

        # convert gff3 transcripts infomation to genome-coordinated bed file
        GENCODE_hg38.write_genomeCoor_bed(genomeCoorFileName="~/hg38.bed", onlyChr=True, pureTransID=True, pureGeneID=True)
        GENCODE_hg19.write_genomeCoor_bed(genomeCoorFileName="~/hg19.bed", onlyChr=True, pureTransID=True, pureGeneID=True)
        GENCODE_mm10.write_genomeCoor_bed(genomeCoorFileName="~/mm10.bed", onlyChr=True, pureTransID=True, pureGeneID=True)
        GENCODE_mm9.write_genomeCoor_bed(genomeCoorFileName="~/mm9.bed", onlyChr=True, pureTransID=True, pureGeneID=True)


        # covert genome-coordinated bed file to transcript-coordinated bed file
        import NCBI_Genome
        NCBI_Genome.NCBI_Genome_Class.genomeCoorBed_To_transCoorBed(os.environ.get("HOME")+"/hg38.bed", os.environ.get("HOME")+"/hg38.transCoor.bed")
        NCBI_Genome.NCBI_Genome_Class.genomeCoorBed_To_transCoorBed(os.environ.get("HOME")+"/hg19.bed", os.environ.get("HOME")+"/hg19.transCoor.bed")
        NCBI_Genome.NCBI_Genome_Class.genomeCoorBed_To_transCoorBed(os.environ.get("HOME")+"/mm10.bed", os.environ.get("HOME")+"/mm10.transCoor.bed")
        NCBI_Genome.NCBI_Genome_Class.genomeCoorBed_To_transCoorBed(os.environ.get("HOME")+"/mm10.bed", os.environ.get("HOME")+"/mm10.transCoor.bed")


        # parse transcriptome from genome fasta file
        # Gencode下载的基因组
        GENCODE_hg38.writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/hg38.fa", os.environ.get("HOME")+"/hg38.fa", pureTransID=True, verbose=True)
        GENCODE_hg19.writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/hg19.fa", os.environ.get("HOME")+"/hg19.fa", pureTransID=True, verbose=True)
        GENCODE_mm10.writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/mm10.fa", os.environ.get("HOME")+"/mm10.fa", pureTransID=True, verbose=True)
        GENCODE_mm9.writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/mm9.fa", os.environ.get("HOME")+"/mm9.fa", pureTransID=True, verbose=True)
    """
    def __init__(self, gtfFileName):
        self.gtfFileName = gtfFileName
        self.gtf_container = self.read_gencode_gtf(gtfFileName)
    def __parse_attributes_gtf(self, attributesString):
        """ 解析attributes列
        测试：
            attri1 = 'gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; havana_gene "OTTHUMG00000000961"; havana_gene_version "2"; transcript_name "DDX11L1-002"; transcript_source "havana"; transcript_biotype "processed_transcript"; havana_transcript "OTTHUMT00000362751"; havana_transcript_version "1"; tag "basic"; transcript_support_level "1";'
            parse_attributes_gtf(attri1)
        """
        attributes = {}
        for Tuple in attributesString.strip(';').split('; '):
            try:
                Key_Value = Tuple.strip(' ').split(' ')
                Key = Key_Value[0]
                Value = ' '.join(Key_Value[1:])
                #(Key, Value) = Tuple.strip(' ').split(' ')
                attributes[ Key ] = Value.strip('"')
            except ValueError:
                print 'Warning: '+Tuple+' cannot be unpacked, skip it'
        return attributes
    def read_gencode_gtf(self, gtfFileName):
        """ 读取 GENCODE gff3 文件: UTR/start_codon/stop_codon都可以通过CDS推测出来，不读取
        测试：
            hg19_gtf_container = read_gencode_gtf('/Share/home/zhangqf8/lipan/DYNAMIC/GTF/GENCODE/GRCh37.p13.gff3')
            hg38_gtf_container = read_gencode_gtf('/Share/home/zhangqf8/lipan/DYNAMIC/GTF/GENCODE/GRCh38.p7.gff3')
            mm10_gtf_container = read_gencode_gtf('/Share/home/zhangqf8/lipan/DYNAMIC/GTF/GENCODE/GRCm38.p4.gff3')
        """
        print '\nRead GENCODE GTF Annotation: ', gtfFileName
        gtf_container = {'gene':{}, 'RNA':{}, 'UTR':{}, 'CDS':{}, 'exon':{}, 'start_codon':{}, 'stop_codon':{}}
        IN = open(gtfFileName)
        count = 0
        line = IN.readline()
        while line:
            count += 1
            if count % 100000 == 0: print '\tRead %d Lines...' % (count,)
            if line[0] == '#':
                line = IN.readline()
                continue
            (Chr, Source, Type, Chr_Start, Chr_End, Score, Strand, Phase, Attributes) = line.strip().split('\t')
            attributes = self.__parse_attributes_gtf(Attributes)
            attributes['chr'] = Chr; attributes['start'] = Chr_Start; attributes['end'] = Chr_End; attributes['strand'] = Strand
            if Type == 'gene':
                if 'gene_name' not in attributes:
                    attributes['gene_name'] = 'NULL'
                gtf_container[ 'gene' ][ attributes['gene_id'] ] = attributes
            if Type == 'transcript':
                if 'gene_name' not in attributes:
                    attributes['gene_name'] = 'NULL'
                if 'transcript_biotype' not in attributes and 'transcript_type' not in attributes:
                    attributes['trans_type'] = 'NULL'
                else:
                    if 'transcript_biotype' in attributes:
                        attributes['trans_type'] = attributes['transcript_biotype']
                        del attributes['transcript_biotype']
                    else:
                        attributes['trans_type'] = attributes['transcript_type']
                        del attributes['transcript_type']
                gtf_container[ 'RNA' ][ attributes['transcript_id'] ] = attributes
            if Type == 'exon':
                gtf_container[ 'exon' ][ attributes['transcript_id'] ] = gtf_container[ 'exon' ].get( attributes['transcript_id'], []) + [attributes]
            if Type == 'CDS':
                gtf_container[ 'CDS' ][ attributes['transcript_id'] ] = gtf_container[ 'CDS' ].get( attributes['transcript_id'], []) + [attributes]
            line = IN.readline()
        return gtf_container
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
    def __formatRNALine(self, rna_ID, pureTransID=False, pureGeneID=False):
        """格式化一个转录本的信息
        测试:
            print __formatRNALine(hg19_gtf_container, 'ENST00000585215.5_1', False, False) # plus strand mRMA
            print __formatRNALine(hg19_gtf_container, 'ENST00000488147.1_1', True, True) # plus strand mRMA
        """
        RNA = self.gtf_container['RNA'][rna_ID]
        "chromosome"
        Chr = RNA['chr']
        "Start, End, Strand, GeneID, TransID, GeneType"
        """
        if 'transcript_biotype' in RNA:
            GeneType = RNA['transcript_biotype']
        elif 'transcript_type' in RNA:
            GeneType = RNA['transcript_type']
        else:
            return ""
        """
        GeneType = RNA['trans_type']
        GeneID = RNA['gene_id'].split('.')[0] if pureGeneID else RNA['gene_id']
        TransID = rna_ID.split('.')[0] if pureTransID else rna_ID
        try:
            (Start, End, Strand, GeneID, TransID, GeneType) = (RNA['start'], RNA['end'], RNA['strand'], RNA['gene_name']+'='+GeneID, TransID, GeneType)
        except KeyError:
            print >>sys.stderr, "Error: gene_name attribute not in "+rna_ID
            raise No_GeneName_Error
        "exon"
        exons = self.gtf_container['exon'][rna_ID]
        exonsList = [ exon['start']+'-'+exon['end'] for exon in exons ]
        "CDS"
        utrString = ''
        if rna_ID in self.gtf_container['CDS']:
            CDSs = self.gtf_container['CDS'][rna_ID]
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
    def write_genomeCoor_bed(self, genomeCoorFileName, onlyChr=False, pureTransID=False, pureGeneID=False):
        """ 把GTF格式的注释转成基因组坐标的bed文件
            chr1    975205  981029  -       ENSG00000187642 ENST00000433179 protein_coding  978881-981029,976499-976624,975205-976269       975205-976171
        测试：
            write_genomeCoor_bed(hg19_gtf_container, genomeCoorFileName="~/hg19.bed", onlyChr=True, pureTransID=True)
            write_genomeCoor_bed(hg38_gtf_container, genomeCoorFileName="~/hg38.bed", onlyChr=True, pureTransID=True)
            write_genomeCoor_bed(mm10_gtf_container, genomeCoorFileName="~/mm10.bed", onlyChr=True, pureTransID=True)
        参数：
            onlyChr: 只输出染色体以chr开头的转录本
            pureTransID: 转录本的ID去除 .version
        """
        import random, os, commands
        tmpFileName = '/tmp/genomeCoor'+str(random.randint(1,1000))+'.bed'
        TMP = open(tmpFileName, 'w')
        for rna_ID in self.gtf_container['RNA']:
            try:
                TranscriptLine = self.__formatRNALine(rna_ID, pureTransID=pureTransID, pureGeneID=pureGeneID)
            except No_GeneName_Error:
                continue
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
    def writeTranscriptome(self, genomeFileName, transcriptomeFileName, pureTransID=False, onlyChr=True, verbose=True):
        """从注释文件个基因组文件中解析出转录组文件
        测试：
            # Gencode下载的基因组
            writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/hg38.fa", os.environ.get("HOME")+"/hg38.fa", hg38_gtf_container, pureTransID=True)
            writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/mm10.fa", os.environ.get("HOME")+"/mm10.fa", mm10_gtf_container, pureTransID=True)
            # GENCODE 下载的基因组
            writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/GENCODE/GRCh38.p7.fa", os.environ.get("HOME")+"/hg38_GENCODE.fa", hg38_gtf_container, pureTransID=True)
            writeTranscriptome(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/GENCODE/GRCm38.p4.fa", os.environ.get("HOME")+"/mm10_GENCODE.fa", mm10_gtf_container, pureTransID=True)
        注意：这个文件可能和GENCODE上直接下载的文件不一样，GENCODE上的mRNA喜欢3'UTR加polyA
        """
        import seq as SeqFunc
        OUT = open(transcriptomeFileName, 'w')
        Seq = SeqFunc.seqClass(genomeFileName)
        count = 0
        for rna_ID in self.gtf_container['RNA']:
            count += 1
            if count % 1000 == 0: print '\tProcessed %.2f%% ...' % (100.0*count/len(self.gtf_container['RNA']))
            "trancript information"
            GeneID = self.gtf_container['RNA'][rna_ID]['gene_id']
            GeneName = self.gtf_container['RNA'][rna_ID]['gene_name']
            GeneType = self.gtf_container['RNA'][rna_ID]['trans_type']
            RNA_info = GeneID+'|'+GeneName+'|'+GeneType
            "strand, chromosome, "
            strand = self.gtf_container['RNA'][rna_ID]['strand']
            Chr = self.gtf_container['RNA'][rna_ID]['chr']
            if onlyChr and not Chr.startswith('chr'): continue
            "Transcript ID"
            TransID = rna_ID.split('.')[0] if pureTransID else rna_ID
            RNA_seq = ''
            exons = self.gtf_container['exon'][rna_ID]
            exonsList = [ [int(exon['start']), int(exon['end'])] for exon in exons ]
            for exon in exonsList:
                try:
                    RNA_seq += Seq.fetch(Chr, exon[0]-1, exon[1], strand)
                except KeyError:
                    if verbose:
                        print 'Warning: KeyError -> %s\t%d\t%d\t%s' % (Chr, exon[0]-1, exon[1], strand)
                    continue
            print >>OUT, '>%s\t%s\n%s' % (TransID, RNA_info, SeqFunc.cutSeq(RNA_seq))
        OUT.close()
    def writeGene(self, genomeFileName, geneFileName, pureGeneID=False, onlyChr=True, verbose=True):
        """从注释文件个基因组文件中解析出基因文件（包含内含子）
        测试：
            # Gencode下载的基因组
            writeGene(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/hg38.fa", os.environ.get("HOME")+"/hg38.gene.fa", hg38_gtf_container, pureGeneID=True)
            writeGene(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/mm10.fa", os.environ.get("HOME")+"/mm10.gene.fa", mm10_gtf_container, pureGeneID=True)
            # GENCODE 下载的基因组
            writeGene(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/GENCODE/GRCh38.p7.fa", os.environ.get("HOME")+"/hg38_GENCODE.gene.fa", hg38_gtf_container, pureGeneID=True)
            writeGene(os.environ.get("HOME")+"/lipan/DYNAMIC/GTF/GENCODE/GRCm38.p4.fa", os.environ.get("HOME")+"/mm10_GENCODE.gene.fa", mm10_gtf_container, pureGeneID=True)
        注意：这个文件可能和GENCODE上直接下载的文件不一样，GENCODE上的mRNA喜欢3'UTR加polyA
        """
        import seq as SeqFunc
        OUT = open(geneFileName, 'w')
        Seq = SeqFunc.seqClass(genomeFileName)
        count = 0
        for gene_ID in self.gtf_container['gene']:
            count += 1
            if count % 1000 == 0: print '\tProcessed %.2f%% ...' % (100.0*count/len(self.gtf_container['gene']))
            "strand, chromosome, "
            strand = self.gtf_container['gene'][gene_ID]['strand']
            Chr = self.gtf_container['gene'][gene_ID]['chr']
            if onlyChr and not Chr.startswith('chr'): continue
            "Transcript ID"
            GeneID = gene_ID.split('.')[0] if pureGeneID else gene_ID
            start = int(self.gtf_container['gene'][gene_ID]['start'])
            end = int(self.gtf_container['gene'][gene_ID]['end'])
            gene_seq = Seq.fetch(Chr, start-1, end, strand)
            gene_name = self.gtf_container['gene'][gene_ID]['gene_name']
            print >>OUT, '>%s\t%s\t%d\t%s:%d-%d:%s\n%s' % (GeneID, gene_name, len(gene_seq), Chr, start, end, strand, SeqFunc.cutSeq(gene_seq))
        OUT.close()














