
### Statistic Ratio of reads cover to all types of RNA

from seq import *
gene_type = anno_Methods.gene_type

def statistic_reads_distribution(input_file_name, output_file_name):
    exon_count = {}
    intron_count = {}
    total = 0
    #IN = open(input_file_name)
    #line = IN.readline()
    lineCount = 0
    for line in open(input_file_name):
        data = line.strip().split()
        if len(data) < 3:
            print "Error", line
            continue
        if data[1] not in ('ENSG00000278996','ENSG00000280441'):
            if len(data) == 4:
                gene_type_list = data[3].split('|')
                strength = 1.0/len(gene_type_list)
                for cur_gene_type in gene_type_list:
                    trans_type = gene_type(cur_gene_type)
                    intron_count[ trans_type ] = intron_count.get(trans_type, 0) + strength
            elif len(data) == 8:
                trans_type = gene_type(data[7])
                exon_count[ trans_type ] = exon_count.get(trans_type, 0) + 1
            else:
                print "Error", line
                continue
            total += 1
        #line = IN.readline()
        lineCount += 1
        if lineCount % 1000000 == 0:
            print "Read ", lineCount, "..."
    OUT = open(output_file_name, "w")
    for trans_type in exon_count:
        print >>OUT, "EXON\t%s\t%s\t%.3f" % (trans_type, exon_count[trans_type], 1.0*exon_count[trans_type]/total)
    for trans_type in intron_count:
        print >>OUT, "INTRON\t%s\t%s\t%.3f" % (trans_type, intron_count[trans_type], 1.0*intron_count[trans_type]/total)
    OUT.close()
    return exon_count, intron_count



# ======》 mES 《=======

exon_count, intron_count = statistic_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/mes_ch_trans.anno", 
    "/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/mes_ch_trans.count")

exon_count, intron_count = statistic_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/mes_np_trans.anno", 
    "/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/mes_np_trans.count")

exon_count, intron_count = statistic_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/mes_cy_trans.anno", 
    "/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/mes_cy_trans.count")



# ======》 HEK293 《=======

exon_count, intron_count = statistic_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/hek_ch_trans.anno", 
    "/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/hek_ch_trans.count")

exon_count, intron_count = statistic_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/hek_np_trans.anno", 
    "/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/hek_np_trans.count")

exon_count, intron_count = statistic_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/hek_cy_trans.anno", 
    "/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/hek_cy_trans.count")



############################
### Update 2018-10-1
############################

# ======》 mES 《=======

exon_count, intron_count = statistic_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/mes_ch_trans.anno", 
    "/tmp/new/mes_ch_trans.count")

exon_count, intron_count = statistic_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/mes_np_trans.anno", 
    "/tmp/new/mes_np_trans.count")

exon_count, intron_count = statistic_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/mes_cy_trans.anno", 
    "/tmp/new/mes_cy_trans.count")



# ======》 HEK293 《=======

exon_count, intron_count = statistic_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/hek_ch_trans.anno", 
    "/tmp/new/hek_ch_trans.count")

exon_count, intron_count = statistic_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/hek_np_trans.anno", 
    "/tmp/new/hek_np_trans.count")

exon_count, intron_count = statistic_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/hek_cy_trans.anno", 
    "/tmp/new/hek_cy_trans.count")















