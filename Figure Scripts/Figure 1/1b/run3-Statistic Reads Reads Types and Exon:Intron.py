

########################################################
#
#   Statistic Reads Reads Types and Exon/Intron
#
########################################################


from ParseTrans import *

mm10_parseTrans = ParseTransClass(genomeCoorBedFile = "/150T/zhangqf/GenomeAnnotation/Gencode/mm10.genomeCoor.bed")
hg38_parseTrans = ParseTransClass(genomeCoorBedFile = "/150T/zhangqf/GenomeAnnotation/Gencode/hg38.genomeCoor.bed")
#hg38_parseTrans = ParseTransClass(genomeCoorBedFile = "/150T/zhangqf/GenomeAnnotation/Gencode/hg38.genomeCoor.bed", removeVersion=True)

def anno_gene_mapped_reads(in_gene_file, in_trans_file, out_file, parseTrans):
    no_annotation_trans = 0
    
    IN_GENE = open(in_gene_file)
    IN_TRANS = open(in_trans_file)
    
    OUT = open(out_file, "w")
    gene_line = IN_GENE.readline()
    trans_line = IN_TRANS.readline()
    #gene_data = gene_line.strip().split()
    #trans_data = trans_line.strip().split()
    while gene_line and trans_line:
        rid_g, gid_gName, g_pos = gene_line.strip().split()
        rid_t, tid, t_pos  = trans_line.strip().split()
        gid, gName = gid_gName.split('|')
        if rid_g == rid_t:
            # Map to Exon
            try:
                gene_type = parseTrans.getTransFeature(tid, verbose=False)['gene_type']
                print >>OUT, "\t".join((rid_g, gid_gName, g_pos))+"\t||\t"+"\t".join((rid_t, tid, t_pos))+"\t"+gene_type
            except KeyError:
                no_annotation_trans += 1
            gene_line = IN_GENE.readline()
            trans_line = IN_TRANS.readline()
        elif rid_g < rid_t:
            # Map to Intron
            try:
                mapped_trans_list = []
                gene_coor = int(g_pos)
                gene_intron = parseTrans.getGeneIntron(gid)
                for trans_id in gene_intron:
                    intron_gene_coor = gene_intron[trans_id][1]
                    for single_intron in intron_gene_coor:
                        if single_intron[0] <= gene_coor < single_intron[1]-10:
                            mapped_trans_list.append( trans_id )
                            break
                total_gene_type = set()
                for trans_id in mapped_trans_list:
                    gene_type = parseTrans.getTransFeature(trans_id, verbose=False)['gene_type']
                    total_gene_type.add(gene_type)
                if total_gene_type:
                    print >>OUT, "\t".join((rid_g, gid_gName, g_pos))+"\t"+"|".join(total_gene_type)
            except KeyError:
                no_annotation_trans += 1
            gene_line = IN_GENE.readline()
        else:
            # cannot explained
            pass
            trans_line = IN_TRANS.readline()
    while gene_line:
        # Map to Intron
        rid_g, gid_gName, g_pos = gene_line.strip().split()
        gid, gName = gid_gName.split('|')
        try:
            mapped_trans_list = []
            gene_coor = int(g_pos)
            gene_intron = parseTrans.getGeneIntron(gid)
            for trans_id in gene_intron:
                intron_gene_coor = gene_intron[trans_id][1]
                for single_intron in intron_gene_coor:
                    if single_intron[0] <= gene_coor < single_intron[1]-10:
                        mapped_trans_list.append( trans_id )
                        break
            total_gene_type = set()
            for trans_id in mapped_trans_list:
                gene_type = parseTrans.getTransFeature(trans_id, verbose=False)['gene_type']
                total_gene_type.add(gene_type)
            if total_gene_type:
                print >>OUT, "\t".join((rid_g, gid_gName, g_pos))+"\t"+"|".join(total_gene_type)
        except KeyError:
            no_annotation_trans += 1
        gene_line = IN_GENE.readline()
    
    print "no_annotation_trans:",no_annotation_trans
    OUT.close()


root = "/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution2/"

######### mES

anno_gene_mapped_reads(root+"mes_ch_gene.tab", 
    root+"mes_ch_trans.tab",
    root+"mes_ch_trans.anno",
    mm10_parseTrans)

anno_gene_mapped_reads(root+"mes_np_gene.tab", 
    root+"mes_np_trans.tab",
    root+"mes_np_trans.anno",
    mm10_parseTrans)

anno_gene_mapped_reads(root+"mes_cy_gene.tab", 
    root+"mes_cy_trans.tab",
    root+"mes_cy_trans.anno",
    mm10_parseTrans)


######### HEK293

anno_gene_mapped_reads(root+"hek_ch_gene.tab", 
    root+"hek_ch_trans.tab",
    root+"hek_ch_trans.anno",
    hg38_parseTrans)

anno_gene_mapped_reads(root+"hek_np_gene.tab", 
    root+"hek_np_trans.tab",
    root+"hek_np_trans.anno",
    hg38_parseTrans)

anno_gene_mapped_reads(root+"hek_cy_gene.tab", 
    root+"hek_cy_trans.tab",
    root+"hek_cy_trans.anno",
    hg38_parseTrans)




