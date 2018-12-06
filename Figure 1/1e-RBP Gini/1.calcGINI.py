
import ParseTrans, tools, numpy, icSHAPE


def parser_intron_exon(geneExons, geneLength, Parser):
    
    # 1. collect
    exons = []
    for tid in geneExons:
        #if Parser.getTransFeature(tid)['gene_type'] in ('mRNA', 'protein_coding'):
        exons += geneExons[tid][1]
    exons.sort(key=lambda x: x[0])
    if len(exons) == 0: return None, None
    
    # 2. combine
    new_exons = [ exons[0] ]
    for i in range(1, len(exons)):
        if new_exons[-1][1] >= exons[i][0]:
            new_exons[-1][1] = max(new_exons[-1][1], exons[i][1])
        else:
            new_exons.append(exons[i])
    
    # 3. get intron
    introns = []
    if new_exons[0][0] != 1:
        print geneExons
        raise Exception("Exons is not started from 1")
    for i in range(1, len(new_exons)):
        introns.append( (new_exons[i-1][1]+1, new_exons[i][0]-1) )
    if new_exons[-1][1] != geneLength:
        print geneExons
        raise Exception("Exons is not ended in geneLength")
    
    # 4. return
    return introns, new_exons

def computeGINI_splice(my_list, valid_cutoff=10, splice=20):
    GINI = []
    idx = 0
    while len(my_list)-idx > splice/2:
        gini = tools.calcGINI(my_list[idx:idx+splice], valid_cutoff)
        if gini != -1:
            GINI.append(gini)
        idx += splice
    return GINI

def calcIntronExonGINI(Shape, Parser):
    geneInfo = Parser.getGeneInfo()
    intron_gini = []
    exon_gini = []
    
    miss_count = 0
    diff_len = 0
    for gid in Shape:
        try:
            geneExons = Parser.getGeneExon(gid)
        except KeyError:
            miss_count += 1
            continue
        if geneInfo[gid]['length'] != len(Shape[gid]):
            diff_len += 1
            continue
        
        if ('mRNA' not in geneInfo[gid]['gene_type']) and ('protein_coding' not in geneInfo[gid]['gene_type']):
            continue
        
        intron_regions, exon_regions = parser_intron_exon(geneExons, geneInfo[gid]['length'], Parser)
        if intron_regions == None: continue
        
        intron_shape_list = []
        for t in intron_regions:
            intron_shape_list += [ float(it) for it in Shape[gid][t[0]:t[1]] if it != 'NULL' ]
        
        exon_shape_list = []
        for t in exon_regions:
            exon_shape_list += [ float(it) for it in Shape[gid][t[0]:t[1]] if it != 'NULL' ]
        
        gini = computeGINI_splice(intron_shape_list)
        if None not in gini:
            intron_gini += gini 
        
        gini = computeGINI_splice(exon_shape_list)
        if None not in gini:
            exon_gini += gini
    
    print "miss_count transcripts: ",miss_count 
    print "diff_len transcripts: ",diff_len 
    
    return intron_gini, exon_gini

def writeFile(fileName, vivo_exon, vivo_intron, vitro_exon, vitro_intron):
    OUT = open(fileName, 'w')
    print >>OUT, "vivo_exon\t"+"\t".join([str(item) for item in vivo_exon])
    print >>OUT, "vivo_intron\t"+"\t".join([str(item) for item in vivo_intron])
    print >>OUT, "vitro_exon\t"+"\t".join([str(item) for item in vitro_exon])
    print >>OUT, "vitro_intron\t"+"\t".join([str(item) for item in vitro_intron])
    OUT.close()


hg38_parser = ParseTrans.ParseTransClass("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed")

input_path = "/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/%s/shape.out"
hg_ch_vivo = tools.loadicSHAPE(input_path % ("hek_gene_ch_vivo", ))
hg_ch_vitro = tools.loadicSHAPE(input_path % ("hek_gene_ch_vitro", ))

hg_vivo_intron, hg_vivo_exon = calcIntronExonGINI(hg_ch_vivo, hg38_parser); print numpy.mean(hg_vivo_intron), numpy.mean(hg_vivo_exon)
hg_vitro_intron, hg_vitro_exon = calcIntronExonGINI(hg_ch_vitro, hg38_parser); print numpy.mean(hg_vitro_intron), numpy.mean(hg_vitro_exon)

writeFile("/tmp/hg_gini.txt", hg_vivo_exon, hg_vivo_intron, hg_vitro_exon, hg_vitro_intron)



