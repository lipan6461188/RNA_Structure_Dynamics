
import ParseTrans, tools, numpy, icSHAPE

hg38_parser = ParseTrans.ParseTransClass("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed")
mm10_parser = ParseTrans.ParseTransClass("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/mm10.genomeCoor.bed")

input_path = "/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/%s/shape.out"
hg_ch_vivo = tools.loadicSHAPE(input_path % ("hek_gene_ch_vivo", ))
hg_ch_vitro = tools.loadicSHAPE(input_path % ("hek_gene_ch_vitro", ))

mm_ch_vivo = tools.loadicSHAPE(input_path % ("mes_gene_ch_vivo", ))
mm_ch_vitro = tools.loadicSHAPE(input_path % ("mes_gene_ch_vitro", ))

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


hg_vivo_intron, hg_vivo_exon = calcIntronExonGINI(hg_ch_vivo, hg38_parser); print numpy.mean(hg_vivo_intron), numpy.mean(hg_vivo_exon)
hg_vitro_intron, hg_vitro_exon = calcIntronExonGINI(hg_ch_vitro, hg38_parser); print numpy.mean(hg_vitro_intron), numpy.mean(hg_vitro_exon)

mm_vivo_intron, mm_vivo_exon = calcIntronExonGINI(mm_ch_vivo, mm10_parser); print numpy.mean(mm_vivo_intron), numpy.mean(mm_vivo_exon)
mm_vitro_intron, mm_vitro_exon = calcIntronExonGINI(mm_ch_vitro, mm10_parser); print numpy.mean(mm_vitro_intron), numpy.mean(mm_vitro_exon)

tools.scipy.stats.mannwhitneyu(hg_vivo_intron, hg_vivo_exon, use_continuity=True, alternative='two-sided')
tools.scipy.stats.mannwhitneyu(hg_vitro_intron, hg_vitro_exon, use_continuity=True, alternative='two-sided')

tools.scipy.stats.mannwhitneyu(mm_vivo_intron, mm_vivo_exon, use_continuity=True, alternative='two-sided')
tools.scipy.stats.mannwhitneyu(mm_vitro_intron, mm_vitro_exon, use_continuity=True, alternative='two-sided')


writeFile("/tmp/hg.txt", hg_vivo_exon, hg_vivo_intron, hg_vitro_exon, hg_vitro_intron)
writeFile("/tmp/mm.txt", mm_vivo_exon, mm_vivo_intron, mm_vitro_exon, mm_vitro_intron)


"""  R to Plot

Mplot <- function(inFile, outFile, Label)
{
    gini_table <- read.csv(inFile, sep="\t", head=FALSE)
    gini_list <- list()

    for(idx in c(1:nrow(gini_table)))
    {
        gini_list[[ as.character(gini_table[idx,1]) ]] <- gini_table[idx,][!is.na(gini_table[idx,])]
    }


    for( item in names(gini_list) )
    {
        gini_list[[item]] <- as.numeric(gini_list[[ item ]])
        gini_list[[item]] <- gini_list[[item]][ !is.na(gini_list[[item]]) ]
    }

    red <- '#e72d27'
    blue <- '#2f52a2'

    # 保存
    pdf(outFile, width=5, height=5)
    par(mar=c(6, 4, 4, 2) + 0.1)
    vioplot(gini_list, col=rep(c(blue, red), 2), names=names(gini_list), las=2, main=Label, ylim=c(0, 1))
    dev.off()   
}

Mplot("/tmp/hg.txt", "~/figs/hg38_exon_intron.pdf", "Exon/Intro Structure")
Mplot("/tmp/mm.txt", "~/figs/mm10_exon_intron.pdf", "Exon/Intro Structure")

"""


############# Remove RBP region


def readGenomeRBPToGeneRBP(inFile, Parser):
    geneRBP = {}
    line_count = 0
    gid_set = set()
    for line in open(inFile):
        line_count += 1
        if line_count % 100000 == 0:
            print "\tlines ", line_count
        data = line.strip().split()
        chr_id,genome_start,genome_end,strand = data[0],int(data[1]),int(data[2]),data[3]
        try:
            geneRegList = Parser.genomeCoor2geneCoor(chr_id, genome_start, genome_end, strand)
        except:
            continue
        for genReg in geneRegList:
            gid,gstart,gend = genReg[3],int(genReg[4]),int(genReg[5])
            if gid not in gid_set:
                gid_set.add(gid)
                geneRBP[gid] = []
            geneRBP[gid].append([gstart, gend])
    return geneRBP


def combineRBPRegions(geneRBP):
    for gid in geneRBP:
        Cur_geneRBP = geneRBP[gid]
        Cur_geneRBP.sort(key=lambda x: x[0])
        idx = 0
        while idx < len(Cur_geneRBP)-1:
            if Cur_geneRBP[idx][0] <= Cur_geneRBP[idx+1][1] and Cur_geneRBP[idx+1][0] <= Cur_geneRBP[idx][1]:
                Cur_geneRBP[idx][1] = max(Cur_geneRBP[idx][1], Cur_geneRBP[idx+1][1])
                del Cur_geneRBP[idx+1]
                continue
            else:
                idx += 1

def calcIntronExonGINI_RBP(Shape, Parser, geneRBP):
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
        if gid in geneRBP:
            for t in intron_regions:
                for gidx in range(t[0], t[1]):
                    if Shape[gid][gidx] == 'NULL': continue
                    in_RBP = False
                    for rbp_reg in geneRBP[gid]:
                        if rbp_reg[0] <= gidx <= rbp_reg[1]:
                            in_RBP = True
                            break
                    if not in_RBP:
                        intron_shape_list.append( float(Shape[gid][gidx]) )
        
        exon_shape_list = []
        if gid in geneRBP:
            for t in exon_regions:
                for gidx in range(t[0], t[1]):
                    if Shape[gid][gidx] == 'NULL': continue
                    in_RBP = False
                    for rbp_reg in geneRBP[gid]:
                        if rbp_reg[0] <= gidx <= rbp_reg[1]:
                            in_RBP = True
                            break
                    if not in_RBP and Shape[gid][gidx] != 'NULL':
                        exon_shape_list.append( float(Shape[gid][gidx]) )
        
        gini = computeGINI_splice(intron_shape_list)
        if None not in gini:
            intron_gini += gini 
        
        gini = computeGINI_splice(exon_shape_list)
        if None not in gini:
            exon_gini += gini
    
    print "miss_count transcripts: ",miss_count 
    print "diff_len transcripts: ",diff_len 
    
    return intron_gini, exon_gini




geneRBP = readGenomeRBPToGeneRBP("/Share/home/zhangqf8/lipan/DYNAMIC/modification/HEK293/m6A/choose_good_protein/clipDB.merged.bed", hg38_parser)

icSHAPE.pbr_count(geneRBP)
combineRBPRegions(geneRBP); print icSHAPE.pbr_count(geneRBP)

hg_vivo_intron, hg_vivo_exon = calcIntronExonGINI_RBP(hg_ch_vivo, hg38_parser, geneRBP); print numpy.mean(hg_vivo_intron), numpy.mean(hg_vivo_exon)
hg_vitro_intron, hg_vitro_exon = calcIntronExonGINI_RBP(hg_ch_vitro, hg38_parser, geneRBP); print numpy.mean(hg_vitro_intron), numpy.mean(hg_vitro_exon)

tools.scipy.stats.mannwhitneyu(hg_vivo_intron, hg_vivo_exon, use_continuity=True, alternative='two-sided')
tools.scipy.stats.mannwhitneyu(hg_vitro_intron, hg_vitro_exon, use_continuity=True, alternative='two-sided')

writeFile("/tmp/hg.txt", hg_vivo_exon, hg_vivo_intron, hg_vitro_exon, hg_vitro_intron)



