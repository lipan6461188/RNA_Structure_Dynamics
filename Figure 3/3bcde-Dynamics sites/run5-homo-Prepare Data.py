
"""
找出human和mouse中同源的碱基区域
"""

from maf import *

def biSearch(item, Set):
    # + 二分查找
    # + Set必需事先从小到大排序
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

def read_icSHAPE_File(fileName, transLenDict):
    " 读取单个的 icSHAPE 文件 "
    shapeDict = {}
    IN = open(fileName)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        shapeDict[arr[0]] = arr[3:]
        if arr[0] not in transLenDict:
            transLenDict[arr[0]] = len( arr[3:] )
        line = IN.readline()
    return { 'keys': sorted(shapeDict.keys()), 'shape': shapeDict }

def read_big_icSHAPE():
    """
    读取目录下所有的icSHAPE文件：
        HEK293: /Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/3-24/HEK293/t100T1
        mES: /Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/3-24/mES_new/t100T1
    """
    HEK293_ROOT = "/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/HEK293/t100T1/"
    mES_ROOT = "/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/mES/t100T1/"
    icSHAPE = {'HEK293': {}, 'mES':{}}
    transLenDict = {}
    for comp_Prefix in ('ch', 'np', 'cy'):
        for repeat in range(1,5):
            HEK293_fileName = HEK293_ROOT + comp_Prefix + str(repeat) + '.out'
            mm10_fileName = mES_ROOT + comp_Prefix + str(repeat) + '.out'
            print 'Now Read '+ HEK293_fileName + '...'
            icSHAPE[ 'HEK293' ][ comp_Prefix+str(repeat) ] = read_icSHAPE_File(HEK293_fileName, transLenDict)
            print 'Now Read '+ mm10_fileName + '...'
            icSHAPE[ 'mES' ][ comp_Prefix+str(repeat) ] = read_icSHAPE_File(mm10_fileName, transLenDict)
    return icSHAPE, transLenDict

def in_icSHAPE(transID, sp_icSHAPE):
    " 检查转录本的 ID 是否在一个icSHAPE的集合中 "
    inComp = []
    for comp_Prefix in ('ch', 'np', 'cy'):
        in_comp = True
        for repeat in ['1', '2', '3', '4']:
            if not biSearch(transID, sp_icSHAPE[comp_Prefix+repeat]['keys']):
                in_comp = False
        inComp.append( in_comp )
    return inComp

def noNULL(List):
    for item in List:
        if item == 'NULL':
            return False
    return True

def getBaseShape_with_sm(homo, icSHAPE, hek_TransID, mes_TransID, comp):
    """
    通过 ctMaf_Homology_Methods.get_slave_base_to_master_base 调用后得到的结果取出对应碱基的icSHAPE值
    """
    big_shapeList = []
    for basePair in homo:
        hek_coor = basePair[0] - 1
        mes_coor = basePair[1] - 1
        align = basePair[2].split('-')
        #if align[0] != align[1]:
        #    continue
        hek_comp1 = icSHAPE['HEK293'][comp+'1']['shape'][hek_TransID][hek_coor]
        hek_comp2 = icSHAPE['HEK293'][comp+'2']['shape'][hek_TransID][hek_coor]
        hek_comp3 = icSHAPE['HEK293'][comp+'3']['shape'][hek_TransID][hek_coor]
        hek_comp4 = icSHAPE['HEK293'][comp+'4']['shape'][hek_TransID][hek_coor]
        mes_comp1 = icSHAPE['mES'][comp+'1']['shape'][mes_TransID][mes_coor]
        mes_comp2 = icSHAPE['mES'][comp+'2']['shape'][mes_TransID][mes_coor]
        mes_comp3 = icSHAPE['mES'][comp+'3']['shape'][mes_TransID][mes_coor]
        mes_comp4 = icSHAPE['mES'][comp+'4']['shape'][mes_TransID][mes_coor]
        shapeList = [hek_comp1, hek_comp2, hek_comp3, hek_comp4, mes_comp1, mes_comp2, mes_comp3, mes_comp4]
        if noNULL( shapeList ):
            big_shapeList.append( [hek_TransID, str(hek_coor+1), mes_TransID, str(mes_coor+1), basePair[2]] + shapeList )
    return big_shapeList

def gene_type(raw_type):
    "Convert Raw Gene Type to Our Defined Gene Type"
    valid_gene_type = ('pseudogene', 'snoRNA', 'snRNA', 'miRNA', 'misc_RNA', 'rRNA')
    #lncRNA_class = ('antisense','lincRNA','processed_transcript','sense_intronic','TEC','sense_overlapping')
    lncRNA_class = ('3prime_overlapping_ncrna','antisense','lincRNA','non_coding','sense_intronic','sense_overlapping','processed_transcript')
    if raw_type in valid_gene_type: return raw_type;
    if re.match('.*pseudogene',raw_type): return 'pseudogene';
    if raw_type == 'protein_coding': return 'mRNA';
    if raw_type in lncRNA_class: return 'lncRNA';
    return 'other'

def loadGeneType(file_name):
    "加载基因的类型"
    typeDict = {}
    IN = open(file_name)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        typeDict[ arr[5] ] = gene_type(arr[6])
        line = IN.readline()
    return typeDict

HEK293_refbed = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Homo_sapiens.trans.bed'
mES_refbed = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Mus_musculus.trans.bed'

HEK293_transType = loadGeneType(HEK293_refbed)
mES_transType = loadGeneType(mES_refbed)

ctMaf_FileName = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/human-mouse/hg38_mm10.ctmaf'
ctMaf_Container = ctMaf_Homology_Methods.read_ctMaf_master_slave(ctMaf_FileName)

icSHAPE, transLenDict = read_big_icSHAPE()


def giveTransID(List, index):
    trueList = []
    falseList = []
    for item in List:
        if item[index] == True:
            trueList.append( item )
        else:
            falseList.append( item )
    trueList.sort(key=lambda x: x[1], reverse=True)
    if len(trueList) >= 1:
        return trueList[0][0]
    else:
        return False




CH_OUT = open('/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/homo/hg38_mm10_ch.txt', 'w')
NP_OUT = open('/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/homo/hg38_mm10_np.txt', 'w')
CY_OUT = open('/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/homo/hg38_mm10_cy.txt', 'w')


for sp1_transID in ctMaf_Container:
    if sp1_transID not in HEK293_transType: continue
    (in_hek_ch, in_hek_np, in_hek_cy) = in_icSHAPE(sp1_transID, icSHAPE['HEK293'])
    if in_hek_ch or in_hek_np or in_hek_cy:
        transLen = transLenDict[sp1_transID]
        sp2_homo = ctMaf_Homology_Methods.get_slave_base_to_master_base(ctMaf_Container, sp1_transID, 1, transLen)
        PairTransSelectList = []
        for sp2_transID in sp2_homo:
            if sp2_transID not in mES_transType: continue
            if HEK293_transType[sp1_transID] == mES_transType[sp2_transID]:
                (in_mes_ch, in_mes_np, in_mes_cy) = in_icSHAPE(sp2_transID, icSHAPE['mES'])
                PairTransSelectList.append( [ sp2_transID, len(sp2_homo[sp2_transID]), in_mes_ch, in_mes_np, in_mes_cy ] )
        " 处理sp_2转录本 "
        ch_sp2_transID = giveTransID(PairTransSelectList, index=2)
        np_sp2_transID = giveTransID(PairTransSelectList, index=3)
        cy_sp2_transID = giveTransID(PairTransSelectList, index=4)
        if ch_sp2_transID and in_hek_ch:
            ch_big_shapeList = getBaseShape_with_sm(sp2_homo[ch_sp2_transID], icSHAPE, sp1_transID, ch_sp2_transID, 'ch')
            for item in ch_big_shapeList:
                print >>CH_OUT, '\t'.join(item)
        if np_sp2_transID and in_hek_np:
            np_big_shapeList = getBaseShape_with_sm(sp2_homo[np_sp2_transID], icSHAPE, sp1_transID, np_sp2_transID, 'np')
            for item in np_big_shapeList:
                print >>NP_OUT, '\t'.join(item)
        if cy_sp2_transID and in_hek_cy:
            cy_big_shapeList = getBaseShape_with_sm(sp2_homo[cy_sp2_transID], icSHAPE, sp1_transID, cy_sp2_transID, 'cy')
            for item in cy_big_shapeList:
                print >>CY_OUT, '\t'.join(item)


CH_OUT.close()
NP_OUT.close()
CY_OUT.close()






























