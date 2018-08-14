

def read_icShape(file_name):
    icDict = {}
    IN = open(file_name)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        trans_id = arr[0]
        try:
            trans_len = int(arr[1])
        except ValueError:
            print arr[:10]
            exit()
        trans_shape = arr[3:]
        icDict[trans_id] = trans_shape
        line = IN.readline()
    return icDict

def subShapeSet(shape, transList):
    subSet = {}
    for trans_id in transList:
        subSet[trans_id] = shape[trans_id]
    return subSet

def intersect(list_1, list_2):
    return [item for item in list_1 if item in list_2]

def noNULL(List):
    for item in List:
        if item == 'NULL':
            return False
    return True

def formatTbt(c1f1, c1f2, c1f3, c1f4, c2f1, c2f2, c2f3, c2f4, outFile):
    big_table = []
    trans_lenList = []
    for trans_id in c1f1:
        trans_lenList.append( [trans_id, len(c1f1[trans_id])] )
        for i in range(len(c1f1[trans_id])):
            c1f1_x = c1f1[trans_id][i]
            c1f2_x = c1f2[trans_id][i]
            c1f3_x = c1f3[trans_id][i]
            c1f4_x = c1f4[trans_id][i]
            c2f1_x = c2f1[trans_id][i]
            c2f2_x = c2f2[trans_id][i]
            c2f3_x = c2f3[trans_id][i]
            c2f4_x = c2f4[trans_id][i]
            if noNULL( (c1f1_x, c1f2_x, c1f3_x, c1f4_x, c2f1_x, c2f2_x, c2f3_x, c2f4_x) ):
                big_table.append( [trans_id, `i+1`, c1f1_x, c1f2_x, c1f3_x, c1f4_x, c2f1_x, c2f2_x, c2f3_x, c2f4_x] )
    trans_lenList.sort(key=lambda x: x[1])
    OUT = open(outFile, 'w')
    for trans_id, trans_len in trans_lenList:
        print >>OUT, '#'+trans_id+'\t'+str(trans_len)
    for arr in big_table:
        print >>OUT, '\t'.join(arr)
    OUT.close()


dataBase = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/mES/t100T1/'
out_dir = '/Share/home/zhangqf8/lipan/DYNAMIC/icSHAPE/10-11/mES/R_input/'

"""
pairs = ( ('ch','np'), ('ch','cy'), ('np','cy'), ('ch','chT'), ('np','npT'), ('cy','cyT'), ('chT','npT'), ('chT','cyT'), ('npT','cyT') )
pairs_other = ( ('cy','cyOT'), ('chT','cyOT'), ('npT','cyOT') )
"""


# ( ('ch','chT'), ('np','npT'), ('cy','cyT'), ('chT','npT'), ('chT','cyT'), ('npT','cyT') )
pairs = ( ('ch','np'), ('ch','cy'), ('np','cy'), ('ch','chT'), ('np','npT'), ('cy','cyT'), ('chT','npT'), ('chT','cyT'), ('npT','cyT') )

for comp1, comp2 in pairs:
    "condition_1"
    comp1_1 = read_icShape(dataBase+comp1+'1.out')
    comp1_2 = read_icShape(dataBase+comp1+'2.out')
    comp1_3 = read_icShape(dataBase+comp1+'3.out')
    comp1_4 = read_icShape(dataBase+comp1+'4.out')
    "condition_2"
    comp2_1 = read_icShape(dataBase+comp2+'1.out')
    comp2_2 = read_icShape(dataBase+comp2+'2.out')
    comp2_3 = read_icShape(dataBase+comp2+'3.out')
    comp2_4 = read_icShape(dataBase+comp2+'4.out')
    "common transcripts"
    common_trans_set = intersect( comp1_1, comp1_2 )
    common_trans_set = intersect( common_trans_set, comp1_3 )
    common_trans_set = intersect( common_trans_set, comp1_4 )
    common_trans_set = intersect( common_trans_set, comp2_1 )
    common_trans_set = intersect( common_trans_set, comp2_2 )
    common_trans_set = intersect( common_trans_set, comp2_3 )
    common_trans_set = intersect( common_trans_set, comp2_4 )
    "sub transcripts"
    fil_comp1_1 = subShapeSet(comp1_1, common_trans_set)
    fil_comp1_2 = subShapeSet(comp1_2, common_trans_set)
    fil_comp1_3 = subShapeSet(comp1_3, common_trans_set)
    fil_comp1_4 = subShapeSet(comp1_4, common_trans_set)
    "sub transcripts"
    fil_comp2_1 = subShapeSet(comp2_1, common_trans_set)
    fil_comp2_2 = subShapeSet(comp2_2, common_trans_set)
    fil_comp2_3 = subShapeSet(comp2_3, common_trans_set)
    fil_comp2_4 = subShapeSet(comp2_4, common_trans_set)
    "output"
    formatTbt(fil_comp1_1, fil_comp1_2, fil_comp1_3, fil_comp1_4, fil_comp2_1, fil_comp2_2, fil_comp2_3, fil_comp2_4, out_dir+comp1+"-"+comp2+".txt")






