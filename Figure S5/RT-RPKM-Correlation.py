
import tools

def readDMSORPKMLen(inFile):
    RPKMLen = {}
    lineCount = 0
    for line in open(inFile):
        if lineCount % 2 == 0:
            data = line.strip().split()
            Len = int(data[1])
            rpkm = tools.numpy.mean([ float(it) for it in data[2].split(',') ])
            RPKMLen[ data[0] ] = rpkm * Len
        lineCount += 1
    return RPKMLen

def readNAIRTTotal(inFile):
    NAITotal = {}
    lineCount = 0
    for line in open(inFile):
        if lineCount % 2 == 1:
            data = line.strip().split()
            Len = int(data[1])
            RTTotal = sum([ float(it) for it in data[3:] ])
            NAITotal[ data[0] ] = RTTotal
        lineCount += 1
    return NAITotal

def buildDF(DMSORPKMLen, NAIRTTotal):
    pairs = []
    for tid in set(DMSORPKMLen)&set(NAIRTTotal):
        rpkm_len = DMSORPKMLen[tid]
        rt_total = NAIRTTotal[tid]
        
        pairs.append( (tools.np.log2(rpkm_len), tools.np.log2(rt_total)) )
    
    df = tools.pd.DataFrame(data=pairs, columns=['log2_RPKMLen', 'log2_RTSum'])
    return df


hek_DMSORPKMLen = readDMSORPKMLen("/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/hek_cy_vivo/rt_D.txt")
hek_NAIRTTotal = readNAIRTTotal("/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/hek_cy_vivo/rt_N.txt")
hek_df = buildDF(hek_DMSORPKMLen, hek_NAIRTTotal)

tools.sns.jointplot(data=hek_df, x='log2_RPKMLen', y='log2_RTSum', color="#8172B2")
tools.plt.savefig("figs/hek_RT_RPKM.png")
tools.plt.show()

mes_DMSORPKMLen = readDMSORPKMLen("/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/mes_cy_vivo/rt_D.txt")
mes_NAIRTTotal = readNAIRTTotal("/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/mes_cy_vivo/rt_N.txt")
mes_df = buildDF(mes_DMSORPKMLen, mes_NAIRTTotal)

tools.sns.jointplot(data=mes_df, x='log2_RPKMLen', y='log2_RTSum', color="#8172B2")
tools.plt.savefig("figs/mes_RT_RPKM.png")
tools.plt.show()








