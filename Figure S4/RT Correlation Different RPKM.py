

import tools
from tools import *

def calc_pearson(RT_1, RT_2, min_RT=2):
    pearson = []
    for tid in set(RT_1) & set(RT_2):
        if tools.numpy.mean(RT_1[tid]) < min_RT or tools.numpy.mean(RT_2[tid]) < min_RT:
            continue
        r,p = tools.scipy.stats.pearsonr(RT_1[tid], RT_2[tid])
        if p<0.01:
            pearson.append(r)
    
    pearson += [-1, 1]
    return pearson

def calc_pearson_dict(RT_1, RT_2, min_RT=2):
    pearson = {}
    for tid in set(RT_1) & set(RT_2):
        if tools.numpy.mean(RT_1[tid]) < min_RT or tools.numpy.mean(RT_2[tid]) < min_RT:
            continue
        r,p = tools.scipy.stats.pearsonr(RT_1[tid], RT_2[tid])
        if p<0.01:
            pearson[tid] = r
    
    return pearson


def read_RT(inFile):
    import numpy
    
    Data = {}
    
    last_tid = ""
    for line in open(inFile):
        if line[0] != '#':
            data = line.strip().split()
            tid, length, rpkm = data[0], data[1], data[2]
            #if min_rpkm > float(rpkm) or float(rpkm) > max_rpkm:
            #    continue
            if tid == last_tid:
                rt_list = [ float(v) for v in data[3:] ]
                Data[tid] = ( float(rpkm), numpy.mean(rt_list), rt_list )
    
            last_tid = tid
    
    return Data

#####################
## Count base with more RT correlation
#####################


############# 1. read RT

IN_ROOT = "/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/"
human_ch_N1 = read_RT(IN_ROOT+"hek_ch_vivo/rt_N1.txt")
human_ch_N2 = read_RT(IN_ROOT+"hek_ch_vivo/rt_N2.txt")
human_np_N1 = read_RT(IN_ROOT+"hek_np_vivo/rt_N1.txt")
human_np_N2 = read_RT(IN_ROOT+"hek_np_vivo/rt_N2.txt")
human_cy_N1 = read_RT(IN_ROOT+"hek_cy_vivo/rt_N1.txt")
human_cy_N2 = read_RT(IN_ROOT+"hek_cy_vivo/rt_N2.txt")

mouse_ch_N1 = read_RT(IN_ROOT+"mes_ch_vivo/rt_N1.txt")
mouse_ch_N2 = read_RT(IN_ROOT+"mes_ch_vivo/rt_N2.txt")
mouse_np_N1 = read_RT(IN_ROOT+"mes_np_vivo/rt_N1.txt")
mouse_np_N2 = read_RT(IN_ROOT+"mes_np_vivo/rt_N2.txt")
mouse_cy_N1 = read_RT(IN_ROOT+"mes_cy_vivo/rt_N1.txt")
mouse_cy_N2 = read_RT(IN_ROOT+"mes_cy_vivo/rt_N2.txt")


#####################
## Filter valid transcript
#####################

def get_valid_trans(icSHAPE):
    valid_trans = []
    for tid in icSHAPE:
        count = len(icSHAPE[tid]) - icSHAPE[tid].count('NULL')
        if count > 20 and 1.0*count/len(icSHAPE[tid])>0.3:
            valid_trans.append(tid)
    return valid_trans

#####################
## Accumulative plot of different RPKM
#####################

meanRT = 2

#### Human

plt.figure(figsize=(10,4))
index = 1
for N1,N2 in [(human_ch_N1, human_ch_N2),(human_np_N1, human_np_N2),(human_cy_N1, human_cy_N2)]:
    
    plt.subplot(1,3,index)
    cutoff = [10, 20, 40, 100, 200]
    
    tids_1 = [ tid for tid in set(N1)&set(N2) if (cutoff[0]<N1[tid][0]) and (N1[tid][1]>=meanRT) and (cutoff[0]<N2[tid][0]) and (N2[tid][1]>=meanRT) ]; print len(tids_1)
    tids_2 = [ tid for tid in set(N1)&set(N2) if (cutoff[1]<N1[tid][0]) and (N1[tid][1]>=meanRT) and (cutoff[1]<N2[tid][0]) and (N2[tid][1]>=meanRT) ]; print len(tids_2)
    tids_3 = [ tid for tid in set(N1)&set(N2) if (cutoff[2]<N1[tid][0]) and (N1[tid][1]>=meanRT) and (cutoff[2]<N2[tid][0]) and (N2[tid][1]>=meanRT) ]; print len(tids_3)
    tids_4 = [ tid for tid in set(N1)&set(N2) if (cutoff[3]<N1[tid][0]) and (N1[tid][1]>=meanRT) and (cutoff[3]<N2[tid][0]) and (N2[tid][1]>=meanRT) ]; print len(tids_4)
    tids_5 = [ tid for tid in set(N1)&set(N2) if (cutoff[4]<N1[tid][0]) and (N1[tid][1]>=meanRT) and (cutoff[4]<N2[tid][0]) and (N2[tid][1]>=meanRT) ]; print len(tids_5)
    
    pearson_1 = calc_pearson({tid:N1[tid][2] for tid in tids_1}, {tid:N2[tid][2] for tid in tids_1}, min_RT=meanRT); print len(pearson_1)
    pearson_2 = calc_pearson({tid:N1[tid][2] for tid in tids_2}, {tid:N2[tid][2] for tid in tids_2}, min_RT=meanRT); print len(pearson_2)
    pearson_3 = calc_pearson({tid:N1[tid][2] for tid in tids_3}, {tid:N2[tid][2] for tid in tids_3}, min_RT=meanRT); print len(pearson_3)
    pearson_4 = calc_pearson({tid:N1[tid][2] for tid in tids_4}, {tid:N2[tid][2] for tid in tids_4}, min_RT=meanRT); print len(pearson_4)
    pearson_5 = calc_pearson({tid:N1[tid][2] for tid in tids_5}, {tid:N2[tid][2] for tid in tids_5}, min_RT=meanRT); print len(pearson_5)
    
    tools.cdf(pearson_1, color='#4C72B0', topdown=True, label=">%s (%s)" % (cutoff[0], len(pearson_1)))
    tools.cdf(pearson_2, color='#55A868', topdown=True, label=">%s (%s)" % (cutoff[1], len(pearson_2)))
    tools.cdf(pearson_3, color='#C44E52', topdown=True, label=">%s (%s)" % (cutoff[2], len(pearson_3)))
    tools.cdf(pearson_4, color='#8172B2', topdown=True, label=">%s (%s)" % (cutoff[3], len(pearson_4)))
    tools.cdf(pearson_5, color='#CCB974', topdown=True, label=">%s (%s)" % (cutoff[4], len(pearson_5)))
        
    plt.legend()
    index += 1

plt.savefig("figs/human_RPKM_2.pdf")
plt.show()



#### Mouse

meanRT = 2

plt.figure(figsize=(10,4))
index = 1
for N1,N2 in [(mouse_ch_N1, mouse_ch_N2),(mouse_np_N1, mouse_np_N2),(mouse_cy_N1, mouse_cy_N2)]:
    
    plt.subplot(1,3,index)
    cutoff = [10, 20, 40, 100, 200]
    
    tids_1 = [ tid for tid in set(N1)&set(N2) if (cutoff[0]<N1[tid][0]) and (N1[tid][1]>=meanRT) and (cutoff[0]<N2[tid][0]) and (N2[tid][1]>=meanRT) ]; print len(tids_1)
    tids_2 = [ tid for tid in set(N1)&set(N2) if (cutoff[1]<N1[tid][0]) and (N1[tid][1]>=meanRT) and (cutoff[1]<N2[tid][0]) and (N2[tid][1]>=meanRT) ]; print len(tids_2)
    tids_3 = [ tid for tid in set(N1)&set(N2) if (cutoff[2]<N1[tid][0]) and (N1[tid][1]>=meanRT) and (cutoff[2]<N2[tid][0]) and (N2[tid][1]>=meanRT) ]; print len(tids_3)
    tids_4 = [ tid for tid in set(N1)&set(N2) if (cutoff[3]<N1[tid][0]) and (N1[tid][1]>=meanRT) and (cutoff[3]<N2[tid][0]) and (N2[tid][1]>=meanRT) ]; print len(tids_4)
    tids_5 = [ tid for tid in set(N1)&set(N2) if (cutoff[4]<N1[tid][0]) and (N1[tid][1]>=meanRT) and (cutoff[4]<N2[tid][0]) and (N2[tid][1]>=meanRT) ]; print len(tids_5)
    
    pearson_1 = calc_pearson({tid:N1[tid][2] for tid in tids_1}, {tid:N2[tid][2] for tid in tids_1}, min_RT=meanRT); print len(pearson_1)
    pearson_2 = calc_pearson({tid:N1[tid][2] for tid in tids_2}, {tid:N2[tid][2] for tid in tids_2}, min_RT=meanRT); print len(pearson_2)
    pearson_3 = calc_pearson({tid:N1[tid][2] for tid in tids_3}, {tid:N2[tid][2] for tid in tids_3}, min_RT=meanRT); print len(pearson_3)
    pearson_4 = calc_pearson({tid:N1[tid][2] for tid in tids_4}, {tid:N2[tid][2] for tid in tids_4}, min_RT=meanRT); print len(pearson_4)
    pearson_5 = calc_pearson({tid:N1[tid][2] for tid in tids_5}, {tid:N2[tid][2] for tid in tids_5}, min_RT=meanRT); print len(pearson_5)
    
    tools.cdf(pearson_1, color='#4C72B0', topdown=True, label=">%s (%s)" % (cutoff[0], len(pearson_1)))
    tools.cdf(pearson_2, color='#55A868', topdown=True, label=">%s (%s)" % (cutoff[1], len(pearson_2)))
    tools.cdf(pearson_3, color='#C44E52', topdown=True, label=">%s (%s)" % (cutoff[2], len(pearson_3)))
    tools.cdf(pearson_4, color='#8172B2', topdown=True, label=">%s (%s)" % (cutoff[3], len(pearson_4)))
    tools.cdf(pearson_5, color='#CCB974', topdown=True, label=">%s (%s)" % (cutoff[4], len(pearson_5)))
    
    plt.legend()
    index += 1

plt.savefig("figs/mouse_RPKM.pdf")
plt.show()




