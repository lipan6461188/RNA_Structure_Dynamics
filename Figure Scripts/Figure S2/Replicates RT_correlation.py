

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


def read_RT(inFile):
    RT = {}
    last_tid = ""
    for line in open(inFile):
        if line[0] != '#':
            data = line.strip().split()
            tid, length, rpkm = data[0], data[1], data[2]
            if tid == last_tid:
                RT[tid] = [ float(v) for v in data[3:] ]
            last_tid = tid
    
    return RT

def cdf(data_list, color='red', topdown=False):
    import re
    import numpy as np
    import matplotlib.pyplot as plt
    
    data_list = np.sort(data_list)
    if topdown:
        p = 1 - 1.0 * np.arange(len(data_list))/(len(data_list) - 1)
    else:
        p = 1.0 * np.arange(len(data_list))/(len(data_list) - 1)
    plt.plot(data_list, p, color=color)


IN_ROOT = "/Share/home/zhangqf8/lipan/DYNAMIC/shape_score/"

####### human ch

ch_D1 = read_RT(IN_ROOT+"hek_ch_vivo/rt_D1.txt")
ch_D2 = read_RT(IN_ROOT+"hek_ch_vivo/rt_D2.txt")

ch_N1 = read_RT(IN_ROOT+"hek_ch_vivo/rt_N1.txt")
ch_N2 = read_RT(IN_ROOT+"hek_ch_vivo/rt_N2.txt")

ch_T1 = read_RT(IN_ROOT+"hek_ch_vitro/rt_N1.txt")
ch_T2 = read_RT(IN_ROOT+"hek_ch_vitro/rt_N2.txt")


ch_D_list = calc_pearson(ch_D1, ch_D2, min_RT=2); print len(ch_D_list)
ch_N_list = calc_pearson(ch_N1, ch_N2, min_RT=2); print len(ch_N_list)
ch_T_list = calc_pearson(ch_T1, ch_T2, min_RT=2); print len(ch_T_list)

cdf(ch_D_list, color='red', topdown=True)
cdf(ch_N_list, color='green', topdown=True)
cdf(ch_T_list, color='blue', topdown=True)

plt.xlim(-1, 1)
plt.savefig("figs/human_ch.pdf")
plt.show()


####### human np

np_D1 = read_RT(IN_ROOT+"hek_np_vivo/rt_D1.txt")
np_D2 = read_RT(IN_ROOT+"hek_np_vivo/rt_D2.txt")

np_N1 = read_RT(IN_ROOT+"hek_np_vivo/rt_N1.txt")
np_N2 = read_RT(IN_ROOT+"hek_np_vivo/rt_N2.txt")

np_T1 = read_RT(IN_ROOT+"hek_np_vitro/rt_N1.txt")
np_T2 = read_RT(IN_ROOT+"hek_np_vitro/rt_N2.txt")

np_D_list = calc_pearson(np_D1, np_D2, min_RT=2); print len(np_D_list)
np_N_list = calc_pearson(np_N1, np_N2, min_RT=2); print len(np_N_list)
np_T_list = calc_pearson(np_T1, np_T2, min_RT=2); print len(np_T_list)

cdf(np_D_list, color='red', topdown=True)
cdf(np_N_list, color='green', topdown=True)
cdf(np_T_list, color='blue', topdown=True)

plt.xlim(-1, 1)
plt.savefig("figs/human_np.pdf")
plt.show()

####### human cy

cy_D1 = read_RT(IN_ROOT+"hek_cy_vivo/rt_D1.txt")
cy_D2 = read_RT(IN_ROOT+"hek_cy_vivo/rt_D2.txt")

cy_N1 = read_RT(IN_ROOT+"hek_cy_vivo/rt_N1.txt")
cy_N2 = read_RT(IN_ROOT+"hek_cy_vivo/rt_N2.txt")

cy_T1 = read_RT(IN_ROOT+"hek_cy_vitro/rt_N1.txt")
cy_T2 = read_RT(IN_ROOT+"hek_cy_vitro/rt_N2.txt")


cy_D_list = calc_pearson(cy_D1, cy_D2, min_RT=2); print len(cy_D_list)
cy_N_list = calc_pearson(cy_N1, cy_N2, min_RT=2); print len(cy_N_list)
cy_T_list = calc_pearson(cy_T1, cy_T2, min_RT=2); print len(cy_T_list)

cdf(cy_D_list, color='red', topdown=True)
cdf(cy_N_list, color='green', topdown=True)
cdf(cy_T_list, color='blue', topdown=True)

plt.xlim(-1, 1)
plt.savefig("figs/human_cy.pdf")
plt.show()







####### mouse ch

ch_D1 = read_RT(IN_ROOT+"mes_ch_vivo/rt_D1.txt")
ch_D2 = read_RT(IN_ROOT+"mes_ch_vivo/rt_D2.txt")

ch_N1 = read_RT(IN_ROOT+"mes_ch_vivo/rt_N1.txt")
ch_N2 = read_RT(IN_ROOT+"mes_ch_vivo/rt_N2.txt")

ch_T1 = read_RT(IN_ROOT+"mes_ch_vitro/rt_N1.txt")
ch_T2 = read_RT(IN_ROOT+"mes_ch_vitro/rt_N2.txt")


ch_D_list = calc_pearson(ch_D1, ch_D2, min_RT=2); print len(ch_D_list)
ch_N_list = calc_pearson(ch_N1, ch_N2, min_RT=2); print len(ch_N_list)
ch_T_list = calc_pearson(ch_T1, ch_T2, min_RT=2); print len(ch_T_list)

cdf(ch_D_list, color='red', topdown=True)
cdf(ch_N_list, color='green', topdown=True)
cdf(ch_T_list, color='blue', topdown=True)

plt.xlim(-1, 1)
plt.savefig("figs/mouse_ch.pdf")
plt.show()


####### mouse np

np_D1 = read_RT(IN_ROOT+"mes_np_vivo/rt_D1.txt")
np_D2 = read_RT(IN_ROOT+"mes_np_vivo/rt_D2.txt")

np_N1 = read_RT(IN_ROOT+"mes_np_vivo/rt_N1.txt")
np_N2 = read_RT(IN_ROOT+"mes_np_vivo/rt_N2.txt")

np_T1 = read_RT(IN_ROOT+"mes_np_vitro/rt_N1.txt")
np_T2 = read_RT(IN_ROOT+"mes_np_vitro/rt_N2.txt")

np_D_list = calc_pearson(np_D1, np_D2, min_RT=2); print len(np_D_list)
np_N_list = calc_pearson(np_N1, np_N2, min_RT=2); print len(np_N_list)
np_T_list = calc_pearson(np_T1, np_T2, min_RT=2); print len(np_T_list)

cdf(np_D_list, color='red', topdown=True)
cdf(np_N_list, color='green', topdown=True)
cdf(np_T_list, color='blue', topdown=True)

plt.xlim(-1, 1)
plt.savefig("figs/mouse_np.pdf")
plt.show()

####### mouse cy

cy_D1 = read_RT(IN_ROOT+"mes_cy_vivo/rt_D1.txt")
cy_D2 = read_RT(IN_ROOT+"mes_cy_vivo/rt_D2.txt")

cy_N1 = read_RT(IN_ROOT+"mes_cy_vivo/rt_N1.txt")
cy_N2 = read_RT(IN_ROOT+"mes_cy_vivo/rt_N2.txt")

cy_T1 = read_RT(IN_ROOT+"mes_cy_vitro/rt_N1.txt")
cy_T2 = read_RT(IN_ROOT+"mes_cy_vitro/rt_N2.txt")








