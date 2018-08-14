
def get_max_min(input_list):
    input_list = sorted(input_list)
    list_len = len(input_list)
    quantile_05 = input_list[ int(0.05*list_len) ]
    quantile_95 = input_list[ int(0.95*list_len) ]
    median = input_list[ int(0.50*list_len) ]
    return quantile_05, median, quantile_95

def get_quantile(input_list, point=0.05):
    input_list = sorted(input_list)
    list_len = len(input_list)
    quantile_05 = input_list[ int(point*list_len) ]
    quantile_95 = input_list[ int((1-point)*list_len) ]
    median = input_list[ int(0.50*list_len) ]
    return quantile_05, median, quantile_95


def smooth(raw_list, winSize, winStep=1):
    s_smooth_v = np.mean(raw_list[:winSize])
    e_smooth_v = np.mean(raw_list[-winSize:])
    start = winSize/2
    end = len(raw_list)-winSize/2
    
    smoothed_v = []
    i = start
    while i<end:
        smoothed_v.append( np.mean(raw_list[i-winSize/2:i+winSize/2]) )
        i += winStep
    
    if winStep == 1:
        smoothed_v = [s_smooth_v]*(winSize/2) + smoothed_v + [e_smooth_v]*(winSize/2)
    
    return smoothed_v

def box_collection(raw_list, winSize, winStep=1):
    s_collect_v = np.mean(raw_list[:winSize])
    e_collect_v = np.mean(raw_list[-winSize:])
    start = winSize/2
    end = len(raw_list)-winSize/2
    
    collect_v = []
    i = start
    while i<end:
        collect_v.append( raw_list[i-winSize/2:i+winSize/2] )
        #for value in raw_list[i-winSize/2:i+winSize/2]:
        #    collect_v.append( (value, i) )
        i += winStep
    
    #collect_df = pd.DataFrame(collect_v, columns=['value', 'type'])
    
    return collect_v

def segment_heatmap(df, col_list, smooth_list=[], yticklabels=False, quantile={}, winSize=10, winStep=1, cbar=True):
    col_num = len(col_list)
    for col_i,col in enumerate(col_list):
        if col_i == 1: yticklabels = False
        plt.subplot(1,col_num+1,col_i+1)
        if col in quantile:
            if len(quantile[col]) == 3:
                quantile_05, median, quantile_95 = quantile[col]
            elif len(quantile[col]) == 2:
                quantile_05,quantile_95 = quantile[col]
                median = np.median(list(df[col]))
            else:
                quantile_05, median, quantile_95 = get_quantile(df[col], point=0.15)
        else:
            quantile_05, median, quantile_95 = get_quantile(df[col], point=0.15) #get_max_min(df[col])
        data_line = pd.DataFrame(df[col])
        if col_i+1 in smooth_list:
            data_line = pd.DataFrame(smooth(list(df[col]), winSize=winSize, winStep=winStep), columns=[col])
        sns.heatmap(data=data_line, cbar=cbar, yticklabels=yticklabels, vmin=quantile_05, center=median, vmax=quantile_95, cmap=sns.color_palette("RdBu_r", 100))#, mask=under_median)
        plt.yticks(rotation=0)
        plt.title("%.2f-%.2f-%.2f" % (quantile_05, median, quantile_95))
        plt.xticks(rotation=90)
    
    plt.tight_layout()


def segment_violin(df, col_list, winSize=10, winStep=1):
    col_num = len(col_list)
    for col_i,col in enumerate(col_list):
        plt.subplot(1,col_num,col_i+1)
        collect_df = box_collection(list(df[col]), winSize=winSize, winStep=winStep)
        
        collect_df.reverse()
        ax = plt.violinplot(dataset=collect_df, showmeans=True, vert=False)
        
        #ax.set_ytickslabel([])
        #ax.set_xtickslabel([])
        #plt.yticks(rotation=0)
        #plt.title("%.2f-%.2f" % (quantile_05, quantile_95))
        plt.xticks(rotation=90)
        plt.xlabel(col)
    
    plt.tight_layout()


def segment_heatmap_deng(data2D, col_list, yticklabels=False):
    sub_data2D = data2D.loc[:,col_list]
    col_num = len(col_list)
    
    plt.subplot(1,col_num+1,1)
    sns.dendrogram( sub_data2D, metric='euclidean', method='average', label=False, axis=0, ax=None, rotate=True, linkage=None)
    
    for col_i,col in enumerate(col_list):
        if col_i == 1: yticklabels = False
        plt.subplot(1,col_num+1,col_i+2)
        quantile_05, median, quantile_95 = get_max_min(data2D[col])
        sns.heatmap(data=pd.DataFrame(data2D[col]), cbar=False, yticklabels=yticklabels, vmin=quantile_05, center=median, vmax=quantile_95, cmap=sns.color_palette("coolwarm", 100))
        plt.yticks(rotation=0)
        plt.title("%.2f-%.2f" % (quantile_05, quantile_95))
        plt.xticks(rotation=90)
    
    plt.tight_layout()




com_trans_list = set(mouse_TS) & set(mouse_TE) & set(mouse_5UTR_gini['ch']) & set(mouse_5UTR_gini['cy'])
print len(com_trans_list) 

mouse_df = tools.init_rect(len(com_trans_list), 4, rowNames=com_trans_list, colNames=('TS' ,'TE', 'ch_5UTR_gini', 'cy_5UTR_gini'))

for tid in com_trans_list:
    mouse_df.loc[tid, 'TS'] = mouse_TS[tid]
    mouse_df.loc[tid, 'TE'] = mouse_TE[tid]
    mouse_df.loc[tid, 'ch_5UTR_gini'] = mouse_5UTR_gini['ch'][tid]
    mouse_df.loc[tid, 'cy_5UTR_gini'] = mouse_5UTR_gini['cy'][tid]


####################
### Heatmap
####################


### Sorted DataSet


quantile = { 'TS':[3.5, 6.5], 'ch_5UTR_gini':[0.50, 0.65], 'cy_5UTR_gini': [0.50, 0.65], 'TE':[-0.40, 0.90] }

sorted_df = mouse_df.sort_values(by='ch_5UTR_gini')
col_list = ('TS', 'ch_5UTR_gini', 'cy_5UTR_gini', 'TE' )
segment_heatmap(sorted_df, col_list, cbar=True, smooth_list=[1,2,3,4], quantile=quantile, winSize=24, winStep=24)
plt.savefig("figs/heatmap.pdf")
plt.show()




