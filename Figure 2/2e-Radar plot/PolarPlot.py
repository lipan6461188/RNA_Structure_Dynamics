
def plot_polar(Examples, colors='red'):
    import numpy as np
    
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    
    if isinstance(colors, list):
        assert len(Examples)<=len(colors)
    else:
        colors = ['red'] * len(Examples)
    
    for ex, color in zip(Examples, colors):
        TE = ex['TE']
        TS = ex['TS']
        ch_GINI = ex['ch_5UTR_gini']
        cy_GINI = ex['cy_5UTR_gini']
        
        dataLenth = 4
        data = np.array([TS, ch_GINI, cy_GINI, TE])
        
        angles = np.linspace(0, 2*np.pi, dataLenth, endpoint=False)
        data = np.concatenate((data, [data[0]])) 
        angles = np.concatenate((angles, [angles[0]])) 
        
        ax.plot(angles, data, 'ro-', linewidth=2, color=color)
    
    ax.set_xticklabels([])
    
    ax.grid(False)
    return ax


def normalized_df(indf, quantile={}):
    import copy
    norm_df = copy.deepcopy(indf)
    
    nr,nl = indf.shape
    min_max = tools.init_rect(nl, 2, rowNames=list(indf.columns), colNames=['min', 'max'])
    for c_index in range(nl):
        col_name = indf.columns[c_index]
        
        if col_name in quantile:
            quantile_05, quantile_95 = quantile[col_name]
            median = np.median(indf.loc[:, col_name])
        else:
            quantile_05, median, quantile_95 = get_quantile(indf.loc[:, col_name], point=0.01)
        min_max.iloc[c_index, 0] = quantile_05
        min_max.iloc[c_index, 1] = quantile_95
        
        large_row = norm_df.iloc[:, c_index] > quantile_95
        small_row = norm_df.iloc[:, c_index] < quantile_05
        
        norm_df.iloc[:, c_index] = 1.0*(indf.iloc[:,c_index]-quantile_05)/(quantile_95-quantile_05)
        
        max_v = 1.0
        min_v = 0.0
        norm_df.loc[large_row, col_name] = max_v
        norm_df.loc[small_row, col_name] = min_v
    
    return norm_df, min_max

def show_Example(id_list, df, Parser):
    for idx in id_list:
        print Parser.getTransFeature( df.index[idx] )['gene_name']
        print df.iloc[idx, :]

def search_example(norm_df, Parser):
    pairs = []
    
    for i in range(len(norm_df)):
        gene_type_1 = Parser.getTransFeature(norm_df.index[i])['gene_type']
        if gene_type_1 not in ('mRNA', 'protein_coding'): continue
        
        ex = dict(norm_df.iloc[i])
        TE1 = ex['TE']
        TS1 = ex['TS']
        ch_GINI1 = ex['ch_5UTR_gini']
        cy_GINI1 = ex['cy_5UTR_gini']
        
        if ch_GINI1==0 or cy_GINI1==0: continue
        for j in range(i+1, len(norm_df)):
            gene_type_2 = Parser.getTransFeature(norm_df.index[j])['gene_type']
            if gene_type_2 not in ('mRNA', 'protein_coding'): continue
            
            ex = dict(norm_df.iloc[j])
            TE2 = ex['TE']
            TS2 = ex['TS']
            ch_GINI2 = ex['ch_5UTR_gini']
            cy_GINI2 = ex['cy_5UTR_gini']
            
            if ch_GINI2==0 or cy_GINI2==0: continue
            
            if TE1<TE2 and TS1<TS2 and ch_GINI1<ch_GINI2 and cy_GINI1<cy_GINI2:
                pairs.append((j, i))
            elif TE1>TE2 and TS1>TS2 and ch_GINI1>ch_GINI2 and cy_GINI1>cy_GINI2:
                pairs.append((i, j))
    
    return pairs

def pair_2_dict(pairs):
    pair_dict = {}
    for l,r in pairs:
        if l not in pair_dict:
            pair_dict[l] = []
        pair_dict[l].append(r)
    return pair_dict

def DFS(G,v,seen=None,path=None):
    if seen is None: seen = []
    if path is None: path = [v]
    
    seen.append(v)
    
    paths = []
    if v not in G: return paths
    for t in G[v]:
        if t not in seen:
            t_path = path + [t]
            paths.append(tuple(t_path))
            paths.extend(DFS(G, t, seen[:], t_path))
    return paths


def find_longest_path(pair_dict, length):
    longest = length-1; length = 0
    for i in range(longest-1, -1, -1):
        path = sorted( DFS(pair_dict, i), key=lambda x: len(x) )
        if path:
            cur_len = len(path[-1])
        if cur_len>length:
            length = cur_len
            longest = i
    return longest, length

def get_longest_path(pair_dict, norm_df, Parser):
    paths = {}
    for row_id in range(len(norm_df)):
        tid = norm_df.index[row_id]
        gene_name = Parser.getTransFeature(tid)['gene_name']
        gene_type = Parser.getTransFeature(tid)['gene_type']
        if gene_type not in ('protein_coding', 'mRNA'): continue
        
        try:
            longest_path = sorted( DFS(pair_dict, row_id), key=lambda x: len(x) )[-1]
        except IndexError:
            continue
        
        gene_name_list = []
        for i in longest_path:
            cur_tid = norm_df.index[i]
            cur_gene_name = Parser.getTransFeature(cur_tid)['gene_name']
            cur_gene_type = Parser.getTransFeature(cur_tid)['gene_type']
            
            if cur_gene_type not in ('protein_coding', 'mRNA'): gene_name_list = []; break
            
            if i == 0 and cur_gene_name != gene_name: print "Error"
            
            gene_name_list.append( cur_gene_name )
        
        if gene_name_list:
            paths[row_id] = gene_name_list
    
    return paths


def draw_trans(trans_list, norm_df, Parser):
    Examples = []
    for i in trans_list:
        Examples.append( dict(norm_df.iloc[i]) )
        print norm_df.index[i]
        print Parser.getTransFeature(norm_df.index[i])['gene_name']
    ax = plot_polar(Examples, colors=sns.color_palette("husl", len(Examples)))
    ax.set_ylim(0,1)
    #plt.show()

def draw_all_trans(trans_list, norm_df, Parser):
    Examples = []
    for i in trans_list:
        Examples.append( dict(norm_df.iloc[i]) )
        print norm_df.index[i]
        print Parser.getTransFeature(norm_df.index[i])['gene_name']
    ax = plot_polar(Examples, colors=['gray']*len(Examples))
    ax.set_ylim(0,1)



com_trans_list = set(mouse_TS) & set(mouse_TE) & set(mouse_5UTR_gini['ch']) & set(mouse_5UTR_gini['cy'])
mouse_df = tools.init_rect(len(com_trans_list), 5, rowNames=com_trans_list, colNames=('TS' ,'TE', 'ch_5UTR_gini', 'cy_5UTR_gini', 'geneName'))

for tid in com_trans_list:
    mouse_df.loc[tid, 'TS'] = mouse_TS[tid]
    mouse_df.loc[tid, 'TE'] = mouse_TE[tid]
    mouse_df.loc[tid, 'ch_5UTR_gini'] = mouse_5UTR_gini['ch'][tid]
    mouse_df.loc[tid, 'cy_5UTR_gini'] = mouse_5UTR_gini['cy'][tid]
    mouse_df.loc[tid, 'geneName'] = mm10_parseTrans.getTransFeature(tid)['gene_name']

sorted_df = mouse_df.sort_values(by='ch_5UTR_gini')
quantile = {'TS':[2.5, 8.0], 'ch_5UTR_gini':[0.50, 0.75], 'cy_5UTR_gini':[0.50, 0.75], 'TE': [-1.5, 2.0]}
norm_df,min_max = normalized_df(sorted_df, quantile=quantile)
norm_df['TS'] = 1-norm_df['TS']
norm_df['TE'] = 1-norm_df['TE']


pairs = search_example(norm_df, mm10_parseTrans); print len(pairs)
pair_dict = pair_2_dict(pairs)

### Find longest path
longest_id = find_longest_path(pair_dict, length=len(norm_df))
trans_list = sorted( DFS(pair_dict, longest_id[0]), key=lambda x: len(x) )[-1]

### Draw Examples
### [TS, ch_GINI, cy_GINI, TE]

trans_list = (423, 290, 253, 225, 218, 216, 164, 99)
draw_trans(trans_list, norm_df, mm10_parseTrans)
plt.savefig("figs/radar.pdf")
plt.show()



