
###### Draw pie and modicate with Illustrater


from seq import *
from icSHAPE import *

def get_key_value_pairs(Dict):
    k_v_pairs = []
    for key in Dict:
        k_v_pairs.append( (key, Dict[key]) )
    k_v_pairs.sort(key=lambda x: x[1], reverse=True)
    keys = [item[0] for item in k_v_pairs]
    values = [item[1] for item in k_v_pairs]
    return keys, values

def pie_plot_reads_distribution(file_name, title):
    Exons = {}
    Introns = {}
    IN = open(file_name)
    line = IN.readline()
    RNA_type = ('pseudogene', 'lncRNA', 'mRNA', 'miRNA', 'misc_RNA', 'snoRNA', 'snRNA')
    while line:
        data = line.strip().split()
        if data[0] == 'EXON':
            if data[1] in RNA_type:
                Exons[ data[1] ] = float(data[2])
        if data[0] == 'INTRON':
            if data[1] in RNA_type:
                Introns[ data[1] ] = float(data[2])
        line = IN.readline()
    IN.close()
    for trans_type in RNA_type:
        if trans_type not in Exons: Exons[trans_type] = 0
        if trans_type not in Introns: Introns[trans_type] = 0
    exon_total = sum(Exons.values())
    intron_total = sum(Introns.values())
    exon_intron_ratio = 1.0 * exon_total / intron_total
    print "exon_intron_ratio: ", exon_intron_ratio
    inner_radius = numpy.sqrt((1+0.04*exon_intron_ratio)/(1+exon_intron_ratio))
    # pie chart
    RNA_Type_Color = [ ['mRNA', '#35a1cd'], ['lncRNA', '#b47a01'], ['misc_RNA', '#b44600'], ['pseudogene', '#ebbb30'], ['snRNA', '#7cba1c'], ['miRNA', '#0f6087'], ['snoRNA', '#ff0000']]
    #RNA_types = ['mRNA', 'lncRNA', 'miscRNA', 'snRNA', 'pseudogene', '']
    colors = [item[1] for item in RNA_Type_Color]
    explode = (0, 0, 0, 0, 0, 0, 0)  # explode 1st slice
    plt.figure(figsize=(11,5))
    plt.subplot(1,2,1)
    #exon_trans_types, exon_trans_ratio = get_key_value_pairs(Exons)
    trans_types = [item[0] for item in RNA_Type_Color]
    exon_trans_ratio = [Exons[item] for item in trans_types]
    intron_trans_ratio = [Introns[item] for item in trans_types]
    plt.pie(exon_trans_ratio, radius=1 ,explode=explode, labels=trans_types, colors=colors, autopct='%1.1f%%', shadow=False, startangle=90)
    plt.title(title+"-Exon")
    plt.subplot(1,2,2)
    #intron_trans_types = []
    #intron_trans_ratio = []
    #for trans_type in exon_trans_types:
    #    intron_trans_types.append(trans_type)
    #    intron_trans_ratio.append(Introns[trans_type])
    #intron_trans_types, intron_trans_ratio = get_key_value_pairs(Introns)
    plt.pie(intron_trans_ratio, radius=inner_radius ,explode=explode, labels=trans_types, colors=colors, autopct='%1.1f%%', shadow=False, startangle=90)
    plt.title(title+"-Intron")
    plt.savefig(HOME+"/figs/"+title+".pdf")
    plt.show()

pie_plot_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/hek_np_trans.count", "hek_np")
pie_plot_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/hek_cy_trans.count", "hek_cy")

pie_plot_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/mes_ch_trans.count", "mes_ch")
pie_plot_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/mes_np_trans.count", "mes_np")
pie_plot_reads_distribution("/Share/home/zhangqf8/lipan/DYNAMIC/statistic_reads_distribution/mes_cy_trans.count", "mes_cy")



