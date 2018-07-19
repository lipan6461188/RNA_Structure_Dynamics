#-*- coding:utf-8 -*-

"""

1. Plot RNA Structure with icSHAPE scores
Plot_RNAStructure_Shape(seq, ss, shape_list, mode='label', correctT=True, highlight_region=[], VARNAProg="/Users/lee/Documents/VARNAv3-93.jar")

2. Plot RNA Structure with Bases
Plot_RNAStructure_Base(seq, ss, mode='fill', correctT=True, highlight_region=[], VARNAProg="/Users/lee/Documents/VARNAv3-93.jar")

3. Plot RNA Structure with given 1-based coordiation
Plot_RNAStructure_highlight(seq, ss, hg_base_list=[], mode='fill', correctT=True, highlight_region=[], VARNAProg="/Users/lee/Documents/VARNAv3-93.jar")

"""


import os, sys, commands, random, platform


def highlight_region_cmd(region_list):
    """
    The input must be [ (1, 10, '#FF0000'), (90, 100, '#00FF00') ]...
    1-based coordination system
    """
    assert( (isinstance(region_list, list) or isinstance(region_list, tuple)) and len(region_list) >= 1 )
    assert( isinstance(region_list[0], list) or isinstance(region_list[0], tuple) )
    CMD = '-highlightRegion "'
    for region in region_list:
        if len(region) == 2:
            CMD += "%s-%s:fill=%s,outline=#FFFFFF,radius=15;" % (region[0], region[1], "00FF00")
        else:
            CMD += "%s-%s:fill=%s,outline=#FFFFFF,radius=15;" % (region[0], region[1], region[2])
    CMD += '"'
    return CMD

def base_color_shape_cmd(shape_list):
    Level1 = ""
    Level2 = ""
    Level3 = ""
    Level4 = ""
    NoData = ""
    for idx in range(len(shape_list)):
        if shape_list[idx] == 'NULL':
            NoData += `idx+1` if NoData == "" else ','+`idx+1`
        elif float(shape_list[idx]) > 0.7:
            Level1 += `idx+1` if Level1 == "" else ','+`idx+1`
        elif float(shape_list[idx]) > 0.5:
            Level2 += `idx+1` if Level2 == "" else ','+`idx+1`
        elif float(shape_list[idx]) > 0.3:
            Level3 += `idx+1` if Level3 == "" else ','+`idx+1`
        else:
            Level4 += `idx+1` if Level4 == "" else ','+`idx+1`
    CMD = ""
    if Level1: CMD += "-applyBasesStyle1on \"%s\" " % (Level1, )
    if Level2: CMD += "-applyBasesStyle2on \"%s\" " % (Level2, )
    if Level3: CMD += "-applyBasesStyle3on \"%s\" " % (Level3, )
    if Level4: CMD += "-applyBasesStyle4on \"%s\" " % (Level4, )
    if NoData: CMD += "-applyBasesStyle5on \"%s\" " % (NoData, )
    return CMD

def base_color_heatmap_cmd(shape_list):
    shape_str = ""
    null_base_idx = []
    for idx,shape in enumerate(shape_list):
        if shape == 'NULL':
            shape_str += '0.0;'
            null_base_idx.append(str(idx+1))
        else:
            shape_str += '%s;' % (shape, )
    CMD = ""
    if null_base_idx:
        CMD += "-basesStyle1 \"label=#828282\" -applyBasesStyle1on \"%s\" " % ( ",".join(null_base_idx), )
    CMD += "-colorMap \"%s\" " % (shape_str[:-1], )
    return CMD

def Plot_RNAStructure_Shape(seq, ss, shape_list, mode='label', correctT=True, highlight_region=[], title="", wait=True, VARNAProg="/Users/lee/Documents/VARNAv3-93.jar"):
    """
    根据icSHAPE值给每一个碱基赋颜色值
    mode must be one of 'label'/'fill'/'heatmap'
    Example: Plot_RNAStructure_Shape("AGCTGGGTTTCCCGATT", "....(((...)))....", [0.8]*4+[0.2]*3+[0.9]*3+[0.2]*3+['NULL']*4 )
    """
    import platform, os
    assert(len(seq) == len(ss) == len(shape_list))
    assert(mode in ('label', 'fill', 'heatmap'))
    if correctT:
        seq = seq.replace('T', 'U')
    CMD = "java -cp "+VARNAProg+" fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN \"%s\" -drawBackbone false " % (seq, ss)
    if mode == 'label':
        CMD += "-basesStyle1 \"label=#B61D22\" " + "-basesStyle2 \"label=#ED9616\" " + "-basesStyle3 \"label=#194399\" " + "-basesStyle4 \"label=#040000\" " + "-basesStyle5 \"label=#828282\" "
        CMD += base_color_shape_cmd(shape_list)
    elif mode == 'fill':
        CMD += "-basesStyle1 \"fill=#B61D22\" " + "-basesStyle2 \"fill=#ED9616\" " + "-basesStyle3 \"fill=#194399\" " + "-basesStyle4 \"fill=#040000\" " + "-basesStyle5 \"fill=#828282\" "
        CMD += base_color_shape_cmd(shape_list)
    else:
        CMD += base_color_heatmap_cmd(shape_list)
    if highlight_region:
        CMD += " " + highlight_region_cmd(highlight_region)
    if title:
        CMD += " -title \"%s\"" % (title, )
    if not wait:
        CMD = "nohup " + CMD + " &"
    print CMD
    if platform.system() == "Darwin":
        os.system(CMD)
    else:
        print "Copy and Paste"

def base_color_base_cmd(seq):
    A = ""
    T = ""
    C = ""
    G = ""
    for idx in range(len(seq)):
        if seq[idx] == 'A':
            A += `idx+1` if A == "" else ','+`idx+1`
        if seq[idx] == 'C':
            C += `idx+1` if C == "" else ','+`idx+1`
        if seq[idx] == 'T' or seq[idx] == 'U':
            T += `idx+1` if T == "" else ','+`idx+1`
        if seq[idx] == 'G':
            G += `idx+1` if G == "" else ','+`idx+1`
    CMD = ""
    if A: CMD += "-applyBasesStyle1on \"%s\" " % (A, )
    if T: CMD += "-applyBasesStyle2on \"%s\" " % (T, )
    if C: CMD += "-applyBasesStyle3on \"%s\" " % (C, )
    if G: CMD += "-applyBasesStyle4on \"%s\" " % (G, )
    return CMD

def Plot_RNAStructure_Base(seq, ss, mode='fill', correctT=True, highlight_region=[], title="", wait=True, VARNAProg="/Users/lee/Documents/VARNAv3-93.jar"):
    """
    根据碱基给每一个碱基赋颜色值
    mode must be one of 'label'/'fill'
    Example: Plot_RNAStructure_Base("AGCTGGGTTTCCCGATT", "....(((...)))....")
    """
    import platform, os
    assert(len(seq) == len(ss))
    assert(mode in ('label', 'fill'))
    if correctT:
        seq = seq.replace('T', 'U')
    CMD = "java -cp "+VARNAProg+" fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN \"%s\" -drawBackbone false " % (seq, ss)
    if mode == 'label':
        CMD += "-basesStyle1 \"label=#CCFFCC\" " + "-basesStyle2 \"label=#CCFFFF\" " + "-basesStyle3 \"label=#FFFFCC\" " + "-basesStyle4 \"label=#FFCCFF\" "
    else:
        CMD += "-basesStyle1 \"fill=#CCFFCC\" " + "-basesStyle2 \"fill=#CCFFFF\" " + "-basesStyle3 \"fill=#FFFFCC\" " + "-basesStyle4 \"fill=#FFCCFF\" "
    CMD += base_color_base_cmd(seq)
    if highlight_region:
        CMD += " " + highlight_region_cmd(highlight_region)
    if title:
        CMD += " -title \"%s\"" % (title, )
    if not wait:
        CMD = "nohup " + CMD + " &"
    print CMD
    if platform.system() == "Darwin":
        os.system(CMD)
    else:
        print "Copy and Paste"


def Plot_RNAStructure_highlight(seq, ss, hg_base_list=[], mode='fill', correctT=True, highlight_region=[], title="", wait=True, VARNAProg="/Users/lee/Documents/VARNAv3-93.jar"):
    """
    高亮其中的某些碱基
    mode must be one of 'label'/'fill'
    Example: Plot_RNAStructure_highlight("AGCTGGGTTTCCCGATT", "....(((...)))....", hg_base_list=[1,2,3,4])
    """
    import platform, os
    assert(len(seq) == len(ss))
    assert(mode in ('label', 'fill'))
    if correctT:
        seq = seq.replace('T', 'U')
    CMD = "java -cp "+VARNAProg+" fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN \"%s\" -drawBackbone false " % (seq, ss)
    if hg_base_list:
        if mode == 'label':
            CMD += "-basesStyle1 \"label=#FF0000\" "
        else:
            CMD += "-basesStyle1 \"fill=#FF0000\" "
        hg_base_list = [`item` for item in hg_base_list]
        CMD += "-applyBasesStyle1on \"%s\" " % (",".join(hg_base_list), )
    if highlight_region:
        CMD += " " + highlight_region_cmd(highlight_region)
    if title:
        CMD += " -title \"%s\"" % (title, )
    if not wait:
        CMD = "nohup " + CMD + " &"
    print CMD
    if platform.system() == "Darwin":
        os.system(CMD)
    else:
        print "Copy and Paste"


def cdf(data, bins=100, color=""):
    import numpy as np
    import matplotlib.pyplot as plt
    counts, bin_edges = np.histogram (data, bins=bins, normed=True)
    cdf = np.cumsum(counts)
    if color:
        plt.plot(bin_edges[1:], cdf/cdf[-1], color=color)
    else:
        plt.plot(bin_edges[1:], cdf/cdf[-1])


