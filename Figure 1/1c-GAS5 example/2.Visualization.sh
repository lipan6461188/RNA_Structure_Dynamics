

##### Visualization

## show SHAPE map with VARNA

import visual, structure, tools

seq_gas5 = tools.readSeq("/Share/home/zhangqf8/lipan/DYNAMIC/human_GAS5/index/Gas5.fa")['ENST00000430245']
shape_gas5_ch = tools.loadicSHAPE("/Share/home/zhangqf8/lipan/DYNAMIC/human_GAS5/shape_score/hek_ch/shape.out")['ENST00000430245']
shape_gas5_np = tools.loadicSHAPE("/Share/home/zhangqf8/lipan/DYNAMIC/human_GAS5/shape_score/hek_np/shape.out")['ENST00000430245']
shape_gas5_cy = tools.loadicSHAPE("/Share/home/zhangqf8/lipan/DYNAMIC/human_GAS5/shape_score/hek_cy/shape.out")['ENST00000430245']

s = seq_gas5.find("CCAGTGGTC")
e = s + 25

shape_ch = shape_gas5_ch[s:e]
shape_np = shape_gas5_np[s:e]
shape_cy = shape_gas5_cy[s:e]

ss = structure.predictStructure(seq_gas5[s:e])

visual.Plot_RNAStructure_Shape(seq_gas5[s:e], ss, shape_ch, mode='heatmap')
visual.Plot_RNAStructure_Shape(seq_gas5[s:e], ss, shape_np, mode='heatmap')
visual.Plot_RNAStructure_Shape(seq_gas5[s:e], ss, shape_cy, mode='heatmap')


