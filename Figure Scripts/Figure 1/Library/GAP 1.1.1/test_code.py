



from ParseTrans import *

Parser = ParseTransClass(genomeCoorBedFile = "/Users/lee/Desktop/test_Ensembl/hg19.genomeCoor.bed")


##### transCoor => genomeCoor +
    # start
Parser.transCoor2genomeCoor("ENST00000456328.2_1", 1)   # ('chr1', 11869, '+')

    # border
Parser.transCoor2genomeCoor("ENST00000456328.2_1", 359) # ('chr1', 12227, '+')
Parser.transCoor2genomeCoor("ENST00000456328.2_1", 360) # ('chr1', 12613, '+')

    # end
Parser.transCoor2genomeCoor("ENST00000456328.2_1", 1657) # ('chr1', 14409, '+')
Parser.transCoor2genomeCoor("ENST00000456328.2_1", 1658) # CoorFunc.out_of_range

##### transCoor => genomeCoor -
    # start
Parser.transCoor2genomeCoor("ENST00000488147.1_1", 1)  #('chr1', 29570, '-')

    # border
Parser.transCoor2genomeCoor("ENST00000488147.1_1", 37)  #('chr1', 29534, '-')
Parser.transCoor2genomeCoor("ENST00000488147.1_1", 38)  #('chr1', 24891, '-')

    # end
Parser.transCoor2genomeCoor("ENST00000488147.1_1", 1351)  #('chr1', 14404, '-')
Parser.transCoor2genomeCoor("ENST00000488147.1_1", 1352)  #CoorFunc.out_of_range




##### genomeCoor => transCoor +
    # start
Parser.genomeCoor2transCoor("chr1", 11869, 11869, '+')      # [['chr1', 11869, 11869, 'ENST00000456328.2_1', 1, 1]]
Parser.genomeCoor2transCoor("chr1", 11869, 11869+10, '+')   # [['chr1', 11869, 11879, 'ENST00000456328.2_1', 1, 11]]

    # Border
Parser.genomeCoor2transCoor("chr1", 12226, 12615, '+')      # [['chr1', 12226, 12227, 'ENST00000456328.2_1', 358, 359], ['chr1', 12613, 12615, 'ENST00000456328.2_1', 360, 362]]

    # end
Parser.genomeCoor2transCoor("chr1", 14409-2, 14409+2, '+')  # [['chr1', 14407, 14409, 'ENST00000456328.2_1', 1655, 1657]]
Parser.genomeCoor2transCoor("chr1", 14409, 14409, '+')      # [['chr1', 14409, 14409, 'ENST00000456328.2_1', 1657, 1657]]


##### genomeCoor => transCoor -
    # start
Parser.genomeCoor2transCoor("chr1", 29570, 29570, '-')      # [['chr1', 29570, 29570, 'ENST00000488147.1_1', 1, 1]]
Parser.genomeCoor2transCoor("chr1", 29570-3, 29570+2, '-')  # [['chr1', 29567, 29570, 'ENST00000488147.1_1', 1, 4]]

    # Border
Parser.genomeCoor2transCoor("chr1", 24891, 29534, '-')      # [['chr1', 29534, 29534, 'ENST00000488147.1_1', 37, 37], ['chr1', 24891, 24891, 'ENST00000488147.1_1', 38, 38]]

    # end
Parser.genomeCoor2transCoor("chr1", 14404, 14404, '-')      # [['chr1', 14404, 14404, 'ENST00000488147.1_1', 1351, 1351]]
Parser.genomeCoor2transCoor("chr1", 14404-3, 14404+2, '-')  # [['chr1', 14404, 14406, 'ENST00000488147.1_1', 1349, 1351]]



##### geneCoor => genomeCoor +
    # start
Parser.geneCoor2genomeCoor("ENSG00000223972.5_2", 1)        # ('chr1', 11869, '+')
Parser.geneCoor2genomeCoor("ENSG00000223972.5_2", 2)        # ('chr1', 11870, '+')
    
    # end
Parser.geneCoor2genomeCoor("ENSG00000223972.5_2", 2541)     # ('chr1', 14409, '+')
Parser.geneCoor2genomeCoor("ENSG00000223972.5_2", 2542)     # CoorFunc.out_of_range

##### geneCoor => genomeCoor -
    # start
Parser.geneCoor2genomeCoor("ENSG00000227232.5_2", 1)        # ('chr1', 29570, '-')
Parser.geneCoor2genomeCoor("ENSG00000227232.5_2", 2)        # ('chr1', 29569, '-')

    # end
Parser.geneCoor2genomeCoor("ENSG00000227232.5_2", 15167)    # ('chr1', 14404, '-')
Parser.geneCoor2genomeCoor("ENSG00000227232.5_2", 15168)    # CoorFunc.out_of_range


##### genomeCoor => geneCoor +
    # start
Parser.genomeCoor2geneCoor("chr1", 11869, 11869+10, '+')    # [['chr1', 11869, 11879, 'ENSG00000223972.5_2', 1, 11]]
Parser.genomeCoor2geneCoor("chr1", 11869, 11869, '+')       # [['chr1', 11869, 11869, 'ENSG00000223972.5_2', 1, 1]]

    # end
Parser.genomeCoor2geneCoor("chr1", 14409-3, 14409+2, '+')   # [['chr1', 14406, 14409, 'ENSG00000223972.5_2', 2538, 2541]]
Parser.genomeCoor2geneCoor("chr1", 14409, 14409, '+')       # [['chr1', 14409, 14409, 'ENSG00000223972.5_2', 2541, 2541]]

##### genomeCoor => geneCoor -
    # start
Parser.genomeCoor2geneCoor("chr1", 29570-3, 29570+2, '-')   # [['chr1', 29567, 29570, 'ENSG00000227232.5_2', 1, 4]]
Parser.genomeCoor2geneCoor("chr1", 29570, 29570, '-')       # [['chr1', 29570, 29570, 'ENSG00000227232.5_2', 1, 1]]

    # end
Parser.genomeCoor2geneCoor("chr1", 14404-3, 14404+2, '-')   # [['chr1', 14404, 14406, 'ENSG00000227232.5_2', 15165, 15167]]
Parser.genomeCoor2geneCoor("chr1", 14404, 14404, '-')       # [['chr1', 14404, 14404, 'ENSG00000227232.5_2', 15167, 15167]]






Parser.geneCoor2genomeCoor("ENSG00000227232.5_2", 10)
Parser.genomeCoor2geneCoor("chr1", 11869, 11869+10, '+')
Parser.genomeCoor2transCoor("chr1", 11869, 11869+10, '+')
Parser.transCoor2genomeCoor("ENST00000456328.2_1", 1)
Parser.geneCoor2transCoor("ENSG00000227232.5_2", 1, 11)
Parser.transCoor2geneCoor("ENST00000456328.2_1", 1, 11)




