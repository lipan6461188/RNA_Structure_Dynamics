#!/usr/bin/env python
#-*- coding:utf-8 -*- 

import re, sys, os, getopt, time, datetime
import NCBI_Genome, GENCODE_Genome

file_information = """
GTF: 
    1. Emsembl: https://asia.ensembl.org/info/data/ftp/index.html
    2. Gencode: https://www.gencodegenes.org/

GFF3:
    1. ftp://ftp.ncbi.nlm.nih.gov: ftp://ftp.ncbi.nlm.nih.gov/genomes/Equus_caballus/GFF/

"""


Usage = """
## --------------------------------------
parse genome coordinate information and transcript coordinate information from Gencode GTF and NCBI GFF3 files
 
Command:
%s -g GFF3/GTF -o output_prefix -s [gencode|ncbi]
 
# what it is:
 -g         genome annotation
 -o         output file prefix
 -s         <gencode/ncbi> data source

# more options [ Strongly Recommend Default Parameters ]
 --removeTransVersion <yes/no>
            ENST00000488147.1 => ENST00000488147 (default: no)
 --removeGeneVersion <yes/no>
            ENSG00000227232.5 => ENSG00000227232 (default: no)
 --removeScaffold <yes/no>
            remove KI270734.1... (default: no)
 --genome   fetch transcriptome from genome file
 --genome_gff_same_chr <yes/no>
            GFF3 file:
                NC_000001.11    RefSeq  region  1       248956422       .       +       .       ID=id0;Dbxref=taxon:9606;Name=1;chromosome=1;gbkey=Src;genome=chromosome;mol_type=genomic DNA
            if genome file:
                    >NC_000001.11
            then yes
            if genome file:
                    >chr1
            then no
            (default: no)

""" % (sys.argv[0], )

def init():
    params = { 'removeTransVersion': False, 'removeGeneVersion': False, 'removeScaffold': False, 'genome_gff_same_chr': False }
    opts, args = getopt.getopt(sys.argv[1:], 'hg:o:s:', ['removeTransVersion=', 'removeGeneVersion=', 'removeScaffold=', 'genome=', 'genome_gff_same_chr='])
    for op, value in opts:
        if op == '-h':
            print Usage;
            sys.exit(-1)
        elif op == '--removeTransVersion':
            params['removeTransVersion'] = True if value == 'yes' else False
            assert value in ('yes', 'no'), 'Error: --removeTransVersion <yes/no>'
        elif op == '--removeGeneVersion':
            params['removeGeneVersion'] = True if value == 'yes' else False
            assert value in ('yes', 'no'), 'Error: --removeGeneVersion <yes/no>'
        elif op == '--removeScaffold':
            params['removeScaffold'] = True if value == 'yes' else False
            assert value in ('yes', 'no'), 'Error: --removeScaffold <yes/no>'
        elif op == '--genome_gff_same_chr':
            params['genome_gff_same_chr'] = True if value == 'yes' else False
            assert value in ('yes', 'no'), 'Error: --removeScaffold <yes/no>'
        elif op == '--genome':
            params['genomeFile'] = value
        elif op == '-g':
            params['input'] = os.path.abspath(value)
        elif op == '-o':
            params['output'] = value
        elif op == '-s':
            params['source'] = value
            assert value in ('ncbi', 'gencode'), 'Error: -t <gencode/ncbi>'
    if 'input' not in params:
        print 'Error: specify -g'
        print Usage
        sys.exit(-1)
    if 'output' not in params:
        print 'Error: specify -o'
        print Usage
        sys.exit(-1)
    if 'source' not in params:
        print 'Error: specify -s'
        print Usage
        sys.exit(-1)
    return params

def main():
    params = init()
    genomeCoorFileName = params['output']+'.genomeCoor.bed'
    transCoorFileName = params['output']+'.transCoor.bed'
    if params['source'] == 'gencode':
        gtf = GENCODE_Genome.GENCODE_Genome_Class(params['input'])
        gtf.write_genomeCoor_bed(genomeCoorFileName, onlyChr=params['removeScaffold'], pureTransID=params['removeTransVersion'], pureGeneID=params['removeGeneVersion'])
        NCBI_Genome.NCBI_Genome_Class.genomeCoorBed_To_transCoorBed(genomeCoorFileName, transCoorFileName)
        if 'genomeFile' in params:
            transcriptomeFileName = params['output']+'_transcriptome.fa'
            gtf.writeTranscriptome(params['genomeFile'], transcriptomeFileName, pureTransID=params['removeTransVersion'], onlyChr=params['removeScaffold'], verbose=True)
    elif params['source'] == 'ncbi':
        gff3 = NCBI_Genome.NCBI_Genome_Class(params['input'])
        gff3.write_genomeCoor_bed(genomeCoorFileName, onlyChr=params['removeScaffold'], pureTransID=params['removeTransVersion'])
        NCBI_Genome.NCBI_Genome_Class.genomeCoorBed_To_transCoorBed(genomeCoorFileName, transCoorFileName)
        if 'genomeFile' in params:
            transcriptomeFileName = params['output']+'_transcriptome.fa'
            genomeChrSym = "chr"
            if params['genome_gff_same_chr']:
                genomeChrSym = ""
            gff3.writeTranscriptome(params['genomeFile'], transcriptomeFileName, genomeChrSym, pureTransID=params['removeTransVersion'], verbose=True)


main()

