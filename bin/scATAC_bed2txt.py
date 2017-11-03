#!/usr/bin/env python

"""
1) convert bed file to txt file 
2) extract reads with selected barcode
Created by Rongxin Fang
"""

import sys
import gzip 

def main():
    from argparse import ArgumentParser
    # parameters 
    parser = ArgumentParser(description='extract reads from sleected cell')
    parser.add_argument('-i', '--input', help='bed file contains read', required=True)
    parser.add_argument('-o', '--outut', help='output txt file', required=True)
    parser.add_argument('-x', '--barocde', help='selected barcode', required=True)
    options = parser.parse_args()
    
    # input parsing
    input_bed = options.input
    output_txt = options.outut
    barocde_txt = options.barocde
    
    # read the barcodes
    sel_barcodes = set([line.split()[0] for line in open(barocde_txt).readlines()])
	
    # open output file
    if input_bed.endswith(".gz"):
        fin = gzip.open(input_bed, "rb")
    else:
        fin = open(input_bed, "r")

    if output_txt.endswith(".gz"):
        fout = gzip.open(output_txt, "wb")
    else:
        fout = open(output_txt, "w")
    
    for line in fin:
        elems = line.split()
        barcode = elems[3].split(":")[0]
        if(barcode in sel_barcodes):
            fout.write("\t".join([elems[0], elems[1], elems[2], barcode])+"\n")
    fin.close()
    fout.close()
    
if __name__ == '__main__':
    main()

