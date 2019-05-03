#!/usr/bin/env python
from __future__ import print_function
import re, sys
import argparse

def gtf_to_bed(input, attributes, cols, output, feature = None):
    
    for index, line in enumerate(input):
            
        if line.startswith("#"):
            continue
        
        fields = line.split("\t")
        if feature is not None:
            if fields[2] != feature:
                continue

        attribute_list = []
        
        for attr in attributes:
            pat = '\\b' + attr + ' \"([^;.]+)\";'
            m = re.search(pat, fields[8])
            if m is None:
                attribute_list.append("NA")
                continue
            
            else:
                attribute_list.append(m.group(1))
        
        if cols == "NA":
            pass
        else:
            for i in cols:   
                try:
                    attribute_list.append(fields[i - 1].rstrip('\n'))
                
                except IndexError:
                    sys.exit("Error column requested does not exist at line %d" % (index + 1))
                
        chrom = fields[0]
        start = int(fields[3]) - 1
        stop = int(fields[4])
        strand = fields[6]
        
        out_attrs = '\t'.join(str(e) for e in attribute_list)
        
        print('%s\t%d\t%d\t.\t.\t%s\t%s' % 
                (chrom, start, stop, strand, out_attrs), file = output) 
                

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description="""convert gtf to bed, while keeping attributes of interest""")
    parser.add_argument('-i','--input_gtf', help="""input gtf file""", required = True)
    parser.add_argument('-a','--attributes' , help="""list of attributes to keep (appended to the 7th to nth columns)""",required = True, nargs='+')
    parser.add_argument('-c','--cols' , help="""additional columns to keep (tab) deliminated""",required = False, nargs='+', type=int)
    parser.add_argument('-o','--output_bed' , help="""output bed file,
            defaults to stdout""",required = False)
    parser.add_argument('-f', '--feature', help="""feature to keep, defaults to all""")
    args = parser.parse_args()
    
    input_gtf = open(args.input_gtf, 'r')
    attributes = args.attributes
    if args.output_bed:
        outfh = open(args.output_bed, 'w')
    else:
        outfh = sys.stdout

    cols=args.cols
    
    if cols:
        gtf_to_bed(input_gtf, attributes, cols, outfh, args.feature)
    else:
        gtf_to_bed(input_gtf, attributes,"NA", outfh, args.feature)
    

    if args.output_bed: 
        outfh.close() 
