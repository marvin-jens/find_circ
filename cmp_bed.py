#!/usr/bin/env python
import sys,os
from collections import defaultdict


def read_to_hash(fname,ds=0,de=0,flank=0,cover=False):
    pos = {}
    for line in file(fname):
        if line.startswith("#"):
            continue
        line = line.strip()
        chrom,start,end,name,score,sense = line.split('\t')[:6]
        start,end = int(start)+ds,int(end)+de

        pos[(chrom,start,end,sense)] = line
        
        if flank:
            for x in xrange(flank):
                pos[(chrom,start-x,end,sense)] = line
                pos[(chrom,start+x,end,sense)] = line
                pos[(chrom,start,end-x,sense)] = line
                pos[(chrom,start,end+x,sense)] = line
        
        #if cover:
            #for x in xrange
    return pos

N = defaultdict(int)
bed1 = read_to_hash(sys.argv[1],flank=0)
N['unique_input1'] = len(bed1)

bed2 = read_to_hash(sys.argv[2])
N['unique_input2'] = len(bed2)

for circ,line in bed2.items():
    if circ in bed1:
        if len(sys.argv) > 3:
            print "%s\t%s" % (bed1[circ].split('\t')[3],line.split('\t')[3])
        else:
            print bed1[circ]

        N['overlap'] += 1        
        del bed1[circ]
    else:
        N['input2_not_in_input1'] += 1
        
for k,l in bed1.items():
    print "MISSING\t%s" % l
    N['input1_not_in_input2'] += 1

for k in sorted(N.keys()):
    sys.stderr.write("%s\t%d\n" % (k,N[k]))
        
found = N['overlap']
detected = N['unique_input2']
total = N['unique_input1']
fp = N['input2_not_in_input1']

if found == total == detected:
    sys.stderr.write("files contain identical splice sites!\n")
        
#print "#sensitivity %d/%d = %.2f %%" % (found,total,float(found)/total*100)
#print "#FDR %d/%d = %.2f %%" % (fp,detected,float(fp)/detected*100)