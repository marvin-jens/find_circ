#!/usr/bin/env python
import sys,os
from collections import defaultdict
from logging import debug,warning,error,info

from optparse import *

usage = """
%prog 1.bed 2.bed [3.bed] [4.bed] [...] > merged.bed

Merge BED or BED-like files on the genomic coordinates. Deals properly
with find_circ.py output and adds a few extra columns. 
"""

parser = OptionParser(usage=usage)
parser.add_option("-f","--flank",dest="flank",type=int,default=0,help="add flanking nucleotides to define more fuzzy overlap (default=0)")
parser.add_option("-s","--stats",dest="stats",default="",help="write statistics to this file (instead of stderr)")
parser.add_option("-c","--old-input",dest="old_input",action="store_true",help="switch on compatibility mode for old input format")
options,args = parser.parse_args()

from numpy import *

def read_to_hash(fname,ds=0,de=0,flank=0,cover=False):
    #print "loading",fname
    pos = {}
    for line in file(fname):
        if line.startswith("#"):
            continue
        line = line.strip()
        chrom,start,end,name,score,sense = line.split('\t')[:6]
        start,end = int(start)+ds,int(end)+de

        #print (chrom,start,end,sense)
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

inputs = [read_to_hash(a,flank=0) for a in args]
names = [os.path.abspath(a) for a in args]

by_name = dict(zip(names,inputs))
merge = {}
support = defaultdict(set)

for name,data in zip(names,inputs):
    N['input_%s'] = len(data)
    merge.update(data)
    for pos in data:
        support[pos].add(name)

from collections import Counter
comb = Counter([tuple(sorted(v)) for v in support.values()])

if options.stats:
    sfile = file(options.stats,"w")
else:
    sfile = sys.stderr

for c in sorted(comb.keys()):
    sfile.write("%s\t%d\n" % ("_AND_".join(c),comb[c]))

def consensus_line(lines,comb):
    samples = []
    counts = defaultdict(int)
    
    def setup_samples(values):
        allsamples = []
        for v in values:
            toadd = v.split(",")
            samples.append(toadd)
            allsamples.extend(toadd)
        return ",".join(sorted(allsamples))
    
    def assign_counts(values):
        for cs,ss in zip(values,samples):
            for samp,count in zip(ss,cs.split(',')):
                counts[samp] += int(count)
            
        res = []
        for k in sorted(counts.keys()):
            res.append(counts[k])
        return ",".join([str(c) for c in res])

    def append_uniq(values):
        v = set()
        for row in values:
            v |= set(row.split(","))

        return ",".join([str(x) for x in sorted(v) if x])

    if options.old_input:
        col_map = {
            3 : lambda values : ",".join(sorted(values)), # combine identifiers
            4 : lambda values : array(values,dtype=int).sum(), # sum n_reads
            6 : lambda values : array(values,dtype=int).sum(), # sum n_uniq
            7 : lambda values : max([int(x) for x in values]), # max of best_uniq_A
            8 : lambda values : max([int(x) for x in values]), # max of best_uniq_B
            9 : lambda values : array(values,dtype=int).sum(), # sum ov_linear_A
            10 : lambda values : array(values,dtype=int).sum(), # sum ov_linear_B
            11 : setup_samples,
            12 : assign_counts,
            13 : lambda values : min([int(x) for x in values]), # min of edits
            14 : lambda values : min([int(x) for x in values]), # min of anchor_overlap
            15 : lambda values : min([int(x) for x in values]), # min of breakpoints
        }
    else:
        col_map = {
            3 : lambda values : ",".join(sorted(values)), # combine identifiers
            4 : lambda values : array(values,dtype=int).sum(), # sum n_reads
            6 : lambda values : array(values,dtype=int).sum(), # sum n_uniq
            7 : lambda values : max([int(x) for x in values]), # max of best_uniq_A
            8 : lambda values : max([int(x) for x in values]), # max of best_uniq_B
            10 : setup_samples,
            11 : assign_counts,
            12 : lambda values : min([int(x) for x in values]), # min of edits
            13 : lambda values : min([int(x) for x in values]), # min of anchor_overlap
            14 : lambda values : min([int(x) for x in values]), # min of breakpoints
        }
    
    from itertools import izip_longest
    parts = []
    for i,column in enumerate(izip_longest(*[l.rstrip().split('\t') for l in lines],fillvalue="")):
        #print i,column
        if i in col_map:
            parts.append(str(col_map[i](column)))
        else:
            parts.append(append_uniq(column))

    return "\t".join(parts)

for pos in merge.keys():
    com = support[pos]
    
    lines = [by_name[name][pos] for name in sorted(com)]

    #for l in lines:
        #print l
    #if len(lines) > 1:
        #print sorted(com)
    print consensus_line(lines,comb)
#for names,inputs in zip(combinations(names


#for circ,line in marv.items():
    #if circ in anna:
        #if len(sys.argv) > 3:
            #print "%s\t%s" % (anna[circ].split('\t')[3],line.split('\t')[3])
        #else:
            #print anna[circ]
        ##print "M",line
        #N['overlap'] += 1        
        #del anna[circ]
    #else:
        #N['input2_not_in_input1'] += 1
    ##print len(anna.keys())
        
#for k,l in anna.items():
    ##if "HEK" in l:
        #print "MISSING\t%s" % l
        #N['input1_not_in_input2'] += 1

#for k in sorted(N.keys()):
    #sys.stderr.write("%s\t%d\n" % (k,N[k]))
        
#found = N['overlap']
#detected = N['unique_input2']
#total = N['unique_input1']
#fp = N['input2_not_in_input1']

#print "#sensitivity %d/%d = %.2f %%" % (found,total,float(found)/total*100)
#print "#FDR %d/%d = %.2f %%" % (fp,detected,float(fp)/detected*100)
