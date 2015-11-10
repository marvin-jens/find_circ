#!/usr/bin/env python
import pysam
import numpy
import os,sys,re
from logging import error, warning, debug
import logging
#logging.basicConfig(level=logging.DEBUG)
import mmap
import optparse
import traceback

COMPLEMENT = {
    'a' : 't',
    't' : 'a',
    'c' : 'g',
    'g' : 'c',
    'k' : 'm',
    'm' : 'k',
    'r' : 'y',
    'y' : 'r',
    's' : 's',
    'w' : 'w',
    'b' : 'v',
    'v' : 'b',
    'h' : 'd',
    'd' : 'h',
    'n' : 'n',
    'A' : 'T',
    'T' : 'A',
    'C' : 'G',
    'G' : 'C',
    'K' : 'M',
    'M' : 'K',
    'R' : 'Y',
    'Y' : 'R',
    'S' : 'S',
    'W' : 'W',
    'B' : 'V',
    'V' : 'B',
    'H' : 'D',
    'D' : 'H',
    'N' : 'N',
}

def complement(s):
    return "".join([COMPLEMENT[x] for x in s])

def rev_comp(seq):
    return complement(seq)[::-1]

class mmap_fasta(object):
    def __init__(self,fname):
        f = file(fname)
        header = f.readline()
        row = f.readline()

        self.ofs = len(header)
        self.lline = len(row)
        self.ldata = len(row.strip())
        self.skip = self.lline-self.ldata
        self.skip_char = row[self.ldata:]
        self.mmap = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)

    def __getslice__(self,start,end):
        l_start = start / self.ldata
        l_end = end / self.ldata
        ofs_start = l_start * self.skip + start + self.ofs
        ofs_end = l_end * self.skip + end + self.ofs
        
        s = self.mmap[ofs_start:ofs_end].replace(self.skip_char,"")
        L = end-start
        if len(s) == L:
            return s
        else:
            return s+"N"*(L-len(s))
        return 

class indexed_fasta(object):
    def __init__(self,fname,split_chrom="",**kwargs):
        self.fname = fname
        self.chrom_stats = {}
        self.split_chrom = split_chrom
        # try to load index
        ipath = fname + '.byo_index'
        if os.access(ipath,os.R_OK):
            self.load_index(ipath)
        else:
            self.index()
            self.store_index(ipath)
                
        f = file(fname)
        self.mmap = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)

    def index(self):
        debug("# indexed_fasta.index('%s')" % self.fname)

        ofs = 0
        f = file(self.fname)
        chrom = "undef"
        chrom_ofs = 0
        size = 0
        nl_char = 0

        for line in f:
            ofs += len(line)
            if line.startswith('>'):
                # store size of previous sequence
                if size: self.chrom_stats[chrom].append(size)

                chrom = line[1:].split()[0].strip()
                if self.split_chrom:
                    # use this to strip garbage from chrom/contig name like chr1:new 
                    # ->split_chrom=':' -> chrom=chr1
                    chrom = chrom.split(self.split_chrom)[0]
                chrom_ofs = ofs
            else:
                if not chrom in self.chrom_stats:
                    # this is the first line of the new chrom
                    size = 0
                    lline = len(line)
                    ldata = len(line.strip())
                    nl_char = lline - ldata
                    self.chrom_stats[chrom] = [chrom_ofs,ldata,nl_char,line[ldata:]]
                size += len(line.strip())

        # store size of previous sequence
        if size: self.chrom_stats[chrom].append(size)
                        
        f.close()

    def store_index(self,ipath):
        debug("# indexed_fasta.store_index('%s')" % ipath)
        
        # write to tmp-file first and in the end rename in order to have this atomic 
        # otherwise parallel building of the same index may screw it up.
        
        import tempfile
        tmp = tempfile.NamedTemporaryFile(mode="w",dir = os.path.dirname(ipath),delete=False)
        for chrom in sorted(self.chrom_stats.keys()):
            ofs,ldata,skip,skipchar,size = self.chrom_stats[chrom]
            tmp.write("%s\t%d\t%d\t%d\t%r\t%d\n" % (chrom,ofs,ldata,skip,skipchar,size))
        
        # make sure everything is on disk
        os.fsync(tmp)
        tmp.close()
        
        # make it accessible to everyone
        import stat
        os.chmod(tmp.name, stat.S_IROTH | stat.S_IRGRP | stat.S_IRUSR)
        
        # this is atomic on POSIX as we have created tmp in the same directory, 
        # therefore same filesystem
        os.rename(tmp.name,ipath)
        

    def load_index(self,ipath):
        debug("# indexed_fasta.load_index('%s')" % ipath)
        self.chrom_stats = {}
        for line in file(ipath):
            chrom,ofs,ldata,skip,skipchar,size = line.rstrip().split('\t')
            self.chrom_stats[chrom] = (int(ofs),int(ldata),int(skip),skipchar[1:-1].decode('string_escape'),int(size))
        
    def get_data(self,chrom,start,end,sense):
        if not self.chrom_stats:
            self.index()

        ofs,ldata,skip,skip_char,size = self.chrom_stats[chrom]
        pad_start = 0
        pad_end = 0
        if start < 0:
            pad_start = -start
            start = 0

        if end > size:
            pad_end = end - size
            end = size

        l_start = start / ldata
        l_end = end / ldata
        ofs_start = l_start * skip + start + ofs
        ofs_end = l_end * skip + end + ofs
        
        s = self.mmap[ofs_start:ofs_end].replace(skip_char,"")
        if pad_start or pad_end:
            s = "N"*pad_start + s + "N"*pad_end

        if sense == "-":
            s = rev_comp(s)
        return s

class Accessor(object):
    supports_write = False

    def __init__(self,path,chrom,sense,sense_specific=False,**kwargs):
        if sense_specific:
            self.covered_strands = [chrom+sense]
        else:
            self.covered_strands = [chrom+'+',chrom+'-']

    def get_data(self,chrom,start,end,sense,**kwargs):
        return []

    def get_oriented(self,chrom,start,end,sense,**kwargs):
        data = self.get_data(chrom,start,end,sense,**kwargs)
        if sense == "-": #  self.sense_specific and
            return data[::-1]
        else:
            return data

    def get_sum(self,chrom,start,end,sense,**kwargs):
        return self.get_data(chrom,start,end,sense,**kwargs).sum()

    def flush(self):
        pass

def to_bool(obj):
    if obj == "False":
        return False
    else:
        return bool(obj)
        
class Track(object):
    """
    Abstraction of chromosome-wide adressable data like sequences, coverage, scores etc.
    Actual access to the data is delegated to accessor objects which are instantiated on-the-fly for
    each chromosome (strand) upon first access and then cached.
    Use of mmap for the accessors is recommended and implemented for sequences and numpy (C-type)
    arrays.

    See io/track_accessors.py for more examples.
    """

    def __init__(self,path,accessor,sense_specific=True,description="unlabeled track",system="hg18",dim=1,auto_flush=False,mode="r",**kwargs):
        self.path = path
        self.mode = mode
        self.acc_cache = {}
        self.accessor = accessor
        self.kwargs = kwargs
        self.sense_specific = to_bool(sense_specific)
        self.dim = int(dim)
        self.description = description
        self.auto_flush = auto_flush
        self.last_chrom = ""
        self.logger = logging.getLogger("Track('%s')" % path)
        self.system = system

        self.logger.debug("Track(auto_flush=%s)" % (str(auto_flush)))
        kwargs['sense_specific'] = self.sense_specific
        kwargs['mode'] = self.mode
        kwargs['system'] = self.system
        kwargs['description'] = self.description
        kwargs['dim'] = self.dim

    def load(self,chrom,sense):
        # automatically flush buffers whenever a new chromosome is seen. reduces memory-footprint for sorted input data
        if self.auto_flush and chrom != self.last_chrom:
            self.logger.debug("Seen new chromosome %s. Flushing accessor caches." % chrom)
            self.flush_all()
        self.last_chrom = chrom

        ID = chrom+sense
        if not ID in self.acc_cache:
            self.logger.debug("Cache miss for %s%s. creating new accessor" % (chrom,sense))
            acc = self.accessor(self.path,chrom,sense,**(self.kwargs))
            
            # register the accessor for as many chromosomes/strands/contigs as it feels responsible for.
            self.logger.debug("Accessor covers strands '%s'" % str(sorted(set(acc.covered_strands))))
            if acc.covered_strands == '*':
                # hint that this accessors covers the complete genome
                self.logger.debug("disabling auto_flush and registering accessor for everything")
                self.acc_cache = DummyCache(acc)
                self.auto_flush = False
            else:
                for ID in acc.covered_strands:
                    self.acc_cache[ID] = acc

        return self.acc_cache[ID]

    def flush(self,chrom,sense):
        ID = self.get_identifier(chrom,sense)
        if ID in self.acc_cache:
            self.logger.warning("Flushing %s%s" % (chrom,sense))
            del self.acc_cache[ID]

    def flush_all(self):
        for a in self.acc_cache.values():
            a.flush()
        self.acc_cache = {}

    def get(self,chrom,start,end,sense,**kwargs):
        acc = self.load(chrom,sense)
        return acc.get_data(chrom,start,end,sense,**kwargs)

    def get_oriented(self,chrom,start,end,sense,**kwargs):
        acc = self.load(chrom,sense)
        return acc.get_oriented(chrom,start,end,sense,**kwargs)

    def get_sum(self,chrom,start,end,sense):
        #print "TRACK.GETSUM"
        acc = self.load(chrom,sense)
        return acc.get_sum(chrom,start,end,sense)
        
    def get_identifier(self,chrom,sense):
        if self.sense_specific:
            return chrom+sense
        else:
            return chrom

        
class GenomeAccessor(Accessor):
    def __init__(self,path,chrom,sense,system='hg19',**kwargs):
        super(GenomeAccessor,self).__init__(path,chrom,sense,system=system,**kwargs)
        debug("# GenomeAccessor mmap: Loading genomic sequence for chromosome %s from '%s'" % (chrom,path))

        self.system = system
        self.data = None
        fname = os.path.join(path)
        try:
            self.data = indexed_fasta(fname)
        except IOError:
            warning("Could not access '%s'. Switching to dummy mode (only Ns)" % fname)
                
            self.get_data = self.get_dummy
            self.get_oriented = self.get_dummy
            self.covered_strands = [chrom+'+',chrom+'-']
        else:
            # register for all chroms/strands
            self.covered_strands = [chrom+'+' for chrom in self.data.chrom_stats.keys()] + [chrom+'-' for chrom in self.data.chrom_stats.keys()]

        # TODO: maybe remove this if not needed
        self.get = self.get_oriented

    def load_indexed(self,path):
        ipath = path+'.index'
        if not os.access(ipath,os.R_OK):
            index = self.build_index(path,ipath)
        else:
            index = self.load_index(ipath)

        self.chrom_ofs = index.chrom_ofs
        
    def get_data(self,chrom,start,end,sense):
        seq = self.data.get_data(chrom,start,end,"+")
            
        if sense == "-":
            seq = complement(seq)

        return seq

    def get_dummy(self,chrom,start,end,sense):
        return "N"*int(end-start)



usage = """

   bowtie2 [mapping options] anchors.fastq.gz | %prog [options] > candidates.bed
"""

parser = optparse.OptionParser(usage=usage)
parser.add_option("-v","--version",dest="version",action="store_true",default=False,help="get version information")
parser.add_option("-S","--system",dest="system",type=str,default="",help="model system database (optional! Requires byo library.)")
parser.add_option("-G","--genome",dest="genome",type=str,default="",help="path to genome (either a folder with chr*.fa or one multichromosome FASTA file)")
parser.add_option("-n","--name",dest="name",default="unknown",help="tissue/sample name to use (default='unknown')")
parser.add_option("-p","--prefix",dest="prefix",default="",help="prefix to prepend to each junction name (default='')")
parser.add_option("-q","--min_uniq_qual",dest="min_uniq_qual",type=int,default=2,help="minimal uniqness for anchor alignments (default=2)")
parser.add_option("-a","--anchor",dest="asize",type=int,default=20,help="anchor size (default=20)")
parser.add_option("-m","--margin",dest="margin",type=int,default=2,help="maximum nts the BP is allowed to reside inside an anchor (default=2)")
parser.add_option("-d","--maxdist",dest="maxdist",type=int,default=2,help="maximum mismatches (no indels) allowed in anchor extensions (default=2)")

parser.add_option("","--noncanonical",dest="noncanonical",default=False,action="store_true",help="relax the GU/AG constraint (will produce many more ambiguous counts)")
parser.add_option("","--randomize",dest="randomize",default=False,action="store_true",help="select randomly from tied, best, ambiguous hits")
parser.add_option("","--allhits",dest="allhits",default=False,action="store_true",help="in case of ambiguities, report each hit")
parser.add_option("","--stranded",dest="stranded",default=False,action="store_true",help="use if the reads are stranded. By default it will be used as control only, use with --strandpref for breakpoint disambiguation.")
parser.add_option("","--strandpref",dest="strandpref",default=False,action="store_true",help="prefer splice sites that match annotated direction of transcription")
parser.add_option("","--halfunique",dest="halfunique",default=False,action="store_true",help="also report junctions where only one anchor aligns uniquely (less likely to be true)")
parser.add_option("","--report_nobridges",dest="report_nobridges",default=False,action="store_true",help="also report junctions lacking at least a single read where both anchors, jointly align uniquely (not recommended. Much less likely to be true.)")

parser.add_option("-R","--reads",dest="reads",default="",help="write spliced reads to this file instead of stderr (RECOMMENDED!)")
parser.add_option("-B","--bam",dest="bam",default="",help="filename to store anchor alignments that were recorded as linear or circular junction candidates")
parser.add_option("-r","--reads2samples",dest="reads2samples",default="",help="path to tab-separated two-column file with read-name prefix -> sample ID mapping")
parser.add_option("-s","--stats",dest="stats",default="runstats.log",help="write numeric statistics on the run to this file")
options,args = parser.parse_args()

if options.version:
    print """find_circ.py version 1.2\n\n(c) Marvin Jens 2012-2015.\nCheck http://www.circbase.org for more information."""
    sys.exit(0)

if options.system:
    import importlib
    system = importlib.import_module("byo.systems.{options.system}".format(**locals()))
    genome = system.genome

if options.genome:
    genome = Track(options.genome,accessor=GenomeAccessor)
    
if not (options.genome or options.system):
    error("need to specify either model system database (-S) or genome (-G).")
    sys.exit(1)

from itertools import izip_longest
def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

if options.reads2samples:
    samples = [line.rstrip().split('\t') for line in file(options.reads2samples)]
else:
    samples = []

if options.reads:
    readfile = file(options.reads,"w")
else:
    readfile = sys.stderr

samples.append(('',options.name))
    
from collections import defaultdict

minmapscore = options.asize * (-2)

class Hit(object):
    def __init__(self):
        self.reads = []
        self.readnames = []
        self.uniq = set()
        self.mapquals_A = []
        self.mapquals_B = []
        self.uniq_bridges = 0
        self.tissues = defaultdict(int)
        self.edits = []
        self.overlaps = []
        self.n_hits = []
        self.signal = "NNNN"
        self.strand_plus = 0
        self.strand_minus = 0
        self.strandmatch = 'NA'
        
    def add(self,read,A,B,dist,ov,strandmatch,signal,n_hits):
        self.signal = signal
        self.strandmatch = strandmatch      
        self.edits.append(dist)
        self.overlaps.append(ov)
        self.n_hits.append(n_hits)

        # by convention have A precede B in the genome.
        if A.pos > B.pos:
            A,B = B,A
            
        # Alignment Score - Secondbest hit score            
        aopt = dict(A.tags)
        bopt = dict(B.tags)
        qA = aopt.get('AS') - aopt.get('XS',minmapscore)
        qB = bopt.get('AS') - bopt.get('XS',minmapscore)

        if qA and qB:
            # both anchors from the *same read* align uniquely
            self.uniq_bridges += 1

        self.mapquals_A.append(qA)
        self.mapquals_B.append(qB)

        # recover the original readname 
        # ('__' is forbidden in input read names!)
        if '__' in A.qname:
            qname = A.qname.split('__')[0][:-2]
        else: # reads have been swapped at some point
            qname = B.qname.split('__')[0][:-2]

        self.readnames.append(qname)
        
        # record the spliced read sequence as it was before mapping
        if A.is_reverse:
            self.strand_minus += 1
            self.reads.append(rev_comp(read))
        else:
            self.strand_plus += 1
            self.reads.append(read)

        # identify the tissue/sample it came from
        sample_name = options.name
        for (prefix,tiss) in samples:
            if qname.startswith(prefix):
                sample_name = tiss
                break

        self.tissues[sample_name] += 1
        
        self.uniq.add((read,sample_name))
        self.uniq.add((rev_comp(read),sample_name))


    def scores(self,chrom,start,end,sense):
        n_reads = len(self.reads)
        n_uniq = len(self.uniq) / 2
        #print sorted(self.mapquals_A,reverse=True)
        #print sorted(self.mapquals_B,reverse=True)
        best_qual_A = sorted(self.mapquals_A,reverse=True)[0]
        best_qual_B = sorted(self.mapquals_B,reverse=True)[0]

        #print self.edits,self.overlaps,self.n_hits
        tissues = sorted(self.tissues.keys())
        tiss_counts = [str(self.tissues[k]) for k in tissues]
        return (n_reads,n_uniq,best_qual_A,best_qual_B,self.uniq_bridges,tissues,tiss_counts,min(self.edits),min(self.overlaps),min(self.n_hits),self.signal,self.strandmatch)
                
        
loci = defaultdict(list)
circs = defaultdict(Hit)
splices = defaultdict(Hit)

N = defaultdict(int)

from numpy import chararray as carray
from numpy import fromstring,byte

def find_breakpoints(A,B,read,chrom,margin=options.margin,maxdist=options.maxdist):

    def mismatches(a,b):
        a,b = fromstring(a,dtype=byte), fromstring(b,dtype=byte)
        
        return (a != b).sum()

    def rnd():
        return numpy.random.random()
    if not options.randomize:
        rnd = lambda : 0

    def strandmatch(ann,sense):
        if ann == sense:
            return "MATCH"
        if ann == "*" or len(ann) > 1:
            return "NA"
        return "MISMATCH"

    L = len(read)
    hits = []
    #print "readlen",L
    #print " "*2+read
    eff_a = options.asize-margin
    internal = read[eff_a:-eff_a].upper()
        
    flank = L - 2*eff_a + 2

    A_flank = genome.get(chrom,A.aend-margin,A.aend-margin + flank,'+').upper()
    B_flank = genome.get(chrom,B.pos - flank+margin,B.pos+margin,'+').upper()

    l = L - 2*eff_a
    for x in range(l+1):
        spliced = A_flank[:x] + B_flank[x+2:]
        dist = mismatches(spliced,internal)        
        
        #bla = A_flank[:x].lower() + B_flank[x+2:]
        #print " "*(eff_a+2)+bla,dist

        ov = 0
        if x < margin:
            ov = margin-x
        if l-x < margin:
            ov = margin-(l-x)
        
        if dist <= maxdist:
            gt = A_flank[x:x+2]
            ag = B_flank[x:x+2]
            gtag = gt+ag
            rc_gtag = rev_comp(gtag)
            
            start,end = B.pos+margin-l+x,A.aend-margin+x+1
            start,end = min(start,end),max(start,end)
            
            strand = "*"
            
            # get strand cues from read if strand-specific sequencing was used.
            if options.stranded:
                if A.is_reverse:
                    strand = '-'
                else:
                    strand = '+'

            if options.noncanonical:
                hits.append((dist,ov,strandmatch(strand,'+'),rnd(),chrom,start,end,gtag,'+'))
                hits.append((dist,ov,strandmatch(strand,'-'),rnd(),chrom,start,end,rc_gtag,'-'))

            else:
                if gtag == 'GTAG':
                    hits.append((dist,ov,strandmatch(strand,'+'),rnd(),chrom,start,end,'GTAG','+'))
                elif gtag == 'CTAC':
                    hits.append((dist,ov,strandmatch(strand,'-'),rnd(),chrom,start,end,'GTAG','-'))

    if len(hits) < 2:
        # unambiguous, return right away
        return hits

    # Hits are sorted, with low edit distance beating low anchor overlap
    hits = sorted(hits)
    # return only hits that are tied with the best candidate by edit distance, anchor overlap 
    # and strand preference (if requested). Ties are either broken by random draw or by selecting 
    # the first hit later, but they are kept here such that they are counted properly.
    best = hits[0]
    if options.strandpref:
        # exploit strand information to break ties only if requested. Otherwise it's only used as control.
        ties = [h for h in hits if (h[0] == best[0]) and (h[1] == best[1]) and (h[2] == best[2])]
    else:
        ties = [h for h in hits if (h[0] == best[0]) and (h[1] == best[1])]
        
    return ties
    
sam = pysam.Samfile('-','r')
if options.bam:
    bam_out = pysam.Samfile(options.bam,'wb',template=sam)
else:
    bam_out = None

try:
    for pair_num,(A,B) in enumerate(grouper(2,sam)):
        #print A
        #print B
        N['total'] += 1
        if A.is_unmapped or B.is_unmapped:
            N['unmapped'] += 1
            continue
        if A.tid != B.tid:
            N['other_chrom'] += 1
            continue
        if A.is_reverse != B.is_reverse:
            N['other_strand'] += 1
            continue

        dist = B.pos - A.pos
        if numpy.abs(dist) < options.asize:
            N['overlapping_anchors'] += 1
            continue
        
        ### anchor pairs that make it up to here are interesting
        if bam_out:
            bam_out.write(A)
            bam_out.write(B)

        #debug("A='%s' B='%s' dist=%d A.is_reverse=%s" % (A,B,dist,A.is_reverse))
        if (A.is_reverse and dist > 0) or (not A.is_reverse and dist < 0):
            # the anchors align in reversed orientation -> circRNA?
            
            read = A.qname.split('__')[1]
            chrom = sam.getrname(A.tid)
            
            if A.is_reverse:
                #print "ISREVERSE"
                A,B = B,A
                read = rev_comp(read)
                            
            bp = find_breakpoints(A,B,read,chrom)
            if not bp:
                N['circ_no_bp'] += 1
            else:
                N['circ_reads'] += 1

            n_hits = len(bp)
            if bp and not options.allhits:
                bp = [bp[0],]

            for h in bp:
                # for some weird reason for circ we need a correction here
                dist,ov,strandmatch,rnd,chrom,start,end,signal,sense = h
                h = (chrom,start+1,end-1,sense)
                circs[h].add(read,A,B,dist,ov,strandmatch,signal,n_hits)

        elif (A.is_reverse and dist < 0) or (not A.is_reverse and dist > 0):
            # the anchors align sequentially -> linear/normal splice junction?
            read = A.qname.split('__')[1]
            chrom = sam.getrname(A.tid)
            
            if A.is_reverse:
                #print "ISREVERSE"
                A,B = B,A
                read = rev_comp(read)
                            
            bp = find_breakpoints(A,B,read,chrom)
            if not bp:
                N['splice_no_bp'] += 1
            else:
                N['spliced_reads'] += 1
            n_hits = len(bp)
            if bp and not options.allhits:
                bp = [bp[0],]
            
            for h in bp:
                #print h
                dist,ov,strandmatch,rnd,chrom,start,end,signal,sense = h
                h = (chrom,start,end,sense)
                splices[h].add(read,A,B,dist,ov,strandmatch,signal,n_hits)
                
                # remember the spliced reads at these sites
                loci[(chrom,start,sense)].append(splices[h])
                loci[(chrom,end,sense)].append(splices[h])
        else:
            N['fallout'] += 1
            warning("unhandled read: A='%s' B='%s'" % (A,B))
            
except KeyboardInterrupt:
    sam_line = pair_num * 2
    fastq_line = pair_num * 8

    logging.warning("KeyboardInterrupt by user while processing alignment pair {pair_num}, on input starting at SAM line {sam_line}, FASTQ line {fastq_line}".format(**locals()))
except:        
    sam_line = pair_num * 2
    fastq_line = pair_num * 8

    logging.error("Unhandled exception raised while processing alignment pair {pair_num}, on input starting at SAM line {sam_line}, FASTQ line {fastq_line}".format(**locals()))
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)

def output(cand,prefix):
    print "#","\t".join(['chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges','best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap','breakpoints','signal','strandmatch','category'])
    n = 1
    for c,hit in cand.items():
        #print c
        chrom,start,end,sense = c
        n_reads,n_uniq,best_qual_A,best_qual_B,uniq_bridges,tissues,tiss_counts,min_edit,min_anchor_ov,n_hits,signal,strandmatch = hit.scores(chrom,start,end,sense)
        
        if options.halfunique:
            if (best_qual_A < options.min_uniq_qual) and (best_qual_B < options.min_uniq_qual):
                N['anchor_not_uniq'] += 1
                continue
        else:
            if (best_qual_A < options.min_uniq_qual) or (best_qual_B < options.min_uniq_qual):
                N['anchor_not_uniq'] += 1
                continue

        if (uniq_bridges == 0) and not options.report_nobridges:
                N['no_uniq_bridges'] += 1
                continue

        name = "%s%s_%06d" % (options.prefix,prefix,n)
        n += 1
        #sys.stderr.write("%s\t%s\n" % (name,"\t".join(sorted(reads))))
        for r_seq,ori_name in zip(hit.reads,hit.readnames):
            readfile.write(">%s %s\n%s\n" % (name,ori_name,r_seq))
      
        categories = []
        if signal == "GTAG":
            categories.append("CANONICAL")
        if strandmatch == "MATCH":
            categories.append("STRANDMATCH")
        if best_qual_A > 0 and best_qual_B > 0 and uniq_bridges > 0:
            categories.append("ANCHOR_UNIQUE")
        if uniq_bridges == 0:
            categories.append("NO_UNIQ_BRIDGES")          
        if n_hits == 1:
            categories.append("UNAMBIGUOUS_BP")
        if min_anchor_ov == 0 and min_edit == 0: 
            categories.append("PERFECT_EXT")
        elif min_anchor_ov <= 1 and min_edit <= 1:
            categories.append("GOOD_EXT")
        elif min_anchor_ov <= 2 and min_edit <= 2:
            categories.append("OK_EXT")

        if not categories:
            categories.append("DUBIOUS")

        if prefix == "circ":
            categories.append("CIRCULAR")
        else:
            categories.append("LINEAR")


        bed = [chrom,start-1,end,name,n_reads,sense,n_uniq,uniq_bridges,best_qual_A,best_qual_B,",".join(tissues),",".join(tiss_counts),min_edit,min_anchor_ov,n_hits,signal,strandmatch,",".join(sorted(categories))]
        print "\t".join([str(b) for b in bed])

stats = file(options.stats,"w")
stats.write(str(N)+"\n")

output(circs,"circ")
output(splices,"norm")
