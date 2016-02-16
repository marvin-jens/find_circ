#!/usr/bin/env python
import os,sys,re
from collections import defaultdict
import numpy as np
from numpy import chararray as carray
from numpy import fromstring,byte
import mmap
import logging
import pysam
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

def k_mers(k,alphabet=['A','C','G','T',]):
    if k == 1:
        prefix = [""]
    else:
        prefix = k_mers(k-1,alphabet)

    for pre in prefix:
        for a in alphabet:
            yield pre+a

# precompute RC of all possible 4-mers for fast
# splice signal check later on.
fast_4mer_RC = {}
for mer in k_mers(4,alphabet=['A','C','G','T','N']):
    fast_4mer_RC[mer] = rev_comp(mer)

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
        self.logger = logging.getLogger("indexed_fasta")
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
        self.logger.info("# indexed_fasta.index('%s')" % self.fname)

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
        self.logger.info("# indexed_fasta.store_index('%s')" % ipath)
        
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
        self.logger.info("# indexed_fasta.load_index('%s')" % ipath)
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
        self.sense_specific = sense_specific
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
        self.logger = logging.getLogger("GenomeAccessor")
        self.logger.info("# mmap: Loading genomic sequence for chromosome %s from '%s'" % (chrom,path))

        self.system = system
        self.data = None
        fname = os.path.join(path)
        try:
            self.data = indexed_fasta(fname)
        except IOError:
            self.logger.warning("Could not access '%s'. Switching to dummy mode (only Ns)" % fname)
                
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
   bwa mem -t<threads> [-p] -A2 -B10 -k 15 -T 1 $GENOME_INDEX reads.fastq.gz | %prog [options]
   
   OR:
   
   %prog [options] <bwa_mem_genome_alignments.bam>
"""

parser = optparse.OptionParser(usage=usage)
parser.add_option("-v","--version",dest="version",action="store_true",default=False,help="get version information")
parser.add_option("-S","--system",dest="system",type=str,default="",help="model system database (optional! Requires byo library.)")
parser.add_option("-G","--genome",dest="genome",type=str,default="",help="path to genome (either a folder with chr*.fa or one multichromosome FASTA file)")
parser.add_option("-o","--output",dest="output",default="find_circ_run",help="where to store output")
parser.add_option("","--stdout",dest="stdout",default=None,choices=['sites','reads','stats','multi'],help="use to direct chosen type of output (sites, reads, stats, multi) to stdout instead of file")
parser.add_option("-n","--name",dest="name",default="unknown",help="tissue/sample name to use (default='unknown')")
#parser.add_option("-p","--prefix",dest="prefix",default="",help="prefix to prepend to each junction name (default='')")
parser.add_option("-q","--min-uniq-qual",dest="min_uniq_qual",type=int,default=1,help="minimal uniqness for anchor alignments (default=2)")
parser.add_option("-a","--anchor",dest="asize",type=int,default=15,help="anchor size (default=15)")
parser.add_option("-m","--margin",dest="margin",type=int,default=2,help="maximum nts the BP is allowed to reside within a segment (default=2)")
parser.add_option("-d","--max-mismatch",dest="maxdist",type=int,default=2,help="maximum mismatches (no indels) allowed in segment extensions (default=2)")
parser.add_option("","--short-threshold",dest="short_threshold",type=int,default=100,help="minimal genomic span [nt] of a circRNA before it is labeled SHORT (default=100)")
parser.add_option("","--huge-threshold",dest="huge_threshold",type=int,default=100000,help="maximal genomic span [nt] of a circRNA before it is labeled HUGE (default=100000)")
parser.add_option("","--debug",dest="debug",default=False,action="store_true",help="Activate LOTS of debug output")
parser.add_option("","--profile",dest="profile",default=False,action="store_true",help="Activate performance profiling")

parser.add_option("","--non-canonical",dest="noncanonical",default=False,action="store_true",help="relax the GU/AG constraint (will produce many more ambiguous counts)")
parser.add_option("","--randomize",dest="randomize",default=False,action="store_true",help="select randomly from tied, best, ambiguous hits")
parser.add_option("","--all-hits",dest="allhits",default=False,action="store_true",help="in case of ambiguities, report each hit")
parser.add_option("","--stranded",dest="stranded",default=False,action="store_true",help="use if the reads are stranded. By default it will be used as control only, use with --strand-pref for breakpoint disambiguation.")
parser.add_option("","--strand-pref",dest="strandpref",default=False,action="store_true",help="prefer splice sites that match annotated direction of transcription")
parser.add_option("","--half-unique",dest="halfunique",default=False,action="store_true",help="also report junctions where only one anchor aligns uniquely (less likely to be true)")
parser.add_option("","--report-nobridges",dest="report_nobridges",default=False,action="store_true",help="also report junctions lacking at least a single read where both anchors, jointly align uniquely (not recommended. Much less likely to be true.)")

#parser.add_option("-R","--reads",dest="reads",default=False,help="write spliced reads to this file instead of stderr (RECOMMENDED!)")
parser.add_option("-B","--bam",dest="bam",default=False,action="store_true",help="store anchor alignments that were recorded as linear or circular junction candidates")

parser.add_option("-t","--throughput",dest="throughput",default=False,action="store_true",help="print information on throughput to stderr (useful for benchmarking)")
parser.add_option("","--chunk-size",dest="chunksize",type=int,default=100000,help="number of reads to be processed in one chunk (default=100000)")
parser.add_option("","--noop",dest="noop",default=False,action="store_true",help="Do not search for any junctions. Only process the alignment stream (useful for benchmarking)")
parser.add_option("","--no-linear",dest="nolinear",default=False,action="store_true",help="Do not investigate linear junctions, unless associated with another backsplice event (saves some time)")
options,args = parser.parse_args()



if options.version:
    print """find_circ.py version 1.99\n\n(c) Marvin Jens 2012-2016.\nCheck http://www.circbase.org for more information."""
    sys.exit(0)

if options.system:
    import importlib
    system = importlib.import_module("byo.systems.{options.system}".format(**locals()))
    genome = system.genome

if options.genome:
    genome = Track(options.genome,accessor=GenomeAccessor)
    
if not (options.genome or options.system):
    error("need to specify either model system database (-S) or genome FASTA file (-G).")
    sys.exit(1)

# prepare logging system
logging.basicConfig(level=logging.INFO,filename=os.path.join(options.output,"run.log"))
logger = logging.getLogger('find_circ')

# prepare output files
if not os.path.isdir(options.output):
    os.mkdir(options.output)

sites_file = file(os.path.join(options.output,"splice_sites.bed"),"w")
reads_file = file(os.path.join(options.output,"spliced_reads.fa"),"w")
stats_file = file(os.path.join(options.output,"stats.log"),"w")
multi_file = file(os.path.join(options.output,"multi_events.tsv"),"w")

# redirect output to stdout, if requested
if options.stdout:
    varname = "{0}_file".format(options.stdout)
    out_file = globals()[varname]
    out_file.write('# redirected to stdout\n')
    logger.info('redirected {0} to stdout'.format(options.stdout))
    globals()[varname] = sys.stdout
    

if args:
    logger.info('reading from {0}'.format(args[0]))
    sam_input = pysam.Samfile(args[0],'rb')
else:
    logger.info('reading from stdin')
    sam_input = pysam.Samfile('-','r')

tid_cache = {}
def fast_chrom_lookup(read):
    tid = read.tid
    if not tid in tid_cache: 
        tid_cache[tid] = sam_input.getrname(tid)

    return tid_cache[tid]

if options.bam:
    bam_filename = os.path.join(options.output,"spliced_alignments.bam")
    bam_out = pysam.Samfile(bam_filename,'wb',template=sam_input)
else:
    bam_out = None

minmapscore = options.asize * (-2)

class Hit(object):
    def __init__(self):
        self.reads = []
        self.counts = 0.
        self.readnames = []
        self.uniq = set()
        self.mapquals_A = []
        self.mapquals_B = []
        self.uniq_bridges = 0.
        self.tissues = defaultdict(int)
        self.edits = []
        self.overlaps = []
        self.n_hits = []
        self.signal = "NNNN"
        self.strand_plus = 0
        self.strand_minus = 0
        self.strandmatch = 'NA'
        self.flags = defaultdict(int)
        self.read_flags = defaultdict(set)
        
    def add_flag(self,flag,frag_name):
        self.flags[flag] += 1
        self.read_flags[frag_name].add(flag)

    def get_flags_counts(self):
        if not self.flags:
            return ["N/A"],[0]
        
        flags = sorted(self.flags)
        counts = [self.flags[f] for f in flags]
        
        return flags,counts
        
    #def add(self,read,A,B,dist,ov,strandmatch,signal,n_hits,weight):
    def add(self,splice):
        self.signal = splice.gtag
        
        # TODO: this belongs to output, not here!
        self.strandmatch = 'N/A'
        if options.stranded:
            if splice.strandmatch:
                self.strandmatch = "MATCH"
            else:
                self.strandmatch = "MISMATCH"

        self.edits.append(splice.dist)
        self.overlaps.append(splice.ov)
        self.n_hits.append(splice.n_hits)
        self.counts += splice.junc_span.weight

        # TODO: Move this logic into JunctionSpan, bc it is closer to the underlying alignments
        # by convention have A precede B in the genome.
        A = splice.junc_span.align_A
        B = splice.junc_span.align_B
        read = splice.junc_span.primary.seq
        weight = splice.junc_span.weight
        if splice.is_backsplice:
            A,B = B,A

        # Alignment Score - Secondbest hit score            
        aopt = dict(A.tags)
        bopt = dict(B.tags)
        qA = aopt.get('AS') - aopt.get('XS',minmapscore)
        qB = bopt.get('AS') - bopt.get('XS',minmapscore)

        if qA and qB:
            # both anchors from the *same read* align uniquely
            self.uniq_bridges += weight

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
            self.strand_minus += weight
            self.reads.append(rev_comp(read))
        else:
            self.strand_plus += weight
            self.reads.append(read)

        sample_name = options.name
        self.tissues[sample_name] += weight
        
        self.uniq.add((read,sample_name))
        self.uniq.add((rev_comp(read),sample_name))


    def scores(self,chrom,start,end,sense):
        n_spanned = len(self.reads)
        n_uniq = len(self.uniq) / 2
        #print sorted(self.mapquals_A,reverse=True)
        #print sorted(self.mapquals_B,reverse=True)
        best_qual_A = sorted(self.mapquals_A,reverse=True)[0]
        best_qual_B = sorted(self.mapquals_B,reverse=True)[0]

        #print self.edits,self.overlaps,self.n_hits
        tissues = sorted(self.tissues.keys())
        tiss_counts = [str(self.tissues[k]) for k in tissues]
        return (n_spanned,self.counts,n_uniq,best_qual_A,best_qual_B,self.uniq_bridges,tissues,tiss_counts,min(self.edits),min(self.overlaps),min(self.n_hits),self.signal,self.strandmatch)
                
circ_splices = defaultdict(Hit)
linear_splices = defaultdict(Hit)

N = defaultdict(float)

class MultiEventStorage(object):
    def __init__(self):
        self.circ_extra_events = defaultdict(list)
    
    def add(self,circ_key,lin_junctions,unspliced):
        extra_lin_events = sorted(lin_junctions)
        extra_unspliced = sorted(unspliced)

        self.circ_extra_events[circ_key].append( (extra_lin_events,extra_unspliced) )
               
    def output(self):
        lines = []
        for circ_key in sorted(self.circ_extra_events.keys()):
            circ_chrom,circ_start,circ_end,circ_sense = circ_key
            circ_cols = [str(c) for c in circ_key]
            
            for extra_lin_events,extra_unspliced in self.circ_extra_events[circ_key]:
           
                if not extra_lin_events:
                    lin_col = ["N/A"]
                else:
                    lin_col = [",".join([ "{0:d}-{1:d}".format(start,end) for chrom,start,end,sense in extra_lin_events])]
                    
                if not extra_unspliced:
                    un_col = ["N/A"]
                else:
                    unlist = []
                    for chrom,start,end,sense in extra_unspliced:
                        if chrom == circ_chrom and (start + options.asize > circ_start) and (end - options.asize < circ_end):
                            unlist.append( "{0:d}-{1:d}".format(start,end) )
                        else:
                            unlist.append( "[{0:s}:{1:d}-{2:d}]".format(chrom,start,end) )
                    un_col = [",".join(unlist)]
                    
                    
                lines.append("\t".join(circ_cols + lin_col + un_col))
                
        return "\n".join(lines)
                
multi_event_storage = MultiEventStorage()


class Splice(object):
    def __init__(self,junc_span,chrom,start,end,strand,dist,ov,gtag):
        self.junc_span = junc_span
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.dist = dist
        self.ov = ov
        self.gtag = gtag
        self.n_hits = 1 # if necessary, this is changed by JunctionSpan.find_breakpoints() after breaking ties

    @property
    def is_canonical(self):
        return self.gtag == 'GTAG'

    @property
    def is_backsplice(self):
        # the underlying alignment segment pair knows
        return self.junc_span.is_backsplice
    
    def __str__(self):
        "python magic to make nice debug output"
        return "Splice({self.chrom}:{self.start}-{self.end}:{self.strand} edits={self.dist} ov={self.ov} score={self.score} backsplice={self.is_backsplice})".format(**locals())

    @property
    def score(self):
        "Used to rank different possible splices. higher is better"
        # canonical beats low extension mismatches, beats low anchor overlap
        s = self.is_canonical*20 - self.dist * 10 - self.ov
        if options.strandpref:
            s += 100 * (self.strand == self.junc_span.strand)

        return s

    @property
    def coord(self):
        return (self.chrom,self.start,self.end,self.strand)

class JunctionSpan(object):
    def __init__(self,align_A,align_B,primary,q_start,q_end,weight):
        self.primary = primary
        self.align_A = align_A
        self.align_B = align_B
        self.q_start = q_start
        self.q_end = q_end
        self.weight = weight
        
        # TODO: fix, depending on mate orientation!
        if primary.is_reverse:
            self.strand = '-'
        else:
            self.strand = '+'

        # segments are ordered by occurence in the reqd query.
        # if this order is the same in the genome it is a linear (normal)
        # splice site. If reversed, it is a backsplice.
        self.dist = align_B.pos - align_A.aend
        full_read = primary.seq
        self.read_part = full_read[q_start:q_end]
    
    @property
    def is_backsplice(self):
        return self.dist < 0
    
    def find_breakpoints(self):
        """
        Extends distal parts of two segment alignments such that the intervening 
        read sequence is completely aligned and that the endpoints of this extension,
        corresponding to the breakpoint inside the read, are flanked by splice signal
        dinucleotides.
        """
        def mismatches(a,b):
            a,b = fromstring(a,dtype=byte), fromstring(b,dtype=byte)
            return (a != b).sum()

        def simple_match(a,b):
            return (a != b)
        
        if options.maxdist == 0:
            # faster match function if no mismatches are allowed anyway!
            mismatches = simple_match

        # short hands
        L = len(self.read_part)
        A = self.align_A
        B = self.align_B
        read = self.read_part
        margin = options.margin
        maxdist = options.maxdist
        is_backsplice = self.is_backsplice
        # effective anchor size is reduced by allowing the breakpoint to reside 
        # options.margin nt's inside the anchor
        eff_a = options.asize-options.margin

        # record all possible breakpoints
        hits = []
        if options.debug:
            print "readlen",L
            print " "*2+read
            print " "*2+A.query[:-options.margin]
            print " "*(2+L-eff_a)+B.query[options.margin:]
            print " "*2+A.seq[:eff_a]
            print " "*(2+L-eff_a)+A.seq[-eff_a:]

        # intervening sequence that needs to be exactly explained by the splicing
        internal = read[eff_a:-eff_a].upper()

        # this much needs to be retrieved from the genome
        chrom = fast_chrom_lookup(A)
        
        flank = L - 2*eff_a + 2
        A_flank = genome.get(chrom,A.pos + eff_a,A.pos + eff_a + flank,'+').upper()
        B_flank = genome.get(chrom,B.aend - eff_a - flank,B.aend - eff_a,'+').upper()

        l = L - 2*eff_a
        # all possible breakpoints
        for x in range(l+1):
            spliced = A_flank[:x] + B_flank[x+2:]
            dist = mismatches(spliced,internal)        
            
            if options.debug:
                bla = A_flank[:x].lower() + B_flank[x+2:]
                if options.debug:
                    print " "*(eff_a+2)+bla,dist
           
            if dist <= maxdist:
                # is the tested breakpoint inside an anchor?
                ov = 0
                if x < margin:
                    ov = margin-x
                if l-x < margin:
                    ov = margin-(l-x)

                gt = A_flank[x:x+2]
                ag = B_flank[x:x+2]
                gtag = gt+ag
                rc_gtag = fast_4mer_RC[gtag]
                
                start,end = B.aend-eff_a-l+x,A.pos+eff_a+x+1
                start,end = min(start,end),max(start,end)
                strand = "*"
                
                # get strand cues from read if strand-specific sequencing was used.
                # TODO: skip if gtag does not match strand!
                if options.stranded:
                    if A.is_reverse:
                        strand = '-'
                    else:
                        strand = '+'

                if is_backsplice:
                    # for circs we need to correct here
                    start += 1
                    end -= 1
                                
                if options.noncanonical:
                    hits.append( Splice(self,chrom,start,end,'+',dist,ov,gtag,) )
                    hits.append( Splice(self,chrom,start,end,'-',dist,ov,rc_gtag) )
                else:
                    if gtag == 'GTAG':
                        hits.append( Splice(self,chrom,start,end,'+',dist,ov,gtag) )
                    elif gtag == 'CTAC':
                        hits.append( Splice(self,chrom,start,end,'-',dist,ov,rc_gtag) )

        if options.debug:
            print hits

        if len(hits) < 2:
            # unambiguous, return right away
            return hits

        # Hits are sorted, with low edit distance beating low anchor overlap
        hits = sorted(hits,key=lambda x : x.score, reverse = True)
        best_score = hits[0].score
        ties = [h for h in hits if (h.score == best_score)]
        n_hits = len(ties)
        
        for h in hits:
            h.n_hits = n_hits

        return ties
    
class MateSegments(object):
    def __init__(self,primary):
        N['total_mates'] += 1
        self.primary = primary
        self.full_seq = primary.seq
        self.segments = [primary]
        self.proper_segs = [primary]
        self.other_chrom_segs = []
        self.other_strand_segs = []
        self.seg_flags = {}

    @property
    def matenum(self):
        if self.primary.is_read2:
            return 2
        else:
            # this also catches unpaired reads
            return 1
        
    @property
    def strand(self):
        if self.primary.is_reverse:
            return '-'
        else:
            return '+'
        
    def __str__(self):
        "Python magic to make nice debug output possible"
        buf = []
        chrom = fast_chrom_lookup(self.primary)

        buf.append("MateSegments({self.primary.qname} mate={self.matenum})".format(**locals()))
        buf.append("\tprimary={chrom}:{self.primary.pos}-{self.primary.aend}:{self.strand}".format(**locals()))
        n_proper = len(self.proper_segs)
        n_oc = len(self.other_chrom_segs)
        n_os = len(self.other_strand_segs)
        flags = ",".join(sorted(self.seg_flags.values()))
        buf.append("\t{n_proper} proper_segs {n_oc} other_chrom {n_os} other_strand. Flags={flags}".format(**locals()) )
        
        return "\n".join(buf)
            
    def is_segment(self,align):
        "tests if align is supplementary to the primary"
        if (align.is_read1 == self.primary.is_read1) and (align.qname == self.primary.qname):
            return True
        else:
            return False

    def is_other_mate(self,align):
        "tests if align is the other mate of the primary alignment."
        if (align.is_read1 != self.primary.is_read1) and (align.qname == self.primary.qname):
            return True
        else:
            return False

    def add_segment(self,align):
        "records align as a new supplementary alignment"
        self.segments.append(align)
        
        if align.tid != self.primary.tid:
            # secondary hit is on another chromosome
            self.other_chrom_segs.append(align)
            self.seg_flags[align] = "SEG_OTHER_CHROM"
            N['seg_other_chrom'] += 1

        elif align.is_reverse != self.primary.is_reverse:
            # secondary hit is on other strand
            self.other_strand_segs.append(align)
            self.seg_flags[align] = "SEG_OTHER_STRAND"
            N['seg_other_strand'] += 1
        else:
            # same chrom, same strand. could be useful
            self.proper_segs.append(align)
            N['seg_proper'] += 1

    def adjacent_segment_pairs(self):
        """
        This generator yields pairs of adjacent segments (meaning adjacent in 
        the reads query sequence) to be processed by the known logic for splice 
        junction determination.
        Input is the recorded list of proper segments as reported by BWA MEM.
        In the simplest case these are two alignments to two exons, with the
        parts belonging to the respective other exon clipped.
        However, for longer reads these may also be split into three or
        more segments.
        """

        if len(self.proper_segs) < 2:
            # unspliced read. Nothing to do
            return

        # the primary reported alignment for a read always contains the
        # full, original read sequence
        full_read = self.primary.seq
        
        # weight to assign to each splice. Now that a read can contain multiple
        # splices (or a circRNA read even the same splice multiple times), we
        # do not want to overcount.
        weight = 1./(len(self.proper_segs)-1.)
        
        def aligned_start_from_cigar(seg):
            """
            Determines the internal start position of this segment,
            relative to the full read sequence.
            """
            start = 0
            for op,count in seg.cigar:
                if op in [4,5]: # soft or hard clip
                    start += count
                elif op == 0: # match: aligned part begins
                    break
            return start

        # TODO: perhaps this can be optimized by making lists in the same order and (i)zipping
        internal_starts = dict([( seg,aligned_start_from_cigar(seg) ) for seg in self.proper_segs])
        internal_ends   = dict([( seg,internal_starts[seg] + len(seg.query) ) for seg in self.proper_segs])
        
        # sort by internal start position
        seg_by_seq = sorted(self.proper_segs,key = lambda seg: internal_starts[seg])
        
        # iterate over consecutive pairs of segments to find the 
        # splices between them
        # [A,B,C,D] -> (A,B),(B,C),(C,D)
        for seg_a,seg_b in zip(seg_by_seq, seg_by_seq[1:]):
            if options.debug:
                print "A",seg_a
                print "B",seg_b

            start_a = internal_starts[seg_a]
            end_a   = internal_ends[seg_a]
            len_a =  end_a - start_a
            
            start_b = internal_starts[seg_b]
            end_b   = internal_ends[seg_b]
            len_b = end_b - start_b
            
            if len_a < options.asize or len_b < options.asize:
                # segment length is below required anchor length
                N['seg_too_short_skip'] += 1
                continue

            r_start = min(start_a,start_b)
            r_end = max(end_a,end_b)
                
            ## anchor pairs that make it up to here are interesting
            if bam_out:
                bam_out.write(seg_a)

            yield JunctionSpan(seg_a,seg_b,self.primary,r_start,r_end,weight)

        if bam_out:
            bam_out.write(seg_b)

def collected_bwa_mem_segments(sam_input):
    """
    Generator that loops over the SAM alignments. It groups alignments 
    belonging to the same original read. 
    """
    from time import time
    t0 = time()
    t_last = t0

    current_mate = None
    other_mate = None
    
    for line_num,align in enumerate(sam_input):
        if align.is_unmapped:
            N['unmapped_reads'] += 1
            continue

        if not current_mate:
            # This only happens for the very first alignment
            current_mate = MateSegments(align)
            continue
            
        if current_mate.is_segment(align):
            #assert align.is_supplementary == True
            current_mate.add_segment(align)
            
        elif current_mate.is_other_mate(align):
            #assert align.is_supplementary == False
            # we are switching to the other mate now:
            other_mate,current_mate = current_mate,MateSegments(align)
        else:
            # we have a completely new read! 
            # Yield what we have so far:
            yield line_num,other_mate,current_mate
            # and start fresh
            #assert align.is_supplementary == False
            other_mate = None
            current_mate = MateSegments(align)

        if options.throughput:
            if line_num and not (line_num % options.chunksize):
                t1 = time()
                M_reads = line_num / 1E6
                mins = (t1 - t0)/60.
                krps = float(options.chunksize)/(t1-t_last)/1000.
                sys.stderr.write("\rprocessed {M_reads:.2f}M reads in {mins:.1f} minutes ({krps:.2f}k reads/second)       \r".format(**locals()))
                t_last = t1

    yield sam_line, other_mate, current_mate

    if options.throughput:
        sys.stderr.write('\n')


def record_hits(frag_name,frag_circ_splices, frag_linear_splices, frag_unspliced):
    """
    This function processes the observed splicing events from a single 
    sequenced fragment (mate pair, or single end read). Long reads, or
    mate pairs can span more than one junction, thus revealing information
    on the internal structure of a circRNA. 
    record_hits() extracts and reports this information.
    """
    global circ_splices, linear_splices
    
    circ_junctions = set()
    lin_junctions = set()
    unspliced = set()
    
    #print "record_hits()"
    for splice,junc_span,mate in frag_circ_splices:
        circ_splices[splice.coord].add(splice)
        circ_junctions.add(splice.coord)
    
    for splice,junc_span,mate in frag_linear_splices:
        linear_splices[splice.coord].add(splice)
        lin_junctions.add(splice.coord)
    
    for primary in frag_unspliced:
        if primary.is_reverse:
            sense = '-'
        else:
            sense = '+'
            
        coord = fast_chrom_lookup(primary),primary.pos,primary.aend,sense
        #print "unspliced",primary.qname,primary.is_read1,primary.is_read2
        #print coord
        unspliced.add(coord)

    if not circ_junctions:
        return

    if len(circ_junctions) > 1:
        for coord in circ_junctions:
            circ_splices[coord].add_flag('WARN_MULTI_BACKSPLICE',frag_name)
        return

    # we are sure now, that there is exactly one circ junction in the fragment.
    circ_coord = circ_junctions.pop()
    circ_chrom,circ_start,circ_end,circ_sense = circ_coord
    circ_hit = circ_splices[circ_coord]

    if options.debug:
        print "record_hits(",circ_coord,") linear: ",sorted(lin_junctions)," unspliced",sorted(unspliced)
    # test if (existing) unspliced reads fall within the circ
    # allow for [anchor-length] nucleotides to lie outside, as these may simply
    # have failed to be soft-clipped/spliced.
    for chrom,start,end,sense in unspliced:
        if circ_chrom == chrom:
            if start + options.asize <= circ_start or end - options.asize >= circ_end:
                circ_hit.add_flag('WARN_OUTSIDE_SEG',frag_name)
            else:
                circ_hit.add_flag('SUPPORT_INSIDE_SEG',frag_name)
        else:
            circ_hit.add_flag('WARN_OTHER_CHROM_SEG',frag_name)
            
    for chrom,start,end,sense in lin_junctions:
        if start <= circ_start or end >= circ_end:
            circ_hit.add_flag('WARN_OUTSIDE_SPLICE_JUNCTION',frag_name)
        else:
            circ_hit.add_flag('SUPPORT_INSIDE_SPLICE_JUNCTION',frag_name)

    if unspliced or lin_junctions:
        # we have a multi-splicing event!
        multi_event_storage.add(circ_coord,lin_junctions, unspliced)


def main():    
    # main loop
    try:
        sam_line = 0
        last_first_mate = ""

        frag_circ_splices = []
        frag_linear_splices = []
        frag_unspliced = []
        
        def process_mate(mate):
            if len(mate.segments) < 2:
                N['unspliced_reads'] += 1
                frag_unspliced.append(mate.primary)
                return
            
            # keep track of which parts of the read are already aligned 
            L = len(mate.full_seq)
            min_s = L
            max_e = 0
            
            #for A,B,r_start,r_end,w in mate.adjacent_segment_pairs():
            for junc_span in mate.adjacent_segment_pairs():
                if options.noop:
                    continue
                
                if junc_span.is_backsplice:
                    if options.debug:
                        print "potential circ"

                    prefix = 'circ'
                    frag_list = frag_circ_splices
                else:
                    # the anchors align sequentially -> linear/normal splice junction?
                    if options.nolinear:
                        N['linear_splice_skipped'] += 1
                        continue

                    prefix = 'linear'
                    frag_list = frag_linear_splices


                # TODO: if --ignore-linear is specified, only do this 
                # if at least one backsplice occurred in at least one mate instead of skipping all linear
                splices = junc_span.find_breakpoints()
                if not splices:
                    N['{0}_no_bp_weighted'.format(prefix)] += junc_span.weight
                    N['{0}_no_bp'.format(prefix)] += 1
                    continue


                N['{0}_spliced_weighted'.format(prefix)] += junc_span.weight
                N['{0}_spliced'.format(prefix)] += 1
                
                min_s = min(min_s,junc_span.q_start)
                max_e = max(max_e,junc_span.q_end)

                if not options.allhits:
                    splices = [splices[0],]

                for splice in splices:
                    frag_list.append( (splice,junc_span,mate) )


            if (max_e < L - options.asize) or (min_s > options.asize):
                # we are still missing a part of the read larger than options.asize
                # that could not be accounted for by the proper_segs
                # treat possible other fragments as if they were an unspliced mate
                frag_unspliced.extend(mate.other_chrom_segs)
                frag_unspliced.extend(mate.other_strand_segs)

            
        for sam_line,mate1,mate2 in collected_bwa_mem_segments(sam_input):
            if options.debug:
                print "one iteration",mate1,mate2

            if mate1: process_mate(mate1)
            if mate2: process_mate(mate2)

            if frag_circ_splices or frag_linear_splices:
                record_hits(last_first_mate,frag_circ_splices,frag_linear_splices,frag_unspliced)
                    
                frag_circ_splices = []
                frag_linear_splices = []
                frag_unspliced = []
            
        if frag_circ_splices or frag_linear_splices:
            record_hits(last_first_mate,frag_circ_splices,frag_linear_splices,frag_unspliced)
                
    except KeyboardInterrupt:
        fastq_line = sam_line * 4

        logging.warning("KeyboardInterrupt by user while processing input starting at SAM line {sam_line}, FASTQ line {fastq_line}".format(**locals()))
    except:        
        fastq_line = sam_line * 4

        logging.error("Unhandled exception raised while processing input starting at SAM line {sam_line}, FASTQ line {fastq_line}".format(**locals()))
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

def output(cand,prefix):
    sites_file.write("#" + "\t".join(['chrom','start','end','name','counts','strand','n_spanned','n_uniq','uniq_bridges','best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap','breakpoints','signal','strandmatch','category']) + "\n")
    n = 1
    for c,hit in cand.items():
        #print c,hit
        chrom,start,end,sense = c
        n_spanned,counts,n_uniq,best_qual_A,best_qual_B,uniq_bridges,tissues,tiss_counts,min_edit,min_anchor_ov,n_hits,signal,strandmatch = hit.scores(chrom,start,end,sense)

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

        name = "%s_%s_%06d" % (options.name,prefix,n)
        n += 1
        
        ## Store the reads associated with this hit. 
        ## TODO: We could also do this right away, when registering the read...
        for r_seq,ori_name in zip(hit.reads,hit.readnames):
            
            flags = hit.read_flags.get(ori_name)
            if flags:
                flag_str = ",".join(sorted(flags))
            else:
                flag_str = "BACKSPLICE"
            reads_file.write(">%s %s %s\n%s\n" % (ori_name,name,r_seq,flag_str))
      
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

        if end-start < options.short_threshold:
            categories.append("SHORT")
        elif end-start > options.huge_threshold:
            categories.append("HUGE")

        flags,flag_counts = hit.get_flags_counts()
        
        bed = [
            chrom,start-1,end,name,counts,sense, # BED6 compatible
            n_spanned,n_uniq,uniq_bridges,best_qual_A,best_qual_B, # mapping reliability
            ",".join(tissues),",".join(tiss_counts), # sample association
            min_edit,min_anchor_ov,n_hits,signal,strandmatch,",".join(sorted(categories) ), # splice site detection reliability
            ",".join(flags),",".join([str(c) for c in flag_counts]), # internal structure support from fragment evidence (BWA MEM)
        ]
        sites_file.write("\t".join([str(b) for b in bed]) + "\n")


if options.profile:
    import cProfile
    cProfile.run('main()')
else:
    main()

stats_file.write(str(N)+"\n")
output(circ_splices,"circ")
output(linear_splices,"norm")
multi_file.write(multi_event_storage.output())