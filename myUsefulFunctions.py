#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 15:36:01 2017

@author: goel
"""

import argparse
import os
import sys


class snvdata:
    """
    Class to store pileup data. Reads the first six mandatory columns and stores them.
    Each line is an object.
    For reading extra column, use the setattr function
    """

    def __init__(self, ls):
        self.chr = ls[0]
        self.pos = int(ls[1])
        self.ref = ls[2]
        self.rc  = int(ls[3])
        self.indelcount, self.bases = [0, []] if self.rc == 0 else self._getbases(ls[4])
        if len(self.bases) != self.rc:
            raise Exception('Number of identified bases if not equals to read count for {}:{}. ReadCount: {}, BaseCount: {}'.format(self.chr, self.pos, self.rc, len(self.bases)))
        self.BQ = [ord(c) - 33 for c in ls[5]]

    def _getbases(self, l):
        from collections import deque
        indelcount = 0
        bases = deque()
        skip = 0
        indel = False
        for c in l:
            if skip > 0 and indel == False:
                skip -= 1
                continue
            if indel == True:
                if c in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                    skip = (skip*10) + int(c)
                    continue
                # skip = int(c)
                else:
                    indel = False
                    skip -= 1
                    continue
            if c == '*':
                # self.indelcount += 1
                # indelcount += 1
                bases.append(c)
                continue
            if c == '$': continue
            if c in ['<', '>']: sys.exit('spliced alignment found')
            if c == '^':
                skip = 1
                continue
            if c in ['+', '-']:
                indel = True
                indelcount += 1
                continue
            bases.append(c)
        return [indelcount, list(bases)]

    def forwardcnt(self):
        try:
            return self.fcnt
        except AttributeError as e:
            self.fcnt = len([1 for i in self.bases if i in ['.', 'A', 'C', 'G', 'T']])
            self.rcnt = self.rc - self.fcnt
            return self.fcnt

    def reversecnt(self):
        try:
            return self.rcnt
        except AttributeError as e:
            self.rcnt = len([1 for i in self.bases if i in [',', 'a', 'c', 'g', 't']])
            self.fcnt = self.rc - self.rcnt
            return self.rcnt

    def basecnt(self, base):
        return len([1 for i in self.bases if i == base])

    def getBQ(self, base):
        return [self.BQ[i] for i in range(len(self.bases)) if self.bases[i] == base]

    def getindelreads(self, bases, reads):
        from collections import deque
        self.indelreads = deque()
        reads = reads.split(",")
        for i in range(len(reads)):
            if bases[0] == '^':
                bases == bases[3:]
                continue
            if bases[0] == '$':
                bases = bases[1:]
            if len(bases) == 1:
                continue
            if bases[1] not in ['+', '-']:
                bases = bases[1:]
            else:
                self.indelreads.append(reads[i])
                bases = bases[2:]
                skip = 0
                while bases[0] in {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}:
                    skip = (skip*10) + int(bases[0])
                    bases = bases[1:]
                bases = bases[skip:]
# END

def cgtpl(cg):
    """
    Takes a cigar string as input and returns a cigar tuple
    """
    for i in "MIDNSHPX=":
        cg = cg.replace(i, ';'+i+',')
    return [i.split(';') for i in cg.split(',')[:-1]]
#end

def cggenlen(cg, gen):
    """
    Takes cigar as input, and return the number of bases covered by it the reference
    or query genome.
    Cigar strings are first converted to cigar tuples.
    """
    if type(cg) == str:
        cg = cgtpl(cg)
    if gen not in ['r', 'q']:
        raise ValueError('gen need to "r" or "q" for reference or query')
        return
    s = set(['M', 'D', 'N', '=', 'X']) if gen == 'r' else set(['M', 'I', 'S', '=', 'X'])
    l = sum([int(i[0]) for i in cg if i[1] in s])
    return l
#end

def mergepdf(fins, fout):
    """
    Merge given PDF files in fins and save the combined file in fout.
    """
    from PyPDF2 import PdfFileMerger, PdfFileReader
    # Call the PdfFileMerger
    mergedObject = PdfFileMerger()
    for fin in fins: mergedObject.append(PdfFileReader(fin, 'rb'))
    # Write all the files into a file which is named as shown below
    mergedObject.write(fout)
    return
#end function

def pminf(array):
    x = 1
    pmin_list = []
    N = len(array)
    for index in range(N):
        if array[index] < x:
            pmin_list.insert(index, array[index])
        else:
            pmin_list.insert(index, x)
    return pmin_list
#end function
 
def cumminf(array):
    cummin = []
    cumulative_min = array[0]
    for p in array:
        if p < cumulative_min:
            cumulative_min = p
        cummin.append(cumulative_min)
    return cummin
#end
 
def cummaxf(array):
    cummax = []
    cumulative_max = array[0]
    for e in array:
        if e > cumulative_max:
            cumulative_max = e
        cummax.append(cumulative_max)
    return cummax
#end
 
def order(*args):
    if len(args) > 1:
        if args[1].lower() == 'false':# if ($string1 eq $string2) {
            return sorted(range(len(args[0])), key=lambda k: args[0][k])
        elif list(args[1].lower()) == list('true'):
            return sorted(range(len(args[0])), key=lambda k: args[0][k], reverse=True)
        else:
            print("{} isn't a recognized parameter".format(args[1]))
            sys.exit()
    elif len(args) == 1:
        return sorted(range(len(args[0])), key=lambda k: args[0][k])
#end
 
def p_adjust(*args):
    """
    copied from 
    https://rosettacode.org/wiki/P-value_correction#Python
    
    """
    method = "bh"
    pvalues = args[0]
    if len(args) > 1:
        methods = {"bh", "fdr", "by", "holm", "hommel", "bonferroni", "hochberg"}
        metharg = args[1].lower()
        if metharg in methods:
            method = metharg
    lp = len(pvalues)
    n = lp
    qvalues = []
 
    if method == 'hochberg':#already all lower case
        o = order(pvalues, 'TRUE')
        cummin_input = []
        for index in range(n):
            cummin_input.insert(index, (index+1)*pvalues[o[index]])
        cummin = cumminf(cummin_input)
        pmin = pminf(cummin)
        ro = order(o)
        qvalues = [pmin[i] for i in ro]
    elif method == 'bh':
        o = order(pvalues, 'TRUE')
        cummin_input = []
        for index in range(n):
            cummin_input.insert(index, (n/(n-index))* pvalues[o[index]])
        ro = order(o)
        cummin = cumminf(cummin_input)
        pmin = pminf(cummin)
        qvalues = [pmin[i] for i in ro]
    elif method == 'by':
        q = 0.0
        o = order(pvalues, 'TRUE')
        ro = order(o)
        for index in range(1, n+1):
            q += 1.0 / index;
        cummin_input = []
        for index in range(n):
            cummin_input.insert(index, q * (n/(n-index)) * pvalues[o[index]])
        cummin = cumminf(cummin_input)
        pmin = pminf(cummin)
        qvalues = [pmin[i] for i in ro]
    elif method == 'bonferroni':
        for index in range(n):
            q = pvalues[index] * n
            if (0 <= q) and (q < 1):
                qvalues.insert(index, q)
            elif q >= 1:
                qvalues.insert(index, 1)
            else:
                print('{} won\'t give a Bonferroni adjusted p'.format(q))
                sys.exit()
    elif method == 'holm':
        o = order(pvalues)
        cummax_input = []
        for index in range(n):
            cummax_input.insert(index, (n - index) * pvalues[o[index]])
        ro = order(o)
        cummax = cummaxf(cummax_input)
        pmin = pminf(cummax)
        qvalues = [pmin[i] for i in ro]
    elif method == 'hommel':
        i = range(1,n+1)
        o = order(pvalues)
        p = [pvalues[index] for index in o]
        ro = order(o)
        pa = []
        q = []
        smin = n*p[0]
        for index in range(n):
            temp = n*p[index] / (index + 1)
            if temp < smin:
                smin = temp
        for index in range(n):
            pa.insert(index, smin)
            q.insert(index, smin)
        for j in range(n-1,1,-1):
            ij = range(1,n-j+2)
            for x in range(len(ij)):
                ij[x] -= 1
            I2_LENGTH = j - 1
            i2 = []
            for index in range(I2_LENGTH+1):
                i2.insert(index, n - j + 2 + index - 1)
            q1 = j * p[i2[0]] / 2.0
            for index in range(1,I2_LENGTH):
                TEMP_Q1 = j * p[i2[index]] / (2.0 + index)
                if TEMP_Q1 < q1:
                    q1 = TEMP_Q1
            for index in range(n - j + 1):
                q[ij[index]] = min(j * p[ij[index]], q1)
            for index in range(I2_LENGTH):
                q[i2[index]] = q[n-j]
            for index in range(n):
                if pa[index] < q[index]:
                    pa[index] = q[index]
            qvalues = [pa[index] for index in ro]
    else:
        print("method {} isn't defined.".format(method))
        sys.exit()
    return qvalues
#end
 
def unlist(nestedList):
    import numpy as np
    """Take a nested-list as input and return a 1d list of all elements in it"""
    outList = []
    for i in nestedList:
        if type(i) in (list, np.ndarray):
            outList.extend(unlist(i))
        else:
            outList.append(i)
    return(outList)


def getValues(l, index):
    """from list l get the values at indices specified by index"""
    return [l[i] for i in index]


def getColors(colorPalette, numOfCol):
	return([colorPalette(i/numOfCol) for i in range(numOfCol)])


def plotdensity(data):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde
    density = gaussian_kde(data)
    xs = np.linspace(min(data), max(data), 1000)
    density.covariance_factor = lambda: .2
    density._compute_covariance()
    plt.plot(xs, density(xs))
    plt.show()
#end

def sublist(lst1, lst2):
    import operator as op
    return(list(map(op.sub,lst1, lst2)))
#end

def intersect(*lists):
    import numpy as np
    return reduce(np.intersect1d,list(lists))


def readblast(f):
    """
    Read tabular format blast output
    """
    import pandas as pd
    return pd.read_table(f, comment="#")


def readfasta(f):
    from gzip import open as gzopen
    from gzip import BadGzipFile
    from collections import deque
    import sys
    out = {}
    chrid = ''
    chrseq = deque()

    # Test if the file is Gzipped or not
    with gzopen(f, 'rb') as fin:
        try:
            fin.read(1)
            isgzip = True
        except BadGzipFile:
            isgzip = False

    try:
        if isgzip:
            with gzopen(f, 'rb') as fin:
                for line in fin:
                    if b'>' in line:
                        if chrid != '':
                            out[chrid] = ''.join(chrseq)
                            chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                            chrseq = deque()
                        else:
                            chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                        if chrid in out.keys():
                            sys.exit(" Duplicate chromosome IDs are not accepted. Chromosome ID {} is duplicated. Provided chromosome with unique IDs".format(chrid))
                    else:
                        chrseq.append(line.strip().decode())
        else:
            with open(f, 'r') as fin:
                for line in fin:
                    if '>' in line:
                        if chrid != '':
                            out[chrid] = ''.join(chrseq)
                            chrid = line.strip().split('>')[1].split(' ')[0]
                            chrseq = deque()
                        else:
                            chrid = line.strip().split('>')[1].split(' ')[0]
                        if chrid in out.keys():
                            sys.exit(" Duplicate chromosome IDs are not accepted. Chromosome ID {} is duplicated. Provided chromosome with unique IDs".format(chrid))
                    else:
                        chrseq.append(line.strip())
    except Exception as e:
        raise Exception(e)

    if chrid != '':
        out[chrid] = ''.join(chrseq)
    # TODO: add check for the validation of input fasta files
    return out


def writefasta(fa, f):
    """
    :param fa: dictionary. Keys are chromosome ids. Values are sequence.
    :param f: Output file
    :return:
    """
    # TODO: Add capability to write fa.gzip files as well
    with open(f, 'w') as fo:
        for k, v in fa.items():
            fo.write('>'+k+'\n')
            for i in range(0, len(v), 60):
                fo.write(v[i:(i+60)]+'\n')
#END


def density_scatter(x, y, ax=None, fig=None, sort=True, bins=20, **kwargs):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import Normalize
    from scipy.interpolate import interpn
    # import seaborn as sns         # Can be used if regression fit is required, but is not working very well with the scatter kwargs
    """
    Scatter plot colored by 2d histogram
    Usage example:
        x = np.random.normal(size=100000)
        y = x * 3 + np.random.normal(size=100000)
        density_scatter( x, y, bins = [30,30] )
    """
    if ax is None:
        fig, ax = plt.subplots()
    data, x_e, y_e = np.histogram2d(x, y, bins=bins, density=True)
    z = interpn((0.5*(x_e[1:] + x_e[:-1]), 0.5*(y_e[1:]+y_e[:-1])), data, np.vstack([x, y]).T, method="splinef2d", bounds_error=False)
    # To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0
    # Sort the points by density, so that the densest points are plotted last
    if sort:
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    # ax = sns.regplot(x=x, y=y, scatter_kws={"s": 0.5}, ax=ax)
    ax.scatter(x, y, c=z, **kwargs)
    norm = Normalize(vmin=np.min(z), vmax=np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm), ax=ax)
    cbar.ax.set_ylabel('Density')
    return ax
# END


# def extractSeq(args):
#     # TODO: Fix this function to work with new parameter style
#     # TODO: Use in-built fasta parser/writer and remove dependency on Bio.SeqIO
#     import pandas as pd
#     from Bio.SeqRecord import SeqRecord
#     from Bio.SeqIO import parse, write
#     filePath = args.fasta.name
#     if args.fin == None:
#        seqID = args.loc[0]
#        querySeq = [fasta for fasta in parse(filePath, 'fasta') if fasta.id == seqID][0]
#
#        start = int(args.s) if args.s is not None else 0
#        end = int(args.e) if args.e is not None else 0
#        end = end if end <= len(querySeq.seq) else len(querySeq.seq)
#
#        querySeq.seq = querySeq.seq[int(start):(int(end)+1)]
#        if args.o != None:
#            write(querySeq, args.o, "fasta")
#        else:
#            print("> "+querySeq.id)
#            print(querySeq.seq)
#     else:
#         fin = pd.read_table(args.fin.name, header=None, delim_whitespace=True)
#         fin.columns = ["chr", "start", "end"]
#         fin[['start', 'end']] = fin[['start', 'end']].astype('int')
#         fin.loc[fin['start'] < 0, 'start'] = 0
#         fin.sort_values(["chr", "start", "end"], inplace=True)
#         outF = deque()
#         for fasta in parse(filePath, 'fasta'):
#             if fasta.id in fin.chr.values:
#                 chrData = fin.loc[fin.chr == fasta.id].copy()
#                 chrData.loc[chrData['end'] > len(fasta.seq), 'end'] = len(fasta.seq)
#                 for row in chrData.itertuples(index=False):
#                     outF.append(SeqRecord(seq=fasta.seq[row.start:row.end], id="_".join(map(str, row)), description=""))
#         if args.o != None:
#             write(outF, args.o.name, "fasta")
#         else:
#            for i in outF:
#                 print("> "+i.id)
#                 print(i.seq)
#END

def extractSeq(args):
    # TODO: Fix this function to work with new parameter style
    print(args, type(args.o.name))
    import pandas as pd
    import warnings
    if args.fin is not None:
        if args.chr is not None or args.start is not None or args.end is not None:
            warnings.warn("Using --fin. Ignoring --chr, -s, -e")
    f = args.fasta.name
    if args.fin is None:
        seqid = args.chr
        q = {c: s for c, s in readfasta(f).items() if c == seqid}
        if len(q) > 1:
            sys.exit("Found multiple chromosomes with same ID. Exiting.")
        start = int(args.start) if args.start is not None else 0
        end = int(args.end) if args.end is not None else len(q[seqid])
        end = end if end <= len(q[seqid]) else len(q[seqid])
        # Output the selected sequence
        if args.o is not None:
            writefasta({seqid: q[seqid][start:(end+1)]}, args.o.name)
        else:
            print("> {}\n{}".format(seqid, q[seqid][start:(end+1)]))
    else:
        fin = pd.read_table(args.fin.name, header=None, delim_whitespace=True)[[0, 1, 2]]
        fin.columns = ["chr", "start", "end"]
        fin[['start', 'end']] = fin[['start', 'end']].astype('int')
        fin.loc[fin['start'] < 0, 'start'] = 0
        fin.sort_values(["chr", "start", "end"], inplace=True)
        out = dict()
        chroms = set(fin.chr.values)
        for c, s in readfasta(f).items():
            print(c, len(s))
            if c in chroms:
                cdf = fin.loc[fin.chr == c].copy()
                cdf.loc[cdf['end'] > len(s), 'end'] = len(s)
                for row in cdf.itertuples(index=False):
                    out['{}_{}_{}'.format(c, row[1], row[2])] = s[row.start:row.end]
        # Output the selected sequence
        if args.o is not None:
            warnings.warn("writing output fasta")
            writefasta(out, args.o.name)
        else:
            for c, s in out.items():
                print("> {}\n{}".format(c, s))
# END

def revcomp(seq):
    assert type(seq) == str
    old = 'ACGTRYKMBDHVacgtrykmbdhv'
    rev = 'TGCAYRMKVHDBtgcayrmkvhdb'
    tab = str.maketrans(old, rev)
    return seq.translate(tab)[::-1]


def subnuc(args):
    from Bio.SeqIO import parse, write
    from Bio.Seq import Seq
    fasta = args.fasta.name
    querySeq = [fasta for fasta in parse(fasta,'fasta')]
    for i in range(len(querySeq)):
        querySeq[i].seq = Seq(str(querySeq[i].seq).replace(args.q, args.t))
#    print(querySeq)
    
    if args.o == None:
        fout = fasta+"_edited"
    else:
        fout = args.o
    with open(fout,"w") as f:
        spacer = ""
        for seq in querySeq:
            f.write(spacer)
            f.write(">"+seq.id+" "+seq.description+"\n")
            f.write('\n'.join(str(seq.seq)[i:i+60] for i in range(0, len(seq.seq), 60)))
            if spacer == "":
                spacer = "\n"


def fileRemove(fName):
    try:
        os.remove(fName)
    except OSError as e:
        if e.errno != 2:    ## 2 is the error number when no such file or directory is present https://docs.python.org/2/library/errno.html
            raise
            

def mergeRanges(ranges):
    """
    Take a 2D numpy array, with each row as a range and return merged ranges
    i.e. ranges which are overlapping would be combined.
    :param ranges:
    :return:
    """
    from collections import deque
    import numpy as np
    if len(ranges) < 2:
        return ranges
    for i in ranges:
        if i[0] > i[1]:
            i[1], i[0] = i[0], i[1]
    ranges = ranges[ranges[:, 0].argsort()]
    min_value = ranges[0, 0]
    max_value = ranges[0, 1]
    out_range = deque()
    for i in ranges[1:]:
        if i[0] > max_value:
            out_range.append([min_value, max_value])
            min_value = i[0]
            max_value = i[1]
        elif i[1] > max_value:
            max_value = i[1]
    out_range.append([min_value, max_value])
    return np.array(out_range)
#END

def subranges(r1, r2):
    """
    Subtract range2 (r2) from range1 (r1) and return non-ooverllaping range.
    Both ranges are considered closed.
    """
    from collections import deque
    import numpy as np
    # Get sorted and merged ranges.
    r1 = mergeRanges(np.array(r1))
    r2 = mergeRanges(np.array(r2))
    i = 0
    j = 0
    lr1 = len(r1)
    lr2 = len(r2)
    cr1 = r1[i]
    outr = deque()
    f = False
    while True:
        if i == lr1: break
        if j == lr2:
            outr.append(cr1)
            for k in range(i+1, lr1):
                outr.append(r1[k])
            f = True
            break
        ## Check no overlap conditions
        if cr1[1] < r2[j][0]:
            outr.append(cr1)
            i += 1
            if i == lr1: break
            cr1 = r1[i]
        elif cr1[0] > r2[j][1]:
            j+=1
        ## Cases where r1-element ends before r2-element
        elif cr1[1] <= r2[j][1]:
            if cr1[0] < r2[j][0]:
                outr.append([cr1[0], r2[j][0]-1])
            i += 1
            if i == lr1: break
            cr1 = r1[i]
        ## Cases where r1-element ends after r2-element
        else:
            if cr1[0] < r2[j][0]:
                outr.append([cr1[0], r2[j][0]-1])
            cr1 = [r2[j][1]+1, cr1[1]]
            j += 1
    return(outr)
#END

def total_size(o, handlers={}, verbose=False):
    try:
        from reprlib import repr
    except ImportError as e :
        sys.exit('missing library' + str(e))
    """ Returns the approximate memory footprint an object and all of its contents.

    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, deque, dict, set and frozenset.
    To search other containers, add handlers to iterate over their contents:

        handlers = {SomeContainerClass: iter,
                    OtherContainerClass: OtherContainerClass.get_elements}

    """
    dict_handler = lambda d: chain.from_iterable(d.items())
    all_handlers = {tuple: iter,
                    list: iter,
                    deque: iter,
                    dict: dict_handler,
                    set: iter,
                    frozenset: iter,
                   }
    all_handlers.update(handlers)     # user handlers take precedence
    seen = set()                      # track which object id's have already been seen
    default_size = getsizeof(0)       # estimate sizeof object without __sizeof__

    def sizeof(o):
        if id(o) in seen:       # do not double count the same object
            return 0
        seen.add(id(o))
        s = getsizeof(o, default_size)

        if verbose:
            print(s, type(o), repr(o), file=stderr)

        for typ, handler in all_handlers.items():
            if isinstance(o, typ):
                s += sum(map(sizeof, handler(o)))
                break
        return s

    return sizeof(o)


def getScaf(args):
    import numpy as np
    from Bio.Alphabet import generic_dna
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqIO import parse, write
    fin = args.fasta.name
    n = args.n
    fout = args.o
    rev = 2 if args.r else 1
    gen = {fasta.id:fasta for fasta in parse(fin,'fasta')}
    chromLen = {chrom.id:len(chrom.seq) for chrom in gen.values()}
    
    chrID = [np.random.choice(list(chromLen.keys()), 1)[0] for i in range(n)]
    cutoff = [np.random.randint(10000,chromLen[i] - 10000, 1)[0] for i in chrID]
    coords = sorted(list(zip(chrID, cutoff)))
    
    cid = ""
    pos = -1
    scafGenome = []
    for i in range(len(coords)):
        coord = coords[i]
        if coord[0] != cid:
            if cid in gen:
                if np.random.randint(0, rev, 1):
                    scafGenome.append(SeqRecord(seq=gen[cid][pos:chromLen[cid]].seq.reverse_complement(), id = cid+"_"+str(chromLen[cid])+"_r", description=""))
                else:
                    scafGenome.append(SeqRecord(seq=gen[cid][pos:chromLen[cid]].seq, id = cid+"_"+str(chromLen[cid]), description=""))  
            cid = coord[0]
            pos = 0
        if np.random.randint(0, rev, 1):
            scafGenome.append(SeqRecord(seq=gen[cid][pos:coord[1]].seq.reverse_complement(), id = cid+"_"+str(coord[1])+"_r", description=""))
        else:
            scafGenome.append(SeqRecord(seq=gen[cid][pos:coord[1]].seq, id = cid+"_"+str(coord[1]), description=""))
        pos = coord[1]
    
    if np.random.randint(0, rev, 1):
        scafGenome.append(SeqRecord(seq=gen[cid][pos:chromLen[cid]].seq.reverse_complement(), id = cid+"_"+str(chromLen[cid])+"_r", description=""))
    else:
        scafGenome.append(SeqRecord(seq=gen[cid][pos:chromLen[cid]].seq, id = cid+"_"+str(chromLen[cid]), description=""))
    
    for key in gen.keys():
        if key not in chrID:
            scafGenome.append(SeqRecord(seq=gen[key].seq, id = key, description=""))
    write(scafGenome,fout,"fasta")

def seqsize(args):
    fin = args.fasta.name
    out = [(chrom, len(seq)) for chrom, seq in readfasta(fin).items()]
    #TODO: Add output file
    for i in out:
        print(i[0], i[1], sep="\t")
    if args.t:
        print("Genome_length", sum([i[1] for i in out]), sep="\t")
#END


def filsize(args):
    """
    Remove molecules
    which are smaller than the specified threshold
    """
    from Bio.SeqIO import parse, write
    fin = args.fasta.name
    size= args.size
    gen = [fasta for fasta in parse(fin,'fasta') if len(fasta.seq) > size]
    write(gen, fin.split(".fna")[0].split("/")[-1]+".filtered.fna", "fasta")
#END

def basrat(args):
    from collections import Counter
    fin = args.fasta.name
    count = {chrom: dict(Counter(seq)) for chrom, seq in readfasta(fin).items()}
    outD = {}
    for v in count.values():
        for k, cnt in v.items():
            if k in outD:
                outD[k] += cnt
            else:
                outD[k] = cnt
    for k, v in outD.items():
        print(k, v, sep="\t")
#END

def getchr(args):
    fin = args.fasta.name
    if args.chrs is None and args.F is not None: chroms = [c.strip() for c in open(args.F.name, 'r').readlines()]
    elif args.chrs is not None and args.F is None: chroms = args.chrs
    else: raise Exception('InvalidValue: Incorrect value for chrs provided')
    out = args.o.name if args.o is not None else fin + ".filtered.fasta"
    with open(out, 'w') as fout:
        for chrom, seq in readfasta(fin).items():
            if (chrom in chroms) != args.v:
                fout.write(">" + chrom + "\n" + "\n".join([seq[i:i+60] for i in range(0, len(seq), 60)]) + '\n')
#END

def genome_ranges(args):
    import os
    import sys
    if args.n < 1:
        sys.exit('Range should be more than 0')
    if os.path.isfile(args.o):
        sys.exit(args.o + ' exists. Cannot overwrite it.')

    from Bio.SeqIO import parse

    start = ''
    with open(args.o, 'w') as fout:
        for fasta in parse(args.fasta, 'fasta'):
            s = len(fasta.seq)
            id = fasta.id
            r1 = list(range(1, s, args.n))
            r2 = list(range(args.n, s, args.n)) + [s]

            for i in range(len(r1)):
                fout.write(start + str(id) + "\t" + str(r1[i]) + "\t" + str(r2[i]))
                start = '\n'
#END

def get_homopoly(args):
    import os
    import sys
    if args.n < 2:
        sys.exit('Minimum allowed homopolymer length = 2')
    if os.path.isfile(args.o):
        sys.exit(args.o + ' exists. Cannot overwrite it.')
    import re
    from collections import OrderedDict
    isbed = 1 if args.b else 0
    start = ''
    with open(args.o, 'w') as fout:
        for chrom, s in readfasta(args.fasta.name).items():
            # id = fasta.id
            outdict = {}
            for b in ['A', 'C', 'G', 'T']:
                hp = b * args.n + '+'
                # hp = b * 4
                r = re.finditer(hp, s)
                for i in r:
                    outdict[i.start() - isbed] = [i.end(), b]
            for i in sorted(outdict.keys()):
                fout.write(start + str(chrom) + "\t" + str(i+1 - args.p) + "\t" + str(outdict[i][0] + args.p) + "\t" + outdict[i][1])
                start = '\n'
#END

def asstat(fin):
    from collections import deque
    fasta_len = deque()
    for chrom, seq in readfasta(fin).items():
        fasta_len.append(len(seq))
    fasta_len = sorted(fasta_len)

    fasta_len_cumsum = deque()
    last = 0
    for i in fasta_len:
        fasta_len_cumsum.append(last + i)
        last = last + i
    half = fasta_len_cumsum[-1]/2
    for i in range(len(fasta_len_cumsum)):
        if fasta_len_cumsum[i] > half:
            return [fasta_len_cumsum[-1], len(fasta_len), fasta_len[-1], fasta_len[0], fasta_len[i], len(fasta_len)-i]
#END

def getasstat(args):
    from collections import deque
    from multiprocessing import Pool
    from functools import partial
    out = args.o.name if args.o is not None else "genomes.n50"
    nc = args.n
    from collections import deque
    fins = deque()
    if args.F is None and args.G is None:
        sys.exit("Provide genome file path in -F or -G")
    elif args.F is not None and args.G is not None:
        sys.exit("Provide genome file path in -F or -G")
    elif args.F is not None:
        with open(args.F.name, 'r') as F:
            for line in F:
                fins.append(line.strip())
    elif args.G is not None:
        fins = args.G
    with Pool(processes=nc) as pool:
        n50_values = pool.map(partial(asstat), fins)
    with open(out, 'w') as fout:
        fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("assemble_id", "assembly_length", "number_of_contig", "longest_contig", "shortest_contig", "N50", "L50"))
        for i in range(len(fins)):
            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(fins[i], n50_values[i][0], n50_values[i][1], n50_values[i][2], n50_values[i][3], n50_values[i][4], n50_values[i][5]))
#END

# def asstat(fin, parse):
#     from numpy import cumsum
#     fasta_len = deque()
#     for fasta in parse(fin, 'fasta'):
#         fasta_len.append(len(str(fasta.seq)))
#     fasta_len = sorted(fasta_len)
#     fasta_len_cumsum = cumsum(fasta_len)
#     half = fasta_len_cumsum[-1]/2
#     for i in range(len(fasta_len_cumsum)):
#         if fasta_len_cumsum[i] > half:
#             return [fasta_len_cumsum[-1], len(fasta_len), fasta_len[-1], fasta_len[0], fasta_len[i], len(fasta_len)-i]
#
#
# def getasstat(args):
#     out = args.o.name if args.o is not None else "genomes.n50"
#     NC = args.n
#     from collections import deque
#     fins = deque()
#     if args.F is None and args.G is None:
#         sys.exit("1Provide genome file path in -F or -G")
#     elif args.F is not None and args.G is not None:
#         sys.exit("Provide genome file path in -F or -G")
#     elif args.F is not None:
#         with open(args.F.name, 'r') as F:
#             for line in F:
#                 fins.append(line.strip())
#     elif args.G is not None:
#         fins = args.G
#     from multiprocessing import Pool
#     from functools import partial
#     from Bio.SeqIO import parse
#     with Pool(processes=NC) as pool:
#         n50_values = pool.map(partial(asstat, parse=parse), fins)
#     with open(out, 'w') as fout:
#         fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("assemble_id", "assembly_length", "number_of_contig", "longest_contig", "shortest_contig", "N50", "L50"))
#         for i in range(len(fins)):
#             fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(fins[i], n50_values[i][0], n50_values[i][1], n50_values[i][2], n50_values[i][3], n50_values[i][4], n50_values[i][5]))


def plthist(args):
    import sys
    import warnings
    if args.f == 'STDIN': fin = sys.stdin
    else:
        try:
            fin = open(args.f, 'r')
        except FileNotFoundError:
            raise FileNotFoundError(' Cannot open file: {}'.format(args.f))
            sys.exit()
  
  # Read data
    data = {}
    for line in fin:
        line = line.strip().split()
        if len(line) > 2: raise warnings.warn('Input has more than 2 columns. Using the first two columns')
        try:
            data[line[1]] = float(line[0])
        except ValueError:
            raise ValueError('First column contains non numerical values')
            sys.exit()
  # Check if keys are numeric
    num = True
    for k in data.keys():
        try:
            a = float(k)
        except ValueError:
            num = False
            break
    # print(num)
    if num:
        if args.xlim is not None:
            MIN, MAX = args.xlim[0], args.xlim[1]
            pltdata = {float(k): v for k, v in data.items() if float(k) <= MAX and float(k) >= MIN}
        else :
            pltdata = {float(k): v for k, v in data.items()}
            MIN = min(list(pltdata.keys()))
            MAX = max(list(pltdata.keys()))
        pltdata2 = {}
        m = MIN
        n = (MAX - MIN)/args.n
        for i in range(args.n):
            pltdata2[(m, m+n)] = 0
            m +=n
        for k, v in pltdata.items():
            for k2 in pltdata2.keys():
                if k >= k2[0] and k <k2[1]:
                    pltdata2[k2] += v
        
        pltdata = {(k[0]+k[1])/2 : v for k, v in pltdata2.items()}
    else:
        pltdata = data

    sort_k = sorted(list(pltdata.keys()))
    # print(pltdata)
    from matplotlib import use as mpuse
    mpuse('agg')
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=[args.W, args.H])
    ax = fig.add_subplot()
    if num:
        ax.bar(sort_k, [int(pltdata[k]) for k in sort_k], width=n)
    else:
        ax.bar(sort_k, [int(pltdata[k]) for k in sort_k])
    ax.set_xlabel(' '.join(args.x))
    ax.set_ylabel(' '.join(args.y))
    if args.t is not None: ax.set_title(args.t)
    if args.xlog: ax.set_xscale('log')
    if args.ylog: ax.set_yscale('log')
    if args.xlim is not None: ax.set_xlim([args.xlim[0], args.xlim[1]])
    if args.ylim is not None: ax.set_ylim([args.ylim[0], args.ylim[1]])
    ax.minorticks_on()
    ax.xaxis.grid(True, which='both', linestyle='--')
    ax.yaxis.grid(True, which='both', linestyle='--')
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig(args.o.name)
    fin.close()


def gfftrans(args):
    print('WARNING: THIS FUNCTION MIGHT HAVE BUGS.')
    from collections import deque
    fasta = readfasta(args.fasta.name)
    trans = {}
    with open(args.gff.name, 'r') as f:
        for line in f:
            line = line.strip().split()
            if line[2] == 'gene':
                if not args.r:
                    trans[line[0] + '_' + line[3] + '_' + line[4]] = fasta[line[0]][int(line[3])-1:int(line[4])-1]
                else:
                    trans[line[0] + '_' + line[3] + '_' + line[4]] = revcomp(fasta[line[0]][int(line[3])-1:int(line[4])-1])
    writefasta(trans, args.o.name)


def vcfdp(args):
    with open(args.vcf.name, 'r') as fin:
        with open(args.o.name, 'w') as fout:
            for line in fin:
                if line[0] == '#': continue
                DP = None
                DP4 = None
                line = line.strip().split()
                for v in line[7].split(';'):
                    v = v.split('=')
                    if v[0] == 'DP': DP = v[1]
                    if v[0] == 'DP4':DP4 = v[1].split(',')
                if DP is None or DP4 is None:
                    raise ValueError("DP and DP4 values are not present for {}".format('\t'.join(line)))
                    sys.exit()
                fout.write('\t'.join([line[0], line[1], line[3], line[4], DP] + DP4) + '\n')
                    
                        
def gffsort(args):
    print('Required format of GFF file:\n')
    print('Chr\tMethod\tType\tStart\tEnd\tScore\tStrand\tPhase\tAttribute(ID=uniq_id;other_info)\n')
    from collections import defaultdict
    header = True
    geneids = defaultdict(dict)
    gffdata = defaultdict(dict)
    first = True
    with open(args.out.name, 'w') as fout:
        with open(args.gff.name, 'r') as fin:
            for line in fin:
                if line.strip() == '': continue
                if line[0] == '#':
                    if header: fout.write(line)
                    continue
                else:
                    header = False
                    line = line.strip().split()
                    if line[2] == 'gene':
                        if first:
                            chr = line[0]
                            id = line[8].split(';')[0].split('=')[1]
                            geneids[line[0]][id] = int(line[3])
                            seq = '\t'.join(line) + '\n'
                            first = False
                        else:
                            gffdata[chr][id] = seq
                            chr = line[0]
                            id = line[8].split(';')[0].split('=')[1]
                            geneids[line[0]][id] = int(line[3])
                            seq = '\t'.join(line) + '\n'
                    else:
                        seq = seq + '\t'.join(line) + '\n'
            gffdata[chr][id] = seq
        
        for chr in sorted(list(geneids.keys())):
            for gid in sorted(geneids[chr].keys(), key = lambda x: geneids[chr][x]):
                fout.write(gffdata[chr][gid])
        

def getcol(args):
    if args.s is None:
        args.s = '\t'
    elif len(args.s) > 1:
        sys.exit('Only single characters are accepted as column separators')
    if len(args.m) == 0 and len(args.c) == 0:
        sys.exit("No columns are selected. Use -m and -c to select columns")

    with open(args.file.name, 'r') as fin:
        with open(args.out.name, 'w') as fout:
            first = True
            for line in fin:
                line = line.strip().split(args.s)
                if first:
                    select = []
                    for m in args.m:
                        for i in range(len(line)):
                            if line[i] == m:
                                select.append(i)
                    for c in args.c:
                        for i in range(len(line)):
                            if c in line[i]:
                                select.append(i)
                    select = sorted(list(set(select)))
                    first = False
                    print('Select columns: {}'.format([line[i] for i in select]))
                fout.write(args.s.join([line[i] for i in select]) + '\n')


def bamcov(args):
    import sys
    import os
    from subprocess import Popen, PIPE
    import warnings

    formatwarning_orig = warnings.formatwarning
    warnings.formatwarning = lambda message, category, filename, lineno, line=None:    formatwarning_orig(message, category, filename, lineno, line='')

    BAM = args.bam.name
    OUT = args.out.name
    BED = args.r.name if args.r is not None else None
    if BED is not None:
        sys.exit("SELECTING REGIONS (-R) HAS NOT BEEN IMPLEMENTED YET. PLEASE DO NOT USE IT.")
    D = args.d
    if not os.path.isfile(BAM):
        sys.exit('BAM file is not found: {}'.format(BAM))

    # if not os.path.isfile(BED):
    #     sys.exit('BED file is not found: '.format(BED))

    # Get chromosome names and lengths
    p = Popen("samtools view -H {}".format(BAM).split(), stdout=PIPE, stderr=PIPE)
    o = p.communicate()
    if o[1] != b'':
        sys.exit("Error in using samtools view to get BAM file header:\n{}".format(o[1].decode()))

    h = o[0].decode().strip().split("\n")
    chrlen = {}
    for line in h:
        line = line.strip().split()
        if line[0] != '@SQ': continue
        if not any(['SN' in x for x in line]):
            sys.exit('SN tag not present in header: {}'.format(' '.join(line)))
        if not any(['LN' in x for x in line]):
            sys.exit('LN tag not present in header: {}'.format(' '.join(line)))
        c = [x for x in line if 'SN' == x[:2]][0].split(':')[1]
        l = [x for x in line if 'LN' == x[:2]][0].split(':')[1]
        chrlen[c] = int(l)


    # Get read-depths
    tname = "TMP_" + os.path.basename(BAM) + ".txt"
    warnings.warn("Creating {} for saving the samtools depth output. Ensure enough storage space is available.".format(tname), stacklevel=-1)
    with open(tname, 'w') as TMP:
        p = Popen("{} {}".format(D, BAM).split(), stdout=TMP, stderr=PIPE)
        o = p.communicate()
        if o[1] != b'':
            sys.exit("Error in using samtools depth:\n{}".format(o[1].decode()))
    chrrd = {c: 0 for c in chrlen.keys()}
    with open(tname, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            try:
                chrrd[line[0]] += int(line[2])
            except KeyError as e:
                sys.exit("Chromosome {} in samtools_depth_output not present in SAM file header. Ensure that all chromosomes and their lengths are present SAM file header.")

    with open(OUT, 'w') as fout:
        for k, v in chrrd.items():
            fout.write("{}\t{}\n".format(k, v/chrlen[k]))
        if args.g:
            fout.write("Genome\t{}\n".format(round(sum(chrrd.values())/sum(chrlen.values()))))
    os.remove(tname)
    # TODO: Add method to filter for BED regions
# END

def run_bam_readcount(tmp_com, outfile):
    from subprocess import Popen, PIPE
    with open(outfile, 'w') as TMP:
        p = Popen(tmp_com.split(), stdout=TMP, stderr=PIPE)
    o = p.communicate()
    # if o[1] != b'Minimum mapping quality is set to 40':
    #     sys.exit("Error in running bam-readcount:\n{}".format(o[1].decode()))


def pbamrc(args):
    # print(args)
    if args.bam is None or args.l is None or args.f is None:
        sys.exit("Require bam, position bed, and refernce fasta files are missing")

    if args.n < 1: sys.exit("Invalid number of cores")
    N = args.n

    from subprocess import Popen, PIPE
    from multiprocessing import Pool
    import pandas as pd
    import numpy as np
    from time import time
    import os

    INDEL = args.I
    bed = pd.read_table(args.l.name, header=None)
    split_N = bed.shape[0] if args.S else N
    splits = np.array_split(range(bed.shape[0]), split_N)
    tmp_df = [bed.iloc[i] for i in splits]
    pre = str(time()).replace(".", "")
    for i in range(split_N):
        tmp_df[i].to_csv(str(pre)+'_'+str(i)+".bed", sep='\t', header=False, index=False)

    command = "bam-readcount"
    if args.q > 0: command += " -q {}".format(args.q)
    if args.b > 0: command += " -b {}".format(args.b)
    if args.d != 10000000: command += " -d {}".format(args.d)
    if args.f is not None: command += " -f {}".format(args.f.name)
    if args.D == 1: command += " -D {}".format(args.D)
    if args.p : command += " -p"
    if args.w != -1 : command += " -w {}".format(args.w)
    if args.i : command += " -i"

    commands = deque()
    for i in range(split_N):
        tname = str(pre)+'_'+str(i)+".rc"
        tmp_com = command + ' -l ' + str(pre)+'_'+str(i)+".bed" + ' ' + args.bam.name
        commands.append([tmp_com, tname])

    with Pool(processes=N) as pool:
        pool.starmap(run_bam_readcount, commands, chunksize=1)
    with open(args.o.name, 'w') as fout:
    # with open('tmp.txt', 'w') as fout:
        for i in range(split_N):
            with open(str(pre)+'_'+str(i)+".rc", 'r') as fin:
                count = 0
                s = deque()
                for line in fin:
                    count += 1
                    if count == 1000000:
                        fout.write("\n".join(s)+"\n")
                        s = deque()
                        count = 0
                    line = line.strip().split()
                    outstr = deque(line[:4])
                    for j in range(5, 10):
                        outstr.append(line[j].split(":")[1])
                    if INDEL:
                        for j in range(10, len(line)):
                            outstr.extend(line[j].split(":")[:2])
                    s.append("\t".join(outstr))
                fout.write("\n".join(s)+"\n")

    for i in range(split_N):
        os.remove(str(pre)+'_'+str(i)+".bed")
        os.remove(str(pre)+'_'+str(i)+".rc")


def run_ppileup(locs, out, bam, pars):
    from subprocess import Popen, PIPE
    with open(out, 'w') as fout:
        for loc in locs:
            # print("/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools mpileup {pars} -r {c}:{s}-{e} {bam}".format(pars=pars, c=loc[0], s=int(loc[1])+1, e=loc[2], bam=bam))
            p = Popen("/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools mpileup {pars} -r {c}:{s}-{e} {bam}".format(pars=pars, c=loc[0], s=int(loc[1])+1, e=loc[2], bam=bam).split(), stdout=fout, stderr=PIPE)
            o = p.communicate()


def ppileup(args):
    # print(args)
    NC = args.n
    BED = args.bed.name
    BAM = args.bam.name
    OUT = args.out.name
    PARAM = args.params if len(args.params) > 0 else ''

    from collections import deque
    import numpy as np
    from multiprocessing import Pool
    from functools import partial
    from time import time
    from subprocess import Popen, PIPE
    import os

    bedpos = deque()
    with open(BED, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            bedpos.append(line)
    splits = np.array_split(range(len(bedpos)), NC)
    # print(splits)
    bedsplit = [[bedpos[i] for i in split] for split in splits]
    # [print(bedsplit[i]) for i in range(NC)]
    pre = int(time())
    outs = [str(pre) + str(i) + '.tmp' for i in range(NC)]
    # outs = ['TMP'+ str(i) + '.tmp' for i in range(NC)]
    with Pool(processes=NC) as pool:
        pool.starmap(partial(run_ppileup, bam=BAM, pars=PARAM), zip(bedsplit, outs))
    with open(OUT, 'w') as fout:
        p = Popen("cat {outs} {OUT}".format(outs=' '.join(outs), OUT=OUT).split(), stdout=fout, stderr=PIPE)
    p.communicate()
    # [os.remove(o) for o in outs]
#END

def splitbam(args):
    BAM = args.bam.name
    TAG = args.tag
    perr = True
    try: import pysam
    except ModuleNotFoundError: print("Pysam not found. Exiting.")
    from collections import deque
    import sys
    rt = 'r' + args.f
    wt = 'w' + args.o
    if args.o == 'b': wex='.bam'
    elif args.o == 's': wex='.sam'
    elif args.o == 'c': wex='.cram'
    else:
        raise ValueError("Incorrect output file type: {}".format(args.o))
        sys.exit()
    b = pysam.AlignmentFile(BAM, rt)
    head = b.header
    # b = pysam.AlignmentFile('/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/tmp2.bam', rt)
    bc = ''
    bcreads = deque()
    for read in b:
        try:
            rbc = read.get_tag(TAG)
        except KeyError:
            if perr:
                print(TAG + " not found in {qname} at position {rname}:{rpos}. All reads without tag would be skipped.".format(qname=read.qname,rname= read.reference_name, rpos=read.reference_start))
                perr = False
            continue
        if rbc == bc: bcreads.append(read)
        elif bc == '':
                bc = rbc
                bcreads.append(read)
        else:
            if len(bcreads) > args.m:
                with pysam.AlignmentFile(bc+wex, wt, header=head) as f:
                    for r in bcreads:
                        f.write(r)
            bc = rbc
            bcreads = deque()
            bcreads.append(read)
    if len(bcreads) > args.m:
        with pysam.AlignmentFile(bc+wex, wt, header=b.head) as f:
            for r in bcreads:
                f.write(r)
#END

def gfatofa(args):
    gfa = args.gfa.name
    fa = args.fa.name
    from collections import OrderedDict
    gendict = OrderedDict()
    with open(gfa, 'r') as fin:
        for line in fin:
            if line[0] == 'S':
                line = line.strip().split()
                gendict[line[1]] = line[2]
    writefasta(gendict, fa)
    print("Finished")
#END

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()
    parser_getchr  = subparsers.add_parser("getchr", help="Get specific chromosomes from the fasta file")
    parser_exseq   = subparsers.add_parser("exseq", help="extract sequence from fasta", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_getscaf = subparsers.add_parser("getscaf", help="generate scaffolds from a given genome")
    parser_seqsize = subparsers.add_parser("seqsize", help="get size of dna sequences in a fasta file")
    parser_filsize = subparsers.add_parser("filsize", help="filter out smaller molecules")
    parser_subnuc  = subparsers.add_parser("subnuc", help="Change character (in all sequences) in the fasta file")
    parser_basrat  = subparsers.add_parser("basrat", help="Calculate the ratio of every base in the genome")
    parser_genome_ranges = subparsers.add_parser("genome_ranges", help="Get a list of genomic ranges of a given size", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_get_homopoly = subparsers.add_parser("get_homopoly", help="Find homopolymeric regions in a given fasta file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_n50 = subparsers.add_parser("asstat", help="Get N50 values for the given list of chromosomes", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_plthist = subparsers.add_parser("plthist", help="Takes frequency output (like from uniq -c) and generates a histogram plot", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_gfftrans = subparsers.add_parser("gfftrans", help="Get transcriptome (gene sequence) for all genes in a gff file. WARNING: THIS FUNCTION MIGHT HAVE BUGS.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_vcfdp = subparsers.add_parser("vcfdp", help="Get DP and DP4 values from a VCF file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_gffsort = subparsers.add_parser("gffsort", help="Sort a GFF file based on the gene start positions", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_getcol = subparsers.add_parser("getcol", help="Select columns from a TSV or CSV file using column names", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_bamcov = subparsers.add_parser("bamcov", help="Get mean read-depth for chromosomes from a BAM file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_pbamrc = subparsers.add_parser("pbamrc", help="Run bam-readcount in a parallel manner by dividing the input bed file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser_parallel = subparsers.add_parser("parallel", help="Run bam-readcount in a parallel manner by dividing the input bed file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_ppileup = subparsers.add_parser("ppileup", help="Currently it is slower than just running mpileup on 1 CPU. Might be possible to optimize later. Run samtools mpileup in parallel when pileup is required for specific positions by dividing the input bed file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_splitbam = subparsers.add_parser("splitbam", help="Split a BAM files based on TAG value. BAM file must be sorted using the TAG.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # TODO: Add functionality for sampling rows
    parser_samplerow = subparsers.add_parser("smprow", help="Select random rows from a text file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #
    parser_gfatofa = subparsers.add_parser("gfatofa", help="Convert a gfa file to a fasta file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()

    # gfatofa
    parser_gfatofa.set_defaults(func=gfatofa)
    parser_gfatofa.add_argument("gfa", help="Input GFA file", type=argparse.FileType('r'))
    parser_gfatofa.add_argument("fa", help="output fasta file name", type=argparse.FileType('w'))

    # Split BAM
    parser_splitbam.set_defaults(func=splitbam)
    parser_splitbam.add_argument("bam", help="BAM file", type=argparse.FileType('r'))
    parser_splitbam.add_argument("tag", help="tag to be used for splitting the BAM file", type=str)
    parser_splitbam.add_argument("-f", help="Alignments format BAM(b)/SAM(s)/CRAM(c) format", type=str, choices=['b', 's', 'c'], default='b')
    parser_splitbam.add_argument("-o", help="Output alignment format BAM(b)/SAM(s)/CRAM(c) format", type=str, choices=['b', 's', 'c'], default='b')
    parser_splitbam.add_argument("-m", help="Minimum number of reads required for a barcode", type=int, default=100)

    parser_ppileup.set_defaults(func=ppileup)
    parser_ppileup.add_argument("n", help="Number of CPU cores to use", type=int, default=1)
    parser_ppileup.add_argument("bed", help="Input bed file with regions to get pileup for", type=argparse.FileType('r'))
    parser_ppileup.add_argument("bam", help="Input bam file", type=argparse.FileType('r'))
    parser_ppileup.add_argument("out", help="Output file name", type=argparse.FileType('w'))
    parser_ppileup.add_argument("params", help="samtools mpileup parameters. Write within quotes", type=str)


    parser_pbamrc.set_defaults(func=pbamrc)
    parser_pbamrc.add_argument("bam", help="Input BAM file. Must have the SQ/SN/ tags for chromosome name and length", type=argparse.FileType('r'))
    parser_pbamrc.add_argument("-q", help="minimum mapping quality of reads used for counting.", type=int, default=0)
    parser_pbamrc.add_argument("-b", help="minimum base quality at a position to use the read for counting.", type=int, default=0)
    parser_pbamrc.add_argument("-d", help="max depth to avoid excessive memory usage.", type=int, default=10000000)
    parser_pbamrc.add_argument("-l", help="file containing a list of regions to report readcounts within.", type=argparse.FileType('r'), required=True)
    parser_pbamrc.add_argument("-f", help="reference sequence in the fasta format.", type=argparse.FileType('r'), required=True)
    parser_pbamrc.add_argument("-D", help="report the mapping qualities as a comma separated list.", default=0, type=int, choices=[0, 1])
    parser_pbamrc.add_argument("-p", help="report results by library.", default=False, action='store_true')
    parser_pbamrc.add_argument("-w", help="maximum number of warnings of each type to emit. -1 gives an unlimited number.", type=int, default=-1)
    parser_pbamrc.add_argument("-i", help="generate indel centric readcounts. Reads containing insertions will not be included in per-base counts", default=False, action='store_true')
    parser_pbamrc.add_argument("-I", help="Output indels as well.", default=False, action='store_true')
    parser_pbamrc.add_argument("-n", help="Number of CPU cores to use", type=int, default=1)
    parser_pbamrc.add_argument("-S", help="Run jobs in sequentially (each line from bed starts a new job). By default, jobs are run in batches (by dividing the BED file)", action='store_true', default=False)
    parser_pbamrc.add_argument("o", help="Output file name", type=argparse.FileType('w'))


    parser_bamcov.set_defaults(func=bamcov)
    parser_bamcov.add_argument("bam", help="Input BAM file. Must have the SQ/SN/ tags for chromosome name and length", type=argparse.FileType('r'))
    parser_bamcov.add_argument("out", help="Output file name", type=argparse.FileType('w'))
    parser_bamcov.add_argument("-r", help="BED file containing genomic regions to use", type=argparse.FileType('r'))
    parser_bamcov.add_argument("-g", help="Also output genomic mean coverage", default=False, action='store_true')
    parser_bamcov.add_argument("-d", help="Samtools depth command to use. Provide within double inverted quotes (\") ", type=str, default="samtools depth -d 0")


    parser_getcol.set_defaults(func=getcol)
    parser_getcol.add_argument("file", help="Input file", type=argparse.FileType('r'))
    parser_getcol.add_argument("out", help="Output file", type=argparse.FileType('w'))
    parser_getcol.add_argument("-s", help="Column separating char if not separated by tab/spaces", type=str)
    parser_getcol.add_argument("-m", help="Select column matching this string. Multiple space separated values can be provided", type=str, nargs='+')
    parser_getcol.add_argument("-c", help="Select column containing this string. Multiple space separated values can be provided", type=str, nargs='+')


    parser_gffsort.set_defaults(func=gffsort)
    parser_gffsort.add_argument("gff", help="Input GFF File", type=argparse.FileType('r'))
    parser_gffsort.add_argument("out", help="Output GFF File", type=argparse.FileType('w'))


    parser_vcfdp.set_defaults(func=vcfdp)
    parser_vcfdp.add_argument("vcf", help="VCF file", type=argparse.FileType('r'))
    parser_vcfdp.add_argument("-o", help="Output file name", type=argparse.FileType('w'), default='vcfdp.txt')
    

    parser_gfftrans.set_defaults(func=gfftrans)
    parser_gfftrans.add_argument("gff", help="Gff file", type=argparse.FileType('r'))
    parser_gfftrans.add_argument("fasta", help="Fasta file", type=argparse.FileType('r'))
    parser_gfftrans.add_argument("-o", help="Output file name", type=argparse.FileType('w'), default='transcriptome.fasta')
    parser_gfftrans.add_argument("-r", help="Reverse complement genes on -ve strand", default=False, action="store_true")


    parser_plthist.set_defaults(func=plthist)
    parser_plthist.add_argument("-f", help="Input file containing the frequency data", type=str, default='STDIN')
    parser_plthist.add_argument("-o", help="Output file name", type=argparse.FileType('w'), default='hist.pdf')
    parser_plthist.add_argument("-W", help="width of the plot (in inches)", type=float, default=4)
    parser_plthist.add_argument("-H", help="height of the plot (in inches)", type=float, default=4)
    parser_plthist.add_argument("-x", help="X-axis label", type=str, default='value', nargs='+')
    parser_plthist.add_argument("-y", help="Y-axis label", type=str, default='counts', nargs='+')
    parser_plthist.add_argument("-xlog", help="X-axis on log scale", default=False, action='store_true')
    parser_plthist.add_argument("-ylog", help="Y-axis on log scale", default=False, action='store_true')
    parser_plthist.add_argument("-xlim", help="Set X-axis limit for numerical data", nargs=2, type=int)
    parser_plthist.add_argument("-ylim", help="Set Y-axis limit for numerical data", nargs=2, type=int)
    parser_plthist.add_argument("-t", help="title of the plot", type=str, default=None)
    parser_plthist.add_argument("-n", help="Number of bins", type=int, default=100)


    parser_n50.set_defaults(func=getasstat)
    parser_n50.add_argument("-G", help="Path to genome/s to calculate n50", nargs='*')
    parser_n50.add_argument("-F", help="File containing path to genomes (one genome per line)", type=argparse.FileType('r'))
    parser_n50.add_argument("-o", help="Output file name", type=argparse.FileType('w'), default=None)
    parser_n50.add_argument("-n", help="number of processors to use", type=int, default=1)


    parser_get_homopoly.set_defaults(func=get_homopoly)
    parser_get_homopoly.add_argument("fasta", help="Input fasta file", type=argparse.FileType('r'))
    parser_get_homopoly.add_argument("-n", help="Minimum homopolymer length", type=int, default=10)
    parser_get_homopoly.add_argument("-o", help="output file name", type=str, default='homopolymer_regions.txt')
    parser_get_homopoly.add_argument("-p", help="size of neighbour padding (output size = n + 2*p)", type=int, default=0)
    parser_get_homopoly.add_argument("-b", help="output in bed file format", default=False, action="store_true")


    parser_genome_ranges.set_defaults(func = genome_ranges)
    parser_genome_ranges.add_argument("fasta", help="Input fasta file", type=argparse.FileType('r'))
    parser_genome_ranges.add_argument("-n", help="Range size", type=int, default=1000000)
    parser_genome_ranges.add_argument("-o", help="output file name", type=str, default='genome_ranges.txt')

    parser_getchr.set_defaults(func=getchr)
    parser_getchr.add_argument("fasta", help="fasta file", type=argparse.FileType('r'))
    parser_getchr.add_argument("-F", help="file containing list of chromosome to select", type=argparse.FileType('r'))
    parser_getchr.add_argument("--chrs", help="list of chomosome IDs to select", nargs='*')
    parser_getchr.add_argument("-v", help="Remove chromosomes provided", default=False, action='store_true')
    parser_getchr.add_argument("-o", help="Output file name", type=argparse.FileType('w'), default=None)
    
    parser_exseq.set_defaults(func=extractSeq)
    parser_exseq.add_argument("fasta", help="fasta file", type=argparse.FileType('r'))
    parser_exseq.add_argument("--chr", help="Chromosome ID", type=str)
    parser_exseq.add_argument("-start", help="Start location", type=int)
    parser_exseq.add_argument("-end", help="End location", type=int)
    parser_exseq.add_argument("--fin", help="File containing locations to extract (BED format)", type=argparse.FileType('r'))
    parser_exseq.add_argument("-o", help="Output file name", type=argparse.FileType('w'), default=None)
    
    parser_getscaf.set_defaults(func=getScaf)
    parser_getscaf.add_argument("fasta", help="genome fasta file", type=argparse.FileType('r'))
    parser_getscaf.add_argument("n", help="number of scaffolds required", type=int)
    parser_getscaf.add_argument("-r", help="reverse complement some scaffolds", default=False, action="store_true")
    parser_getscaf.add_argument("-o", help="output file name", default="scaf.fasta")

    parser_seqsize.set_defaults(func=seqsize)
    parser_seqsize.add_argument("fasta", help="genome fasta file", type=argparse.FileType('r'))
    parser_seqsize.add_argument("-t", help="also output total genome size", default=False, action="store_true")

    parser_filsize.set_defaults(func=filsize)
    parser_filsize.add_argument("fasta", help="genome fasta file", type=argparse.FileType('r'))
    parser_filsize.add_argument("size", help="molecule cut-off in bp, all smaller molecules will be filtered out", type=int, default=0)
    
    parser_subnuc.set_defaults(func=subnuc)
    parser_subnuc.add_argument("fasta", help="genome fasta file", type=argparse.FileType('r'))
    parser_subnuc.add_argument("q", help="character to change", default="", type=str)
    parser_subnuc.add_argument("t", help="character to change to", default="", type=str)
    parser_subnuc.add_argument("-o", help="output file name", type=argparse.FileType('w'))
    
    parser_basrat.set_defaults(func=basrat)
    parser_basrat.add_argument("fasta", help="genome fasta file", type=argparse.FileType('r'))
    
    args = parser.parse_args()

    from functools import reduce
    from itertools import cycle
    from sys import getsizeof, stderr
    from itertools import chain
    from collections import deque


#    print(args)
    args.func(args)

