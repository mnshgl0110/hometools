#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 15:36:01 2017

@author: goel
"""

import argparse
import os
import sys

import sys
 
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
            return sorted(range(len(args[0])), key = lambda k: args[0][k])
        elif list(args[1].lower()) == list('true'):
            return sorted(range(len(args[0])), key = lambda k: args[0][k], reverse = True)
        else:
            print("{} isn't a recognized parameter".format(args[1]))
            sys.exit()
    elif len(args) == 1:
        return sorted(range(len(args[0])), key = lambda k: args[0][k])
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


def plotDensity(data):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde
    density = gaussian_kde(data)
    xs = np.linspace(min(data), max(data),1000)
    density.covariance_factor = lambda : .2
    density._compute_covariance()
    plt.plot(xs,density(xs))
    plt.show()


def subList(lst1, lst2):
    import operator as op
    return(list(map(op.sub,lst1, lst2)))


def intersect(*lists):
    import numpy as np
    return reduce(np.intersect1d,list(lists))


def readfasta(f):
    from gzip import open as gzopen
    from collections import deque
    out = {}
    chrid = ''
    chrseq = deque()
    try:
        with gzopen(f, 'rb') as fin:
            for line in fin:
                break
                if b'>' in line:
                    if chrid != '':
                        out[chrid] = ''.join(chrseq)
                        chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                        chrseq = deque()
                    else:
                        chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                else:
                    chrseq.append(line.strip().decode())
    except OSError:
        try:
            with open(f, 'r') as fin:
                for line in fin:
                    if '>' in line:
                        if chrid != '':
                            out[chrid] = ''.join(chrseq)
                            chrid = line.strip().split('>')[1].split(' ')[0]
                            chrseq = deque()
                        else:
                            chrid = line.strip().split('>')[1].split(' ')[0]
                    else:
                        chrseq.append(line.strip())
        except Exception as e:
            raise Exception(e)
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


def extractSeq(args):
    import pandas as pd
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqIO import parse, write
    filePath = args.fasta.name
    if args.fin == None:
       seqID = args.loc[0]
       querySeq = [fasta for fasta in parse(filePath, 'fasta') if fasta.id == seqID][0]

       start = int(args.loc[1]) if int(args.loc[1]) >= 0 else 0
       end = int(args.loc[2]) if int(args.loc[2]) <= len(querySeq.seq) else len(querySeq.seq)

       querySeq.seq = querySeq.seq[int(start):(int(end)+1)]
       if args.o != None:
           write(querySeq, args.o, "fasta")
       else:
           print("> "+querySeq.id)
           print(querySeq.seq)
    else:
        fin = pd.read_table(args.fin.name, header=None, delim_whitespace=True)
        fin.columns = ["chr", "start", "end"]
        fin[['start', 'end']] = fin[['start', 'end']].astype('int')
        fin.loc[fin['start'] < 0, 'start'] = 0
        fin.sort_values(["chr", "start", "end"], inplace=True)
        outF = deque()
        for fasta in parse(filePath, 'fasta'):
            if fasta.id in fin.chr.values:
                chrData = fin.loc[fin.chr == fasta.id].copy()
                chrData.loc[chrData['end'] > len(fasta.seq), 'end'] = len(fasta.seq)
                for row in chrData.itertuples(index=False):
                    outF.append(SeqRecord(seq=fasta.seq[row.start:row.end], id="_".join(map(str, row)), description=""))
        if args.o != None:
            write(outF, args.o.name, "fasta")
        else:
           for i in outF:
                print("> "+i.id)
                print(i.seq)
                

def revcomp(seq):
    assert type(seq) == str
    old = 'ACGTacgt'
    rev = 'TGCAtgca'
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
    Take a 2D numpy array, with each row as a range and return merged ranges i.e. ranges which are overlapping would be
    combined as well.
    :param ranges:
    :return:
    """
    from collections import deque
    import numpy as np
    if len(ranges) < 2:
        return ranges
    ranges = ranges[ranges[:, 0].argsort()]
    for i in ranges:
        if i[0] > i[1]:
            garb = i[0]
            i[0] = i[1]
            i[1] = garb
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
    from Bio.SeqIO import parse, write
    fin = args.fasta.name
    out = [(fasta.id, len(fasta.seq)) for fasta in parse(fin,'fasta')]
    for i in out:
        print(i[0],i[1],sep="\t")
    if args.t:
        print("Genome_length",sum([i[1] for i in out]), sep="\t")


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


def basrat(args):
    from Bio.SeqIO import parse
    from collections import Counter 
    fin = args.fasta.name
    count = {fasta.id:dict(Counter(fasta.seq)) for fasta in parse(fin, "fasta")}
    outD = {}
    for v in count.values():
        for k, cnt in v.items():
            if k in outD:
                outD[k] += cnt
            else:
                outD[k] = cnt
    for k,v in outD.items():
        print(k, v, sep = "\t")
    

def getchr(args):
    from Bio.SeqIO import parse
    from gzip import open as gzopen
    fin = args.fasta.name
    if args.chrs is None and args.F is not None: chroms = [c.strip() for c in open(args.F.name, 'r').readlines()]
    elif args.chrs is not None and args.F is None: chroms = args.chrs
    else: raise Exception('InvalidValue: Incorrect value for chrs provided')
    out = args.o.name if args.o is not None else fin + ".filtered.fasta"
    with open(out, 'w') as fout:
        try:
            finh = gzopen(fin, 'rt')
            for fasta in parse(finh, 'fasta'):
                if (fasta.id in chroms) != args.v:
                    fout.write(">" + fasta.id + "\n" + "\n".join([str(fasta.seq[i:i+60]) for i in range(0, len(fasta.seq), 60)]) + '\n')
        except OSError:
            finh = open(fin, 'rt')
            for fasta in parse(finh, 'fasta'):
                if (fasta.id in chroms) != args.v:
                    fout.write(">" + fasta.id + "\n" + "\n".join([str(fasta.seq[i:i+60]) for i in range(0, len(fasta.seq), 60)]) + '\n')


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


def get_homopoly(args):
    import os
    import sys

    if args.n < 2:
        sys.exit('Minimum allowed homopolymer length = 2')
    if os.path.isfile(args.o):
        sys.exit(args.o + ' exists. Cannot overwrite it.')

    import re
    from collections import OrderedDict
    from Bio.SeqIO import parse
    isbed = 1 if args.b else 0

    start = ''
    with open(args.o, 'w') as fout:
        for fasta in parse(args.fasta, 'fasta'):
            s = str(fasta.seq)
            id = fasta.id
            outdict = {}
            for b in ['A', 'C', 'G', 'T']:
                hp = b * args.n + '+'
                # hp = b * 4
                r = re.finditer(hp, s)
                for i in r:
                    outdict[i.start() - isbed] = [i.end(), b]
            for i in sorted(outdict.keys()):
                fout.write(start + str(id) + "\t" + str(i+1 - args.p) + "\t" + str(outdict[i][0] + args.p) + "\t" + outdict[i][1])
                start = '\n'


def asstat(fin, parse):
    from numpy import cumsum
    fasta_len = deque()
    for fasta in parse(fin, 'fasta'):
        fasta_len.append(len(str(fasta.seq)))
    fasta_len = sorted(fasta_len)
    fasta_len_cumsum = cumsum(fasta_len)
    half = fasta_len_cumsum[-1]/2
    for i in range(len(fasta_len_cumsum)):
        if fasta_len_cumsum[i] > half:
            return [fasta_len_cumsum[-1], len(fasta_len), fasta_len[-1], fasta_len[0], fasta_len[i], len(fasta_len)-i]


def getasstat(args):
    out = args.o.name if args.o is not None else "genomes.n50"
    NC = args.n

    from collections import deque
    fins = deque()
    if args.F is None and args.G is None:
        sys.exit("1Provide genome file path in -F or -G")
    elif args.F is not None and args.G is not None:
        sys.exit("Provide genome file path in -F or -G")
    elif args.F is not None:
        with open(args.F.name, 'r') as F:
            for line in F:
                fins.append(line.strip())
    elif args.G is not None:
        fins = args.G

    from multiprocessing import Pool
    from functools import partial
    from Bio.SeqIO import parse
    with Pool(processes=NC) as pool:
        n50_values = pool.map(partial(asstat, parse=parse), fins)

    with open(out, 'w') as fout:
        fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("assemble_id", "assembly_length", "number_of_contig", "longest_contig", "shortest_contig", "N50", "L50"))
        for i in range(len(fins)):
            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(fins[i], n50_values[i][0], n50_values[i][1], n50_values[i][2], n50_values[i][3], n50_values[i][4], n50_values[i][5]))


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

    from matplotlib import pyplot as plt
    fig = plt.figure(figsize = [args.W, args.H])
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



if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()
    parser_getchr  = subparsers.add_parser("getchr", help="Get specific chromosomes from the fasta file")
    parser_exseq   = subparsers.add_parser("exseq", help="extract sequence from fasta")
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

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()


    parser_getcol.set_defaults(func=getcol)
    parser_getcol.add_argument("file", help="Input file File", type=argparse.FileType('r'))
    parser_getcol.add_argument("out", help="Output file File", type=argparse.FileType('w'))
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
    parser_plthist.add_argument("-t", help="title of the plot", type=str, default=None)
    parser_plthist.add_argument("-n", help="Number of bins", type=int, default=100)


    parser_n50.set_defaults(func=getasstat)
    parser_n50.add_argument("-F", help="File containing path to genomes (one genome per line)", type=argparse.FileType('r'))
    parser_n50.add_argument("-G", help="Path to genome/s to calculate n50", nargs='*')
    parser_n50.add_argument("-o", help="Output file name", type=argparse.FileType('w'), default=None)
    parser_n50.add_argument("-n", help="number of processors to use", type=int, default=1)


    parser_get_homopoly.set_defaults(func=get_homopoly)
    parser_get_homopoly.add_argument("fasta", help="Input fasta file", type=argparse.FileType('r'))
    parser_get_homopoly.add_argument("-n", help="Minimum homopolymer length", type=int, default=10)
    parser_get_homopoly.add_argument("-o", help="output file name", type=str, default='homopolymer_regions.txt')
    parser_get_homopoly.add_argument("-p", help="size of neighbour padding (output size = n + 2*p)", type=int, default=0)
    parser_get_homopoly.add_argument("-b", help="output in bed file format", default=False, action="store_true")


    parser_genome_ranges.set_defaults(func=genome_ranges)
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
    group = parser_exseq.add_mutually_exclusive_group(required=True)
    group.add_argument("--loc", help="Location to extract: chr start end", nargs=3)
    group.add_argument("--fin", help="File containing locations to extract", type=argparse.FileType('r'))
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

