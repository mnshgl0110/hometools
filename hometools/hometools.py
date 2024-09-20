#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 15:36:01 2017

@author: Manish Goel
"""

import argparse
import os
import sys


class Namespace:
    """
    Use this to create args object for parsing to functions
    args = Namespace(a=1, b='c')
    print(args.a, args.b)
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
# END
inttr = lambda x: [int(x[0]), x[1]]


############################# Set logger #######################################

def setlogconfig(lg):
    """
    :param lg: Log-level
    :return:
    """
    import logging.config
    logging.config.dictConfig({
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
            'log_file': {
                'format': "%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s",
            },
            'stdout': {
                'format': "%(name)s - %(levelname)s - %(message)s",
            },
        },
        'handlers': {
            'stdout': {
                'class': 'logging.StreamHandler',
                'formatter': 'stdout',
                'level': 'WARNING',
            },
        },
        'loggers': {
            '': {
                'level': lg,
                'handlers': ['stdout'],
                # 'handlers': ['stdout', 'log_file'],
            },
        },
    })
# END


def mylogger(logname):
    from hometools.classes import CustomFormatter
    import logging
    logger = logging.getLogger(logname)
    handler = logging.StreamHandler()
    handler.setFormatter(CustomFormatter())
    logger.addHandler(handler)
    logging.basicConfig(level=logging.INFO)
    logger.propagate = False
    return logger
# END


logger = mylogger(__name__)

############################# FASTA ############################################
def readfasta(f):
    # TODO: This takes too long when used with getchr for large genomes. Try to optimise FASTA reading when the entire genome is not needed.
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
# END


def writefasta(fa, f):
    """
    :param fa: dictionary. Keys are chromosome ids. Values are sequence.
    :param f: Output file
    :return:
    Can output .gz file if output file name ends with .gz
    """
    # TODO: write bgzip file
    # from pysam import tabix_compress as gzopen
    from gzip import open as gzopen
    logger = mylogger("writefasta")
    isgzip = f.rsplit('.', 1)[-1] == 'gz'
    # b, nl >> Begin character, new line
    op, openstr, b, nl = (gzopen, 'wb', b'>', b'\n') if isgzip else (open, 'w', '>', '\n')
    with op(f, openstr) as fo:
        for k, v in fa.items():
            if isgzip:
                k = k.encode()
                v = v.encode()
            fo.write(b+k+nl)
            fo.write(nl.join([v[i:i+60] for i in range(0, len(v), 60)]) + nl)
# END


def readfasta_iterator(fin, isgzip=False):
    # from gzip import open as gzopen
    from collections import deque
    """
    :param fin: File handle (already opened using normal open or gzip.open)
    :param isgzip: True when reading gzipped files, false for text files
    """
    seq = deque()
    b, nl, empty, sep = (b'>', b'\n', b'', b' ') if isgzip else ('>', '\n', '', ' ')
    title = empty
    out = lambda x: (x[0].decode(), x[1].decode()) if isgzip else x
    for line in fin:
        if line.strip() == empty: continue
        if line[0] == b[0]:
            if title == empty:
                title = line[1:].rstrip().split(sep)[0]
            else:
                # yield title.decode(), empty.join(seq).decode()
                yield out((title, empty.join(seq)))
                seq = deque()
                title = line[1:].rstrip().split(sep)[0]
        else:
            seq.append(line.strip())
    yield out((title, empty.join(seq)))
# END


def readfaidxbed(f):
    """
    Reads .faidx file from a genome assembly and returns a BED file consisting
    for entire chromosomes
    """
    from collections import deque
    import pybedtools as bt
    fabed = deque()
    with open(f, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            fabed.append([line[0], 0, int(line[1])])
    return list(fabed)
# END


############################# CIGAR ############################################

def cgtpl(cg, to_int=False):
    """
    Takes a cigar string as input and returns a cigar tuple
    """
    for i in "MIDNSHPX=":
        cg = cg.replace(i, ';'+i+',')
    if to_int:
        return [inttr(i.split(';')) for i in cg.split(',')[:-1]]
    else:
        return [i.split(';') for i in cg.split(',')[:-1]]
# END


def cgstr(cg):
    """
    Converts a CIGAR tuple into a CIGAR string
    """
    return "".join(["".join([str(i[0]), i[1]]) for i in cg])
# END


def cggenlen(cg, gen, full=False):
    """
    Takes cigar as input, and return the number of bases covered by it the reference
    or query genome.
    Cigar strings are first converted to cigar tuples.
    Use full=True to get the query length of the entire sequence and not just the aligned region
    """
    if type(cg) == str:
        cg = cgtpl(cg)
    try:
        assert(gen in ['r', 'q'])
    except AssertionError:
        raise ValueError('gen need to "r" or "q" for reference or query')
    if gen == 'r':
        s = {'M', 'D', 'N', '=', 'X'}
    elif full:
        s = {'M', 'I', 'S', '=', 'X', 'H'}
    else:
        s = {'M', 'I', 'S', '=', 'X'}
    l = sum([int(i[0]) for i in cg if i[1] in s])
    return l
# END


def cgwalk(cg, n, ngen='r'):
    """
    Returns how many bases would be travelled in the non-focal genome after moving n bases in the focal (ngen) genome.
    """
    cnt = 0                   # Count in the focal genome
    ocnt = 0                  # Count in the other genome
    nset = {'M', 'D', 'N', '=', 'X'} if ngen == 'r' else {'M', 'I', '=', 'X'}
    oset = {'M', 'I', '=', 'X'} if ngen == 'r' else {'M', 'D', 'N', '=', 'X'}
    for c in cg:
        if c[1] in {'H', 'P', 'S'}: continue
        if c[1] in nset:
            if n <= cnt + c[0]:
                if c[1] in {'M', '=', 'X'}:
                    return ocnt + (n - cnt)
                elif c[1] in {'I', 'D', 'S', 'N'}:
                    return ocnt
            else:
                cnt += c[0]
        if c[1] in oset:
            ocnt += c[0]
    raise ValueError("n is larger than the number of bases covered by cigartuple")
# END

############################# logging ##########################################

############################# UTIL #############################################

def isgzip(f):
    '''
    Checks if a given file is gzipped. Returns true if the is gzipped.
    '''
    from gzip import open as gzopen
    from gzip import BadGzipFile
    # Test if the file is Gzipped or not
    with gzopen(f, 'rb') as fin:
        try:
            fin.read(1)
            isgz = True
        except BadGzipFile:
            isgz = False
    return isgz
# END


def randomstring(l):
    """
    l = length of the string
    returns a random string of size l
    """
    from random import choices
    from string import ascii_letters as letters
    return ''.join(choices(letters, k=l))
# END


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
# END


################################ RANGES FUNCTIONS ###############################
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
# END


def sumranges(r):
    """
    Takes a 2D numpy array of non-overlapping ranges and output total length of the ranges
    """
    return sum(r[:, 1] - r[:, 0] + 1)
# END


def subranges(r1, r2):
    """
    Subtract range2 (r2) from range1 (r1) and return non-overlaping range.
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
# END


def findoverlaps(r1, r2):
    '''
    Find r2 elements that overlap with r1 elements
    '''
    # pybedtools seems to be working quite well.
    # bed1 = BedTool.from_dataframe(df1)
    # bed2 = BedTool.from_dataframe(df2)
    # for b in bed1.window(bed2, w=0):
    #   print(b)
    # df columns: [chr, start, end, ....]
    pass
# END


############################# Other ############################################


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
# END


def cumminf(array):
    cummin = []
    cumulative_min = array[0]
    for p in array:
        if p < cumulative_min:
            cumulative_min = p
        cummin.append(cumulative_min)
    return cummin
# END


def cummaxf(array):
    cummax = []
    cumulative_max = array[0]
    for e in array:
        if e > cumulative_max:
            cumulative_max = e
        cummax.append(cumulative_max)
    return cummax
# END


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
# END


def p_adjust(pvalues, method='bh'):
    """
    Takes list of values and do multiple hypothesis adjustment
    copied from
    https://rosettacode.org/wiki/P-value_correction#Python

    """
    # method = "bh"
    # pvalues = args[0]
    # if len(args) > 1:
    methods = {"bh", "fdr", "by", "holm", "hommel", "bonferroni", "hochberg"}
        # metharg = args[1].lower()
    if method not in methods:
            logger.error(f"Incorrect method provided. Available methods: {methods}")

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
            cummin_input.insert(index, (n/(n-index)) * pvalues[o[index]])
        ro = order(o)
        cummin = cumminf(cummin_input)
        pmin = pminf(cummin)
        qvalues = [pmin[i] for i in ro]
    elif method == 'by':
        q = 0.0
        o = order(pvalues, 'TRUE')
        ro = order(o)
        for index in range(1, n+1):
            q += 1.0 / index
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
        # sys.exit()
    return qvalues
# END


def summary(arr):
    '''
    Takes a list or np.array and print summary statistics similar to the summary()
    function in R
    '''
    import numpy as np
    arr = np.array(arr)
    print(f'Minimum: {round(np.min(arr), 6)}\n1st Qu.: {round(np.quantile(arr, 0.25), 6)}\nMedian: {round(np.median(arr), 6)}\nMean: {round(np.mean(arr), 6)}\n3rd Qu.: {round(np.quantile(arr, 0.75), 6)}\nMax: {round(np.max(arr), 6)}')
    return
# END


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
# END


def undict(input):
    """
    Takes a nested dict and converts it into a list of list. `lvl` corresponds to
    the number of iterations to be done.

    Example:
    ```
    d = {'a': {'a1': {'a2': 1}}, 'b': {'b1': {'b2': 2}}}

    d1 = undict(d)
    print(d1)
    [['a', 'a1', 'a2', 1], ['b', 'b1', 'b2', 2]]
    ```
    """
    from collections import deque
    from collections.abc import Iterable
    def getlist(var):
        """
        takes an object as input and creates a list with the object as its element
        """
        if isinstance(var, Iterable):
            if not isinstance(var, str):
                return list(var)
            else:
                return [var]
        else:
            return [var]
    # END
    if isinstance(input, dict):
        outdict = deque()
        for k in input.keys():
            childdict, wasdict = undict(input[k])
            if wasdict:
                outdict.extend([getlist(k) + i for i in childdict])
            else:
                outdict.extend([getlist(k) + childdict])
        return outdict, True
    else:
        return getlist(input), False

# END


def getvalues(l, index):
    """from list l get the values at indices specified by index"""
    return [l[i] for i in index]
# END


def sublist(lst1, lst2):
    import operator as op
    return(list(map(op.sub,lst1, lst2)))
# END


def intersect(*lists):
    import numpy as np
    from functools import reduce
    return reduce(np.intersect1d, list(lists))
# END





############################## API #############################################
## Functions callable from python scripts


def revcomp(seq):
    # end = end if end <= len(q[seqid]) else len(q[seqid])
    # # Output the selected sequence
    assert(type(seq) == str)
    old = 'ACGTRYKMBDHVacgtrykmbdhv'
    rev = 'TGCAYRMKVHDBtgcayrmkvhdb'
    tab = str.maketrans(old, rev)
    return seq.translate(tab)[::-1]
# END


def canonical(seq):
    """
    For an input string (seq, k-mer), return the canonical K-mer which is defined as
    the minimum string between seq and revcomp(seq)
    """
    return min(seq, revcomp(seq))
# END


def canonical_kmers(size):
    """
    For a given size (n), return all canonical k-mers of size (n)
    """
    from itertools import product
    return sorted(set(map(canonical, map(''.join, product(*['ATGC']*size)))))
# END


def fileRemove(fName):
    try:
        os.remove(fName)
    except OSError as e:
        if e.errno != 2:    ## 2 is the error number when no such file or directory is present https://docs.python.org/2/library/errno.html
            raise
# END


def view(d, n=5):
    """
    For some objects, like dict, it is difficult to view just a part of it (like head).
    This function provides a `head` like utility.

    NOTE: Should not be used to iterators
    """
    # TODO: Make it work properly
    try:
        iter = d.__iter__()
        # TODO: Add iterable for rows of dataframe
        # if type(g).__name__ == 'DataFrame':

    except AttributeError as e:
        print(f"Object of type: {type(d)} is not iterable. Exiting")
        return
    # dicttypes = {'dict', 'defaultdict', }
    count = 0
    isdict = True if 'dict' in str(type(d)).lower() else False
    for i in iter:
        print(f"{i}: {d[i]}" if isdict else i)
        count += 1
        if count == n: break
    return
# END


def printdf(df, n=5):
    """
    Takes a DataFrame as input and print all columns
    """
    logger.info(f'printing first {n} rows')
    print(df.head(n).to_string())
    return
# END


def rtigercos(bed):
    """
    Reads RTIGER output BED file and return the CO table
    :param bed:
    :return:
    """
    import pandas as pd
    from collections import deque
    bed = pd.read_table(bed, header=None)
    cos = deque()
    for g in bed.groupby([0]):
        for i in range(g[1].shape[0] - 1):
            cos.append([g[0], g[1].iat[i, 2], g[1].iat[i+1, 1], g[1].iat[i, 3], g[1].iat[i+1, 3]])
    return pd.DataFrame(cos)
# END


def total_size(o, handlers={}, verbose=False):
    from collections import deque
    from sys import getsizeof, stderr
    from itertools import chain
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
# END


def readblast(f):
    """
    Read tabular format blast output
    """
    import pandas as pd
    return pd.read_table(f, comment="#")
# END


def readsyriout(f, annos=['SYN', 'SYNAL', 'INV', 'TRANS', 'INVTR', 'DUP', 'INVDP']):
    """
    There is a version in reganno.py as well
    """
    from collections import OrderedDict, deque
    import logging
    import pandas as pd
    import numpy as np

    logger = logging.getLogger("readsyriout")
    # Reads syri.out. Select: achr, astart, aend, bchr, bstart, bend, srtype
    # assert(all([v in ['SYN', 'SYNAL', 'INV', 'TRANS', 'INVTR', 'DUP', 'INVDP', 'NOTAL'] for v in annos]))
    syri_regs = deque()
    with open(f, 'r') as fin:
        for line in fin:
            l = line.strip().split()
            if l[10] in annos:
                syri_regs.append(l)

    try:
        df = pd.DataFrame(list(syri_regs))[[0, 1, 2, 5, 6, 7, 10, 3, 4]]
        df = df.loc[df[0] != '-']
    except KeyError:
        raise ImportError("Incomplete input file {}, syri.out file should have 11 columns.".format(f))
    df[[0, 5, 10]] = df[[0, 5, 10]].astype(str)
    try:
        if 'NOTAL' in annos:
            df[[1, 2]] = df[[1, 2]].astype(int)
        else:
            df[[1, 2, 6, 7]] = df[[1, 2, 6, 7]].astype(int)
    except ValueError:
        raise ValueError("Non-numerical values used as genome coordinates in {}. Exiting".format(f))

    df.sort_values([0, 1, 2], inplace=True)
    # chr ID map: DO NOT DELETE AS WILL BE USED IN SOME FUNCTIONS
    # chrid = []
    # chrid_dict = OrderedDict()
    # for i in np.unique(df[0]):
    #     chrid.append((i, np.unique(df.loc[(df[0] == i) & (df[10] == 'SYN'), 5])[0]))
    #     chrid_dict[i] = np.unique(df.loc[(df[0] == i) & (df[10] == 'SYN'), 5])[0]
    # df.columns = ['achr', 'astart', 'aend', 'bchr', 'bstart', 'bend',  'type']

    return df
# END


def samtocoords(f):
    '''
    Reads a SAM file and converts it to a align coords file
    '''
    # TODO: adjust this function to work similarly to bam2coords
    from pandas import DataFrame
    from collections import deque
    # logger = logging.getLogger('SAM reader')
    logger = mylogger("SAM reader")
    rc = {}        # Referece chromosomes
    rcs = {}        # Selected chromosomes
    al = deque()    # Individual alignment
    try:
        with open(f, 'r') as fin:
            for l in fin:
                if l[:3] == '@SQ':
                    c, s = 0, 0
                    for h in l.strip().split()[1:]:
                        h = h.split(':')
                        if h[0] == 'SN': c = h[1]
                        if h[0] == 'LN': s = int(h[1])
                    rcs[c] = s
                    continue
                elif l[0] == '@': continue

                l = l.split('\t')[:6]
                # if l[1] == '2064': break
                if l[2] == '*':
                    logger.warning(l[0]+ ' do not align with any reference sequence and cannot be analysed. Remove all unplaced scaffolds and contigs from the assemblies.')  # Skip rows corresponding to non-mapping sequences (contigs/scaffolds)
                    continue

                if 'M' in l[5]:
                    logger.error('Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: ' + l[5])
                    sys.exit()
                cgt = [[int(j[0]), j[1]] for j in [i.split(';') for i in l[5].replace('S', ';S,').replace('H', ';H,').replace('=', ';=,').replace('X', ';X,').replace('I', ';I,').replace('D', ';D,').split(',')[:-1]]]
                if len(cgt) > 2:
                    if True in [True if i[1] in ['S', 'H'] else False for i in cgt[1:-1]]:
                        logger.error("Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + l[5])
                        sys.exit()

                bf = '{:012b}'.format(int(l[1]))

                rs = int(l[3])
                re = rs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'D']])

                if bf[7] == '0':    # forward alignment
                    if cgt[0][1] == '=':
                        qs = 1
                    elif cgt[0][1] in ['S', 'H']:
                        qs = cgt[0][0] + 1
                    else:
                        print('ERROR: CIGAR string starting with non-matching base')
                    qe = qs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
                elif bf[7] == '1':  # inverted alignment
                    if cgt[-1][1] == '=':
                        qs = 1
                    elif cgt[-1][1] in ['S', 'H']:
                        qs = cgt[-1][0] + 1
                    else:
                        print('ERROR: CIGAR string starting with non-matching base')
                    qe = qs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
                    qs, qe = qe, qs

                al.append([
                    rs,
                    re,
                    qs,
                    qe,
                    abs(re-rs) + 1,
                    abs(qs-qe) + 1,
                    format((sum([i[0] for i in cgt if i[1] == '=']) / sum(
                        [i[0] for i in cgt if i[1] in ['=', 'X', 'I', 'D']])) * 100, '.2f'),
                    1,
                    1 if bf[7] == '0' else -1,
                    l[2],
                    l[0],
                    "".join([str(i[0])+i[1] for i in cgt if i[1] in ['=', 'X', 'I', 'D']])
                ])
                rcs[l[2]] = 1
            rcs = list(rcs.keys())
            for k in list(rc.keys()):
                if k not in rcs: logger.warning(l[0]+ ' do not align with any query sequence and cannot be analysed. Remove all unplaced scaffolds and contigs from the assemblies.')
    except Exception as e:
        logger.error('Error in reading SAM file: ' + str(e))
        sys.exit()
    al = DataFrame(list(al))
    al[6] = al[6].astype('float')
    al.sort_values([9,0,1,2,3,10], inplace = True, ascending=True)
    al.index = range(len(al.index))
    return al
# END


def extractseq(args):
    # TODO: Fix this function to work with new parameter style
    """
    parser_exseq.set_defaults(func=extractseq)
    parser_exseq.add_argument("fasta", help="fasta file", type=argparse.FileType('r'))
    parser_exseq.add_argument("-c", "--chr", help="Chromosome ID", type=str)
    parser_exseq.add_argument("-s", "--start", help="Start location (BED format, 0-base, end included)", type=int)
    parser_exseq.add_argument("-e", "--end", help="End location (BED format, 0-base, end excluded)", type=int)
    parser_exseq.add_argument("--fin", help="File containing locations to extract (BED format)", type=argparse.FileType('r'))
    parser_exseq.add_argument("-o", help="Output file name", type=argparse.FileType('w'), default=None)
    parser_exseq.add_argument("-r", help="Reverse complement the sequence", default=False, action='store_true')

    """
    import pandas as pd
    import warnings
    if args.fin is not None:
        if args.chr is not None or args.start is not None or args.end is not None:
            warnings.warn("Using --fin. Ignoring --chr, -s, -e")
    fasta = args.fasta.name
    rev = args.r
    output = args.o.name if args.o is not None else args.o

    if args.fin is None:
        chrom = args.chr
        start = int(args.start) if args.start is not None else 0
        end = int(args.end) if args.end is not None else args.end
        for c, s in readfasta_iterator(open(fasta, 'r'), isgzip(fasta)):
            if c == chrom:
                end = end if end is not None else len(s)
                outseq = revcomp(s[start:end]) if rev else s[start:end]
                break
        # if len(q) > 1:
        #     sys.exit("Found multiple chromosomes with same ID. Exiting.")
        # start = int(args.start) if args.start is not None else 0
        # end = int(args.end) if args.end is not None else len(q[chrom])
        # outseq = revcomp(q[chrom][start:end]) if r else q[chrom][start:end]
        if output is not None:
            writefasta({chrom: outseq}, output)
        else:
            print(">{}\n{}".format(chrom, outseq))
    else:
        fin = pd.read_table(args.fin.name, header=None, delim_whitespace=True)[[0, 1, 2]]
        fin.columns = ["chr", "start", "end"]
        fin[['start', 'end']] = fin[['start', 'end']].astype('int')
        fin.loc[fin['start'] < 0, 'start'] = 0
        fin.sort_values(["chr", "start", "end"], inplace=True)
        out = dict()
        chroms = set(fin.chr.values)
        for c, s in readfasta_iterator(open(fasta, 'r'), isgzip(fasta)):
            if c in chroms:
                cdf = fin.loc[fin.chr == c].copy()
                cdf.loc[cdf['end'] > len(s), 'end'] = len(s)
                for row in cdf.itertuples(index=False):
                    out['{}_{}_{}'.format(c, row[1], row[2])] = revcomp(s[row.start:row.end]) if rev else s[row.start:row.end]
        # Output the selected sequence
        if output is not None:
            writefasta(out, output)
        else:
            for c, s in out.items():
                print(">{}\n{}".format(c, s))
    return
# END


def subnuc(args):
    from Bio.SeqIO import parse, write
    from Bio.Seq import Seq
    fasta = args.fasta.name
    querySeq = [fasta for fasta in parse(fasta, 'fasta')]
    for i in range(len(querySeq)):
        querySeq[i].seq = Seq(str(querySeq[i].seq).replace(args.q, args.t))
#    print(querySeq)

    if args.o == None:
        fout = fasta+"_edited"
    else:
        fout = args.o
    with open(fout, "w") as f:
        spacer = ""
        for seq in querySeq:
            f.write(spacer)
            f.write(">"+seq.id+" "+seq.description+"\n")
            f.write('\n'.join(str(seq.seq)[i:i+60] for i in range(0, len(seq.seq), 60)))
            if spacer == "":
                spacer = "\n"
# END

## Plotting

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
    plt.xticks(rotation=args.xltilt)
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig(args.o.name)
    fin.close()
# END


def pltbar(args):
    '''
    Generate a bar plot for the input. First column needs to be features, second
    column values
    '''
    import sys
    from collections import OrderedDict
    logger = mylogger("pltbar")
    if args.f == 'STDIN':
        fin = sys.stdin

    else:
        try:
            fin = open(args.f, 'r')
        except FileNotFoundError:
            logger.error('Cannot open file: {}'.format(args.f))
            sys.exit()
    # Read data
    data = OrderedDict()
    for line in fin:
        line = line.strip().split()
        if len(line) > 2:
            logger.warning('Input has more than 2 columns. Using the first two columns')
        try:
            data[line[0]] = float(line[1])
        except ValueError:
            logger.error('Second column contains non numerical values. Exiting.')
            sys.exit()
    pltdata = data
    if args.sx:
        keys = sorted(list(pltdata.keys()))
    elif args.sy:
        keys = sorted(list(pltdata.keys()), key=lambda x: pltdata[x])
    else:
        keys = list(pltdata.keys())
    from matplotlib import use as mpuse
    mpuse('agg')
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=[args.W, args.H])
    ax = fig.add_subplot()
    ax.bar(keys, [int(pltdata[k]) for k in keys])
    ax.set_xlabel(' '.join(args.x))
    ax.set_ylabel(' '.join(args.y))
    if args.t is not None:
        ax.set_title(args.t)
    if args.ylog:
        ax.set_yscale('log')
    if args.ylim is not None:
        ax.set_ylim([args.ylim[0], args.ylim[1]])
    ax.minorticks_on()
    ax.xaxis.grid(True, which='both', linestyle='--')
    ax.yaxis.grid(True, which='both', linestyle='--')
    plt.xticks(rotation=args.xltilt)
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig(args.o.name)
    fin.close()
# END


def plotal(args):
    """
    Input file format:
    genome1:chr1:start-end  genome2:chr2:start-end  Colour
    1:1:2-15        2:1:1-14        #006c66
    2:1:1-14        3:1:3-16        #006c66
    3:1:3-16        4:1:1-14        #006c66

    :param args:
    :return:
    """
    from collections import deque
    import pandas as pd
    from matplotlib import pyplot as plt
    from matplotlib.pyplot import get_cmap
    import matplotlib
    from hometools.plot import bezierpath
    # Parse arguments
    finname = args.align.name
    out = args.out.name
    DPI = args.D
    # Read alignments
    als = deque()
    with open(finname, 'r') as fin:
        for line in fin:
            if line[0] == '#':
                continue
            line = line.strip().split()
            if len(line) == 0:
                continue
            s = line[0].split(':')
            s += s[2].split('-')
            s.pop(2)
            e = line[1].split(':')
            e += e[2].split('-')
            e.pop(2)
            als.append(s + e + [line[2]])
    als = pd.DataFrame(als)
    als[[2, 3, 6, 7]] = als[[2, 3, 6, 7]].astype(int)
    ngen = len(set(pd.concat([als[0], als[4]])))
    gdict = dict(zip(pd.unique(pd.concat([als[0], als[4]])), range(ngen-1, -1, -1)))
    maxp = als[[2, 3, 6, 7]].max().max()
    fig = plt.figure(figsize=[5, 4])
    ax = fig.add_subplot()
    ax.set_ylim([-0.1, ngen - 1 + 0.1])
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)


    if ngen <= 10:
        CHRCOLS = [matplotlib.colors.to_hex(get_cmap('tab10')(i)) for i in range(ngen)]
    else:
        CHRCOLS = [matplotlib.colors.to_hex(get_cmap('gist_rainbow')(int(255/ngen) * i)) for i in range(0, ngen)]
        if ngen % 2 != 0:
            m = CHRCOLS[int((ngen/2) + 1)]
            CHRCOLS = [j for i in range(int(ngen/2)) for j in [CHRCOLS[i]] + [CHRCOLS[int(i + ngen/2)]]] + [m]
        else:
            CHRCOLS = [j for i in range(int(ngen/2)) for j in [CHRCOLS[i]] + [CHRCOLS[int(i + ngen/2)]]]

    for i in range(ngen):
        ax.hlines(ngen - 1 - i, 1, maxp,
                  color=CHRCOLS[i],
                  linewidth=5,
                  zorder=2)
    # print(gdict)
    for al in als.itertuples(index=False):
        # zorder = 0 if abs(gdict[al[0]] - gdict[al[4]]) == 1 else 3
        ax.add_patch(bezierpath(al[2], al[3], al[6], al[7], gdict[al[0]], gdict[al[4]], False, al[8], 1))
    plt.tight_layout(pad=0, h_pad=0, w_pad=0)
    plt.savefig(out, dpi=DPI)
    plt.close()
    return
# END


def getscaf(args):
    import numpy as np
    logger = mylogger("Get Scaffolds")
    fin = args.fasta.name
    n = args.n
    fout = args.o
    prefix = args.p + '_' if args.p is not None else ''
    rev = 2 if args.r else 1
    logger.info('Reading input fasta')
    gen = readfasta(fin)
    chromlen = {chrom: len(seq) for chrom, seq in gen.items()}
    logger.info(f'Selecting {n-len(gen)} breakpoints to generate {n} scaffolds')
    chrid = [np.random.choice(list(chromlen.keys()), 1)[0] for i in range(n-len(gen))]
    cutoff = [np.random.randint(10000, chromlen[i] - 10000, 1)[0] for i in chrid]
    coords = sorted(list(zip(chrid, cutoff)))
    cid = ""
    pos = -1
    logger.info(f'Generating scaffolds')
    scafgenome = {}
    for i in range(len(coords)):
        coord = coords[i]
        if coord[0] != cid:
            if cid in gen:
                if np.random.randint(0, rev, 1):
                    scafgenome[f'{prefix}{cid}_{chromlen[cid]}_r'] = revcomp(gen[cid][pos:chromlen[cid]])
                else:
                    scafgenome[f'{prefix}{cid}_{chromlen[cid]}'] = gen[cid][pos:chromlen[cid]]
            cid = coord[0]
            pos = 0
        if np.random.randint(0, rev, 1):
            scafgenome[f'{prefix}{cid}_{coord[1]}_r'] = revcomp(gen[cid][pos:coord[1]])
        else:
            scafgenome[f'{prefix}{cid}_{coord[1]}'] = gen[cid][pos:coord[1]]
        pos = coord[1]
    if np.random.randint(0, rev, 1):
        scafgenome[f'{prefix}{cid}_{chromlen[cid]}_r'] = revcomp(gen[cid][pos:chromlen[cid]])
    else:
        scafgenome[f'{prefix}{cid}_{chromlen[cid]}'] = gen[cid][pos:chromlen[cid]]
    for key in gen.keys():
        if key not in chrid:
            scafgenome[f'{prefix}{key}'] = gen[key]
    logger.info(f'Writing scaffolds to {fout}')
    writefasta(scafgenome, fout)
    logger.info(f'Finished')
    return
# END


def seqsize(args):
    from gzip import open as gzopen
    fin = args.fasta.name
    f = gzopen(fin, 'r') if isgzip(fin) else open(fin, 'r')
    # out = [(chrom, len(seq)) for chrom, seq in readfasta(fin).items()]
    out = [(chrom, len(s)) for chrom, s in readfasta_iterator(f, isgzip(fin))]
    #TODO: Add output file
    for i in out:
        print(i[0], i[1], sep="\t")
    if args.t:
        print("Genome_length", sum([i[1] for i in out]), sep="\t")
# END


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
# END


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
# END


def getchr(args):
    fin = args.fasta.name
    if args.chrs is None and args.F is not None: chroms = [c.strip().split(' ')[0] for c in open(args.F.name, 'r').readlines()]
    elif args.chrs is not None and args.F is None: chroms = args.chrs
    else: raise Exception('InvalidValue: Incorrect value for chrs provided')
    out = args.o.name if args.o is not None else fin + ".filtered.fasta"
    fa = readfasta(fin)
    chrs = list(fa.keys())
    for c in chrs:
        if (c not in chroms) != args.v:
            fa.pop(c)
    writefasta(fa, out)
# END


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
# END


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
# END


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
# END


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
# END


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
# END


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
# END


def gffsort(args):
    logger.info('Required format of GFF file:\n')
    logger.info('Chr\tMethod\tType\tStart\tEnd\tScore\tStrand\tPhase\tAttribute(ID=uniq_id;other_info)\n')
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
            for gid in sorted(geneids[chr].keys(), key=lambda x: geneids[chr][x]):
                fout.write(gffdata[chr][gid])
# END


def getcol(args):
    """
    Takes an input file and selects columns from it
    """
    logger = mylogger("getcol")
    fin = args.file.name
    fout = args.out.name
    s = '\t' if args.s is None else args.s
    ms = [] if args.m is None else args.m
    cs = [] if args.c is None else args.c
    reorder = args.reorder
    if len(s) > 1:
        logger.error('Only single characters are accepted as column separators')
        sys.exit()
    if len(ms) == 0 and len(cs) == 0:
        logger.error("No columns are selected. Use -m and -c to select columns")
        sys.exit()
    select = []
    first = True
    logger.info(f'Reading from: {fin}')
    with open(fin, 'r') as f:
        logger.info(f'Writing to: {fout}')
        with open(fout, 'w') as out:
            for line in f:
                if line == '':
                    continue
                if line is None:
                    continue
                line = line.strip().split(s)
                if first:
                    for m in ms:
                        for i, v in enumerate(line):
                            if v == m:
                                select.append(i)
                    for c in cs:
                        for i, v in enumerate(line):
                            if c in v:
                                select.append(i)
                    if not reorder:
                        select = sorted(list(set(select)))
                    first = False
                    logger.info(f'Selected columns: {",".join([line[i] for i in select])}')
                out.write(s.join([line[i] for i in select]) + '\n')
    logger.info('Finished')
    return
# END


def bamcov(args):
    import sys
    import os
    from subprocess import Popen, PIPE
    import warnings
    from collections import deque
    import numpy as np
    formatwarning_orig = warnings.formatwarning
    warnings.formatwarning = lambda message, category, filename, lineno, line=None:    formatwarning_orig(message, category, filename, lineno, line='')
    print(args)
    BAM = args.bam.name
    OUT = args.out.name
    BED = None if args.r is None else args.r.name
    # if BED is not None:
    #     sys.exit("SELECTING REGIONS (-R) HAS NOT BEEN IMPLEMENTED YET. PLEASE DO NOT USE IT.")
    D = args.d
    if not os.path.isfile(BAM):
        sys.exit('BAM file is not found: {}'.format(BAM))
    try :
        if not os.path.isfile(BED):
            sys.exit('BED file is not found: '.format(BED))
        else:
            # Update samtools depth command to include the BED file as well
            D += f" -b {BED}"
    except TypeError:
        pass

    # Get chromosome names and lengths
    p = Popen("samtools view -H {}".format(BAM).split(), stdout=PIPE, stderr=PIPE)
    o = p.communicate()
    if o[1] != b'':
        sys.exit("Error in using samtools view to get BAM file header:\n{}".format(o[1].decode()))
    h = o[0].decode().strip().split("\n")
    chrlen = {}
    if not BED:
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
    else:
        chrrange = deque()
        chrom = ''
        with open(BED, 'r') as fin:
            for line in fin:
                line = line.strip().split()
                if line[0] != chrom:
                    if chrom != '':
                        chrrange = mergeRanges(np.array(chrrange))
                        chrlen[chrom] = np.sum(chrrange[:, 1] - chrrange[:, 0] + 1)
                        chrom = line[0]
                        chrrange = deque([[int(line[1]), int(line[2])]])
                    else:
                        chrom = line[0]
                        chrrange = deque([[int(line[1]), int(line[2])]])
                else:
                    chrrange.append([int(line[1]), int(line[2])])
            chrrange = mergeRanges(np.array(chrrange))
            chrlen[chrom] = np.sum(chrrange[:, 1] - chrrange[:, 0] + 1)


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
            fout.write("Genome\t{}\n".format(round(sum(chrrd.values())/sum(chrlen.values()), ndigits=2)))
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
    return
# END


def pbamrc(args):
    # print(args)
    if args.bam is None or args.l is None or args.f is None:
        sys.exit("Require bam, position bed, and refernce fasta files are missing")

    if args.n < 1: sys.exit("Invalid number of cores")
    N = args.n

    # from subprocess import Popen, PIPE
    from multiprocessing import Pool
    from collections import deque
    import pandas as pd
    import numpy as np
    from time import time
    import os

    INDEL = args.I
    bed = pd.read_table(args.l.name, header=None)
    split_N = bed.shape[0] if args.S else N
    splits = np.array_split(range(bed.shape[0]), split_N)
    tmp_df = [bed.iloc[i] for i in splits]
    pre = str(time()).replace(".", "") + str(np.random.randint(1000000000))
    for i in range(split_N):
        tmp_df[i].to_csv(str(pre)+'_'+str(i)+".bed", sep='\t', header=False, index=False)

    command = "bam-readcount" if args.bamrcpath is None else args.bamrcpath.name
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
# END


def bamrc2af(args):
    """
    Reads the output of pbamrc and a corresponding VCF file and returns the allele frequencies of the alt alleles.
    Currently, working for SNPs only
    """
    from gzip import open as gzopen
    logger = mylogger("bamrc2af")
    rcfin = args.bamrc.name
    vcffin = args.vcf.name
    outfin = 'bamrc_af.txt' if args.out is None else args.out.name
    minrc = args.min_rc

    logger.info('Reading VCF')
    posdict = dict()
    op = gzopen if isgzip(vcffin) else open
    with op(vcffin, 'r') as vcf:
        for line in vcf:
            # break
            if line[0] == 35: continue
            line = line.decode()
            line = line.strip().split()
            if line[4].upper() not in 'ACGT' : continue
            posdict[tuple(line[:2])] = line[3], line[4]

    logger.info('Reading bamrc')
    # Get AF from bam readcount
    basedict = dict(zip('ACGT', range(4, 8)))
    with open(rcfin, 'r') as rc, open(outfin, 'w') as out:
        for line in rc:
            line = line.strip().split()
            # if line[3] == '0': continue
            if int(line[3]) < minrc: continue
            try:
                ref, alt = posdict[(line[0], line[1])]
            except KeyError:
                logger.warning(f'Position {line[0]}:{line[1]} not found in VCF. Skipping it.')
                continue
            refi = basedict[ref]
            alti = basedict[alt]
            out.write(f'{line[0]}\t{line[1]}\t{ref}\t{alt}\t{round(int(line[refi])/int(line[3]) , 2)}\t{round(int(line[alti])/int(line[3]), 2)}\n')
    logger.info('Finished')
# END


def run_ppileup(locs, out, bam, pars):
    from subprocess import Popen, PIPE
    with open(out, 'w') as fout:
        for loc in locs:
            # print("/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools mpileup {pars} -r {c}:{s}-{e} {bam}".format(pars=pars, c=loc[0], s=int(loc[1])+1, e=loc[2], bam=bam))
            p = Popen("/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools mpileup {pars} -r {c}:{s}-{e} {bam}".format(pars=pars, c=loc[0], s=int(loc[1])+1, e=loc[2], bam=bam).split(), stdout=fout, stderr=PIPE)
            o = p.communicate()
# END


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
# END


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
                print(TAG + " not found in {qname} at position {rname}:{rpos}. All reads without tag would be skipped.".format(qname=read.qname, rname=read.reference_name, rpos=read.reference_start))
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
# END


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
# END


def sampfa(args):
    """
    Selects random regions from a fasta file and saves them to a new fasta file.

    :param args: An argparse object containing the following attributes:

                - n: number of regions to select

                - s: size of each region to select

                - fa: input fasta file object

                - chr: list of chromosome IDs to select regions from

                - v: if True, selects regions from chromosomes not in the provided list of chromosomes to select from

                - o: output file object (default: "selected_regions.fa")

    :type args: argparse.Namespace
    :raises ValueError: If n or s are less than 1
    :return: None
    """
    n = args.n
    s = int(args.s/2)
    fin = args.fa.name
    inchrs = args.chr
    v = args.v
    fout = args.o.name if args.o is not None else "selected_regions.fa"
    if n < 1:
        raise ValueError("n is too small")
    if args.s < 1:
        raise ValueError("s is too small")
    import random
    from collections import deque
    fasta = readfasta(fin)
    if inchrs is not None:
        if not v:
            select = inchrs
        else:
            select = list(set(fasta.keys()) - set(inchrs))
        fasta = {k: v for k, v in fasta.items() if k in select}
    chrsize = {k: len(v) for k, v in fasta.items()}
    chrs = list(fasta.keys())
    pos = deque()
    for i in range(n):
        c = random.choice(chrs)
        p = random.randint(s, chrsize[c] - s - 1)
        pos.append([c, p])
    # For each position select the region to create fasta_dict
    writefasta({f"{p[0]}_{p[1]-s}_{p[1]+s+1}": fasta[p[0]][(p[1]-s):(p[1]+s+1)] for p in pos}, fout)
    return
# END


def faline(args):
    from collections import deque
    input = args.input.name
    if args.o is not None and args.i is True:
        raise ValueError("Use either -o or -i. Exiting.")
    if args.o is None and args.i is False:
        raise ValueError("Use either -o or -i. Exiting.")
    if args.o == args.input.name:
        raise ValueError("Same filename for input and output file. Use -i instead. Exiting.")
    if args.o is not None:
        out = args.o
    else:
        out = randomstring(10) + ".fa"
    s = args.s
    # TODO: See if this can work with zipped files as well
    seq = deque()
    with open(input, 'r') as fin:
        with open(out, 'w') as fout:
            for line in fin:
                if '>' == line[0]:
                    if len(seq) == 0:
                        fout.write(line)
                    else:
                        seq = ''.join(seq)
                        if s == -1:
                            fout.write(seq + '\n')
                        else:
                            fout.write('\n'.join([seq[i:i+s] for i in range(0, len(seq), s)]) + '\n')
                        fout.write(line)
                        seq = deque()
                else:
                    seq.append(line.strip())
            seq = ''.join(seq)
            if s == -1:
                fout.write(seq + '\n')
            else:
                fout.write('\n'.join([seq[i:i+s] for i in range(0, len(seq), s)]) + '\n')
    if args.i:
        os.rename(out, input)
# END


def shannon_seq(seq):
    """
    For a given sequence, return its Shannon-index
    """
    from collections import deque, Counter
    from math import log
    freq = Counter(seq)
    return round(sum([-1*(freq[b]/len(seq))*(log(freq[b]/len(seq))) for b in {'A', 'C', 'G', 'T'} if b in freq]), 2)
# END


def shannon(args):
    fasta = args.fasta.name
    output = args.output.name
    window = args.w
    size = args.s
    t = args.t
    if size > window:
        raise ValueError("Step size cannot be more than window size. Exiting.")
    from multiprocessing import Pool
    import os
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    import numpy as np
    shannon_values = {}
    for chrom, s in readfasta_iterator(open(fasta, 'r')):
        print(chrom)
        l = len(s)
        seq = [s[i:(i+window)] for i in range(0, l, size)]
        with Pool(processes=t) as pool:
            shannons = pool.map(shannon_seq, seq)
        shannon_values[chrom] = dict(zip(range(0, l, size), shannons))
    with open(output, 'w') as fout:
        for k, v in shannon_values.items():
            for p, sv in v.items():
                fout.write(f"{k}\t{p}\t{p+window}\t{sv}\n")
    return f"Finished calculting shannon index. Output saved in {output}."
# END


def asmreads(args):
    """
    Print reads that constitute the assembly graph at the given position
    """
    try:
        c = args.pos.split(':')[0]
        s, e = map(int, args.pos.split(':')[1].split('-'))
    except Exception as e:
        raise Exception
    gfa = args.gfa.name
    agp = args.agp.name if args.agp is not None else None
    def getcontigfromagp(agp, c, s, e):
        """
        Reads AGP and select contig overlapping the given position
        """
        with open(agp, 'r') as fin:
            for line in fin:
                if line[0] == '#': continue
                line = line.strip().split()
                if line[0] != c: continue
                if line[4] != 'W': continue
                if int(line[1]) <= s and e <= int(line[2]):
                    if line[8] != '-':
                        return line[5], s - int(line[1]) + 1, e - int(line[1]) + 1
                    else:
                        return line[5], int(line[7]) - (e - int(line[1])), int(line[7]) - (s - int(line[1]))
        return None
    # END
    def getreadsfromgfa(gfa, c, s, e):
        from collections import deque
        reads = deque()
        with open(gfa, 'r') as fin:
            cfnd = False
            sfnd = False
            for line in fin:
                if line[0] == 'S': continue
                line = line.strip().split()
                if line[1] != c:
                    if cfnd: break
                    continue
                cfnd = True
                if (s > int(line[2]) + int(line[6])) or (e < int(line[2])):
                    if sfnd: break
                    continue
                sfnd = True
                reads.append(line[4])
        return reads
    # END

    if agp is not None:
        try:
            contig, start, end = getcontigfromagp(agp, c, s, e)
        except TypeError:
            return
    if contig is None:
        raise ValueError("Input genomic coordinates overlap more than 1 contig or overlap a gap. This case cannot be handled currently. Provide smaller input coordinate.")
    for r in getreadsfromgfa(gfa, contig, start, end):
        print(r)
# END


def reg_str_to_list(regstr):
    """
    Converts region string ("contig:start-end") to region tuple. Follow same
    standard as pysam: https://pysam.readthedocs.io/en/latest/glossary.html#term-region

    Here, regions in string format are 1-based closed regions, whereas regions in
    tuple format would be 0-based half open.

    Example: Chr1:15001-20000 would become (Chr1, 15000, 20000).
    Explanation:
        * 15001 becomes 15000 because 1-based get converted to 0-based.
        * 20000 would become 19999 for 0-based and since it is open at the end, next position (i.e. 20000) would be used in the region
    """
    try:
        c, p = regstr.split(':')
        s, e = map(int, p.split('-'))
    except Exception as e:
        raise Exception(e)
    if s < 1 or e < 1:
        raise ValueError('start and end should be more than 0')
    if s > e:
        raise ValueError('start cannot be more than end')
    return [c, s-1, e]
# END


def mapbp(sfin, mapfin, d, posstr, **kwargs):
    """
    Outputs mapping positions for the given reference genome coordinate
    """
    import pysam
    from collections import defaultdict, deque

    def getsyripos(sout, pos):
        """
        :param: sout = syri output file (sorted and tabix indexed) handler (type: pysam.libctabix.TabixFile)
        :param: pos = reference genome coordinate in (chrom, start, end) format
        """
        return [b.split('\t') for b in sout.fetch(*pos)]
    # END
    QUERY_REG = None
    if 'QUERY_REG' in kwargs:
        QUERY_REG = kwargs['QUERY_REG']
    outd = deque()
    pos = reg_str_to_list(posstr)
    pos = pos[:2] + [pos[1] + 1]
    logger.info(f"Getting mapping coordinate for position {pos[0]}:{pos[1]+1}")
    # if syri output is provided, select alignments selected by syri
    if sfin is not None:
        logger.info(f'Syri output is provided. reading file: {sfin}')
        qryreg = defaultdict(deque)
        sout = pysam.TabixFile(sfin)
        # possyri = pos.copy()
        poso = getsyripos(sout, pos)    # Positions overlapping
        if len(poso) == 1:
            if poso[0][6] == 'NOTAL':
                logger.info('No alignment found')
                return
        for p in poso:
            if 'AL' in p[6]:
                if QUERY_REG is not None:
                    if (int(p[4]) == QUERY_REG[0]) and (int(p[5]) == QUERY_REG[1]):
                        qryreg[p[3]].append([int(p[4]), int(p[5])])
                else:
                    qryreg[p[3]].append([int(p[4]), int(p[5])])
    # TODO: Consider adding reader for PAF alignments as well.
    # print(qryreg)
    logger.info(f"Reading BAM file: {mapfin}")
    bam = pysam.AlignmentFile(mapfin)
    for al in bam.fetch(*pos):
        dircg = al.cigartuples #if al.is_forward else al.cigartuples[::-1]
        qs = al.qstart + (dircg[0][1] if dircg[0][0] == 5 else 0) + 1
        qe = qs + al.qlen - 1
        if sfin is not None:
            if al.query_name not in qryreg:
                continue
            if [qs, qe] not in qryreg[al.query_name]:
                continue
        n = pos[1] + 1 - al.reference_start
        p = qs + cgwalk(cgtpl(al.cigarstring, to_int=True), n) - 1
        if not al.is_forward:
            qlen = cggenlen(al.cigarstring, 'q', full=True)
            p = qlen - p + 1
        out = f"{al.query_name}:{p}-{p}"
        if d:
            out = out + '\t+' if al.is_forward else out + '\t-'
        outd.append(out)
    return outd

# END


def mapbp_cli(args):
    sfin = args.anno.name if args.anno is not None else None           # input syri.out file
    mapfin = args.map.name
    d = args.d
    pos = args.pos
    outd = mapbp(sfin, mapfin, d, pos)
    if outd is not None:
        print('\n'.join(outd))
    logger.info('Finished')
    return
# END

def fachrid(args):
    fa = args.fa.name
    out = args.out.name
    names = args.names.name
    iddict = {}
    with open(names, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            iddict[line[0]] = line[1]
    fasta = readfasta(fa)
    chrs = list(fasta.keys())
    outfasta = {}
    for c in chrs:
        if c in iddict:
            outfasta[iddict[c]] = fasta[c]
    writefasta(outfasta, out)
# END


def runsryi(args):
    from subprocess import Popen
    ref = args.ref.name
    qry = args.qry.name
    refi = args.rid.name if args.rid is not None else None
    N = args.n
    prefix = args.p
    altype = args.alignment

    alfile = f'{prefix}.{altype}'
    r = refi if refi is not None else ref
    # align the genomes
    if altype == 'paf':
        command = f'minimap2 -cx asm5 -t {N} --eqx -o {alfile} {r} {qry}'
        proc = Popen(command.split())
        proc.wait()
    else:
        command = f'minimap2 -ax asm5 -t {N} --eqx {r} {qry} | samtools sort -@ {N} -O {altype.upper()} -o {alfile} - '
        # Use shell=True so make the pipe work
        proc = Popen(command, shell=True)
        proc.wait()
        if altype == 'bam':
            proc = Popen(f'samtools index -@ {N} {alfile}'.split())
            proc.wait()
    # Run syri
    proc = Popen(f'syri -c {alfile} -r {ref} -q {qry} -F {altype[0].upper()} --prefix {prefix} --nc {N}'.split())
    proc.wait()
    return
# END


def bam2coords(args):
    '''
    Converts alignment SAM/BAM to alignment coords file
    Strip-down and simplified version. Does not have many checks that are used
    in the convertor used in syri

    The output coords file can then be bgzip and tabix indexed:
    bgzip alignment.coords
    tabix -p bed alignment.coords.gz
    '''

    # def readSAMBAM(fin, type='B'):
    import pysam
    import logging
    import numpy as np
    import pandas as pd

    logger = logging.getLogger('Reading BAM/SAM file')

    fin = args.fin.name
    t = args.t
    try:
        if t == 'B':
            findata = pysam.AlignmentFile(fin, 'rb')
        elif t == 'S':
            return samtocoords(fin)
        else:
            raise ValueError("Wrong parameter")
    except ValueError as e:
        logger.error("Error in opening BAM/SAM file. " + str(e))
        sys.exit()
    except OSError as e:
        logger.error("Error in reading input file." + str(e))
        sys.exit()
    except Exception as e:
        logger.error("Unexpected error in opening BAM/SAM file. " + str(e))
        sys.exit()

    try:
        cgdict = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X'}
        coords = {}
        index = 0
        # findiden = True
        for aln in findata:
            index += 1

            ## Pass non-alinging chromosomes
            if aln.cigarstring is None:
                logger.warning(aln.query_name + ' do not align with any reference chromosome and cannot be analysed')
                continue

            ## Check CIGAR:
            # Commented the check below as that is not required in general
            # if False in [False if i[0] not in [1,2,4,5,7,8] else True for i in aln.cigartuples]:
            #     logger.error("Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: " + str(aln.cigarstring))
            #     sys.exit()
            if len(aln.cigartuples) > 2:
                if True in [True if i[0] in [4, 5] else False for i in aln.cigartuples[1:-1]]:
                    logger.error("Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + aln.cigarstring)
                    sys.exit()

            ## Parse information from the aln object
            astart = aln.reference_start+1
            aend = aln.reference_end
            is_inv = True if np.binary_repr(aln.flag, 12)[7] == '1' else False
            if not is_inv:
                if aln.cigartuples[0][0] in [4, 5]:
                    bstart = aln.cigartuples[0][1]+1
                else:
                    bstart = 1
                bend = bstart + aln.query_alignment_length - 1
            else:
                if aln.cigartuples[-1][0] in [4, 5]:
                    bend = aln.cigartuples[-1][1]+1
                else:
                    bend = 1
                bstart = bend + aln.query_alignment_length - 1
            alen = abs(aend - astart) + 1
            blen = abs(bend - bstart) + 1

            # if findiden:
            #     iden = format((sum([i[1] for i in aln.cigartuples if i[0] == 7])/sum([i[1] for i in aln.cigartuples if i[0] in [1, 2, 7, 8]]))*100, '.2f')
            adir = 1
            bdir = -1 if is_inv else 1
            achr = aln.reference_name
            bchr = aln.query_name
            cg = "".join([str(i[1]) + cgdict[i[0]] for i in aln.cigartuples])
            # coords[index] = [astart, aend, bstart, bend, alen, blen, iden, adir, bdir, achr, bchr, cg]
            coords[index] = [achr, astart, aend, bchr, bstart, bend, alen, blen, adir, bdir, cg]

        ## Return alignments
        coords = pd.DataFrame.from_dict(coords, orient='index')
        coords.sort_values([0, 1, 2, 4, 5, 3], inplace=True, ascending=True)
        coords.index = range(len(coords.index))
        # coords[6] = coords[6].astype('float')
        print(coords.to_csv(index=False, header=False, sep='\t'), end='')
    except Exception as e:
        logger.error("Error in reading BAM/SAM file. " + str(e))
        sys.exit()
# END


def syriidx(args):
    """
    Currently, outputs reference coordinates in the BED format and the query coordinate in the REGION format.
    """
    # syriidx
    from subprocess import Popen, PIPE
    logger = mylogger("syriidx")

    fin = args.syriout.name
    notal = args.notal
    filter = args.f
    annos = ['SYN', 'SYNAL', 'INV', 'TRANS', 'INVTR', 'DUP', 'INVDP']
    if not filter:
        annos += ['CPG', 'CPL', 'DEL', 'DUPAL', 'HDR', 'INS', 'INVAL', 'INVDPAL', 'INVTRAL', 'NOTAL', 'SNP', 'TDM', 'TRANSAL']
    else:
        logger.info("syri annotations will be filtered")
    if notal:
        annos += ['NOTAL']
    annos = set(annos)
    logger.info("reading annotations")
    df = readsyriout(fin, annos)
    df[1] -= 1          # Change coordinate to BED format
    outfin = f'{fin}.bed'
    df.to_csv(outfin, index=False, header=False, sep='\t')
    logger.info("Compressing annotations")
    p = Popen(f"bgzip -f {outfin}".split(), stdout=PIPE, stderr=PIPE)
    o = p.communicate()
    if o[1] != b'':
        sys.exit("Error in bgzip:\n{}".format(o[1].decode()))
    logger.info("Indexing annotations")

    p = Popen(f"tabix -fp bed {outfin}.gz".split(), stdout=PIPE, stderr=PIPE)
    o = p.communicate()
    if o[1] != b'':
        sys.exit("Error in tabix:\n{}".format(o[1].decode()))
    return
# END


def hyellow(s):
    """
    Returns the input string s formatted in yellow color.

    :param s: input string
    :type s: str
    :return: input string s formatted in yellow color
    :rtype: str
    """
    return('\033[93m' + s + '\033[39m')
# END


def revcompseq(args):
    """
    Reverses and complements sequences in an input fasta file and saves the output to a file.

    :param args: An argparse object containing the following attributes:

                - fasta: input fasta file

                - o: output file name (default: 'revcomp.fasta')

                - chr: list of chromosome IDs to reverse complement

    :type args: argparse.Namespace
    :return: None
    """
    import sys
    # Set inputs
    logger = mylogger('revcompseq')
    fasta = args.fasta.name
    if args.o is None:
        out = 'revcomp.fasta'
        logger.warning('Setting output file to revcomp.fasta')
    else:
        out = args.o.name
    # Check for chromosome IDs
    outfasta = dict()
    if args.chr is None:
        logger.error('No sequence ID provided. Exiting.')
        sys.exit()
    else:
        chrids = set(args.chr)
    # Get reversed complemented fasta
    for chrom, s in readfasta_iterator(open(fasta, 'r'), isgzip(fasta)):
        if chrom in chrids:
            logger.info(f'Reverse complementing chromosome {chrom}')
            outfasta[chrom] = revcomp(s)
        else:
            outfasta[chrom] = s
    writefasta(outfasta, out)
    logger.info('Finished')
# END


def syri2bed(args):
    '''
    Take output of syri and converts it into bedpe format
    '''
    input = args.input.name
    output = args.output.name
    logger = mylogger("syri2bedpe")
    f = args.f
    # TODO: Complete this function

    # with open(input, 'r') as fin, open(output, 'w') as fout:
    #     for line in fin:
    #         line = line.strip().split()
    logger.warning('This function is not working yet.')
    return
# END


def splitfa(args):
    """
    Takes an input fasta file and splits its chromosomes/sequences. Each chromosome is saved in a separate file.

    Args:
    args: An argparse object containing the following attributes:
        - fasta: input fasta file
        - prefix: prefix for output file names

    Returns:
    None
    """
    logger = mylogger('splifa')
    fasta = args.fasta.name
    prefix = args.prefix
    logger.info(f"Reading {fasta}")
    for chrom, s in readfasta_iterator(open(fasta, 'r'), isgzip(fasta)):
        outname = f'{prefix}{chrom}.fa'
        logger.info(f'Writing {outname}')
        writefasta({chrom: s}, outname)
    logger.info('Finished')
    return
# END


def cntkmer(args):
    """
    Takes an input fasta file and counts the number of times an input kmer is present.

    Args:
    args: An argparse object containing the following attributes:

        - fasta: input fasta file
        - kmer: kmer to count
        - canonical: Count both the kmer and its reverse complement

    Returns:
    None
    """
    def getloc(seq, k, loc):
        import re
        location = []
        if loc:
            location = [r.start() for r in re.finditer(f'(?=({k}))', seq)]
            count = len(location)
        else:
            count = len(re.findall(f'(?=({k}))', seq))
        return location, count
    # END

    from collections import deque
    import pandas as pd
    # logger = mylogger('cntkmer')
    fasta = args.fasta.name
    kmer = args.kmer.upper()
    cano = args.canonical
    loc = args.loc
    if loc:
        locfin = args.locfin.name if args.locfin.name is not None else f"kmer_{kmer}_locations.bed"
    count = 0
    locations = deque()
    # kmerlower = kmer.lower()
    kmerrc = revcomp(kmer) if cano else ''
    # kmerrclower = revcomp(kmer).lower() if cano else ''
    if cano:
        kmers = [kmer, kmerrc]
    else:
        kmers = [kmer]
    # selectedkmer = [s for s in [kmer, kmerlower, kmerrc, kmerrclower] if s != '']
    logger.info(f'Kmer strings to search: {kmers}')
    logger.info(f"Reading {fasta}")
    lenk = len(kmer)
    for chrom, seq in readfasta_iterator(open(fasta, 'r'), isgzip(fasta)):
        seq = seq.upper()
        for k in kmers:
            location, cnt = getloc(seq, k, loc)
            if loc:
                locations.extend([(chrom, l, l+lenk, k) for l in location])
            count += cnt
    if loc:
        if len(locations) > 0:
            locations = pd.DataFrame(locations)
            locations.sort_values([0, 1], inplace=True)
            locations.to_csv(locfin, sep='\t', index=False, header=False)
    print(f'Number of occurence of K-mer {kmer}: {count}')
    logger.info('Finished')
    logger.propagate = False
    return count
# END


def xls2csv(args):
    """
    Converts one or more sheets from an excel file into one txt file
    """
    import pandas as pd
    xls = args.xls.name
    fout = args.output.name
    sheets = 0 if args.s is None else None if args.s == ['all'] else args.s
    addname = args.n
    logger.info(f"Reading {xls}.")
    data = pd.read_excel(xls, sheet_name=sheets)
    if isinstance(data, pd.core.frame.DataFrame):
        outdict = data
    elif isinstance(data, dict):
        outdict = pd.DataFrame()
        for k, v in data.items():
            df = v.copy()
            if addname:
                df['sheet'] = k
            outdict = pd.concat([outdict, df])
    logger.info(f"Writing to {fout}.")
    outdict.to_csv(fout, index=False, header=True, sep='\t' if fout[-4:] == ".tsv" else ',')
    logger.info('Finished')
    return
# END


def gz2bgz(args):
    import numpy as np
    from gzip import open as gzopen
    from Bio import bgzf
    import os
    from time import time
    logger = mylogger('gz2bgz')
    fin = args.fin.name
    out = args.out.name if args.out is not None else None
    try:
        assert isgzip(fin)
    except AssertionError:
        logger.error('Input file is not gzip compressed. Exiting')
        sys.exit()
    if out is None:
        usetmp = True
        fout = f'{str(time()).replace(".", "")}{np.random.randint(1000000000)}.gz'
        logger.info(f'Saving to temporary file {fout}')
    else:
        usetmp = False
        fout = out
        logger.info(f'Saving to output file {fout}')
    with gzopen(fin, 'rb') as f, bgzf.open(fout, 'wb') as fo:
        for line in f:
            fo.write(line)
    if usetmp:
        logger.info(f'Ranaming temporary file {fout} to input file {fin}')
        os.rename(fout, fin)
    logger.info('Finished')
    return
# END


def topmsa(args):
    from hometools.topological_msa import topological_sort
    logger = mylogger('topmsa')
    topological_sort(args)
    logger.info('Finished')
# END


def main(cmd):
    parser = argparse.ArgumentParser("Collections of command-line functions to perform common pre-processing and analysis functions.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()

    # <editor-fold desc="FASTA Commands">
    parser_getchr = subparsers.add_parser("getchr", help=hyellow("FASTA: Get specific chromosomes from the fasta file"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sampfa  = subparsers.add_parser("sampfa", help=hyellow("FASTA: Sample random sequences from a fasta file"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_exseq = subparsers.add_parser("exseq", help= hyellow("FASTA: extract sequence from fasta"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_getscaf = subparsers.add_parser("getscaf", help=hyellow("FASTA: generate scaffolds from a given genome"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_seqsize = subparsers.add_parser("seqsize", help=hyellow("FASTA: get size of dna sequences in a fasta file"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_filsize = subparsers.add_parser("filsize", help=hyellow("FASTA: filter out smaller molecules"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_subnuc = subparsers.add_parser("subnuc", help=hyellow("FASTA: Change character (in all sequences) in the fasta file"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_basrat = subparsers.add_parser("basrat", help=hyellow("FASTA: Calculate the ratio of every base in the genome"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_genome_ranges = subparsers.add_parser("genome_ranges", help=hyellow("FASTA: Get a list of genomic ranges of a given size"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_get_homopoly = subparsers.add_parser("get_homopoly", help=hyellow("FASTA: Find homopolymeric regions in a given fasta file"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_n50 = subparsers.add_parser("asstat", help=hyellow("FASTA: Get N50 values for the given list of chromosomes"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_shannon = subparsers.add_parser("shannon", help=hyellow("FASTA: Get Shanon entropy across the length of the chromosomes using sliding windows"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_fachrid = subparsers.add_parser("fachrid", help=hyellow("FASTA: Change chromosome IDs"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_faline = subparsers.add_parser("faline", help=hyellow("FASTA: Convert fasta file from single line to multi line or vice-versa"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_revcompseq = subparsers.add_parser("revcompseq", help=hyellow("FASTA: Reverse complement specific chromosomes in a fasta file"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_splitfa = subparsers.add_parser("splitfa", help=hyellow("FASTA: Split fasta files to individual sequences. Each sequence is saved in a separate file."), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_cntkmer = subparsers.add_parser("countkmer", help=hyellow("FASTA: Count the number of occurence of a given kmer in a fasta file."), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # </editor-fold>

    # <editor-fold desc="BAM Commands">
    parser_bamcov = subparsers.add_parser("bamcov", help="BAM: Get mean read-depth for chromosomes from a BAM file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_pbamrc = subparsers.add_parser("pbamrc", help="BAM: Run bam-readcount in a parallel manner by dividing the input bed file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_bamrc2af = subparsers.add_parser("bamrc2af", help="BAM: Reads the output of pbamrc and a corresponding VCF file and saves the allele frequencies of the ref/alt alleles.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_splitbam = subparsers.add_parser("splitbam", help="BAM: Split a BAM files based on TAG value. BAM file must be sorted using the TAG.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_mapbp = subparsers.add_parser("mapbp", help="BAM: For a given reference coordinate get the corresponding base and position in the reads/segments mapping the reference position", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_bam2coords = subparsers.add_parser("bam2coords", help="BAM: Convert BAM/SAM file to alignment coords", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_ppileup = subparsers.add_parser("ppileup", help="BAM: Currently it is slower than just running mpileup on 1 CPU. Might be possible to optimize later. Run samtools mpileup in parallel when pileup is required for specific positions by dividing the input bed file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # </editor-fold>

    # <editor-fold desc="syri CLI">
    parser_runsyri = subparsers.add_parser("runsyri", help=hyellow("syri: Parser to align and run syri on two genomes"),
                                           formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_syriidx = subparsers.add_parser("syriidx", help=hyellow(
        "syri: Generates index for syri.out. Filters non-SR annotations, then bgzip, then tabix index"),
                                           formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_syri2bed = subparsers.add_parser("syri2bed", help=hyellow("syri: Converts syri output to bedpe format"),
                                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # </editor-fold>

    # <editor-fold desc="Plotting">
    parser_plthist = subparsers.add_parser("plthist", help="Plot: Takes frequency output (like from uniq -c) and generates a histogram plot", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_plotal = subparsers.add_parser("plotal", help="Plot: Visualise pairwise-whole genome alignments between multiple genomes", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_plotbar = subparsers.add_parser("pltbar", help="Plot: Generate barplot. Input: a two column file with first column as features and second column as values", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # </editor-fold>

    # <editor-fold desc='Assembly graphs'>
    parser_asmreads = subparsers.add_parser("asmreads", help=hyellow("GFA: For a given genomic region, get reads that constitute the corresponding assembly graph"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_gfatofa = subparsers.add_parser("gfatofa", help=hyellow("GFA: Convert a gfa file to a fasta file"), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # </editor-fold>

    # <editor-fold desc='GFF'>
    parser_gfftrans = subparsers.add_parser("gfftrans", help="GFF: Get transcriptome (gene sequence) for all genes in a gff file. WARNING: THIS FUNCTION MIGHT HAVE BUGS.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_gffsort = subparsers.add_parser("gffsort", help="GFF: Sort a GFF file based on the gene start positions", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # </editor-fold>

    # <editor-fold desc='VCF'>
    parser_vcfdp = subparsers.add_parser("vcfdp", help=hyellow("VCF: Get DP and DP4 values from a VCF file."), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # </editor-fold>

    # <editor-fold desc='Tables'>
    parser_getcol = subparsers.add_parser("getcol", help="Table:Select columns from a TSV or CSV file using column names", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # TODO: Add functionality for sampling rows
    parser_samplerow = subparsers.add_parser("smprow", help="Table:Select random rows from a text file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_xls2tsv = subparsers.add_parser("xls2csv", help="Table:Convert excel tables to .tsv/.csv", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # </editor-fold>

    # <editor-fold desc='Misc'>
    parser_gz2bgz = subparsers.add_parser("gz2bgz",
                                          help="Misc:converts a gz to bgzip",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_topmsa = subparsers.add_parser("topmsa",
                                          help="Misc:Create topological MSA for non-duplicated nodes (genomic regions).",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # </editor-fold>

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()

    # topmsa
    parser_topmsa.set_defaults(func=topmsa)
    parser_topmsa.add_argument('input', help="Input file containing sequences to align. One sequence per line. Sequence values are separated by commas", type=argparse.FileType('r'))
    parser_topmsa.add_argument('output', help="Output file containing MSA.", type=argparse.FileType('w'))

    # gz2bgz
    parser_gz2bgz.set_defaults(func=gz2bgz)
    parser_gz2bgz.add_argument('fin', help="Input gzip file", type=argparse.FileType('r'))
    parser_gz2bgz.add_argument('--out', help='Output file name', type=argparse.FileType('w'))

    # bamrc2af
    parser_bamrc2af.set_defaults(func=bamrc2af)
    parser_bamrc2af.add_argument("bamrc", help="BAM readcount file generated using bamrc", type=argparse.FileType('r'))
    parser_bamrc2af.add_argument("vcf", help="VCF file", type=argparse.FileType('r'))
    parser_bamrc2af.add_argument("out", help="Output file", type=argparse.FileType('w'))
    parser_bamrc2af.add_argument("--min_rc", help="Minimum required read count. Position with lower number of reads would be filtered out", type=int, default=1)

    # xls2csv
    parser_xls2tsv.set_defaults(func=xls2csv)
    parser_xls2tsv.add_argument("xls", help="Input excel file", type=argparse.FileType('r'))
    parser_xls2tsv.add_argument("output", help="Output file name", type=argparse.FileType('w'))
    parser_xls2tsv.add_argument("-s", help="Sheet names or 0-based index for selecting multiple sheets. Default selects first sheet (index 0). 'all' selects all sheets.", type=str, nargs='+')
    parser_xls2tsv.add_argument("-n", help="Add sheet name as extra column. Used when multiple sheets are selected.", default=False, action='store_true')


    # cntkmer
    parser_cntkmer.set_defaults(func=cntkmer)
    parser_cntkmer.add_argument("fasta", help="Input fasta file", type=argparse.FileType('r'))
    parser_cntkmer.add_argument("kmer", help="Kmer to find in fasta", type=str)
    parser_cntkmer.add_argument("--loc", help="Output locations of kmer to BED format", default=False, action='store_true')
    parser_cntkmer.add_argument("--locfin", help="BED file name", type=argparse.FileType('w'))
    parser_cntkmer.add_argument("--canonical", help="Count both the kmer and its reverse complement", default=False, action='store_true')


    # splitfa
    parser_splitfa.set_defaults(func=splitfa)
    parser_splitfa.add_argument("fasta", help="Input fasta file", type=argparse.FileType('r'))
    parser_splitfa.add_argument("--prefix", help="Prefix to add before file names", type=str, default='')

    # revcompseq
    parser_revcompseq.set_defaults(func=revcompseq)
    parser_revcompseq.add_argument("fasta", help="Input fasta file", type=argparse.FileType('r'))
    parser_revcompseq.add_argument("--chr", help="Sequence ID to reverse complement", type=str, action='append')
    parser_revcompseq.add_argument("-o", help="Output file name", type=argparse.FileType('w'))


    # plotbar
    parser_plotbar.set_defaults(func=pltbar)
    parser_plotbar.add_argument("-f", help="Input file containing the frequency data", type=str, default='STDIN')
    parser_plotbar.add_argument("-o", help="Output file name", type=argparse.FileType('w'), default='bar.pdf')
    parser_plotbar.add_argument("-W", help="width of the plot (in inches)", type=float, default=4)
    parser_plotbar.add_argument("-H", help="height of the plot (in inches)", type=float, default=4)
    parser_plotbar.add_argument("-x", help="X-axis label", type=str, default='feature', nargs='+')
    parser_plotbar.add_argument("-y", help="Y-axis label", type=str, default='value', nargs='+')
    parser_plotbar_sort = parser_plotbar.add_mutually_exclusive_group()
    parser_plotbar_sort.add_argument("--sx", help="Sort features", default=False, action='store_true')
    parser_plotbar_sort.add_argument("--sy", help="Sort values", default=False, action='store_true')
    parser_plotbar.add_argument("--xltilt", help="X-axis label tilt (range: 0-90)", default=0, type=int)
    parser_plotbar.add_argument("--ylog", help="Y-axis on log scale", default=False, action='store_true')
    # parser_plotbar.add_argument("-xlim", help="Set X-axis limit for numerical data", nargs=2, type=int)
    parser_plotbar.add_argument("--ylim", help="Set Y-axis limit for numerical data", nargs=2, type=int)
    parser_plotbar.add_argument("-t", help="title of the plot", type=str, default=None)


    # bam2coords
    parser_bam2coords.set_defaults(func=bam2coords)
    parser_bam2coords.add_argument("fin", help='Input BAM/SAM file', type=argparse.FileType('r'))
    parser_bam2coords.add_argument("t", help='File type (B: BAM, S: SAM)', type=str, choices=['B', 'S'])


    # runsyri
    parser_runsyri.set_defaults(func=runsryi)
    parser_runsyri.add_argument("ref", help='Reference genome', type=argparse.FileType('r'))
    parser_runsyri.add_argument("qry", help='Query genome', type=argparse.FileType('r'))
    parser_runsyri.add_argument("--rid", help='Reference index generated using minima2', type=argparse.FileType('r'))
    parser_runsyri.add_argument("-n", help='Number of cores to use', type=int, default=1)
    parser_runsyri.add_argument("-p", help='prefix', type=str, default='out')
    parser_runsyri.add_argument("-alignment", help='Output alignment type', choices=['sam', 'bam', 'paf'], default='paf', type=str)

    # syriidx
    parser_syriidx.set_defaults(func=syriidx)
    parser_syriidx.add_argument("syriout", help='syri output file', type=argparse.FileType('r'))
    parser_syriidx.add_argument("-f", help='Only output syntenic and SR regions.', default=False, action='store_true')
    parser_syriidx.add_argument("--notal", help='Also include reference NOTAL regions', default=False, action='store_true')

    # syri2bed
    parser_syri2bed.set_defaults(func=syri2bed)


    # fachrid
    parser_fachrid.set_defaults(func=fachrid)
    # TODO : Clarify input alignment file format
    parser_fachrid.add_argument("fa", help='Input fasta file', type=argparse.FileType('r'))
    parser_fachrid.add_argument("out", help='Output fasta file', type=argparse.FileType('w'))
    parser_fachrid.add_argument("names", help='Table listing old and new names', type=argparse.FileType('r'))


    # plotal
    parser_plotal.set_defaults(func=plotal)
    # TODO : Clarify input alignment file format
    parser_plotal.add_argument("align", help='Input alignment file', type=argparse.FileType('r'))
    parser_plotal.add_argument("out", help='Output file name', type=argparse.FileType('w'))
    parser_plotal.add_argument("-D", help='DPI', type=int, default=300)


    # mapbp
    parser_mapbp.set_defaults(func=mapbp_cli)
    parser_mapbp.add_argument("pos", help='Genome position in \'chr:start-end\' format.', type=str)
    # parser_mapbp.add_argument("map", help='Alignment file in BAM/PAF format.', type=argparse.FileType('r'))
    parser_mapbp.add_argument("map", help='Alignment file in BAM format.', type=argparse.FileType('r'))
    parser_mapbp.add_argument("--anno", help='Syri annotation file. Only alignments present in the syri output would be selected. Need syri.out to be sorted and indexed with tabix. Use: hometools syridx.', type=argparse.FileType('r'))
    parser_mapbp.add_argument("-d", help='Output alignment strand (sense(+)/antisense(-))', default=False, action='store_true')


    # asmreads
    parser_asmreads.set_defaults(func=asmreads)
    parser_asmreads.add_argument("pos", help='Genome position in \'chr:start-end\' format.', type=str)
    parser_asmreads.add_argument("gfa", help='GFA file consisting of contigs and reads (https://github.com/GFA-spec/GFA-spec/blob/master/GFA2.md)', type=argparse.FileType('r'))
    parser_asmreads.add_argument("--agp", help='For scaffolded assembly, AGP file containing information on contig order and orientation (https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/).', type=argparse.FileType('r'))

    # shanon
    parser_shannon.set_defaults(func=shannon)
    parser_shannon.add_argument("fasta", help='Input fasta file', type=argparse.FileType('r'))
    parser_shannon.add_argument("output", help='Output file', type=argparse.FileType('w'))
    parser_shannon.add_argument("-w", help='Window size', type=int, default=10000)
    parser_shannon.add_argument("-s", help='Step size', type=int, default=1000)
    # parser_shannon.add_argument("-p", help='Use canonical probabilities (0.25) for the occurence of a base. Default is to calculate probability from the sequence substring.', action='strore_true')
    parser_shannon.add_argument("-t", help='Number of threads to use.', type=int, default=1)

    # faline
    parser_faline.set_defaults(func=faline)
    parser_faline.add_argument("input", help='Input fasta file', type=argparse.FileType('r'))
    parser_faline.add_argument("-o", help='Output fasta file', type=str)
    parser_faline.add_argument("-i", help='Change input fasta file inplace', action="store_true")
    parser_faline.add_argument("-s", help='Length of the sequence line. Use -1 for single line fasta.', type=int, default=60)

    # sampfa
    parser_sampfa.set_defaults(func=sampfa)
    parser_sampfa.add_argument("fa", help='Input fasta file', type=argparse.FileType('r'))
    parser_sampfa.add_argument("-n", help='Number of regions to select', type=int, default=1)
    parser_sampfa.add_argument("-s", help='Size of region to select', type=int, default=100)
    parser_sampfa.add_argument("--chr", help='Chromosome from which the regions should be selected. Multiple chromosomes can be selected.', type=str, action='append')
    parser_sampfa.add_argument("-v", help='Invert chromosome selection', default=False, action='store_true')
    parser_sampfa.add_argument("-o", help='Output file name', type=argparse.FileType('w'))

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
    parser_pbamrc.add_argument("--bamrcpath", help="Location of bam-readcount executable", type=argparse.FileType('r'))


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
    parser_getcol.add_argument("--reorder", help="Change output column order based on the order of the selected columns in -m and -c", default=False, action='store_true')



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
    parser_plthist.add_argument("--xltilt", help="X-axis label tilt (range: 0-90)", default=0, type=int)
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

    parser_exseq.set_defaults(func=extractseq)
    parser_exseq.add_argument("fasta", help="fasta file", type=argparse.FileType('r'))
    parser_exseq.add_argument("-c", "--chr", help="Chromosome ID", type=str)
    parser_exseq.add_argument("-s", "--start", help="Start location (BED format, 0-base, start included)", type=int)
    parser_exseq.add_argument("-e", "--end", help="End location (BED format, 1-base, end excluded)", type=int)
    parser_exseq.add_argument("--fin", help="File containing locations to extract (BED format)", type=argparse.FileType('r'))
    parser_exseq.add_argument("-o", help="Output file name", type=argparse.FileType('w'), default=None)
    parser_exseq.add_argument("-r", help="Reverse complement the sequence", default=False, action='store_true')

    parser_getscaf.set_defaults(func=getscaf)
    parser_getscaf.add_argument("fasta", help="genome fasta file", type=argparse.FileType('r'))
    parser_getscaf.add_argument("n", help="number of scaffolds required", type=int)
    parser_getscaf.add_argument("-r", help="reverse complement some scaffolds", default=False, action="store_true")
    parser_getscaf.add_argument("-o", help="output file name", default="scaf.fasta")
    parser_getscaf.add_argument("-p", help="Chromosome ID prefix", type=str)

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
    # print(args)
    from itertools import cycle
    from sys import getsizeof, stderr
    from itertools import chain
    from collections import deque
    args.func(args)

