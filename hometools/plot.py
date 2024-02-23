"""
This module contains utility functions for working with plots
"""


def getcolors(cp, count):
    """
    Output `count` colours form the `cp`.
    """
    return [cp(i/count) for i in range(count)]
# END


def plotref(refidx, varpos, syrireg, figout, bw=100000):
    """
    Plots distribution of bed file on a genome annotated with SRs
    :param refidx: path to ref.faidx
    :param syrireg: path to syri regions in BEDPE format
    :param varpos: variant positions in BED
    :param figout: output figure path
    :return:
    """
    from matplotlib import pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.collections import PatchCollection
    import pybedtools as bt     # TODO: remove this dependency
    from collections import deque, defaultdict
    import numpy as np
    import logging
    from hometools.hometools import mergeRanges, readfaidxbed
    import sys
    logger = logging.getLogger('plotref')

    def _readbed(vp, refbed):
        _chrs = set([r[0] for r in refbed])
        bincnt = defaultdict(deque)
        skipchrs = []
        curchr = ''
        pos = deque()
        added_chrs = list()
        with open(vp, 'r') as fin:
            for line in fin:
                line = line.strip().split()
                if len(line) < 3:
                    logger.warning("Incomplete information in BED file at line: {}. Skipping it.".format("\t".join(line)))
                    continue
                if line[0] not in _chrs:
                    if line[0] not in skipchrs:
                        logger.info("Chromosome in BED is not present in FASTA or not selected for plotting. Skipping it. BED line: {}".format("\t".join(line)))
                        skipchrs.append(line[0])
                    continue
                if curchr == '':
                    curchr = line[0]
                    pos.append([int(line[1]), int(line[2])])
                elif curchr == line[0]:
                    pos.append([int(line[1]), int(line[2])])
                else:
                    if line[0] in added_chrs:
                        logger.error("BED file: {} is not sorted. For plotting tracks, sorted bed file is required. Exiting.".format(vp))
                        sys.exit()
                    if len(pos) > 1:
                        rngs = mergeRanges(np.array(pos))
                    else:
                        rngs = pos
                    chrpos = np.array(list(set([i for r in rngs for i in range(r[0], r[1])])))
                    # Get bin breakpoints for the chromosome
                    bins = np.concatenate((np.arange(0, [r[2] for r in refbed if r[0] == curchr][0], bw), np.array([r[2] for r in refbed if r[0] == curchr])), axis=0)
                    binval = np.histogram(chrpos, bins)[0]
                    bincnt[curchr] = deque([((bins[i] + bins[i+1])/2, binval[i]/bw) for i in range(len(binval))])
                    added_chrs.append(curchr)
                    # Set the new chromosome
                    curchr = line[0]
                    pos = deque([[int(line[1]), int(line[2])]])
            if len(pos) > 1:
                rngs = mergeRanges(np.array(pos))
            else:
                rngs = pos
            chrpos = np.array(list(set([i for r in rngs for i in range(r[0], r[1])])))
            # Get bin breakpoints for the chromosome
            bins = np.concatenate((np.arange(0, [r[2] for r in refbed if r[0] == curchr][0], bw), np.array([r[2] for r in refbed if r[0] == curchr])), axis=0)
            binval = np.histogram(chrpos, bins)[0]
            bincnt[curchr] = deque([((bins[i] + bins[i+1])/2, binval[i]/bw) for i in range(len(binval))])
        return bincnt
    # END
    chr_height = 0.3
    th = 0.6    # Track height for plotting bedfile data
    refbed = readfaidxbed(refidx)
    # TODO: remove dependency on pybedtools
    # syrioutbed = bt.BedTool('/srv/netscratch/dep_mercier/grp_schneeberger/projects/read_vs_assembly/results/human/hg002/reads15kb/svcalls/syri/hg002.reads15kb.winnowmap.hap1syri.bedpe').sort()
    syrioutbed = bt.BedTool(syrireg).sort()
    bedbin = _readbed(varpos, refbed)
    fig = plt.figure(figsize=[12, 10])
    ax = fig.add_subplot()
    ax.set_ylim([0, len(refbed)+1])
    ax.set_xlim([0, max([r[2] for r in refbed])])
    colors = {'SYN': 'lightgrey', 'INV': '#FFA500', 'TRANS': '#9ACD32', 'DUP': '#00BBFF'}
    peakpos = defaultdict()
    for i in range(len(refbed)):
        patches = deque()
        r = refbed[i]
        y = len(refbed) - i
        ax.add_patch(Rectangle((r[1], y), r[2], chr_height, fill=False, linewidth=0.5))
        # ax.add_patch(Rectangle((r[1], y), r[2], chr_height, linewidth=0.5, color='black'))
        bed = syrioutbed.filter(lambda b: b.chrom == r[0]).saveas()
        for b in bed:
            patches.append(Rectangle((b.start, y), b.stop-b.start, chr_height, color=colors[b[6]], linewidth=0))
        ax.add_collection(PatchCollection(patches, match_original=True))
        chrpos = [k[0] for k in bedbin[r[0]]]
        tpos = [k[1] for k in bedbin[r[0]]]
        tposmax = max(tpos)
        y0 = len(refbed) - i + chr_height
        ypos = [(t*th/tposmax)+y0 for t in tpos]
        ax.fill_between(chrpos, ypos, y0, color='blue', lw=0.1, zorder=2)
        y_mid = y0 + (th/2)
        peakpos[refbed[i][0]] = [(chrpos[j], tpos[j]) for j in range(len(ypos)) if ypos[j] > y_mid]
    plt.tight_layout()
    fig.savefig(figout)
    plt.close()
    return peakpos
# END


def bezierpath(rs, re, qs, qe, ry, qy, v, col, alpha, label='', lw=0, zorder=0):
    import matplotlib.patches as patches
    from matplotlib.path import Path
    smid = (qs-rs)/2    # Start increment
    emid = (qe-re)/2    # End increment
    hmid = (qy-ry)/2    # Height increment
    if not v:
        verts = [(rs, ry),
                 (rs, ry+hmid),
                 (rs+2*smid, ry+hmid),
                 (rs+2*smid, ry+2*hmid),
                 (qe, qy),
                 (qe, qy-hmid),
                 (qe-2*emid, qy-hmid),
                 (qe-2*emid, qy-2*hmid),
                 (rs, ry),
                 ]
    else:
        verts = [(ry, rs),
                 (ry+hmid, rs),
                 (ry+hmid, rs+2*smid),
                 (ry+2*hmid, rs+2*smid),
                 (qy, qe),
                 (qy-hmid, qe),
                 (qy-hmid, qe-2*emid),
                 (qy-2*hmid, qe-2*emid),
                 (ry, rs),
                 ]
    codes = [
        Path.MOVETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.LINETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor=col, lw=lw, alpha=alpha, label=label, edgecolor=col, zorder=zorder)
    return patch
# END


def bezierpath2(verts, col, alpha, label='', lw=0, zorder=0):
    """
    General purpose function for drawing bezier curve along the given points
    """
    import matplotlib.patches as patches
    from matplotlib.path import Path
    verts = verts.tolist() + [verts[0].tolist()]
    codes = [Path.MOVETO] + [Path.CURVE4]*(len(verts)-2) + [Path.CLOSEPOLY]
    print(verts, codes)
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor=col, lw=lw, alpha=alpha, label=label, edgecolor=col, zorder=zorder)
    return patch
# END


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
# END


def cleanax(ax):
    """
    Remove the top and right spine and set plot to have tight layout
    """
    from matplotlib import pyplot as plt
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # ax.legend(bbox_to_anchor=(1.01, 1))
    plt.tight_layout(pad=0.1)
    return ax
# END


def density_scatter(x, y, ax=None, fig=None, sort=True, bins=20, **kwargs):
    """
    This ta
    """
    # TODO: See if this can be converted to a command-line API
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


def plot2emf(fin, fout):
    """
    Converts a PDF/SVG file to emf using inkscape. Other input file format may
    also work depending on inkscapes capabilities.
    """
    from subprocess import call
    return call(f'inkscape {fin} -M {fout}', shell=True)
# END


def inlineplotsr(ax, args):
    """
    Method for plotting plotsr plots on a given axis. Takes axis and plotsr arguments
    as keywords
    Takes diction
    """
    import logging
    from pandas import concat as pdconcat
    from pandas import unique
    from plotsr.func import setlogconfig, readbasecfg, readsyriout, readbedout, filterinput, validalign2fasta, selectchrom, selectregion, createribbon, drawax, pltchrom, pltsv, drawmarkers, readtrack, drawtracks
    from collections import deque, OrderedDict
    import os
    from math import ceil
    import matplotlib
    from matplotlib import pyplot as plt
    import sys

    ## Define loggers
    setlogconfig(args.log)
    # filehandler = getfilehandler(args.logfin.name, args.log)
    # global getlogger
    # getlogger = definelogger(filehandler)
    logger = logging.getLogger("Plotsr")

    ###################################################################
    # Check python and pandas version. Add other checks (if required!!)
    ###################################################################
    # logger.debug('checking arguments')
    # try:
    #     assert sys.version_info.major == 3
    #     assert sys.version_info.minor >= 8
    # except AssertionError:
    #     logger.warning('\nPlotsr is tested for Python >=3.8. Currently using Python {}.{}. This may result in errors.'.format(sys.version_info.major, sys.version_info.minor))
    # except KeyboardInterrupt:
    #     raise()
    # except Exception as E:
    #     sys.exit(E)

    ## Validate input
    if args.sr is None and args.bp is None:
        logger.error("No structural annotations provided. Use --sr or -bp to provide path to input files")
        sys.exit()

    if args.sr is not None and args.bp is not None:
        logger.error("Both --sr and --bp cannot be used. Use single file type for all input structural annotations files. User converter to reformat BEDPE/syri.out files")
        sys.exit()

    # Check if both --chr and --reg are defined
    if args.chr is not None and args.reg is not None:
        logger.error("Both --chr and --reg are provided. Only one parameter can be provided at a time. Exiting.")
        sys.exit()

    # Check if both --chr and --chrord are defined
    if args.chr is not None and args.chrord is not None:
        logger.error("Both --chr and --chrord are provided. Only one parameter can be provided at a time. Exiting.")
        sys.exit()

    # Check if --rtr is used without --reg
    if args.rtr and args.reg is None:
        logger.error("Cannot use --rtr without --reg. Exiting.")
        sys.exit()


    ###################################################################
    # Declare variable using argument values
    ###################################################################

    # Set Figure height and width. Change later based on chromosome number and size
    logger.info('Starting')
    FS = args.f             # Font size
    H = args.H              # Height
    W = args.W              # Width
    O = args.o              # Output file name
    D = args.d              # Output file DPI
    R = args.R              # Create ribbons
    V = args.v              # Vertical chromosomes
    S = args.S              # Space between homologous chromosomes
    B = None if args.markers is None else args.markers.name              # Annotation bed file
    TRACKS = None if args.tracks is None else args.tracks.name
    REG = None if args.reg is None else args.reg.strip().split(":")
    RTR = args.rtr
    CHRS = args.chr
    ITX = args.itx
    CHRNAME= args.chrname.name if args.chrname is not None else None
    print(plt.rcParams['font.size'])
    plt.rcParams['font.family'] = 'Arial'
    # AC = args.aligncolour

    # print(ITX)
    # sys.exit()

    ## Get config
    cfg = readbasecfg('', V) if args.cfg is None else readbasecfg(args.cfg.name, V)
    if S < 0.1 or S > 0.75:
        logger.warning('Value for S outside of normal range 0.1-0.75.')

    ## Check output file extension
    if len(O.split('.')) == 1:
        logger.warning("Output filename has no extension. Plot would be saved as a pdf")
        O = O + ".pdf"
    elif O.split('.')[-1] not in ['pdf', 'png', 'svg']:
        logger.warning("Output file extension is not in {'pdf','png', 'svg'}. Plot would be saved as a pdf")
        O = O.rsplit(".", 1)[0] + ".pdf"

    ## Set matplotlib backend
    # try :
    #     matplotlib.use(args.b)
    #     # matplotlib.use('Qt5Agg')    # TODO: Delete this line
    # except :
    #     sys.exit('Matplotlib backend cannot be selected')

    # fins = ['col_lersyri.out', 'ler_cvisyri.out', 'cvi_erisyri.out', 'eri_shasyri.out', 'sha_kyosyri.out', 'kyo_an1syri.out', 'an1_c24syri.out'] #TODO: Delete this line
    # Read alignment coords
    alignments = deque()
    chrids = deque()
    if args.sr is not None:
        for f in args.sr:
            fin = f.name
            # for fin in fins: #TODO: Delete this line
            al, cid = readsyriout(fin)
            alignments.append([os.path.basename(fin), al])
            chrids.append((os.path.basename(fin), cid))
    elif args.bp is not None:
        for f in args.bp:
            fin = f.name
            al, cid = readbedout(fin)
            alignments.append([os.path.basename(fin), al])
            chrids.append((os.path.basename(fin), cid))

    # Get groups of homologous chromosomes. Use the order from the user if provided.
    cs = set(unique(alignments[0][1]['achr']))
    if args.chrord is None:
        chrs = [k for k in chrids[0][1].keys() if k in alignments[0][1]['achr'].unique()]
    else:
        chrs = deque()
        with open(args.chrord.name, 'r') as fin:
            for line in fin:
                c = line.strip()
                if c not in cs:
                    logger.error("Chromosome {} in {} is not a chromosome in alignment file {}. Exiting.".format(c, args.chrord.name, alignments[0][0]))
                    sys.exit()
                chrs.append(c)
        chrs = list(chrs)
        # Check that the chrorder file contains all chromosomes
        if len(chrs) != len(cs):
            logger.error("Number of chromsomes in {} is less than the number of chromsomes in the alignment file {}. Either list the order of all chromosomes or use --chr if chromosome selection is requires. Exiting.".format(args.chrord.name, alignments[0][0]))
            sys.exit()

    chrgrps = OrderedDict()
    for c in chrs:
        cg = deque([c])
        cur = c
        for i in range(len(chrids)):
            n = chrids[i][1][cur]
            cg.append(n)
            cur = n
        chrgrps[c] = cg

    # Filter alignments to select long alignments between homologous chromosomes
    for i in range(len(alignments)):
        alignments[i][1] = filterinput(args, alignments[i][1], chrids[i][1], ITX)

        # Check chromsome IDs and sizes
    chrlengths, genomes = validalign2fasta(alignments, args.genomes.name)
    # chrlengths, genomes = validalign2fasta(alignments, 'genomes.txt') # TODO: Delete this line


    # Select only chromosomes selected by --chr
    if CHRS is not None:
        alignments, chrs, chrgrps, chrlengths = selectchrom(CHRS, cs, chrgrps, alignments, chrlengths, chrids)


    if REG is not None:
        alignments, chrs, chrgrps = selectregion(REG, RTR, chrlengths, alignments, chrids)

    # Combine Ribbon is selected than combine rows
    if R:
        for i in range(len(alignments)):
            alignments[i][1] = createribbon(alignments[i][1])

    # invert coord for inverted query genome
    for i in range(len(alignments)):
        df = alignments[i][1].copy()
        invindex = ['INV' in i for i in df['type']]
        g = set(df.loc[invindex, 'bstart'] < df.loc[invindex, 'bend'])
        if len(g) == 2:
            logger.error("Inconsistent coordinates in input file {}. For INV, INVTR, INVDUP annotations, either bstart < bend for all annotations or bstart > bend for all annotations. Mixing is not permitted.".format(alignments[i][0]))
            sys.exit()
        elif False in g:
            continue
        df.loc[invindex, 'bstart'] = df.loc[invindex, 'bstart'] + df.loc[invindex, 'bend']
        df.loc[invindex, 'bend'] = df.loc[invindex, 'bstart'] - df.loc[invindex, 'bend']
        df.loc[invindex, 'bstart'] = df.loc[invindex, 'bstart'] - df.loc[invindex, 'bend']
        alignments[i][1] = df.copy()


    # from matplotlib import pyplot as plt
    # plt = matplotlib.pyplot
    # plt.rcParams['font.size'] = FS
    # try:
    #     if H is None and W is None:
    #         H = len(chrs)
    #         W = 3
    #         fig = plt.figure(figsize=[W, H])
    #     elif H is not None and W is None:
    #         fig = plt.figure(figsize=[H, H])
    #     elif H is None and W is not None:
    #         fig = plt.figure(figsize=[W, W])
    #     else:
    #         fig = plt.figure(figsize=[W, H])
    # except Exception as e:
    #     logger.error("Error in initiliazing figure. Try using a different backend.\n{}".format(e.with_traceback()))
    #     sys.exit()
    # ax = fig.add_subplot(111, frameon=False)

    allal = pdconcat([alignments[i][1] for i in range(len(alignments))])
    if ITX:
        minl = 0
        MCHR = 0.01     # TODO : read spacing between neighbouring chromosome from config file
        maxchr = max([sum(chrlengths[i][1].values()) for i in range(len(chrlengths))])
        maxl = int(maxchr/(MCHR + 1 - (MCHR*len(chrgrps))))
    elif REG is None:
        minl, maxl = 0, -1
    else:
        minl = min(allal[['astart', 'bstart']].apply(min))
        maxl = max(allal[['aend', 'bend']].apply(max))
    labelcnt = 0
    if 'SYN' in allal['type'].array:
        labelcnt += 1
    if 'INV' in allal['type'].array:
        labelcnt += 1
    if 'TRA' in allal['type'].array or 'INVTR' in allal['type'].array:
        labelcnt += 1
    if 'DUP' in allal['type'].array or 'INVDP' in allal['type'].array:
        labelcnt += 1

    ## Draw Axes
    ax = drawax(ax, chrgrps, chrlengths, V, S, cfg, ITX, minl=minl, maxl=maxl, chrname=CHRNAME)

    ## Draw Chromosomes
    ax, indents, chrlabels = pltchrom(ax, chrs, chrgrps, chrlengths, V, S, genomes, cfg, ITX, minl=minl, maxl=maxl)

    if cfg['genlegcol'] < 1:
        ncol = ceil(len(chrlengths)/labelcnt)
    else:
        ncol = int(cfg['genlegcol'])

    # Get Genome legend
    if cfg['legend']:
        bbox_to_anchor = cfg['bbox']
        if not ITX:
            l1 = plt.legend(handles=chrlabels, loc='lower left', bbox_to_anchor=bbox_to_anchor, ncol=ncol, mode=None, borderaxespad=0., frameon=False, title='Genomes')
            l1._legend_box.align = "left"
            plt.gca().add_artist(l1)

    # Plot structural annotations
    ax, svlabels = pltsv(ax, alignments, chrs, V, chrgrps, chrlengths, indents, S, cfg, ITX, maxl)

    if cfg['legend']:
        bbox_to_anchor[0] += cfg['bboxmar']
        plt.legend(handles=svlabels, loc='lower left', bbox_to_anchor=bbox_to_anchor, ncol=1, mode='expand', borderaxespad=0., frameon=False, title='Annotations')._legend_box.align = "left"

    # Plot markers
    if B is not None:
        ax = drawmarkers(ax, B, V, chrlengths, indents, chrs, chrgrps, S, cfg, ITX, minl=minl, maxl=maxl)

    # Draw tracks
    if TRACKS is not None:
        tracks = readtrack(TRACKS, chrlengths)
        # tracks = readtrack(f, chrlengths) #TODO: delete this
        ax = drawtracks(ax, tracks, S, chrgrps, chrlengths, V, ITX, cfg, minl=minl, maxl=maxl)
    #
    # # Save the plot
    # try:
    #     fig.savefig(O, dpi=D, bbox_inches='tight', pad_inches=0.01)
    #     logger.info("Plot {O} generated.".format(O=O))
    # except Exception as e:
    #     sys.exit('Error in saving the figure. Try using a different backend.' + '\n' + e.with_traceback())
    # logger.info('Finished')

    return ax


def loghist(x, ax, bins=10):
    """
    Generate the plot on a logarithmic histogram
    """
    import numpy as np
    hist, bins = np.histogram(x, bins=bins)
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    ax.hist(x, bins=logbins)
    ax.set_xscale('log')
    return ax
# END