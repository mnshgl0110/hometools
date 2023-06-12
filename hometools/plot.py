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

